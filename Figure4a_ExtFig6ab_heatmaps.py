
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pyBigWig

# =========================================================
# REGION
# =========================================================
region_chrom = "chr7"
region_start= 118935704
region_end=118966473
bw_bin_size = 10
bw_dir = "/Volumes/jo5/rep123atmiatacrerun/bwa/merged_library/bigwig"
gRNA_start = 118_951_249
gRNA_end   = 118_951_268
gRNA_center = (gRNA_start + gRNA_end) / 2


# BIGWIG FILES
# =========================================================

bw_files = {
    "DMSO_sgCTRL": [
        "DMSO_sgCTRL_REP1.mLb.clN.bigWig",
        "DMSO_sgCTRL_REP2.mLb.clN.bigWig",
        "DMSO_sgCTRL_REP3.mLb.clN.bigWig",
    ],
    "DMSO_sgLAD": [
        "DMSO_sgLAD_REP1.mLb.clN.bigWig",
        "DMSO_sgLAD_REP2.mLb.clN.bigWig",
        "DMSO_sgLAD_REP3.mLb.clN.bigWig",
    ],
    "ATMi_sgCTRL": [
        "ATMi_sgCTRL_REP1.mLb.clN.bigWig",
        "ATMi_sgCTRL_REP2.mLb.clN.bigWig",
        "ATMi_sgCTRL_REP3.mLb.clN.bigWig",
    ],
    "ATMi_sgLAD": [
        "ATMi_sgLAD_REP1.mLb.clN.bigWig",
        "ATMi_sgLAD_REP2.mLb.clN.bigWig",
        "ATMi_sgLAD_REP3.mLb.clN.bigWig",
    ],
}

# Add full paths
bw_files = {
    group: [os.path.join(bw_dir, f) for f in files]
    for group, files in bw_files.items()
}

# FUNCTIONS
# =========================================================

def extract_bigwig_region_binned(bw_path, chrom, start, end, bin_size=10):
    bw = pyBigWig.open(bw_path)

    starts = np.arange(start, end, bin_size)
    ends = np.minimum(starts + bin_size, end)

    values = []
    for s, e in zip(starts, ends):
        v = bw.stats(chrom, int(s), int(e), type="mean")[0]
        values.append(0.0 if v is None else v)

    bw.close()

    return pd.DataFrame({
        "chrom": chrom,
        "start": starts,
        "end": ends,
        "center": (starts + ends) / 2,
        "signal": values
    })


def extract_group_replicates(group_name, bw_paths):
    dfs = []

    for i, bw_path in enumerate(bw_paths, start=1):
        df = extract_bigwig_region_binned(
            bw_path,
            region_chrom,
            region_start,
            region_end,
            bin_size=bw_bin_size
        )

        df = df.rename(columns={"signal": f"REP{i}"})
        dfs.append(df)

    merged = dfs[0]

    for df in dfs[1:]:
        merged = merged.merge(
            df[["chrom", "start", "end", "center", df.columns[-1]]],
            on=["chrom", "start", "end", "center"],
            how="inner"
        )

    rep_cols = [c for c in merged.columns if c.startswith("REP")]
    merged["mean"] = merged[rep_cols].mean(axis=1)

    return merged

# LOAD ALL GROUPS
# =========================================================

track_dfs = {
    group: extract_group_replicates(group, paths)
    for group, paths in bw_files.items()
}

plot_order = [
    "DMSO_sgCTRL",
    "DMSO_sgLAD",
    "ATMi_sgCTRL",
    "ATMi_sgLAD",
]

colors = {
    "DMSO_sgCTRL": "grey",
    "DMSO_sgLAD": "black",
    "ATMi_sgCTRL": "tomato",
    "ATMi_sgLAD": "firebrick",
}
##################BASE PLOT



fig, ax0 = plt.subplots(figsize=(20, 3), dpi=200)

for group in plot_order:
    df = track_dfs[group]
    rep_cols = [c for c in df.columns if c.startswith("REP")]

    for rep in rep_cols:
        ax0.plot(
            df["center"],
            df[rep],
            color=colors[group],
            linewidth=0.8,
            alpha=0.18
        )

    ax0.plot(
        df["center"],
        df["mean"],
        color=colors[group],
        linewidth=1.8,
        alpha=0.95,
        label=group
    )

ax0.axvline(
    gRNA_center,
    color="black",
    linestyle=":",
    linewidth=1.4,
    alpha=0.5
)


ax0.spines["top"].set_visible(False)
ax0.spines["right"].set_visible(False)
ax0.legend(frameon=False, loc="upper right", fontsize=9)

ax0.set_xlabel(f"{region_chrom}:{region_start:,}-{region_end:,}")
ax0.set_ylabel("ATAC signal")

plt.tight_layout()
plt.show()


######################################################################################################################################################
#above was what to run to check that the data is loaded in correctly, compare with IGV etc...if all good, proceed!
#Below here is the script to plot a heatmap
########################################################################################################################################################

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

# =========================================================
# gRNA POSITION
# =========================================================

gRNA_start = 118_951_249
gRNA_end   = 118_951_268
gRNA_center = (gRNA_start + gRNA_end) / 2


# =========================================================
# PLOT ORDER
# =========================================================

plot_order = [
    "DMSO_sgCTRL",
    "DMSO_sgLAD",
    "ATMi_sgCTRL",
    "ATMi_sgLAD",
]

# Optional: nicer labels for plotting
plot_labels = {
    "DMSO_sgCTRL": "DMSO sgCTRL",
    "DMSO_sgLAD": "DMSO sgLAD",
    "ATMi_sgCTRL": "ATMi sgCTRL",
    "ATMi_sgLAD": "ATMi sgLAD",
}


# =========================================================
# BUILD SHARED COLOUR SCALE
# =========================================================

all_values = []

for group in plot_order:
    df = track_dfs[group]
    rep_cols = [c for c in df.columns if c.startswith("REP")]
    all_values.append(df[rep_cols].to_numpy())

all_values = np.concatenate(all_values, axis=None)

vmin = 0
#vmax = np.nanpercentile(all_values, 99)   # avoids one extreme bin dominating
vmax= 0.5

# =========================================================
# HEATMAP PLOT
# =========================================================

fig = plt.figure(figsize=(20, 7), dpi=200)

gs = fig.add_gridspec(
    nrows=len(plot_order),
    ncols=2,
    width_ratios=[40, 1],
    height_ratios=[1, 1, 1, 1],
    hspace=0.18,
    wspace=0.03
)

heatmap_axes = []

for i, group in enumerate(plot_order):
    ax = fig.add_subplot(gs[i, 0])
    heatmap_axes.append(ax)

    df = track_dfs[group]
    rep_cols = [c for c in df.columns if c.startswith("REP")]

    # rows = replicates, columns = genomic bins
    mat = df[rep_cols].to_numpy().T

    im = ax.imshow(
        mat,
        aspect="auto",
        interpolation="nearest",
        cmap="magma",
        vmin=vmin,
        vmax=vmax,
        extent=[
            region_start,
            region_end,
            len(rep_cols) + 0.5,
            0.5
        ]
    )

    # gRNA centre dotted line
    ax.axvline(
        gRNA_center,
        color="white",
        linestyle=":",
        linewidth=1.6,
        alpha=0.5
    )

    ax.set_yticks([1, 2, 3])
    ax.set_yticklabels(["R1", "R2", "R3"], fontsize=10)

    ax.set_ylabel(
        plot_labels[group],
        rotation=0,
        ha="right",
        va="center",
        fontsize=11
    )

    ax.tick_params(axis="x", labelbottom=False)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)

# only bottom panel gets x-axis labels
heatmap_axes[-1].tick_params(axis="x", labelbottom=True)
heatmap_axes[-1].set_xlabel(
    f"{region_chrom}:{region_start:,}-{region_end:,}",
    fontsize=11
)

# colourbar
cax = fig.add_subplot(gs[:, 1])
cbar = fig.colorbar(im, cax=cax)
cbar.set_label("Signal", rotation=270, labelpad=15)

plt.show()

#%%%
##########################################################################################################################

#s139a
###########################################################################################################################



region_chrom = "chr7"
region_start = 118943204
region_end   = 118958973

#region_start = 118950983
#region_end   = 118951568

bw_bin_size = 10

# Number of 10 bp bins used for smoothing
# 25 bins = 250 bp smoothing window
smoothing_window_bins = 25

bw_dir = "/Volumes/zoe/s139a_out_swapcorrected/bwa/merged_library/bigwig/"


gRNA_start = 118_951_249
gRNA_end   = 118_951_268
gRNA_center = (gRNA_start + gRNA_end) / 2


# =========================================================
# BIGWIG FILES
# =========================================================

bw_files = {
    "WT_sgCTRL": [
        "WT_sgCTRL_REP1.mLb.clN.bigWig",
        "WT_sgCTRL_REP2.mLb.clN.bigWig",
        "WT_sgCTRL_REP3.mLb.clN.bigWig",
    ],
    "WT_sgLAD": [
        "WT_sgLAD_REP1.mLb.clN.bigWig",
        "WT_sgLAD_REP2.mLb.clN.bigWig",
        "WT_sgLAD_REP3.mLb.clN.bigWig",
    ],
    "SA_sgCTRL": [
        "SA_sgCTRL_REP1.mLb.clN.bigWig",
        "SA_sgCTRL_REP2.mLb.clN.bigWig",
        "SA_sgCTRL_REP3.mLb.clN.bigWig",
    ],
    "SA_sgLAD": [
        "SA_sgLAD_REP1.mLb.clN.bigWig",
        "SA_sgLAD_REP2.mLb.clN.bigWig",
        "SA_sgLAD_REP3.mLb.clN.bigWig",
    ],
}


# Add full paths
bw_files = {
    group: [os.path.join(bw_dir, f) for f in files]
    for group, files in bw_files.items()
}


# =========================================================
# FUNCTIONS
# =========================================================

def extract_bigwig_region_binned(
    bw_path,
    chrom,
    start,
    end,
    bin_size=10
):
    bw = pyBigWig.open(bw_path)

    starts = np.arange(start, end, bin_size)
    ends = np.minimum(starts + bin_size, end)

    values = []

    for s, e in zip(starts, ends):
        value = bw.stats(
            chrom,
            int(s),
            int(e),
            type="mean"
        )[0]

        values.append(
            0.0 if value is None else value
        )

    bw.close()

    return pd.DataFrame({
        "chrom": chrom,
        "start": starts,
        "end": ends,
        "center": (starts + ends) / 2,
        "signal": values
    })


def smooth_signal(
    values,
    window_bins=20,
    center=True
):
    """
    Smooth a one-dimensional signal using a rolling mean.

    Parameters
    ----------
    values : array-like
        Signal values to smooth.

    window_bins : int
        Number of genomic bins in the smoothing window.

    center : bool
        Centre the rolling window on each genomic position.

    Returns
    -------
    pandas.Series
        Smoothed signal.
    """

    return (
        pd.Series(values, dtype=float)
        .rolling(
            window=window_bins,
            center=center,
            min_periods=1
        )
        .mean()
    )


def extract_group_replicates(group_name, bw_paths):
    dfs = []

    for i, bw_path in enumerate(bw_paths, start=1):
        df = extract_bigwig_region_binned(
            bw_path,
            region_chrom,
            region_start,
            region_end,
            bin_size=bw_bin_size
        )

        df = df.rename(
            columns={"signal": f"REP{i}"}
        )

        dfs.append(df)

    merged = dfs[0]

    for df in dfs[1:]:
        signal_column = [
            column
            for column in df.columns
            if column.startswith("REP")
        ][0]

        merged = merged.merge(
            df[
                [
                    "chrom",
                    "start",
                    "end",
                    "center",
                    signal_column
                ]
            ],
            on=[
                "chrom",
                "start",
                "end",
                "center"
            ],
            how="inner"
        )

    rep_cols = [
        column
        for column in merged.columns
        if column.startswith("REP")
    ]

    # Raw mean across replicates
    merged["mean"] = merged[rep_cols].mean(axis=1)

    # Smooth each replicate separately
    for rep in rep_cols:
        merged[f"{rep}_smooth"] = smooth_signal(
            merged[rep],
            window_bins=smoothing_window_bins
        )

    smooth_rep_cols = [
        f"{rep}_smooth"
        for rep in rep_cols
    ]

    # Mean of the smoothed replicate traces
    merged["mean_smooth"] = merged[
        smooth_rep_cols
    ].mean(axis=1)

    return merged


# =========================================================
# LOAD ALL GROUPS
# =========================================================

track_dfs = {
    group: extract_group_replicates(group, paths)
    for group, paths in bw_files.items()
}


gRNA_start = 118_951_249
gRNA_end   = 118_951_268
gRNA_center = (gRNA_start + gRNA_end) / 2


# =========================================================
# PLOT ORDER
# =========================================================

plot_order = [
    "WT_sgCTRL",
    "WT_sgLAD",
    "SA_sgCTRL",
    "SA_sgLAD",
]

plot_labels = {
    "WT_sgCTRL": "WT sgCTRL",
    "WT_sgLAD": "WT sgLAD",
    "SA_sgCTRL": "H2AXS139A sgCTRL",
    "SA_sgLAD": "H2AXS139A sgLAD",
}

replicate_order = ["REP1", "REP2", "REP3"]

replicate_labels = {
    "REP1": "R1",
    "REP2": "R2",
    "REP3": "R3",
}


# =========================================================
# HEATMAP PLOT
# =========================================================

fig = plt.figure(
    figsize=(20, 9),
    dpi=200
)

gs = fig.add_gridspec(
    nrows=3,
    ncols=2,
    width_ratios=[40, 1],
    height_ratios=[1, 1, 1],
    hspace=0.22,
    wspace=0.03
)

heatmap_axes = []

for i, rep in enumerate(replicate_order):

    ax = fig.add_subplot(gs[i, 0])
    cax = fig.add_subplot(gs[i, 1])

    heatmap_axes.append(ax)

    # Each row is one sample for the current replicate
    mat = np.vstack([
        track_dfs[group][rep].to_numpy()
        for group in plot_order
    ])

    # Independent colour scale for this replicate
    rep_values = mat[np.isfinite(mat)]

    vmin = 0
    vmax = np.nanpercentile(rep_values, 99.5)

    im = ax.imshow(
        mat,
        aspect="auto",
        interpolation="nearest",
        cmap="magma",
        vmin=vmin,
        vmax=vmax,
        extent=[
            region_start,
            region_end,
            len(plot_order) + 0.5,
            0.5
        ]
    )

    # Individual colourbar for this replicate
    cbar = fig.colorbar(
        im,
        cax=cax
    )

    cbar.set_label(
        "Signal",
        rotation=270,
        labelpad=15
    )

    cbar.ax.tick_params(labelsize=9)

    # gRNA centre dotted line
    ax.axvline(
        gRNA_center,
        color="white",
        linestyle=":",
        linewidth=1.6,
        alpha=0.9
    )

    ax.set_yticks(
        np.arange(1, len(plot_order) + 1)
    )

    ax.set_yticklabels(
        [plot_labels[group] for group in plot_order],
        fontsize=10
    )

    ax.set_ylabel(
        replicate_labels[rep],
        rotation=0,
        ha="right",
        va="center",
        fontsize=12,
        fontweight="bold",
        labelpad=25
    )

    ax.tick_params(
        axis="x",
        labelbottom=False
    )

    ax.tick_params(
        axis="y",
        length=0,
        pad=6
    )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)


# Only bottom panel gets x-axis labels
heatmap_axes[-1].tick_params(
    axis="x",
    labelbottom=True
)

heatmap_axes[-1].set_xlabel(
    f"{region_chrom}:{region_start:,}-{region_end:,}",
    fontsize=11
)

plt.show()

##########################################################################################################################

#53BP1ko
###########################################################################################################################



region_chrom = "chr7"
region_start = 118943204
region_end   = 118958973

#region_start = 118950983
#region_end   = 118951568
bw_bin_size = 10

bw_dir = '/Volumes/zoe/s139a_out_swapcorrected/bwa/merged_library/bigwig/'


gRNA_start = 118_951_249
gRNA_end   = 118_951_268
gRNA_center = (gRNA_start + gRNA_end) / 2
# =========================================================
# BIGWIG FILES
# =========================================================

bw_files = {
    "WT_sgCTRL": [
        "/Volumes/zoe/s139a_out_swapcorrected/bwa/merged_library/bigwig/WT_sgCTRL_REP3.mLb.clN.bigWig",
        "/Volumes/zoe/53bp1r4r5_OUT/bwa/merged_library/bigwig/WT_sgCTRL_REP1.mLb.clN.bigWig",
        "/Volumes/zoe/53bp1r4r5_OUT/bwa/merged_library/bigwig/WT_sgCTRL_REP2.mLb.clN.bigWig",
    ],
    "WT_sgLAD": [
        "/Volumes/zoe/s139a_out_swapcorrected/bwa/merged_library/bigwig/WT_sgLAD_REP3.mLb.clN.bigWig",
        "/Volumes/zoe/53bp1r4r5_OUT/bwa/merged_library/bigwig/WT_sgLAD_REP1.mLb.clN.bigWig",
        "/Volumes/zoe/53bp1r4r5_OUT/bwa/merged_library/bigwig/WT_sgLAD_REP2.mLb.clN.bigWig",
    ],
    "53BP1ko_sgCTRL": [
        "/Volumes/zoe/s139a_out_swapcorrected/bwa/merged_library/bigwig/53ko_sgCTRL_REP1.mLb.clN.bigWig",
        "/Volumes/zoe/53bp1r4r5_OUT/bwa/merged_library/bigwig/53ko_sgCTRL_REP1.mLb.clN.bigWig",
        "/Volumes/zoe/53bp1r4r5_OUT/bwa/merged_library/bigwig/53ko_sgLAD_REP2.mLb.clN.bigWig",
    ],
    "53BP1ko_sgLAD": [
        "/Volumes/zoe/s139a_out_swapcorrected/bwa/merged_library/bigwig/53ko_sgLAD_REP1.mLb.clN.bigWig",
        "/Volumes/zoe/53bp1r4r5_OUT/bwa/merged_library/bigwig/53ko_sgLAD_REP1.mLb.clN.bigWig",
        "/Volumes/zoe/53bp1r4r5_OUT/bwa/merged_library/bigwig/53ko_sgCTRL_REP2.mLb.clN.bigWig",
    ],
}

# Add full paths
bw_files = {
    group: [os.path.join(bw_dir, f) for f in files]
    for group, files in bw_files.items()
}

# =========================================================
# FUNCTIONS
# =========================================================

def extract_bigwig_region_binned(bw_path, chrom, start, end, bin_size=10):
    bw = pyBigWig.open(bw_path)

    starts = np.arange(start, end, bin_size)
    ends = np.minimum(starts + bin_size, end)

    values = []
    for s, e in zip(starts, ends):
        v = bw.stats(chrom, int(s), int(e), type="mean")[0]
        values.append(0.0 if v is None else v)

    bw.close()

    return pd.DataFrame({
        "chrom": chrom,
        "start": starts,
        "end": ends,
        "center": (starts + ends) / 2,
        "signal": values
    })


def extract_group_replicates(group_name, bw_paths):
    dfs = []

    for i, bw_path in enumerate(bw_paths, start=1):
        df = extract_bigwig_region_binned(
            bw_path,
            region_chrom,
            region_start,
            region_end,
            bin_size=bw_bin_size
        )

        df = df.rename(columns={"signal": f"REP{i}"})
        dfs.append(df)

    merged = dfs[0]

    for df in dfs[1:]:
        merged = merged.merge(
            df[["chrom", "start", "end", "center", df.columns[-1]]],
            on=["chrom", "start", "end", "center"],
            how="inner"
        )

    rep_cols = [c for c in merged.columns if c.startswith("REP")]
    merged["mean"] = merged[rep_cols].mean(axis=1)

    return merged


# =========================================================
# LOAD ALL GROUPS
# =========================================================

track_dfs = {
    group: extract_group_replicates(group, paths)
    for group, paths in bw_files.items()
}




# =========================================================
# PLOT ORDER
# =========================================================

plot_order = [
    "WT_sgCTRL",
    "WT_sgLAD",
    "53BP1ko_sgCTRL",
    "53BP1ko_sgLAD",
]

plot_labels = {
    "WT_sgCTRL": "WT sgCTRL",
    "WT_sgLAD": "WT sgLAD",
    "53BP1ko_sgCTRL": "53BP1ko sgCTRL",
    "53BP1ko_sgLAD": "53BP1ko sgLAD",
}

replicate_order = ["REP1", "REP2", "REP3"]

replicate_labels = {
    "REP1": "R1",
    "REP2": "R2",
    "REP3": "R3",
}


# =========================================================
# HEATMAP PLOT
# =========================================================

fig = plt.figure(
    figsize=(20, 9),
    dpi=200
)

gs = fig.add_gridspec(
    nrows=3,
    ncols=2,
    width_ratios=[40, 1],
    height_ratios=[1, 1, 1],
    hspace=0.22,
    wspace=0.03
)

heatmap_axes = []

for i, rep in enumerate(replicate_order):

    ax = fig.add_subplot(gs[i, 0])
    cax = fig.add_subplot(gs[i, 1])

    heatmap_axes.append(ax)

    # Each row is one sample for the current replicate
    mat = np.vstack([
        track_dfs[group][rep].to_numpy()
        for group in plot_order
    ])

    # Independent colour scale for this replicate
    rep_values = mat[np.isfinite(mat)]

    vmin = 0
    vmax = np.nanpercentile(rep_values, 99.5)

    im = ax.imshow(
        mat,
        aspect="auto",
        interpolation="nearest",
        cmap="magma",
        vmin=vmin,
        vmax=vmax,
        extent=[
            region_start,
            region_end,
            len(plot_order) + 0.5,
            0.5
        ]
    )

    # Individual colourbar for this replicate
    cbar = fig.colorbar(
        im,
        cax=cax
    )

    cbar.set_label(
        "Signal",
        rotation=270,
        labelpad=15
    )

    cbar.ax.tick_params(labelsize=9)

    # gRNA centre dotted line
    ax.axvline(
        gRNA_center,
        color="white",
        linestyle=":",
        linewidth=1.6,
        alpha=0.9
    )

    ax.set_yticks(
        np.arange(1, len(plot_order) + 1)
    )

    ax.set_yticklabels(
        [plot_labels[group] for group in plot_order],
        fontsize=10
    )

    ax.set_ylabel(
        replicate_labels[rep],
        rotation=0,
        ha="right",
        va="center",
        fontsize=12,
        fontweight="bold",
        labelpad=25
    )

    ax.tick_params(
        axis="x",
        labelbottom=False
    )

    ax.tick_params(
        axis="y",
        length=0,
        pad=6
    )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)


# Only bottom panel gets x-axis labels
heatmap_axes[-1].tick_params(
    axis="x",
    labelbottom=True
)

heatmap_axes[-1].set_xlabel(
    f"{region_chrom}:{region_start:,}-{region_end:,}",
    fontsize=11
)

plt.show()




