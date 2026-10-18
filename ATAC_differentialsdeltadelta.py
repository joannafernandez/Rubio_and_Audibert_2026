#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 15:05:29 2026

@author: patricfernandez
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyBigWig


# =========================================================
# REGION SETTINGS
# =========================================================

region_chrom = "chr7"
region_start = 118943204
region_end   = 118958973

bw_bin_size = 10

gRNA_start = 118_951_249
gRNA_end   = 118_951_268
gRNA_center = (gRNA_start + gRNA_end) / 2


# =========================================================
# SMOOTHING SETTINGS
# =========================================================

smoothing_bins = 20

# With bw_bin_size = 10:
# 20 bins × 10 bp = 200 bp smoothing


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


def extract_group_replicates(
    group_name,
    bw_paths
):
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

        rep_column = [
            c for c in df.columns
            if c.startswith("REP")
        ][0]

        merged = merged.merge(
            df[
                [
                    "chrom",
                    "start",
                    "end",
                    "center",
                    rep_column
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
        c for c in merged.columns
        if c.startswith("REP")
    ]

    merged["mean"] = merged[rep_cols].mean(axis=1)

    return merged


def load_track_set(bw_files):
    """
    Load all groups in one experimental dataset.
    """

    return {
        group: extract_group_replicates(group, paths)
        for group, paths in bw_files.items()
    }


def calculate_difference_track(
    track_dfs,
    lad_group,
    ctrl_group,
    smoothing_bins=20
):
    """
    Calculate sgLAD - sgCTRL using the mean signal
    across replicates.
    """

    lad_df = track_dfs[lad_group][
        [
            "chrom",
            "start",
            "end",
            "center",
            "mean"
        ]
    ].rename(
        columns={"mean": "sgLAD_mean"}
    )

    ctrl_df = track_dfs[ctrl_group][
        [
            "chrom",
            "start",
            "end",
            "center",
            "mean"
        ]
    ].rename(
        columns={"mean": "sgCTRL_mean"}
    )

    diff_df = lad_df.merge(
        ctrl_df,
        on=[
            "chrom",
            "start",
            "end",
            "center"
        ],
        how="inner"
    )

    diff_df["mean_difference"] = (
        diff_df["sgLAD_mean"]
        - diff_df["sgCTRL_mean"]
    )

    diff_df["smoothed_difference"] = (
        diff_df["mean_difference"]
        .rolling(
            window=smoothing_bins,
            center=True,
            min_periods=1
        )
        .mean()
    )

    return diff_df


# =========================================================
# DATASET 1: DMSO / ATMi
# =========================================================

atmi_bw_dir = (
    "/Volumes/jo5/rep123atmiatacrerun/bwa/merged_library/bigwig"
)

atmi_bw_files = {
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

atmi_bw_files = {
    group: [
        os.path.join(atmi_bw_dir, filename)
        for filename in filenames
    ]
    for group, filenames in atmi_bw_files.items()
}

atmi_track_dfs = load_track_set(atmi_bw_files)


# =========================================================
# DATASET 2: WT / H2AX S139A
# =========================================================

s139a_bw_dir = (
    "/Volumes/zoe/s139a_out_swapcorrected/bwa/merged_library/bigwig/"
)

s139a_bw_files = {
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

s139a_bw_files = {
    group: [
        os.path.join(s139a_bw_dir, filename)
        for filename in filenames
    ]
    for group, filenames in s139a_bw_files.items()
}

s139a_track_dfs = load_track_set(s139a_bw_files)


# =========================================================
# DATASET 3: WT / 53BP1 KO
# =========================================================

bp53_bw_files = {
    "WT_sgCTRL": [
        "/Volumes/zoe/s139a_out_swapcorrected/bwa/merged_library/bigwig/"
        "WT_sgCTRL_REP3.mLb.clN.bigWig",

        "/Volumes/zoe/53bp1r4r5_OUT/bwa/merged_library/bigwig/"
        "WT_sgCTRL_REP1.mLb.clN.bigWig",

        "/Volumes/zoe/53bp1r4r5_OUT/bwa/merged_library/bigwig/"
        "WT_sgCTRL_REP2.mLb.clN.bigWig",
    ],

    "WT_sgLAD": [
        "/Volumes/zoe/s139a_out_swapcorrected/bwa/merged_library/bigwig/"
        "WT_sgLAD_REP3.mLb.clN.bigWig",

        "/Volumes/zoe/53bp1r4r5_OUT/bwa/merged_library/bigwig/"
        "WT_sgLAD_REP1.mLb.clN.bigWig",

        "/Volumes/zoe/53bp1r4r5_OUT/bwa/merged_library/bigwig/"
        "WT_sgLAD_REP2.mLb.clN.bigWig",
    ],

    "53BP1ko_sgCTRL": [
        "/Volumes/zoe/s139a_out_swapcorrected/bwa/merged_library/bigwig/"
        "53ko_sgCTRL_REP1.mLb.clN.bigWig",

        "/Volumes/zoe/53bp1r4r5_OUT/bwa/merged_library/bigwig/"
        "53ko_sgCTRL_REP1.mLb.clN.bigWig",

        "/Volumes/zoe/53bp1r4r5_OUT/bwa/merged_library/bigwig/"
        "53ko_sgLAD_REP2.mLb.clN.bigWig",
    ],

    "53BP1ko_sgLAD": [
        "/Volumes/zoe/s139a_out_swapcorrected/bwa/merged_library/bigwig/"
        "53ko_sgLAD_REP1.mLb.clN.bigWig",

        "/Volumes/zoe/53bp1r4r5_OUT/bwa/merged_library/bigwig/"
        "53ko_sgLAD_REP1.mLb.clN.bigWig",

        "/Volumes/zoe/53bp1r4r5_OUT/bwa/merged_library/bigwig/"
        "53ko_sgCTRL_REP2.mLb.clN.bigWig",
    ],
}

bp53_track_dfs = load_track_set(bp53_bw_files)


# =========================================================
# CALCULATE DIFFERENCE TRACKS
# =========================================================

# Axis 1: ATMi
dmso_difference = calculate_difference_track(
    track_dfs=atmi_track_dfs,
    lad_group="DMSO_sgLAD",
    ctrl_group="DMSO_sgCTRL",
    smoothing_bins=smoothing_bins
)
atmi_difference = calculate_difference_track(
    track_dfs=atmi_track_dfs,
    lad_group="ATMi_sgLAD",
    ctrl_group="ATMi_sgCTRL",
    smoothing_bins=smoothing_bins
)


# Axis 2: WT and H2AX S139A
s139a_wt_difference = calculate_difference_track(
    track_dfs=s139a_track_dfs,
    lad_group="WT_sgLAD",
    ctrl_group="WT_sgCTRL",
    smoothing_bins=smoothing_bins
)

s139a_difference = calculate_difference_track(
    track_dfs=s139a_track_dfs,
    lad_group="SA_sgLAD",
    ctrl_group="SA_sgCTRL",
    smoothing_bins=smoothing_bins
)


# Axis 3: WT and 53BP1 KO
bp53_wt_difference = calculate_difference_track(
    track_dfs=bp53_track_dfs,
    lad_group="WT_sgLAD",
    ctrl_group="WT_sgCTRL",
    smoothing_bins=smoothing_bins
)

bp53ko_difference = calculate_difference_track(
    track_dfs=bp53_track_dfs,
    lad_group="53BP1ko_sgLAD",
    ctrl_group="53BP1ko_sgCTRL",
    smoothing_bins=smoothing_bins
)


# =========================================================
# THREE-AXIS DIFFERENTIAL PLOT
# =========================================================

fig, axes = plt.subplots(
    nrows=3,
    ncols=1,
    figsize=(20, 9),
    dpi=200,
    sharex=True
)

ax1, ax2, ax3 = axes


# ---------------------------------------------------------
# AXIS 1: ATMi
# ---------------------------------------------------------

ax1.plot(
    dmso_difference["center"],
    dmso_difference["smoothed_difference"],
    color="#75B1CE",
    linewidth=1.7,
    label="DMSO"
)

ax1.plot(
    atmi_difference["center"],
    atmi_difference["smoothed_difference"],
    color="#526C94",
    linewidth=1.7,
    label="ATMi"
)

ax1.set_ylabel(
    "ATMi\nsgLAD − sgCTRL",
    rotation=0,
    ha="right",
    va="center",
    fontsize=11,
    labelpad=20
)

ax1.legend(
    frameon=False,
    loc="upper right",
    fontsize=10
)


# ---------------------------------------------------------
# AXIS 2: WT versus H2AX S139A
# ---------------------------------------------------------

ax2.plot(
    s139a_wt_difference["center"],
    s139a_wt_difference["smoothed_difference"],
    color="#75B1CE",
    linewidth=1.7,
    label="WT"
)

ax2.plot(
    s139a_difference["center"],
    s139a_difference["smoothed_difference"],
    color="#526C94",
    linewidth=1.7,
    label="H2AX S139A"
)

ax2.set_ylabel(
    "H2AX genotype\nsgLAD − sgCTRL",
    rotation=0,
    ha="right",
    va="center",
    fontsize=11,
    labelpad=20
)

ax2.legend(
    frameon=False,
    loc="upper right",
    fontsize=10
)


# ---------------------------------------------------------
# AXIS 3: WT versus 53BP1 KO
# ---------------------------------------------------------

ax3.plot(
    bp53_wt_difference["center"],
    bp53_wt_difference["smoothed_difference"],
    color="#75B1CE",
    linewidth=1.7,
    label="WT"
)

ax3.plot(
    bp53ko_difference["center"],
    bp53ko_difference["smoothed_difference"],
    color="#526C94",
    linewidth=1.7,
    label="53BP1 KO"
)

ax3.set_ylabel(
    "53BP1 genotype\nsgLAD − sgCTRL",
    rotation=0,
    ha="right",
    va="center",
    fontsize=11,
    labelpad=20
)

ax3.legend(
    frameon=False,
    loc="upper right",
    fontsize=10
)


# =========================================================
# COMMON FORMATTING
# =========================================================

for ax in axes:

    # Zero reference line
    ax.axhline(
        0,
        color="black",
        linestyle="--",
        linewidth=1,
        alpha=0.65
    )

    # gRNA centre
    ax.axvline(
        gRNA_center,
        color="black",
        linestyle=":",
        linewidth=1.5,
        alpha=0.8
    )

    ax.set_xlim(
        region_start,
        region_end
    )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


ax3.set_xlabel(
    f"{region_chrom}:{region_start:,}-{region_end:,}",
    fontsize=11
)

plt.tight_layout()
plt.show()


#%%%

# =========================================================
# CALCULATE DIFFERENCE-OF-DIFFERENCE TRACKS
# =========================================================

def calculate_delta_of_differences(
    reference_df,
    comparison_df,
    reference_name,
    comparison_name
):
    """
    Calculate:

        reference (sgLAD - sgCTRL)
        -
        comparison (sgLAD - sgCTRL)

    Both input DataFrames should have been generated by
    calculate_difference_track().
    """

    ref = reference_df[
        [
            "chrom",
            "start",
            "end",
            "center",
            "smoothed_difference"
        ]
    ].rename(
        columns={
            "smoothed_difference": f"{reference_name}_difference"
        }
    )

    comp = comparison_df[
        [
            "chrom",
            "start",
            "end",
            "center",
            "smoothed_difference"
        ]
    ].rename(
        columns={
            "smoothed_difference": f"{comparison_name}_difference"
        }
    )

    delta_df = ref.merge(
        comp,
        on=[
            "chrom",
            "start",
            "end",
            "center"
        ],
        how="inner"
    )

    delta_df["delta_difference"] = (
        delta_df[f"{reference_name}_difference"]
        - delta_df[f"{comparison_name}_difference"]
    )

    return delta_df


# =========================================================
# DIFFERENCE-OF-DIFFERENCES
# =========================================================

# ATMi(sgLAD - sgCTRL) - DMSO(sgLAD - sgCTRL)
dmso_vs_atmi_delta = calculate_delta_of_differences(
    reference_df=atmi_difference,
    comparison_df=dmso_difference,
    reference_name="ATMi",
    comparison_name="DMSO"
)


# S139A(sgLAD - sgCTRL) - WT(sgLAD - sgCTRL)
wt_vs_s139a_delta = calculate_delta_of_differences(
    reference_df=s139a_difference,
    comparison_df=s139a_wt_difference,
    reference_name="S139A",
    comparison_name="WT"
)


# 53BP1ko(sgLAD - sgCTRL) - WT(sgLAD - sgCTRL)
wt_vs_53bp1ko_delta = calculate_delta_of_differences(
    reference_df=bp53ko_difference,
    comparison_df=bp53_wt_difference,
    reference_name="53BP1KO",
    comparison_name="WT"
)


# =========================================================
# 3 ROW × 2 COLUMN PLOT
# =========================================================

fig, axes = plt.subplots(
    nrows=3,
    ncols=2,
    figsize=(20, 10),
    dpi=200,
    sharex=True,
    gridspec_kw={
        "width_ratios": [1, 1]
    }
)

ax1, ax1_delta = axes[0]
ax2, ax2_delta = axes[1]
ax3, ax3_delta = axes[2]


# =========================================================
# ROW 1: DMSO versus ATMi
# =========================================================

# Left: individual sgLAD - sgCTRL tracks
ax1.plot(
    dmso_difference["center"],
    dmso_difference["smoothed_difference"],
    color="#75B1CE",
    linewidth=1.5,
    label="DMSO"
)

ax1.plot(
    atmi_difference["center"],
    atmi_difference["smoothed_difference"],
    color="#526C94",
    linewidth=1.5,
    label="ATMi"
)

ax1.set_ylabel(
    "ATM inhibition\nsgLAD − sgCTRL",
    rotation=0,
    ha="right",
    va="center",
    fontsize=3,
    labelpad=25
)

ax1.legend(
    frameon=False,
    loc="upper right",
    fontsize=3
)


# Right: DMSO difference - ATMi difference
ax1_delta.plot(
    dmso_vs_atmi_delta["center"],
    dmso_vs_atmi_delta["delta_difference"],
    color="#DC6B83",
    linewidth=1.5
)

ax1_delta.set_ylabel(
    "ATMi - DMSO",
    fontsize=3
)

ax1_delta.set_title(
    "Difference of differences",
    fontsize=3
)


# =========================================================
# ROW 2: WT versus H2AX S139A
# =========================================================

# Left: individual sgLAD - sgCTRL tracks
ax2.plot(
    s139a_wt_difference["center"],
    s139a_wt_difference["smoothed_difference"],
    color="#75B1CE",
    linewidth=1.5,
    label="WT"
)

ax2.plot(
    s139a_difference["center"],
    s139a_difference["smoothed_difference"],
    color="#526C94",
    linewidth=1.5,
    label="H2AX S139A"
)

ax2.set_ylabel(
    "H2AX genotype\nsgLAD − sgCTRL",
    rotation=0,
    ha="right",
    va="center",
    fontsize=3,
    labelpad=3
)

ax2.legend(
    frameon=False,
    loc="upper right",
    fontsize=3
)


# Right: WT difference - S139A difference
ax2_delta.plot(
    wt_vs_s139a_delta["center"],
    wt_vs_s139a_delta["delta_difference"],
    color="#DC6B83",
    linewidth=1.5
)

ax2_delta.set_ylabel(
    "H2AX S139A - WT",
    fontsize=3
)


# =========================================================
# ROW 3: WT versus 53BP1 KO
# =========================================================

# Left: individual sgLAD - sgCTRL tracks
ax3.plot(
    bp53_wt_difference["center"],
    bp53_wt_difference["smoothed_difference"],
    color="#75B1CE",
    linewidth=1.5,
    label="WT"
)

ax3.plot(
    bp53ko_difference["center"],
    bp53ko_difference["smoothed_difference"],
    color="#526C94",
    linewidth=1.5,
    label="53BP1 KO"
)

ax3.set_ylabel(
    "53BP1 genotype\nsgLAD − sgCTRL",
    rotation=0,
    ha="right",
    va="center",
    fontsize=3,
    labelpad=3
)

ax3.legend(
    frameon=False,
    loc="upper right",
    fontsize=3
)


# Right: WT difference - 53BP1 KO difference
ax3_delta.plot(
    wt_vs_53bp1ko_delta["center"],
    wt_vs_53bp1ko_delta["delta_difference"],
    color="#DC6B83",
    linewidth=1.5
)

ax3_delta.set_ylabel(
    "53BP1 KO - WT",
    fontsize=3
)

# =========================================================
# CENTRE ZERO ON DIFFERENCE-OF-DIFFERENCES AXES
# =========================================================

delta_axes_and_dfs = [
    (ax1_delta, dmso_vs_atmi_delta),
    (ax2_delta, wt_vs_s139a_delta),
    (ax3_delta, wt_vs_53bp1ko_delta),
]

for ax, delta_df in delta_axes_and_dfs:

    # Largest absolute value in this difference track
    max_abs = np.nanmax(
        np.abs(delta_df["delta_difference"].to_numpy())
    )

    # Add a little vertical padding
    y_limit = max_abs * 1.10

    # Equal negative and positive limits place zero in the centre
    ax.set_ylim(-y_limit, y_limit)
# =========================================================
# COMMON FORMATTING
# =========================================================

for ax in axes.flat:

    # Zero reference line
    ax.axhline(
        0,
        color="black",
        linestyle="--",
        linewidth=1,
        alpha=0.65
    )

    # gRNA centre
    ax.axvline(
        gRNA_center,
        color="black",
        linestyle=":",
        linewidth=1.5,
        alpha=0.8
    )

    ax.set_xlim(
        region_start,
        region_end
    )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


# Column titles
ax1.set_title(
    "Within-condition difference",
    fontsize=3
)


# Bottom x-axis labels
ax3.set_xlabel(
    f"{region_chrom}:{region_start:,}-{region_end:,}",
    fontsize=3
)

ax3_delta.set_xlabel(
    f"{region_chrom}:{region_start:,}-{region_end:,}",
    fontsize=3
)


plt.tight_layout()
plt.show()