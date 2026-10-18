#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 10:12:54 2026

@author: patricfernandez
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pyBigWig

path = "/Volumes/jo4/laminnospike/03_peak_calling/03_bed_to_bigwig/derived_bigwigs/DMSO_5kb_rawcounts.tsv"

df = pd.read_csv(path, sep="\t")
df.columns = [c.replace("'", "").lstrip("#") for c in df.columns]
df = df.rename(columns={"chr": "chrom"})
df = df.fillna(0)

df.head()

df = df.fillna(0)

df["delta_R1"] = df["LB1_DMSO_LAD_R1.adj_log2.bw"] - df["LB1_DMSO_CTRL_R1.adj_log2.bw"]
df["delta_R2"] = df["LB1_DMSO_LAD_R2.adj_log2.bw"] - df["LB1_DMSO_CTRL_R2.adj_log2.bw"]
df["delta_R3"] = df["LB1_DMSO_LAD_R3.adj_log2.bw"] - df["LB1_DMSO_CTRL_R3.adj_log2.bw"]

df["mean_LAD"] = df[
    ["LB1_DMSO_LAD_R1.adj_log2.bw",
     "LB1_DMSO_LAD_R2.adj_log2.bw",
     "LB1_DMSO_LAD_R3.adj_log2.bw"]
].mean(axis=1)

df["mean_CTRL"] = df[
    ["LB1_DMSO_CTRL_R1.adj_log2.bw",
     "LB1_DMSO_CTRL_R2.adj_log2.bw",
     "LB1_DMSO_CTRL_R3.adj_log2.bw"]
].mean(axis=1)

df["mean_delta"] = df[["delta_R1", "delta_R2", "delta_R3"]].mean(axis=1)


#%% grna loc

region_chrom = "chr7"
region_start = 114000000
region_end   = 124000000

#%%

regions = {
    "chr7 locus": {
        "chrom": "chr7",
        "start": 115_000_000,
        "end":   123_000_000,
        "gRNA_center": (118_951_249 + 118_951_268) / 2
    },

    "chr9 locus": {
        "chrom": "chr9",
        "start": 32_200_000,
        "end":   40_200_000,
        "gRNA_center": None
    }
}


# Both regions are 10 Mb:
# 10 Mb / 2000 bins = approximately 5 kb per extracted value
n_bins = 2000
# 50 bins × 5 kb = approximately 250 kb smoothing
window_bins = 50


# read in avaeraged 53BP1 bigwigs

bp53_mean_files = {
    "DMSO_CTRL": Path(
        "/Volumes/jo4a/53bp1nospike/03_peak_calling/"
        "03_bed_to_bigwig/53BP1_DMSO_CTRL_averaged_signal.bigWig"
    ),

    "DMSO_LAD": Path(
        "/Volumes/jo4a/53bp1nospike/03_peak_calling/"
        "03_bed_to_bigwig/53BP1_DMSO_LAD_averaged_signal.bigWig"
    )
}



# read in LAD annotations

lad_bed = "/Users/patricfernandez/Downloads/4DNFIABD95PM.bed"

lad_df = pd.read_csv(
    lad_bed,
    sep="\t",
    header=None,
    comment="#"
)

lad_df = lad_df.iloc[:, :3].copy()
lad_df.columns = ["chrom", "start", "end"]

lad_df["start"] = pd.to_numeric(lad_df["start"])
lad_df["end"] = pd.to_numeric(lad_df["end"])


#%%
# FUNCTION TO EXTRACT BIGWIG SIGNAL

def extract_bigwig_region(
    bw_path,
    chrom,
    start,
    end,
    n_bins=2000,
    stat_type="mean"
):
    """
    Extract binned signal directly from a bigWig file.
    """

    bw_path = Path(bw_path)

    if not bw_path.exists():
        raise FileNotFoundError(
            f"Could not find bigWig:\n{bw_path}"
        )

    with pyBigWig.open(str(bw_path)) as bw:

        chromosome_sizes = bw.chroms()

        if chrom not in chromosome_sizes:
            raise ValueError(
                f"{chrom} is not present in {bw_path.name}"
            )

        chrom_length = chromosome_sizes[chrom]

        if start < 0 or end > chrom_length:
            raise ValueError(
                f"Requested region {chrom}:{start:,}-{end:,} "
                f"is outside chromosome bounds 0-{chrom_length:,}"
            )

        values = bw.stats(
            chrom,
            start,
            end,
            nBins=n_bins,
            type=stat_type
        )

    values = np.array(
        [
            np.nan if value is None else float(value)
            for value in values
        ],
        dtype=float
    )

    bin_edges = np.linspace(start, end, n_bins + 1)
    positions = (bin_edges[:-1] + bin_edges[1:]) / 2

    return positions, values



# SMOOTHING FUNCTION

def smooth_signal(values, window=50):
    """
    Apply a centred rolling mean while ignoring missing values.
    """

    return (
        pd.Series(values)
        .rolling(
            window=window,
            center=True,
            min_periods=1
        )
        .mean()
        .to_numpy()
    )

# EXTRACT AVERAGED 53BP1 SIGNAL

bp53_signal = {}

for region_name, region in regions.items():

    bp53_signal[region_name] = {}

    for track_name, bw_path in bp53_mean_files.items():

        positions, values = extract_bigwig_region(
            bw_path=bw_path,
            chrom=region["chrom"],
            start=region["start"],
            end=region["end"],
            n_bins=n_bins
        )

        bp53_signal[region_name][track_name] = {
            "positions": positions,
            "raw": values,
            "smooth": smooth_signal(
                values,
                window=window_bins
            )
        }


# PREPARE MEAN LB1 SIGNAL FROM df

lb1_signal = {}

for region_name, region in regions.items():

    region_df = df[
        (df["chrom"] == region["chrom"]) &
        (df["end"] > region["start"]) &
        (df["start"] < region["end"])
    ].copy()

    region_df = region_df.sort_values("start")

    region_df["center"] = (
        region_df["start"] + region_df["end"]
    ) / 2

    region_df["mean_CTRL_smooth"] = (
        region_df["mean_CTRL"]
        .rolling(
            window=window_bins,
            center=True,
            min_periods=1
        )
        .mean()
    )

    region_df["mean_LAD_smooth"] = (
        region_df["mean_LAD"]
        .rolling(
            window=window_bins,
            center=True,
            min_periods=1
        )
        .mean()
    )

    lb1_signal[region_name] = region_df

#%% plotting 


ctrl_colour = "#3b3b3b"
lad_colour  = "#bb382c"

# DMSO = solid
# ATMi = dashed
dmso_linestyle = "-"



region_names = list(regions.keys())

fig, axes = plt.subplots(
    nrows=2,
    ncols=2,
    figsize=(22, 9),
    sharex=False,
    sharey=False
)


#av 53BP1 signal
for col, region_name in enumerate(region_names):

    ax = axes[0, col]
    region = regions[region_name]

    # DMSO CTRL
    track = bp53_signal[region_name]["DMSO_CTRL"]

    ax.plot(
        track["positions"],
        track["smooth"],
        color=ctrl_colour,
        linestyle=dmso_linestyle,
        linewidth=2.2,
        label="DMSO CTRL"
    )

    # DMSO LAD
    track = bp53_signal[region_name]["DMSO_LAD"]

    ax.plot(
        track["positions"],
        track["smooth"],
        color=lad_colour,
        linestyle=dmso_linestyle,
        linewidth=2.2,
        label="DMSO LAD"
    )


    if region["gRNA_center"] is not None:
        ax.axvline(
            region["gRNA_center"],
            color="black",
            linestyle=":",
            linewidth=1.5
        )

    ax.set_xlim(
        region["start"],
        region["end"]
    )

#    ax.set_title(
 #       (
  #          f"{region_name}\n"
   #         f"{region['chrom']}:"
    #        f"{region['start']:,}–{region['end']:,}"
     #   ),
      #  fontsize=15,
       # fontweight="bold")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.tick_params(
        axis="both",
        labelsize=18
    )

    # Hide x ticks and labels on the upper row
    ax.tick_params(
        axis="x",
        which="both",
        bottom=False,
        labelbottom=False
    )

    ax.legend(
        frameon=False,
        fontsize=10,
        loc="upper right",
        ncol=2
    )

axes[0, 0].set_ylabel(
    "Mean smoothed 53BP1 signal",
    fontsize=12
)


#LaminB1 averaged signal 

for col, region_name in enumerate(region_names):

    ax = axes[1, col]
    region = regions[region_name]
    region_df = lb1_signal[region_name]

    # Shade annotated LADs
    lad_region = lad_df[
        (lad_df["chrom"] == region["chrom"]) &
        (lad_df["end"] > region["start"]) &
        (lad_df["start"] < region["end"])
    ]

    for _, lad_row in lad_region.iterrows():

        shade_start = max(
            lad_row["start"],
            region["start"]
        )

        shade_end = min(
            lad_row["end"],
            region["end"]
        )

        ax.axvspan(
            shade_start,
            shade_end,
            color="#A5CDFF",
            alpha=0.28,
            zorder=0
        )

    ax.plot(
        region_df["center"],
        region_df["mean_CTRL_smooth"],
        color="grey",
        linewidth=2.2,
        label="Mean CTRL",
        zorder=2
    )

    ax.plot(
        region_df["center"],
        region_df["mean_LAD_smooth"],
        color="tomato",
        linewidth=2.2,
        label="Mean LAD",
        zorder=2
    )

    if region["gRNA_center"] is not None:
        ax.axvline(
            region["gRNA_center"],
            color="black",
            linestyle=":",
            linewidth=1.5,
            zorder=3
        )

    ax.set_xlim(
        region["start"],
        region["end"]
    )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.tick_params(
        axis="both",
        labelsize=18
    )

    ax.legend(
        frameon=False,
        fontsize=10,
        loc="upper right"
    )

axes[1, 0].set_ylabel(
    "Mean smoothed LB1 signal",
    fontsize=12
)



# MATCH Y-AXES ACROSS EACH ROW

for row in range(2):

    left_ymin, left_ymax = axes[row, 0].get_ylim()
    right_ymin, right_ymax = axes[row, 1].get_ylim()

    shared_ymin = min(left_ymin, right_ymin)
    shared_ymax = max(left_ymax, right_ymax)

    axes[row, 0].set_ylim(
        shared_ymin,
        shared_ymax
    )

    axes[row, 1].set_ylim(
        shared_ymin,
        shared_ymax
    )


#axis labels
axes[1, 0].set_xlabel(
    "Genomic position on chr7",
    fontsize=13
)

axes[1, 1].set_xlabel(
    "Genomic position on chr9",
    fontsize=13
)


#fig.suptitle(
 #   "Mean 53BP1 and Lamin B1 signal across chr7 and chr9 loci",
  #  fontsize=18,
  #  fontweight="bold",
  #  y=0.995
#)

plt.tight_layout()
plt.show()
