#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 15:41:56 2026

@author: joannaernandez
Written with assistance from GPT 5.6
"""
#run this cell when plotting ATMi condition
#atmiatmiatmiatmiatmiatmiatmiatmiatmiatmiatmiatmiatmiatmiatmiatmiatmiatmiatmiatmiatmiatmiatmiatmiatmiatmiatmiatmiatmi

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pyBigWig


# LOAD BINNED BIGWIG TABLE
path = "/Volumes/jo4/laminnospike/03_peak_calling/03_bed_to_bigwig/derived_bigwigs/ATMi_5kb_rawcounts.tsv"

df = pd.read_csv(path, sep="\t")
df.columns = [c.replace("'", "").lstrip("#") for c in df.columns]
df = df.rename(columns={"chr": "chrom"})
df = df.fillna(0)

df["delta_R1"] = df["LB1_ATMi_LAD_R1.adj_log2.bw"] - df["LB1_ATMi_CTRL_R1.adj_log2.bw"]
df["delta_R3"] = df["LB1_ATMi_LAD_R3.adj_log2.bw"] - df["LB1_ATMi_CTRL_R3.adj_log2.bw"]

df["mean_LAD"] = df[
    ["LB1_ATMi_LAD_R1.adj_log2.bw",
     "LB1_ATMi_LAD_R3.adj_log2.bw"]
].mean(axis=1)

df["mean_CTRL"] = df[
    ["LB1_ATMi_CTRL_R1.adj_log2.bw",
     "LB1_ATMi_CTRL_R3.adj_log2.bw"]
].mean(axis=1)

df["mean_delta"] = df[["delta_R1", "delta_R3"]].mean(axis=1)



#%%
#run this cell when plotting dmso condition
#dmsodmsodmsodmsodmsodmsodmsodmsodmsodmsodmsodmsodmsodmsodmsodmsodmsodmsodmsodmsodmsodmsodmsodmsodmsodmsodmsodmsodmsodmso

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


#%%


#%% domainogram structuring
bin_size = 5000

window_sizes_bp = [25_000, 35_000, 50_000, 70_000, 100_000, 140_000,
                   200_000, 280_000, 400_000, 560_000, 800_000, 1_000_000]

window_sizes_bins = [w // bin_size for w in window_sizes_bp]

chrom_order = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
df["chrom"] = pd.Categorical(df["chrom"], categories=chrom_order, ordered=True)
df = df.sort_values(["chrom", "start"]).reset_index(drop=True)



def sliding_mean(arr, window):
    """
    Compute sliding mean for a 1D numpy array using a fixed window size.
    Returns an array of length len(arr) - window + 1.
    """
    arr = np.asarray(arr, dtype=float)
    csum = np.cumsum(np.insert(arr, 0, 0.0))
    return (csum[window:] - csum[:-window]) / window


def compute_windowed_scores(df, window_bins, bin_size=5000):
    """
    For one window size, compute sliding-window means chromosome by chromosome.
    Uses window CENTER for plotting, which gives the more classic triangular
    domainogram appearance.
    """
    results = []

    for chrom, sub in df.groupby("chrom", observed=True):
        sub = sub.sort_values("start").reset_index(drop=True)

        if len(sub) < window_bins:
            continue

        delta_vals = sub["mean_delta"].to_numpy()
        lad_vals   = sub["mean_LAD"].to_numpy()
        ctrl_vals  = sub["mean_CTRL"].to_numpy()

        delta_win = sliding_mean(delta_vals, window_bins)
        lad_win   = sliding_mean(lad_vals, window_bins)
        ctrl_win  = sliding_mean(ctrl_vals, window_bins)

        starts = sub["start"].to_numpy()[:len(delta_win)]
        ends   = starts + window_bins * bin_size
        center = starts + (window_bins * bin_size) / 2

        out = pd.DataFrame({
            "chrom": chrom,
            "start": starts,
            "end": ends,
            "center": center,
            "window_bins": window_bins,
            "window_bp": window_bins * bin_size,
            "mean_delta_window": delta_win,
            "mean_LAD_window": lad_win,
            "mean_CTRL_window": ctrl_win
        })

        results.append(out)

    return pd.concat(results, ignore_index=True)


def process_one_window_size(df, window_bins, bin_size=5000, baseline_q=0.30):
    """
    For one window size:
    1. compute sliding-window means
    2. baseline-filter windows
    3. rank retained windows
    4. assign top/bottom 5% signal
    """
    win = compute_windowed_scores(df, window_bins=window_bins, bin_size=bin_size).copy()

    # baseline thresholds for this window size
    lad_thresh = win["mean_LAD_window"].quantile(baseline_q)
    ctrl_thresh = win["mean_CTRL_window"].quantile(baseline_q)

    # keep windows unless both LAD and CTRL are baseline-low
    win["keep"] = ~(
        (win["mean_LAD_window"] <= lad_thresh) &
        (win["mean_CTRL_window"] <= ctrl_thresh)
    )

    # rank only the kept windows
    win_kept = win.loc[win["keep"]].copy()
    win_kept["percentile"] = win_kept["mean_delta_window"].rank(method="average", pct=True)

    # binary paper-style coloring
    win_kept["signal"] = 0
    win_kept.loc[win_kept["percentile"] <= 0.05, "signal"] = -1
    win_kept.loc[win_kept["percentile"] >= 0.95, "signal"] = 1

    return win_kept


# BUILD DOMAINOGRAM DF

domainogram_layers = []

for w in window_sizes_bins:
    print(f"Processing {w} bins ({w*5000/1000:.0f} kb)")
    win_kept = process_one_window_size(df, window_bins=w, bin_size=5000, baseline_q=0.30)
    domainogram_layers.append(win_kept)

domainogram_df = pd.concat(domainogram_layers, ignore_index=True)

#%% plotting LAD
region_chrom = "chr7"
region_start = 114000000
region_end   = 124000000

#unsilence to plot iLAD
#region_chrom = "chr9"
#region_start = 31200000
#region_end   = 41200000


plot_df = domainogram_df[
    (domainogram_df["chrom"] == region_chrom) &
    (domainogram_df["end"] > region_start) &
    (domainogram_df["start"] < region_end)
].copy()

plot_df["window_kb"] = plot_df["window_bp"] / 1000


# pivot using CENTER, not midpoint
matrix_binary = plot_df.pivot_table(
    index="window_kb",
    columns="center",
    values="signal",
    aggfunc="first",
    fill_value=0
).sort_index(ascending=False)

#%%

domainogram_df = domainogram_df.copy()

domainogram_df["tail_score"] = 0.0

low_mask = domainogram_df["percentile"] <= 0.05
domainogram_df.loc[low_mask, "tail_score"] = (
    (domainogram_df.loc[low_mask, "percentile"] / 0.05) - 1.0
)

high_mask = domainogram_df["percentile"] >= 0.95
domainogram_df.loc[high_mask, "tail_score"] = (
    (domainogram_df.loc[high_mask, "percentile"] - 0.95) / 0.05
)

plot_df = domainogram_df[
    (domainogram_df["chrom"] == region_chrom) &
    (domainogram_df["end"] > region_start) &
    (domainogram_df["start"] < region_end)
].copy()

plot_df["window_kb"] = plot_df["window_bp"] / 1000

matrix_gradual = plot_df.pivot_table(
    index="window_kb",
    columns="center",
    values="tail_score",
    aggfunc="first",
    fill_value=0
).sort_index(ascending=True)


#%%

raw_region = df[
    (df["chrom"] == region_chrom) &
    (df["end"] > region_start) &
    (df["start"] < region_end)
].copy()

raw_region["center"] = (raw_region["start"] + raw_region["end"]) / 2
raw_region["LAD_smooth"] = raw_region["mean_LAD"].rolling(20, center=True).mean()
raw_region["CTRL_smooth"] = raw_region["mean_CTRL"].rolling(20, center=True).mean()

#%%


#%%

#nospikeversion
#dmsodmsodmsodmsodmsodmsodmsodmsodmsodmsodmsodmsodmsodmsodmsodmsodmsodmsodmsodmsodmsodmsodmsodmsodmsodmsodmsodmsodmsodmso

bw_53bp1_ctrl = "/Volumes/jo4a/53bp1nospike/03_peak_calling/03_bed_to_bigwig/53BP1_DMSO_CTRL_averaged_signal.bigWig"
bw_53bp1_lad  = "/Volumes/jo4a/53bp1nospike/03_peak_calling/03_bed_to_bigwig/53BP1_DMSO_LAD_averaged_signal.bigWig"
# display resolution for the 53BP1 track

bw_bin_size = 5000


#%%


bw_53bp1_ctrl = "/Volumes/jo4a/53bp1nospike/03_peak_calling/03_bed_to_bigwig/53BP1_ATMi_CTRL_averaged_signal.bigWig"
bw_53bp1_lad  = "/Volumes/jo4a/53bp1nospike/03_peak_calling/03_bed_to_bigwig/53BP1_ATMi_LAD_averaged_signal.bigWig"
# display resolution for the 53BP1 track

bw_bin_size = 5000
#%%

#chr7
gRNA_start = 118_951_249
gRNA_end   = 118_951_268
gRNA_center = (gRNA_start + gRNA_end) / 2

#unsilence for ilad
#chr9
#gRNA_start = 36210899
#gRNA_end   = 36210918
#gRNA_center = (gRNA_start + gRNA_end) / 2

#%%load in LAD annotations


lad_bed = '/Users/patricfernandez/Downloads/4DNFIABD95PM.bed'

# assume standard BED: chrom, start, end, ...
lad_df = pd.read_csv(
    lad_bed,
    sep="\t",
    header=None,
    comment="#"
)

# keep first three BED columns
lad_df = lad_df.iloc[:, :3].copy()
lad_df.columns = ["chrom", "start", "end"]

# subset to plotted region
lad_region = lad_df[
    (lad_df["chrom"] == region_chrom) &
    (lad_df["end"] > region_start) &
    (lad_df["start"] < region_end)
].copy()

lad_region.head()


#%%
# =========================================================
# HELPER: read bigwig into fixed bins across the chosen region
# =========================================================

def extract_bigwig_region_binned(bw_path, chrom, start, end, bin_size=5000):
    bw = pyBigWig.open(bw_path)

    starts = np.arange(start, end, bin_size)
    ends = np.minimum(starts + bin_size, end)

    values = []
    for s, e in zip(starts, ends):
        v = bw.stats(chrom, int(s), int(e), type="mean")[0]
        values.append(0.0 if v is None else v)

    bw.close()

    out = pd.DataFrame({
        "chrom": chrom,
        "start": starts,
        "end": ends,
        "center": (starts + ends) / 2,
        "signal": values
    })
    return out

# =========================================================
# 53BP1 TRACKS
# =========================================================

bp1_ctrl = extract_bigwig_region_binned(
    bw_53bp1_ctrl, region_chrom, region_start, region_end, bin_size=bw_bin_size
)
bp1_lad = extract_bigwig_region_binned(
    bw_53bp1_lad, region_chrom, region_start, region_end, bin_size=bw_bin_size
)

bp1_df = bp1_ctrl.rename(columns={"signal": "BP1_CTRL"}).merge(
    bp1_lad[["chrom", "start", "end", "center", "signal"]].rename(columns={"signal": "BP1_LAD"}),
    on=["chrom", "start", "end", "center"],
    how="inner"
)

# optional light smoothing for display only
bp1_df["BP1_CTRL_smooth"] = bp1_df["BP1_CTRL"].rolling(5, center=True).mean()
bp1_df["BP1_LAD_smooth"]  = bp1_df["BP1_LAD"].rolling(5, center=True).mean()

# =========================================================
# LAMIN RAW TRACKS (from your existing df)
# =========================================================

raw_region = df[
    (df["chrom"] == region_chrom) &
    (df["end"] > region_start) &
    (df["start"] < region_end)
].copy()

raw_region["center"] = (raw_region["start"] + raw_region["end"]) / 2
raw_region["LAD_smooth"] = raw_region["mean_LAD"].rolling(20, center=True).mean()
raw_region["CTRL_smooth"] = raw_region["mean_CTRL"].rolling(20, center=True).mean()

# =========================================================
# DOMAINOGRAM MATRIX (from your existing domainogram_df)
# assumes you already created domainogram_df and tail_score
# =========================================================

plot_df = domainogram_df[
    (domainogram_df["chrom"] == region_chrom) &
    (domainogram_df["end"] > region_start) &
    (domainogram_df["start"] < region_end)
].copy()

plot_df["window_kb"] = plot_df["window_bp"] / 1000

# ---------------------------------------------------------
# Important:
# Window centres can fall either on exact 5 kb positions
# or halfway between them, depending on whether the window
# has an even or odd number of bins.
#
# Therefore use a 2.5 kb coordinate grid, not a 5 kb grid.
# ---------------------------------------------------------

half_bin = bin_size / 2  # 2500 bp

plot_df["center_binned"] = (
    np.round(plot_df["center"] / half_bin) * half_bin
)

x_grid = np.arange(
    region_start,
    region_end + half_bin,
    half_bin
)

plot_df = plot_df[
    (plot_df["center_binned"] >= x_grid.min()) &
    (plot_df["center_binned"] <= x_grid.max())
].copy()

# Pivot WITHOUT fill_value first.
# This lets us distinguish true masked 0s from missing coordinates.
matrix_gradual = plot_df.pivot_table(
    index="window_kb",
    columns="center_binned",
    values="tail_score",
    aggfunc="mean"
).sort_index(ascending=True)

# Reindex to the common 2.5 kb grid
matrix_gradual = matrix_gradual.reindex(columns=x_grid)

# ---------------------------------------------------------
# Fill only tiny one-column gaps caused by the odd/even
# window-centre grid mismatch.
#
# This is not broad interpolation across the region.
# It only prevents barcode-like missing-column artefacts.
# ---------------------------------------------------------

matrix_gradual = (
    matrix_gradual
    .ffill(axis=1, limit=1)
    .bfill(axis=1, limit=1)
)

# Remaining missing values are true background/masked values
matrix_gradual = matrix_gradual.fillna(0)






#%% plotting 

fig = plt.figure(figsize=(20, 8), dpi=200)
gs = fig.add_gridspec(
    nrows=3, ncols=2,
    width_ratios=[40, 1],
    height_ratios=[1, 1, 3],
    hspace=0.08, wspace=0.05
)

ax0 = fig.add_subplot(gs[0, 0])              # 53BP1
ax1 = fig.add_subplot(gs[1, 0], sharex=ax0)  # Lamin
ax2 = fig.add_subplot(gs[2, 0], sharex=ax0)  # Domainogram
cax = fig.add_subplot(gs[2, 1])              # Colorbar

# TOP PANEL: 53BP1

ax0.plot(
    bp1_df["center"],
    bp1_df["BP1_CTRL_smooth"],
    color="grey",
    linewidth=1.2,
    label="53BP1 CTRL",alpha=0.7
)

ax0.plot(
    bp1_df["center"],
    bp1_df["BP1_LAD_smooth"],
    color="tomato",
    linewidth=1.2,
    label="53BP1 LAD", alpha=0.7
)

ax0.set_ylim(top=2)

#ax0.set_ylabel("53BP1\nsignal")
#ax0.set_title("53BP1 averaged signal", loc="left")
#ax0.legend(frameon=False, loc="upper right", fontsize=8)
ax0.spines["top"].set_visible(False)
ax0.spines["right"].set_visible(False)
ax0.tick_params(axis="x", labelbottom=False)


# MIDDLE PANEL: LaminB1

ax1.plot(
    raw_region["center"],
    raw_region["CTRL_smooth"],
    color="grey",
    linewidth=1.2,
    label="LB1 CTRL"
)

ax1.plot(
    raw_region["center"],
    raw_region["LAD_smooth"],
    color="tomato",
    linewidth=1.2,
    label="LB1 LAD", alpha =0.7
)

################3
# OVERLAY LADs ON ax1

for _, row in lad_region.iterrows():
    ax1.axvspan(
        row["start"],
        row["end"],
        color="#A5CDFF",
        alpha=0.28,
        zorder=0
    )
    
########  

ax1.axvline(
    gRNA_center,
    color="black",
    linestyle="--",
    linewidth=0.8,
    zorder=10
)
ax0.axvline(
    gRNA_center,
    color="black",
    linestyle="--",
    linewidth=0.8,
    zorder=10
)
#ax1.set_ylabel("LaminB1\nsignal")
#ax1.set_title("LaminB1 signal", loc="left")
#ax1.legend(frameon=False, loc="upper right", fontsize=8)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)
ax1.tick_params(axis="x", labelbottom=False)


# BOTTOM PANEL: DOMAINOGRAM HEATMAP

n_rows = matrix_gradual.shape[0]
y_positions = np.arange(n_rows)
y_labels = matrix_gradual.index.to_numpy(dtype=float)

im = ax2.imshow(
    matrix_gradual.values,
    aspect="auto",
    interpolation="bilinear",
    cmap="RdBu_r",
    vmin=-1,
    vmax=1,
    extent=[
        region_start,
        region_end,
        -0.5,
        n_rows - 0.5
    ],
    origin="lower"
)

ax2.axvline(
    gRNA_center,
    color="black",
    linestyle="--",
    linewidth=0.8,
    zorder=10
)

# x-axis
x_ticks = np.linspace(region_start, region_end, 6)
ax2.set_xticks(x_ticks)
ax2.set_xticklabels([f"{x / 1e6:.1f}" for x in x_ticks])

# y-axis as categorical window-size rows
ax2.set_yticks(y_positions)
ax2.set_yticklabels([f"{int(y)} kb" for y in y_labels])
ax2.tick_params(axis="y", labelsize=8)

ax2.set_xlabel(f"{region_chrom} position (Mb)")
ax2.set_ylabel("Window size")

ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)

# FORCE IDENTICAL X LIMITS

for a in [ax0, ax1, ax2]:
    a.set_xlim(region_start, region_end)

# COLORBAR

cbar = fig.colorbar(im, cax=cax)
cbar.set_label("Tail-scaled percentile score")
cbar.set_ticks([-1, -0.5, 0, 0.5, 1])
cbar.set_ticklabels(["0%", "2.5%", "masked", "97.5%", "100%"])

plt.show()

