#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 10:53:52 2026

@author: joannafernandez
Written with assistance from GPT 5.6 Thinking
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os 

#import data file and name cols 
def Find(file, rep):
    genes = pd.read_csv(file, delimiter=",")
    genes['replicate'] = rep
    genes_all = genes[['condition', 'nuclei_no', 'yH2AX', 'pATM', 'pATMcoloc', 'replicate']]
    genes_large = genes[['condition', 'nuclei_no', 'large_yH2AX','large_pATM', 'large_coloc','replicate']]
    return genes_all, genes_large

####replicates
rep5all, rep5large = Find('/Users/patricfernandez/rep5_foci.csv', 2)
rep4all, rep4large = Find('/Users/patricfernandez/rep4_foci.csv', 1)
rep6all, rep6large = Find('/Users/patricfernandez/rep6_foci.csv', 3)


def processing_raw (df4,df5,df6):
    
    dfy = pd.concat([df4, df5, df6])
    dfy = dfy.dropna()
    cols_to_int = ["nuclei_no", "yH2AX", "pATM", "pATMcoloc", "replicate"]

    dfy[cols_to_int] = dfy[cols_to_int].astype(int)
    dfy["condition"] = dfy["condition"].astype(str)
    
    #now, calculate per replicate, how many cells have yH2AX (non-zero) also by condition
    
    
    return dfy

allrepsall = processing_raw(rep4all, rep5all, rep6all)

#%% where the "unit" is an individual focus

df_foci = allrepsall.copy()

foci_summary = (
    df_foci
    .groupby(["condition", "replicate"], as_index=False)
    .agg(
        total_nuclei=("nuclei_no", "size"),
        total_yH2AX_foci=("yH2AX", "sum"),
        total_pATM_foci=("pATM", "sum"),
        total_pATMcoloc_foci=("pATMcoloc", "sum")
    )
)

# Mean number of foci per nucleus
foci_summary["yH2AX_foci_per_nucleus"] = (
    foci_summary["total_yH2AX_foci"]
    / foci_summary["total_nuclei"]
)

foci_summary["pATM_foci_per_nucleus"] = (
    foci_summary["total_pATM_foci"]
    / foci_summary["total_nuclei"]
)

foci_summary["pATMcoloc_foci_per_nucleus"] = (
    foci_summary["total_pATMcoloc_foci"]
    / foci_summary["total_nuclei"]
)

# Colocalised pATM foci relative to γH2AX foci
foci_summary["percent_coloc_relative_to_yH2AX"] = np.where(
    foci_summary["total_yH2AX_foci"] > 0,
    (
        foci_summary["total_pATMcoloc_foci"]
        / foci_summary["total_yH2AX_foci"]
    ) * 100,
    np.nan
)

#%% plotting 

rep_palette = {
    1: "#ed3658",
    2: "#36bbac",
    3: "#FFFFFF"
}

condition_order = [
    "WT_sgCTRL",
    "WT_sgLAD",
    "53ko_sgCTRL",
    "53ko_sgLAD"
]

condition_colours = {
    "WT_sgCTRL": "#B3B3B3",
    "WT_sgLAD": "#FB0D44",
    "53ko_sgCTRL": "#B3B3B3",
    "53ko_sgLAD": "#FB0D44"
}

bar_colours = [
    condition_colours[condition]
    for condition in condition_order
]

# Choose the column to plot
value_col = "percent_coloc_relative_to_yH2AX"

# Set condition order
foci_summary["condition"] = pd.Categorical(
    foci_summary["condition"],
    categories=condition_order,
    ordered=True
)

foci_summary = foci_summary.sort_values(
    ["condition", "replicate"]
)

# Mean and SEM across biological replicates
plot_summary = (
    foci_summary
    .groupby("condition", observed=False)
    .agg(
        mean_value=(value_col, "mean"),
        sem_value=(value_col, "sem")
    )
    .reindex(condition_order)
    .reset_index()
)

# Plot
fig, ax = plt.subplots(figsize=(8, 7))

x = np.arange(len(condition_order))
bar_width = 0.65

ax.bar(
    x,
    plot_summary["mean_value"],
    width=bar_width,
    color=bar_colours,
    edgecolor="black",
    linewidth=1.2,
    zorder=1
)

# SEM across the three biological replicates
ax.errorbar(
    x,
    plot_summary["mean_value"],
    yerr=plot_summary["sem_value"],
    fmt="none",
    ecolor="black",
    elinewidth=1.5,
    capsize=5,
    capthick=1.5,
    zorder=4
)

# Fixed offsets for the three replicate dots
replicate_offsets = {
    1: -0.10,
    2:  0.00,
    3:  0.10
}

for i, condition in enumerate(condition_order):

    condition_data = foci_summary[
        foci_summary["condition"] == condition
    ]

    for _, row in condition_data.iterrows():

        offset = replicate_offsets.get(row["replicate"], 0)

        rep = row["replicate"]

        ax.scatter(
            i + offset,
            row[value_col],
            s=140,
            facecolor=rep_palette.get(rep, "white"),
            edgecolor="black",
            linewidth=1.3,
            zorder=5,
            label=f"Rep {rep}"
        )

ax.set_xticks(x)
ax.set_xticklabels(
    condition_order,
    rotation=30,
    ha="right"
)

#ax.axvline(
 #   x=1.5,
  #  color="grey",
   # linestyle=":",
    #linewidth=1.5,
   # zorder=0)

ax.set_ylabel(
    "Colocalised pATM foci relative to\nγH2AX foci (%)",
    fontsize=14
)

ax.set_xlabel("")
ax.tick_params(axis="y", labelsize=18)

# Leave enough space above the highest point/error bar
upper_limit = np.nanmax([
    foci_summary[value_col].max(),
    (
        plot_summary["mean_value"]
        + plot_summary["sem_value"]
    ).max()
])

ax.set_ylim(0, upper_limit * 1.15)


# Clean replicate legend
handles, labels = ax.get_legend_handles_labels()

unique = {}
for h, l in zip(handles, labels):
    if l.startswith("Rep") and l not in unique:
        unique[l] = h

ax.legend(
    unique.values(),
    unique.keys(),
    title="Replicate",
    frameon=False,
    fontsize=12,
    title_fontsize=13,
    loc="upper right"
)


sns.despine()
plt.tight_layout()
plt.show()


#%% stats

import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

value_col = "percent_coloc_relative_to_yH2AX"
stats_df = foci_summary[
    ["condition", "replicate", value_col]
].copy()

stats_df["genotype"] = stats_df["condition"].map({
    "WT_sgCTRL": "WT",
    "WT_sgLAD": "WT",
    "53ko_sgCTRL": "53ko",
    "53ko_sgLAD": "53ko"})

stats_df["guide"] = stats_df["condition"].map({
    "WT_sgCTRL": "sgCTRL",
    "WT_sgLAD": "sgLAD",
    "53ko_sgCTRL": "sgCTRL",
    "53ko_sgLAD": "sgLAD"})

stats_df["replicate"] = stats_df["replicate"].astype("category")

model = smf.ols(
    formula=(
        f'{value_col} ~ '
        'C(genotype, Sum) * C(guide, Sum) + '
        'C(replicate, Sum)'
    ),
    data=stats_df
).fit()

anova_tbl = sm.stats.anova_lm(
    model,
    typ=3
)

print(anova_tbl)






