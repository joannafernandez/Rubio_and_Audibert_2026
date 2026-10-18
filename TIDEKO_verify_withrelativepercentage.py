#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 17:44:20 2026

@author: joannafernandez
"""

#%% import and processing raw data 

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pyBigWig

#import data file and name cols 
def Find(file):
    genes = pd.read_csv(file, delimiter=",")
    genes.rename(columns={'Unnamed: 0': 'replicate'}, inplace=True)
    genes.rename(columns={'Unnamed: 1': 'indels'}, inplace=True)

    return genes


R0perc = Find('/Users/patricfernandez/koverify_perc.csv')
R0perc= R0perc.dropna() #remove empty rows
R0pval = Find('/Users/patricfernandez/koverify_pval.csv')
R0pval= R0pval.dropna() #remove empty rows

R0perc.loc[:, ~R0perc.columns.isin(['replicate', 'indels'])] = R0perc.loc[:, ~R0perc.columns.isin(['replicate', 'indels'])].mask(R0pval > 0.001, 0)

#formating for plotting
melted_df = R0perc.melt(id_vars=['replicate', 'indels'],
                               var_name='condition',
                               value_vars= ['sgMDC1', 'sgRNF168', 'sgLaminB1', 'sgLBR'],
                               value_name='percentage')
melted_df['percentage'] = melted_df['percentage'].astype(int)

# assign frameshift type
def classify_frameshift(indel):
    if indel == 0:
        return "uncut"
    elif indel % 3 == 0:
        return "in-frame"
    else:
        return "frameshift"
    
melted_df["class"] = melted_df["indels"].apply(classify_frameshift)




#%%

#alternative
melted_df['class'] = melted_df['indels'].apply(classify_frameshift)
melted_df = melted_df[melted_df['class'].notna()]

#find total % of sequences that are either NHEJ or MMEJ. (uncut/wt sequences are only ever "0")
summed_df = (
    melted_df
    .groupby(['condition', 'replicate', 'class', 'indels'], as_index=False)
    .agg(summed_percentage=('percentage', 'sum'))
)


subset = summed_df.copy()

#Calculate total NHEJ + uncut per condition/replicate
total = subset.groupby(['condition', 'replicate'])['summed_percentage'].transform('sum')
#Calculate percentage within that total
subset['relative_percentage'] = (subset['summed_percentage'] / total) * 100
# Split 'condition' into two new columns: 'condition' and 'time'
#subset[['condition', 'time']] = subset['condition'].str.split('_', n=1, expand=True)

melted_df = subset.copy()

#%% seperate for plots

sgMDC1_df = melted_df[melted_df["condition"] == "sgMDC1"].copy()
sgRNF168_df = melted_df[melted_df["condition"] == "sgRNF168"].copy()
sgLaminB1_df = melted_df[melted_df["condition"] == "sgLaminB1"].copy()
sgLBR_df = melted_df[melted_df["condition"] == "sgLBR"].copy()

#%% plotting sgMDC1
custom_palette = {
    "in-frame": "#A1C379",
    "uncut": "#F3C773",
    "frameshift": "#DC6B83"
}

# restrict x-axis to indels between -15 and 5
plot_df = sgMDC1_df[
    (sgMDC1_df["indels"] >= -15) &
    (sgMDC1_df["indels"] <= 5)
].copy()

plot_dots = (
    plot_df
    .groupby(["indels", "class"], group_keys=False, observed=True)
    .filter(lambda x: (x["relative_percentage"] != 0).any())
    .copy()
)

indel_order = sorted(plot_df["indels"].unique())
class_order = ["in-frame", "frameshift", "uncut"]

# Use this for the manual scatter section
sgMDC1_dots = plot_dots.copy()




plt.figure(figsize=(14, 5))
sns.set_context("talk")

ax = sns.barplot(
    data=plot_df,
    x="indels",
    y="relative_percentage",
    hue="class",
    order=indel_order,
    hue_order=class_order,
    palette=custom_palette,
    estimator="mean",
    errorbar=None,
    edgecolor="black",
    linewidth=1.4,
    alpha=0.75
)

# make x positions manually so dots can be coloured by replicate
x_lookup = {indel: i for i, indel in enumerate(indel_order)}

hue_offsets = {
    "in-frame": -0.25,
    "frameshift": 0,
    "uncut": 0.25
}

# replicate colours
rep_palette = {
    1: "#ed3658",
    2: "#36bbac",
    3: "#FFFFFF"}

# rep_palette = {
#     "1": "#4D4D4D",
#     "2": "#969696",
#     "3": "#D9D9D9"}

rng = np.random.default_rng(42)

for cls in class_order:
    sub_cls = sgMDC1_dots[sgMDC1_dots["class"] == cls].copy()

    # plot replicate dots, coloured by replicate
    for rep, sub_rep in sub_cls.groupby("replicate"):
        x = sub_rep["indels"].map(x_lookup).astype(float)
        x = x + hue_offsets[cls]
        x = x + rng.normal(0, 0.035, size=len(sub_rep))

        ax.scatter(
            x,
            sub_rep["relative_percentage"],
            s=120,
            facecolors=rep_palette[rep],
            edgecolors="black",
            linewidths=1.8,
            zorder=5,
            label=f"Rep {rep}"
        )

    # add mean ± SEM error bars for each indel within this class
    for indel in indel_order:
        sub = sgMDC1_df[
            (sgMDC1_df["indels"] == indel) &
            (sgMDC1_df["class"] == cls)
        ]

        vals = sub["relative_percentage"].to_numpy()

        # skip if there are no values or all values are zero
        if vals.size == 0 or np.all(vals == 0):
            continue

        mean_val = vals.mean()

        if vals.size > 1:
            sem = vals.std(ddof=1) / np.sqrt(vals.size)
        else:
            sem = 0

        xpos = x_lookup[indel] + hue_offsets[cls]

        ax.errorbar(
            xpos,
            mean_val,
            yerr=sem,
            fmt="none",
            ecolor="black",
            elinewidth=1.6,
            capsize=5,
            capthick=1.6,
            zorder=6)
        

handles, labels = ax.get_legend_handles_labels()
n_classes = len(class_order)

ax.legend(
    handles[:n_classes],
    labels[:n_classes],
    title="Class",
    frameon=False,
    fontsize=12,
    title_fontsize=13)

ax.set_xlabel("Indel size", fontsize=3)
ax.set_ylabel("Percentage", fontsize=3)

ax.tick_params(axis="x", labelsize=22, width=1.6, length=4, rotation=90)
ax.tick_params(axis="y", labelsize=22, width=1.6, length=4)

ax.spines["left"].set_linewidth(2.0)
ax.spines["bottom"].set_linewidth(2.0)
ax.spines["left"].set_color("black")
ax.spines["bottom"].set_color("black")

sns.despine()
plt.tight_layout()
plt.show()


#%% plotting sgLBR

# restrict x-axis to indels between -15 and 5
plot_df = sgLBR_df[
    (sgLBR_df["indels"] >= -15) &
    (sgLBR_df["indels"] <= 5)
].copy()

plot_dots = (
    plot_df
    .groupby(["indels", "class"], group_keys=False, observed=True)
    .filter(lambda x: (x["relative_percentage"] != 0).any())
    .copy()
)

indel_order = sorted(plot_df["indels"].unique())
class_order = ["in-frame", "frameshift", "uncut"]

# Use this for the manual scatter section
sgLBR_dots = plot_dots.copy()

plt.figure(figsize=(14, 5))
sns.set_context("talk")

ax = sns.barplot(
    data=sgLBR_df,
    x="indels",
    y="relative_percentage",
    hue="class",
    order=indel_order,
    hue_order=class_order,
    palette=custom_palette,
    estimator="mean",
    errorbar=None,
    edgecolor="black",
    linewidth=1.4,
    alpha=0.75)

# make x positions manually so dots can be coloured by replicate
x_lookup = {indel: i for i, indel in enumerate(indel_order)}

hue_offsets = {
    "in-frame": -0.25,
    "frameshift": 0,
    "uncut": 0.25
}


rng = np.random.default_rng(42)

for cls in class_order:
    sub_cls = sgLBR_dots[sgLBR_dots["class"] == cls].copy()

    # plot replicate dots, coloured by replicate
    for rep, sub_rep in sub_cls.groupby("replicate"):
        x = sub_rep["indels"].map(x_lookup).astype(float)
        x = x + hue_offsets[cls]
        x = x + rng.normal(0, 0.035, size=len(sub_rep))

        ax.scatter(
            x,
            sub_rep["relative_percentage"],
            s=120,
            facecolors=rep_palette[rep],
            edgecolors="black",
            linewidths=1.8,
            zorder=5,
            label=f"Rep {rep}"
        )

    # add mean ± SEM error bars for each indel within this class
    for indel in indel_order:
        sub = sgLBR_df[
            (sgLBR_df["indels"] == indel) &
            (sgLBR_df["class"] == cls)
        ]

        vals = sub["relative_percentage"].to_numpy()

        # skip if there are no values or all values are zero
        if vals.size == 0 or np.all(vals == 0):
            continue

        mean_val = vals.mean()

        if vals.size > 1:
            sem = vals.std(ddof=1) / np.sqrt(vals.size)
        else:
            sem = 0

        xpos = x_lookup[indel] + hue_offsets[cls]

        ax.errorbar(
            xpos,
            mean_val,
            yerr=sem,
            fmt="none",
            ecolor="black",
            elinewidth=1.6,
            capsize=5,
            capthick=1.6,
            zorder=6)


handles, labels = ax.get_legend_handles_labels()
n_classes = len(class_order)

ax.legend(
    handles[:n_classes],
    labels[:n_classes],
    title="Class",
    frameon=False,
    fontsize=12,
    title_fontsize=13
)

ax.set_xlabel("Indel size", fontsize=3)
ax.set_ylabel("Percentage", fontsize=3)

ax.tick_params(axis="x", labelsize=22, width=1.6, length=4, rotation=90)
ax.tick_params(axis="y", labelsize=22, width=1.6, length=4)

ax.spines["left"].set_linewidth(2.0)
ax.spines["bottom"].set_linewidth(2.0)
ax.spines["left"].set_color("black")
ax.spines["bottom"].set_color("black")

sns.despine()
plt.tight_layout()
plt.show()

#%% plotting sgRNF168

# restrict x-axis to indels between -15 and 5
plot_df = sgRNF168_df[
    (sgRNF168_df["indels"] >= -15) &
    (sgRNF168_df["indels"] <= 5)
].copy()

plot_dots = (
    plot_df
    .groupby(["indels", "class"], group_keys=False, observed=True)
    .filter(lambda x: (x["relative_percentage"] != 0).any())
    .copy()
)

indel_order = sorted(plot_df["indels"].unique())
class_order = ["in-frame", "frameshift", "uncut"]

# Use this for the manual scatter section
sgRNF168_dots = plot_dots.copy()



plt.figure(figsize=(14, 5))
sns.set_context("talk")

ax = sns.barplot(
    data=plot_df,
    x="indels",
    y="relative_percentage",
    hue="class",
    order=indel_order,
    hue_order=class_order,
    palette=custom_palette,
    estimator="mean",
    errorbar=None,
    edgecolor="black",
    linewidth=1.4,
    alpha=0.75)

# make x positions manually so dots can be coloured by replicate
x_lookup = {indel: i for i, indel in enumerate(indel_order)}

hue_offsets = {
    "in-frame": -0.25,
    "frameshift": 0,
    "uncut": 0.25}

rng = np.random.default_rng(42)

for cls in class_order:
    sub_cls = sgRNF168_dots[sgRNF168_dots["class"] == cls].copy()

    # plot replicate dots, coloured by replicate
    for rep, sub_rep in sub_cls.groupby("replicate"):
        x = sub_rep["indels"].map(x_lookup).astype(float)
        x = x + hue_offsets[cls]
        x = x + rng.normal(0, 0.035, size=len(sub_rep))

        ax.scatter(
            x,
            sub_rep["relative_percentage"],
            s=120,
            facecolors=rep_palette[rep],
            edgecolors="black",
            linewidths=1.8,
            zorder=5,
            label=f"Rep {rep}"
        )

    # add mean ± SEM error bars for each indel within this class
    for indel in indel_order:
        sub = sgRNF168_df[
            (sgRNF168_df["indels"] == indel) &
            (sgRNF168_df["class"] == cls)
        ]

        vals = sub["relative_percentage"].to_numpy()

        # skip if there are no values or all values are zero
        if vals.size == 0 or np.all(vals == 0):
            continue

        mean_val = vals.mean()

        if vals.size > 1:
            sem = vals.std(ddof=1) / np.sqrt(vals.size)
        else:
            sem = 0

        xpos = x_lookup[indel] + hue_offsets[cls]

        ax.errorbar(
            xpos,
            mean_val,
            yerr=sem,
            fmt="none",
            ecolor="black",
            elinewidth=1.6,
            capsize=5,
            capthick=1.6,
            zorder=6)


handles, labels = ax.get_legend_handles_labels()
n_classes = len(class_order)

ax.legend(
    handles[:n_classes],
    labels[:n_classes],
    title="Class",
    frameon=False,
    fontsize=12,
    title_fontsize=13)

ax.set_xlabel("Indel size", fontsize=3)
ax.set_ylabel("Percentage", fontsize=3)

ax.tick_params(axis="x", labelsize=22, width=1.6, length=4, rotation=90)
ax.tick_params(axis="y", labelsize=22, width=1.6, length=4)

ax.spines["left"].set_linewidth(2.0)
ax.spines["bottom"].set_linewidth(2.0)
ax.spines["left"].set_color("black")
ax.spines["bottom"].set_color("black")

sns.despine()
plt.tight_layout()
plt.show()


#%% plotting sgLaminB1


# restrict x-axis to indels between -15 and 5
plot_df = sgLaminB1_df[
    (sgLaminB1_df["indels"] >= -15) &
    (sgLaminB1_df["indels"] <= 5)
].copy()

plot_dots = (
    plot_df
    .groupby(["indels", "class"], group_keys=False, observed=True)
    .filter(lambda x: (x["relative_percentage"] != 0).any())
    .copy()
)

indel_order = sorted(plot_df["indels"].unique())
class_order = ["in-frame", "frameshift", "uncut"]

# Use this for the manual scatter section
sgLB1_dots = plot_dots.copy()



plt.figure(figsize=(14, 5))
sns.set_context("talk")

ax = sns.barplot(
    data=plot_df,
    x="indels",
    y="relative_percentage",
    hue="class",
    order=indel_order,
    hue_order=class_order,
    palette=custom_palette,
    estimator="mean",
    errorbar=None,
    edgecolor="black",
    linewidth=1.4,
    alpha=0.75)


x_lookup = {indel: i for i, indel in enumerate(indel_order)}
hue_offsets = {
    "in-frame": -0.25,
    "frameshift": 0,
    "uncut": 0.25}


rng = np.random.default_rng(42)

for cls in class_order:
    sub_cls = sgLB1_dots[sgLB1_dots["class"] == cls].copy()

    # plot replicate dots, coloured by replicate
    for rep, sub_rep in sub_cls.groupby("replicate"):
        x = sub_rep["indels"].map(x_lookup).astype(float)
        x = x + hue_offsets[cls]
        x = x + rng.normal(0, 0.035, size=len(sub_rep))

        ax.scatter(
            x,
            sub_rep["relative_percentage"],
            s=120,
            facecolors=rep_palette[rep],
            edgecolors="black",
            linewidths=1.8,
            zorder=5,
            label=f"Rep {rep}"
        )

    # add mean ± SEM error bars for each indel within this class
    for indel in indel_order:
        sub = sgLaminB1_df[
            (sgLaminB1_df["indels"] == indel) &
            (sgLaminB1_df["class"] == cls)
        ]

        vals = sub["relative_percentage"].to_numpy()

        # skip if there are no values or all values are zero
        if vals.size == 0 or np.all(vals == 0):
            continue

        mean_val = vals.mean()

        if vals.size > 1:
            sem = vals.std(ddof=1) / np.sqrt(vals.size)
        else:
            sem = 0

        xpos = x_lookup[indel] + hue_offsets[cls]

        ax.errorbar(
            xpos,
            mean_val,
            yerr=sem,
            fmt="none",
            ecolor="black",
            elinewidth=1.6,
            capsize=5,
            capthick=1.6,
            zorder=6)

handles, labels = ax.get_legend_handles_labels()
n_classes = len(class_order)

ax.legend(
    handles[:n_classes],
    labels[:n_classes],
    title="Class",
    frameon=False,
    fontsize=12,
    title_fontsize=13)

ax.set_xlabel("Indel size", fontsize=3)
ax.set_ylabel("Percentage", fontsize=3)

ax.tick_params(axis="x", labelsize=22, width=1.6, length=4, rotation=90)
ax.tick_params(axis="y", labelsize=22, width=1.6, length=4)

ax.spines["left"].set_linewidth(2.0)
ax.spines["bottom"].set_linewidth(2.0)
ax.spines["left"].set_color("black")
ax.spines["bottom"].set_color("black")

sns.despine()
plt.tight_layout()
plt.show()
