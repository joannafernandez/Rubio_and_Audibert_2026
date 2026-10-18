#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 11:36:16 2026

@author: joannafernandez

written with assistance from GPT 5.5
"""

#%%
#run this code in batches. This is a file with plottung and paths to ALL TIDE sgLAD cut experiments (with treatment)

#%%S139A
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os 

#import data file and name cols 
def Find(file):
    genes = pd.read_csv(file, delimiter=",")
    genes.rename(columns={'Unnamed: 0': 'replicate'}, inplace=True)
    genes.rename(columns={'Unnamed: 1': 'indels'}, inplace=True)

    return genes

####s139a only
R0perc = Find('/Users/patricfernandez/Downloads/s139a_only_percentages.csv')
R0perc= R0perc.dropna() #remove empty rows
R0pval = Find('/Users/patricfernandez/Downloads/s139a_only_pvals.csv')
R0pval= R0pval.dropna() #remove empty rows

R0perc.loc[:, ~R0perc.columns.isin(['replicate', 'indels'])] = R0perc.loc[:, ~R0perc.columns.isin(['replicate', 'indels'])].mask(R0pval > 0.001, 0)


melted_df = R0perc.melt(id_vars=['replicate', 'indels'],
                               var_name='condition',
                               value_vars= [ 'WT_24hr', 'WT_48hr', 'S139A_24hr', 'S139A_48hr'],
                               value_name='percentage')
melted_df['percentage'] = melted_df['percentage'].astype(int)


#assign indel type 
def classify_indel(indel):
    if indel == 0:
        return "uncut"
    elif -4 <= indel <= 2:
        return "NHEJ"
    elif -20 <= indel <= -3:
        return "MMEJ"
    else:
        return None


melted_df['class'] = melted_df['indels'].apply(classify_indel)
melted_df = melted_df[melted_df['class'].notna()]

#find total % of sequences that are either NHEJ or MMEJ. (uncut/wt sequences are only ever "0")
summed_df = (
    melted_df
    .groupby(['condition', 'replicate', 'class'], as_index=False)
    .agg(summed_percentage=('percentage', 'sum'))
)


subset = summed_df.copy()

#Calculate total NHEJ + uncut per condition/replicate
total = subset.groupby(['condition', 'replicate'])['summed_percentage'].transform('sum')
#Calculate percentage within that total
subset['relative_percentage'] = (subset['summed_percentage'] / total) * 100
# Split 'condition' into two new columns: 'condition' and 'time'
subset[['condition', 'time']] = subset['condition'].str.split('_', n=1, expand=True)


#find total % of sequences that are either NHEJ or MMEJ. (uncut/wt sequences are only ever "0")
summed_df = (
    melted_df
    .groupby(['condition', 'replicate', 'class'], as_index=False)
    .agg(summed_percentage=('percentage', 'sum'))
)

subset = summed_df[summed_df['class'].isin(['MMEJ', 'NHEJ', 'uncut'])].copy()

#Calculate total NHEJ + uncut per condition/replicate
total = subset.groupby(['condition', 'replicate'])['summed_percentage'].transform('sum')
#Calculate percentage within that total
subset['relative_percentage'] = (subset['summed_percentage'] / total) * 100
# Split 'condition' into two new columns: 'condition' and 'time'
subset[['condition', 'time']] = subset['condition'].str.split('_', n=1, expand=True)

import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests


# --- Collapse MMEJ + NHEJ into "mutation" ---
mutation_df = (
    subset[subset['class'].isin(['NHEJ','MMEJ'])]
        .groupby(['time','condition','replicate'], as_index=False)['relative_percentage']
        .sum()
)

# Rename for clarity
mutation_df = mutation_df.rename(columns={'relative_percentage': 'mutation_percentage'})

df = mutation_df.copy()

df['group'] = df['condition'].astype(str) + '_' + df['time'].astype(str)


df['group'] = pd.Categorical(
    df['group'],
    categories=['WT_24hr', 'S139A_24hr', 'WT_48hr', 'S139A_48hr'],
    ordered=True
)

model = smf.ols(
    'mutation_percentage ~ C(group)',
    data=df
).fit()

anova_tbl = sm.stats.anova_lm(model, typ=2)
print(anova_tbl)

comparisons = {
    '24hr_S139A_vs_WT': 'C(group)[T.S139A_24hr] = 0',
    '48hr_S139A_vs_WT': 'C(group)[T.S139A_48hr] - C(group)[T.WT_48hr] = 0'
}

rows = []

for name, contrast in comparisons.items():
    tt = model.t_test(contrast)
    ci_low, ci_high = tt.conf_int().ravel()

    rows.append({
        'comparison': name,
        'difference': float(tt.effect),
        'SE': float(tt.sd),
        't': float(tt.tvalue),
        'p_raw': float(tt.pvalue),
        'CI95_low': ci_low,
        'CI95_high': ci_high
    })

posthoc = pd.DataFrame(rows)
posthoc['p_sidak'] = multipletests(posthoc['p_raw'], method='sidak')[1]

print(posthoc)

#########plotting


# 1) Collapse to mutation vs uncut

df0 = subset.copy()

df_mut = (
    df0[df0["class"].isin(["NHEJ", "MMEJ"])]
        .groupby(["time", "condition", "replicate"], as_index=False)["relative_percentage"]
        .sum()
        .assign(**{"class": "mutation"})
)

df_uncut = (
    df0[df0["class"] == "uncut"]
        .loc[:, ["time", "condition", "replicate", "relative_percentage"]]
        .assign(**{"class": "uncut"})
)

df = pd.concat([df_mut, df_uncut], ignore_index=True)


# 2) Categorical order + palette

time_order = ["24hr", "48hr"]          # or your existing time_order
condition_order = ["WT", "S139A"]        # or your existing condition_order
class_order = ["mutation", "uncut"]    # NEW

df["time"] = pd.Categorical(df["time"], categories=time_order, ordered=True)
df["condition"] = pd.Categorical(df["condition"], categories=condition_order, ordered=True)
df["class"] = pd.Categorical(df["class"], categories=class_order, ordered=True)

# palette (use yours if you want)
custom_palette = {
    "mutation": "#BBBBBB",   # change if desired
    "uncut": "#EFB54F"
}


# 3) Means for stacked bars

means = (
    df.groupby(["time", "condition", "class"], observed=True)["relative_percentage"]
      .mean()
      .unstack("class")
      .reindex(index=pd.MultiIndex.from_product([time_order, condition_order],
                                                names=["time", "condition"]),
               columns=class_order)
      .fillna(0.0)
)


#okay so seperat ebut with colours 


rep_palette = {
    1: "#ed3658",
    2: "#36bbac",
    3: "#FFFFFF",
    4: "#EFD45B"
}

fig, axes = plt.subplots(
    nrows=1,
    ncols=2,
    figsize=(7.2, 4.8),
    sharey=False,
    gridspec_kw={"wspace": 0.30}
)

x_pos = {
    "WT": 0,
    "S139A": 1
}

bar_width = 0.65

# Reproducible jitter
rng = np.random.default_rng(42)

for ax, t in zip(axes, time_order):

    for cond in condition_order:

        x = x_pos[cond]
        hrow = means.loc[(t, cond)]
        bottom = 0.0

        # Stacked mutation + uncut bars
        for cls in class_order:

            h = float(hrow.get(cls, 0.0))

            ax.bar(
                x,
                h,
                width=bar_width,
                bottom=bottom,
                color=custom_palette[cls],
                edgecolor="black",
                linewidth=1.5,
                alpha=0.8,
                zorder=2
            )

            bottom += h

        # Replicate points for mutation
        dft = df[
            (df["time"] == t) &
            (df["condition"] == cond) &
            (df["class"] == "mutation")
        ].copy()

        if not dft.empty:

            dft["replicate"] = dft["replicate"].astype(int)

            jitter = rng.uniform(
                -bar_width * 0.12,
                bar_width * 0.12,
                size=len(dft)
            )

            for i, (_, row) in enumerate(dft.iterrows()):

                ax.scatter(
                    x + jitter[i],
                    row["relative_percentage"],
                    s=180,
                    facecolors=rep_palette[row["replicate"]],
                    edgecolors="black",
                    linewidths=1.8,
                    zorder=4,
                    clip_on=False
                )

            vals = dft["relative_percentage"].to_numpy(dtype=float)

            if vals.size > 1:

                mean_val = vals.mean()
                sem = vals.std(ddof=1) / np.sqrt(vals.size)

                ax.errorbar(
                    x,
                    mean_val,
                    yerr=sem,
                    fmt="none",
                    ecolor="black",
                    elinewidth=1.8,
                    capsize=6,
                    capthick=1.8,
                    zorder=5
                )

    # Axis formatting
    ax.set_xticks([x_pos["WT"], x_pos["S139A"]])
    ax.set_xticklabels(["WT", "S139A"], fontsize=11)

    ax.set_title(
        t.replace("hr", " hr"),
        fontsize=14,
        pad=10
    )

    ax.set_xlim(-0.55, 1.55)
    ax.set_ylim(0, 102)
    ax.set_yticks([0, 20, 40, 60, 80, 100])

    ax.axhline(
        50,
        linestyle="--",
        color="grey",
        linewidth=1,
        zorder=1
    )

    ax.spines["left"].set_linewidth(2.0)
    ax.spines["bottom"].set_linewidth(2.0)

    ax.tick_params(
        axis="both",
        width=1.6,
        length=4
    )

    ax.tick_params(
        axis="y",
        labelsize=14
    )

    ax.set_xlabel("")

axes[0].set_ylabel("Relative percentage", fontsize=14)
#axes[1].set_ylabel("Relative percentage", fontsize=14)

sns.despine(fig=fig)

plt.tight_layout()
plt.show()


#%%


import data file and name cols 
def Find(file):
    genes = pd.read_csv(file, delimiter=",")
    genes.rename(columns={'Unnamed: 0': 'replicate'}, inplace=True)
    genes.rename(columns={'Unnamed: 1': 'indels'}, inplace=True)

    return genes


R0perc = Find('/Users/patricfernandez/Downloads/s139a53bp1_perc.csv')
R0perc= R0perc.dropna() #remove empty rows
R0pval = Find('/Users/patricfernandez/Downloads/s193a53bp1_pvals.csv')
R0pval= R0pval.dropna() #remove empty rows
#where p>0.001 (non-significant indel) remove from consideration
R0perc.loc[:, ~R0perc.columns.isin(['replicate', 'indels'])] = R0perc.loc[:, ~R0perc.columns.isin(['replicate', 'indels'])].mask(R0pval > 0.001, 0)

#formating for plotting
melted_df = R0perc.melt(id_vars=['replicate', 'indels'],
                               var_name='condition',
                               value_vars= [ 'WT_24hr', 'WT_48hr', 'S139A_24hr', 'S139A_48hr','53BP1KO_24hr', '53BP1KO_48hr'],
                               value_name='percentage')
melted_df['percentage'] = melted_df['percentage'].astype(int)


#assign indel type 
def classify_indel(indel):
    if indel == 0:
        return "uncut"
    elif -4 <= indel <= 2:
        return "NHEJ"
    elif -20 <= indel <= -3:
        return "MMEJ"
    else:
        return None


melted_df['class'] = melted_df['indels'].apply(classify_indel)
melted_df = melted_df[melted_df['class'].notna()]

#find total % of sequences that are either NHEJ or MMEJ. (uncut/wt sequences are only ever "0")
summed_df = (
    melted_df
    .groupby(['condition', 'replicate', 'class'], as_index=False)
    .agg(summed_percentage=('percentage', 'sum'))
)


subset = summed_df.copy()

#Calculate total NHEJ + uncut per condition/replicate
total = subset.groupby(['condition', 'replicate'])['summed_percentage'].transform('sum')
#Calculate percentage within that total
subset['relative_percentage'] = (subset['summed_percentage'] / total) * 100
# Split 'condition' into two new columns: 'condition' and 'time'
subset[['condition', 'time']] = subset['condition'].str.split('_', n=1, expand=True)


#find total % of sequences that are either NHEJ or MMEJ. (uncut/wt sequences are only ever "0")
summed_df = (
    melted_df
    .groupby(['condition', 'replicate', 'class'], as_index=False)
    .agg(summed_percentage=('percentage', 'sum'))
)

subset = summed_df[summed_df['class'].isin(['MMEJ', 'NHEJ', 'uncut'])].copy()

#Calculate total NHEJ + uncut per condition/replicate
total = subset.groupby(['condition', 'replicate'])['summed_percentage'].transform('sum')
#Calculate percentage within that total
subset['relative_percentage'] = (subset['summed_percentage'] / total) * 100
# Split 'condition' into two new columns: 'condition' and 'time'
subset[['condition', 'time']] = subset['condition'].str.split('_', n=1, expand=True)


# --- Collapse MMEJ + NHEJ into "mutation" ---
mutation_df = (
    subset[subset['class'].isin(['NHEJ','MMEJ'])]
        .groupby(['time','condition','replicate'], as_index=False)['relative_percentage']
        .sum()
)

# Rename for clarity
mutation_df = mutation_df.rename(columns={'relative_percentage': 'mutation_percentage'})

df = mutation_df.copy()

df['group'] = df['condition'].astype(str) + '_' + df['time'].astype(str)


########
df['group'] = pd.Categorical(
    df['group'],
    categories=['WT_24hr', '53BP1KO_24hr', 'WT_48hr', '53BP1KO_48hr'],
    ordered=True
)


model = smf.ols(
    'mutation_percentage ~ C(group)',
    data=df
).fit()

anova_tbl = sm.stats.anova_lm(model, typ=2)
print(anova_tbl)

comparisons = {
    '24hr_53BP1KO_vs_WT': 'C(group)[T.53BP1KO_24hr] = 0',
    '48hr_53BP1KO_vs_WT': 'C(group)[T.53BP1KO_48hr] - C(group)[T.WT_48hr] = 0'
}

rows = []

for name, contrast in comparisons.items():
    tt = model.t_test(contrast)
    ci_low, ci_high = tt.conf_int().ravel()

    rows.append({
        'comparison': name,
        'difference': float(tt.effect),
        'SE': float(tt.sd),
        't': float(tt.tvalue),
        'p_raw': float(tt.pvalue),
        'CI95_low': ci_low,
        'CI95_high': ci_high
    })

posthoc = pd.DataFrame(rows)
posthoc['p_sidak'] = multipletests(posthoc['p_raw'], method='sidak')[1]

print(posthoc)

# 1) Collapse to mutation vs uncut

df0 = subset.copy()

df_mut = (
    df0[df0["class"].isin(["NHEJ", "MMEJ"])]
        .groupby(["time", "condition", "replicate"], as_index=False)["relative_percentage"]
        .sum()
        .assign(**{"class": "mutation"})
)

df_uncut = (
    df0[df0["class"] == "uncut"]
        .loc[:, ["time", "condition", "replicate", "relative_percentage"]]
        .assign(**{"class": "uncut"})
)

df = pd.concat([df_mut, df_uncut], ignore_index=True)

# 2) Categorical order + palette

time_order = ["24hr", "48hr"]          # or your existing time_order
condition_order = ["WT", "53BP1KO"]        # or your existing condition_order
class_order = ["mutation", "uncut"]    # NEW

df["time"] = pd.Categorical(df["time"], categories=time_order, ordered=True)
df["condition"] = pd.Categorical(df["condition"], categories=condition_order, ordered=True)
df["class"] = pd.Categorical(df["class"], categories=class_order, ordered=True)

# palette (use yours if you want)
custom_palette = {
    "mutation": "#BBBBBB",   # change if desired
    "uncut": "#EFB54F"
}


# 3) Means for stacked bars

means = (
    df.groupby(["time", "condition", "class"], observed=True)["relative_percentage"]
      .mean()
      .unstack("class")
      .reindex(index=pd.MultiIndex.from_product([time_order, condition_order],
                                                names=["time", "condition"]),
               columns=class_order)
      .fillna(0.0)
)


rep_palette = {
    1: "#ed3658",
    2: "#36bbac",
    3: "#FFFFFF",
    4: "#EFD45B"
}

fig, axes = plt.subplots(
    nrows=1,
    ncols=2,
    figsize=(7.2, 4.8),
    sharey=False,
    gridspec_kw={"wspace": 0.30}
)

x_pos = {
    "WT": 0,
    "53BP1KO": 1
}

bar_width = 0.65

# Reproducible jitter
rng = np.random.default_rng(42)

for ax, t in zip(axes, time_order):

    for cond in condition_order:

        x = x_pos[cond]
        hrow = means.loc[(t, cond)]
        bottom = 0.0

        # Stacked mutation + uncut bars
        for cls in class_order:

            h = float(hrow.get(cls, 0.0))

            ax.bar(
                x,
                h,
                width=bar_width,
                bottom=bottom,
                color=custom_palette[cls],
                edgecolor="black",
                linewidth=1.5,
                alpha=0.8,
                zorder=2
            )

            bottom += h

        # Replicate points for mutation
        dft = df[
            (df["time"] == t) &
            (df["condition"] == cond) &
            (df["class"] == "mutation")
        ].copy()

        if not dft.empty:

            dft["replicate"] = dft["replicate"].astype(int)

            jitter = rng.uniform(
                -bar_width * 0.12,
                bar_width * 0.12,
                size=len(dft)
            )

            for i, (_, row) in enumerate(dft.iterrows()):

                ax.scatter(
                    x + jitter[i],
                    row["relative_percentage"],
                    s=180,
                    facecolors=rep_palette[row["replicate"]],
                    edgecolors="black",
                    linewidths=1.8,
                    zorder=4,
                    clip_on=False
                )

            vals = dft["relative_percentage"].to_numpy(dtype=float)

            if vals.size > 1:

                mean_val = vals.mean()
                sem = vals.std(ddof=1) / np.sqrt(vals.size)

                ax.errorbar(
                    x,
                    mean_val,
                    yerr=sem,
                    fmt="none",
                    ecolor="black",
                    elinewidth=1.8,
                    capsize=6,
                    capthick=1.8,
                    zorder=5
                )

    # Axis formatting
    ax.set_xticks([x_pos["WT"], x_pos["53BP1KO"]])
    ax.set_xticklabels(["WT", "53BP1KO"], fontsize=11)

    ax.set_title(
        t.replace("hr", " hr"),
        fontsize=14,
        pad=10
    )

    ax.set_xlim(-0.55, 1.55)
    ax.set_ylim(0, 102)
    ax.set_yticks([0, 20, 40, 60, 80, 100])

    ax.axhline(
        50,
        linestyle="--",
        color="grey",
        linewidth=1,
        zorder=1
    )

    ax.spines["left"].set_linewidth(2.0)
    ax.spines["bottom"].set_linewidth(2.0)

    ax.tick_params(
        axis="both",
        width=1.6,
        length=4
    )

    ax.tick_params(
        axis="y",
        labelsize=14
    )

    ax.set_xlabel("")

axes[0].set_ylabel("Relative percentage", fontsize=14)
#axes[1].set_ylabel("Relative percentage", fontsize=14)

sns.despine(fig=fig)

plt.tight_layout()
plt.show()




#%%

######################################################
#ATMi sgLAD


#import data file and name cols 
def Find(file):
    genes = pd.read_csv(file, delimiter=",")
    genes.rename(columns={'Unnamed: 0': 'replicate'}, inplace=True)
    genes.rename(columns={'Unnamed: 1': 'indels'}, inplace=True)

    return genes

####s139a only
R0perc = Find("/Users/patricfernandez/Downloads/ATMi_percentages.csv")
R0perc= R0perc.dropna() #remove empty rows
R0pval = Find("/Users/patricfernandez/Downloads/ATMi_pvals.csv")
R0pval= R0pval.dropna() #remove empty rows

R0perc.loc[:, ~R0perc.columns.isin(['replicate', 'indels'])] = R0perc.loc[:, ~R0perc.columns.isin(['replicate', 'indels'])].mask(R0pval > 0.001, 0)


melted_df = R0perc.melt(id_vars=['replicate', 'indels'],
                               var_name='condition',
                               value_vars= ['DMSO_6hr', 'DMSO_12hr', 'DMSO_24hr',
                                      'DMSO_48hr', 'ATMi_6hr', 'ATMi_12hr', 'ATMi_24hr', 'ATMi_48hr'],
                               value_name='percentage')
melted_df['percentage'] = melted_df['percentage'].astype(int)


#assign indel type 
def classify_indel(indel):
    if indel == 0:
        return "uncut"
    elif -4 <= indel <= 2:
        return "NHEJ"
    elif -20 <= indel <= -3:
        return "MMEJ"
    else:
        return None


melted_df['class'] = melted_df['indels'].apply(classify_indel)
melted_df = melted_df[melted_df['class'].notna()]

#find total % of sequences that are either NHEJ or MMEJ. (uncut/wt sequences are only ever "0")
summed_df = (
    melted_df
    .groupby(['condition', 'replicate', 'class'], as_index=False)
    .agg(summed_percentage=('percentage', 'sum'))
)


subset = summed_df.copy()

#Calculate total NHEJ + uncut per condition/replicate
total = subset.groupby(['condition', 'replicate'])['summed_percentage'].transform('sum')
#Calculate percentage within that total
subset['relative_percentage'] = (subset['summed_percentage'] / total) * 100
# Split 'condition' into two new columns: 'condition' and 'time'
subset[['condition', 'time']] = subset['condition'].str.split('_', n=1, expand=True)


#find total % of sequences that are either NHEJ or MMEJ. (uncut/wt sequences are only ever "0")
summed_df = (
    melted_df
    .groupby(['condition', 'replicate', 'class'], as_index=False)
    .agg(summed_percentage=('percentage', 'sum'))
)

subset = summed_df[summed_df['class'].isin(['MMEJ', 'NHEJ', 'uncut'])].copy()

#Calculate total NHEJ + uncut per condition/replicate
total = subset.groupby(['condition', 'replicate'])['summed_percentage'].transform('sum')
#Calculate percentage within that total
subset['relative_percentage'] = (subset['summed_percentage'] / total) * 100
# Split 'condition' into two new columns: 'condition' and 'time'
subset[['condition', 'time']] = subset['condition'].str.split('_', n=1, expand=True)




mutation_df = (
    subset[subset["class"].isin(["NHEJ", "MMEJ"])]
        .groupby(["time", "condition", "replicate"], as_index=False)["relative_percentage"]
        .sum()
)

mutation_df = mutation_df.rename(
    columns={"relative_percentage": "mutation_percentage"}
)

df = mutation_df.copy()


# 2) Make combined group variable


time_order = ["6hr", "12hr", "24hr", "48hr"]
condition_order = ["DMSO", "ATMi"]

df["time"] = pd.Categorical(
    df["time"],
    categories=time_order,
    ordered=True
)

df["condition"] = pd.Categorical(
    df["condition"],
    categories=condition_order,
    ordered=True
)

df["group"] = df["condition"].astype(str) + "_" + df["time"].astype(str)

group_order = [
    "DMSO_6hr", "ATMi_6hr",
    "DMSO_12hr", "ATMi_12hr",
    "DMSO_24hr", "ATMi_24hr",
    "DMSO_48hr", "ATMi_48hr"
]

df["group"] = pd.Categorical(
    df["group"],
    categories=group_order,
    ordered=True
)

# 3) One-way ANOVA


model = smf.ols(
    "mutation_percentage ~ C(group)",
    data=df
).fit()

anova_tbl = sm.stats.anova_lm(model, typ=2)
print(anova_tbl)


# 4) Planned comparisons:
#    ATMi vs DMSO at each timepoint


comparisons = {
    "6hr_ATMi_vs_DMSO":  "C(group)[T.ATMi_6hr] = 0",
    "12hr_ATMi_vs_DMSO": "C(group)[T.ATMi_12hr] - C(group)[T.DMSO_12hr] = 0",
    "24hr_ATMi_vs_DMSO": "C(group)[T.ATMi_24hr] - C(group)[T.DMSO_24hr] = 0",
    "48hr_ATMi_vs_DMSO": "C(group)[T.ATMi_48hr] - C(group)[T.DMSO_48hr] = 0"
}

rows = []

for name, contrast in comparisons.items():
    tt = model.t_test(contrast)
    ci_low, ci_high = tt.conf_int().ravel()

    rows.append({
        "comparison": name,
        "diff_ATMi_minus_DMSO": float(tt.effect),
        "SE": float(tt.sd),
        "t": float(tt.tvalue),
        "p_raw": float(tt.pvalue),
        "CI95_low": ci_low,
        "CI95_high": ci_high
    })

posthoc = pd.DataFrame(rows)

# Šidák correction for the 4 planned comparisons
posthoc["p_sidak"] = multipletests(
    posthoc["p_raw"],
    method="sidak"
)[1]

print(posthoc)



#%%



df0 = subset.copy()

df_mut = (
    df0[df0["class"].isin(["NHEJ", "MMEJ"])]
        .groupby(["time", "condition", "replicate"], as_index=False)["relative_percentage"]
        .sum()
        .assign(**{"class": "mutation"})
)

df_uncut = (
    df0[df0["class"] == "uncut"]
        .loc[:, ["time", "condition", "replicate", "relative_percentage"]]
        .assign(**{"class": "uncut"})
)

df = pd.concat([df_mut, df_uncut], ignore_index=True)

# 2) Categorical order + palette

time_order = ["6hr", "12hr", "24hr", "48hr"]          # or your existing time_order
condition_order = ["DMSO", "ATMi"]        # or your existing condition_order
class_order = ["mutation", "uncut"]    # NEW

df["time"] = pd.Categorical(df["time"], categories=time_order, ordered=True)
df["condition"] = pd.Categorical(df["condition"], categories=condition_order, ordered=True)
df["class"] = pd.Categorical(df["class"], categories=class_order, ordered=True)

# palette (use yours if you want)
custom_palette = {
    "mutation": "#BBBBBB",   # change if desired
    "uncut": "#EFB54F"
}


# 3) Means for stacked bars

means = (
    df.groupby(["time", "condition", "class"], observed=True)["relative_percentage"]
      .mean()
      .unstack("class")
      .reindex(index=pd.MultiIndex.from_product([time_order, condition_order],
                                                names=["time", "condition"]),
               columns=class_order)
      .fillna(0.0)
)
rep_palette = {
    1: "#ed3658",
    2: "#36bbac",
    5: "#FFFFFF"
}

fig, axes = plt.subplots(
    nrows=1,
    ncols=4,
    figsize=(7.2, 4.8),
    sharey=False,
    gridspec_kw={"wspace": 0.30}
)

x_pos = {
    "DMSO": 0,
    "ATMi": 1
}

bar_width = 0.65

# Reproducible jitter
rng = np.random.default_rng(42)

for ax, t in zip(axes, time_order):

    for cond in condition_order:

        x = x_pos[cond]
        hrow = means.loc[(t, cond)]
        bottom = 0.0

        # Stacked mutation + uncut bars
        for cls in class_order:

            h = float(hrow.get(cls, 0.0))

            ax.bar(
                x,
                h,
                width=bar_width,
                bottom=bottom,
                color=custom_palette[cls],
                edgecolor="black",
                linewidth=1.5,
                alpha=0.8,
                zorder=2
            )

            bottom += h

        # Replicate points for mutation
        dft = df[
            (df["time"] == t) &
            (df["condition"] == cond) &
            (df["class"] == "mutation")
        ].copy()

        if not dft.empty:

            dft["replicate"] = dft["replicate"].astype(int)

            jitter = rng.uniform(
                -bar_width * 0.12,
                bar_width * 0.12,
                size=len(dft)
            )

            for i, (_, row) in enumerate(dft.iterrows()):

                ax.scatter(
                    x + jitter[i],
                    row["relative_percentage"],
                    s=180,
                    facecolors=rep_palette[row["replicate"]],
                    edgecolors="black",
                    linewidths=1.8,
                    zorder=4,
                    clip_on=False
                )

            vals = dft["relative_percentage"].to_numpy(dtype=float)

            if vals.size > 1:

                mean_val = vals.mean()
                sem = vals.std(ddof=1) / np.sqrt(vals.size)

                ax.errorbar(
                    x,
                    mean_val,
                    yerr=sem,
                    fmt="none",
                    ecolor="black",
                    elinewidth=1.8,
                    capsize=6,
                    capthick=1.8,
                    zorder=5
                )

    # Axis formatting
    ax.set_xticks([x_pos["DMSO"], x_pos["ATMi"]])
    ax.set_xticklabels(["DMSO", "ATMi"], fontsize=11)

    ax.set_title(
        t.replace("hr", " hr"),
        fontsize=14,
        pad=10
    )

    ax.set_xlim(-0.55, 1.55)
    ax.set_ylim(0, 102)
    ax.set_yticks([0, 20, 40, 60, 80, 100])

    ax.axhline(
        50,
        linestyle="--",
        color="grey",
        linewidth=1,
        zorder=1
    )

    ax.spines["left"].set_linewidth(2.0)
    ax.spines["bottom"].set_linewidth(2.0)

    ax.tick_params(
        axis="both",
        width=1.6,
        length=4
    )

    ax.tick_params(
        axis="y",
        labelsize=14
    )

    ax.set_xlabel("")

axes[0].set_ylabel("Relative percentage", fontsize=14)
#axes[1].set_ylabel("Relative percentage", fontsize=14)

sns.despine(fig=fig)

plt.tight_layout()
plt.show()


#%%

#%%


###ATMi sgiLAD

atmiperc = Find("/Users/patricfernandez/ATMi_iLAD_perc.csv")
atmiperc= atmiperc.dropna() #remove empty rows
atmipval = Find("/Users/patricfernandez/ATMi_iLAD_pval.csv")
atmipval= atmipval.dropna() #remove empty rows

atmiperc.loc[:, ~atmiperc.columns.isin(['replicate', 'indels'])] = atmiperc.loc[:, ~atmiperc.columns.isin(['replicate', 'indels'])].mask(atmipval > 0.001, 0)


melted_df = atmiperc.melt(id_vars=['replicate', 'indels'],
                               var_name='condition',
                               value_vars= ['DMSO_6hr', 'DMSO_12hr', 'DMSO_24hr',
                                      'DMSO_48hr', 'ATMi_6hr', 'ATMi_12hr', 'ATMi_24hr', 'ATMi_48hr'],
                               value_name='percentage')
melted_df['percentage'] = melted_df['percentage'].astype(int)



#assign indel type 
def classify_indel(indel):
    if indel == 0:
        return "uncut"
    elif -4 <= indel <= 2:
        return "NHEJ"
    elif -20 <= indel <= -3:
        return "MMEJ"
    else:
        return None


melted_df['class'] = melted_df['indels'].apply(classify_indel)
melted_df = melted_df[melted_df['class'].notna()]

#find total % of sequences that are either NHEJ or MMEJ. (uncut/wt sequences are only ever "0")
summed_df = (
    melted_df
    .groupby(['condition', 'replicate', 'class'], as_index=False)
    .agg(summed_percentage=('percentage', 'sum'))
)


subset = summed_df.copy()

#Calculate total NHEJ + uncut per condition/replicate
total = subset.groupby(['condition', 'replicate'])['summed_percentage'].transform('sum')
#Calculate percentage within that total
subset['relative_percentage'] = (subset['summed_percentage'] / total) * 100
# Split 'condition' into two new columns: 'condition' and 'time'
subset[['condition', 'time']] = subset['condition'].str.split('_', n=1, expand=True)


#find total % of sequences that are either NHEJ or MMEJ. (uncut/wt sequences are only ever "0")
summed_df = (
    melted_df
    .groupby(['condition', 'replicate', 'class'], as_index=False)
    .agg(summed_percentage=('percentage', 'sum'))
)

subset = summed_df[summed_df['class'].isin(['MMEJ', 'NHEJ', 'uncut'])].copy()

#Calculate total NHEJ + uncut per condition/replicate
total = subset.groupby(['condition', 'replicate'])['summed_percentage'].transform('sum')
#Calculate percentage within that total
subset['relative_percentage'] = (subset['summed_percentage'] / total) * 100
# Split 'condition' into two new columns: 'condition' and 'time'
subset[['condition', 'time']] = subset['condition'].str.split('_', n=1, expand=True)




mutation_df = (
    subset[subset["class"].isin(["NHEJ", "MMEJ"])]
        .groupby(["time", "condition", "replicate"], as_index=False)["relative_percentage"]
        .sum()
)

mutation_df = mutation_df.rename(
    columns={"relative_percentage": "mutation_percentage"}
)

df = mutation_df.copy()



# 2) Make combined group variable
time_order = ["6hr", "12hr", "24hr", "48hr"]
condition_order = ["DMSO", "ATMi"]

df["time"] = pd.Categorical(
    df["time"],
    categories=time_order,
    ordered=True
)

df["condition"] = pd.Categorical(
    df["condition"],
    categories=condition_order,
    ordered=True
)

df["group"] = df["condition"].astype(str) + "_" + df["time"].astype(str)

group_order = [
    "DMSO_6hr", "ATMi_6hr",
    "DMSO_12hr", "ATMi_12hr",
    "DMSO_24hr", "ATMi_24hr",
    "DMSO_48hr", "ATMi_48hr"
]

df["group"] = pd.Categorical(
    df["group"],
    categories=group_order,
    ordered=True
)

# 3) One-way ANOVA


model = smf.ols(
    "mutation_percentage ~ C(group)",
    data=df
).fit()

anova_tbl = sm.stats.anova_lm(model, typ=2)
print(anova_tbl)

# 4) Planned comparisons:
#    ATMi vs DMSO at each timepoint

comparisons = {
    "6hr_ATMi_vs_DMSO":  "C(group)[T.ATMi_6hr] = 0",
    "12hr_ATMi_vs_DMSO": "C(group)[T.ATMi_12hr] - C(group)[T.DMSO_12hr] = 0",
    "24hr_ATMi_vs_DMSO": "C(group)[T.ATMi_24hr] - C(group)[T.DMSO_24hr] = 0",
    "48hr_ATMi_vs_DMSO": "C(group)[T.ATMi_48hr] - C(group)[T.DMSO_48hr] = 0"
}

rows = []

for name, contrast in comparisons.items():
    tt = model.t_test(contrast)
    ci_low, ci_high = tt.conf_int().ravel()

    rows.append({
        "comparison": name,
        "diff_ATMi_minus_DMSO": float(tt.effect),
        "SE": float(tt.sd),
        "t": float(tt.tvalue),
        "p_raw": float(tt.pvalue),
        "CI95_low": ci_low,
        "CI95_high": ci_high
    })

posthoc = pd.DataFrame(rows)

# Šidák correction for the 4 planned comparisons
posthoc["p_sidak"] = multipletests(
    posthoc["p_raw"],
    method="sidak"
)[1]

print(posthoc)


#%%


df0 = subset.copy()

df_mut = (
    df0[df0["class"].isin(["NHEJ", "MMEJ"])]
        .groupby(["time", "condition", "replicate"], as_index=False)["relative_percentage"]
        .sum()
        .assign(**{"class": "mutation"})
)

df_uncut = (
    df0[df0["class"] == "uncut"]
        .loc[:, ["time", "condition", "replicate", "relative_percentage"]]
        .assign(**{"class": "uncut"})
)

df = pd.concat([df_mut, df_uncut], ignore_index=True)

# 2) Categorical order + palette

time_order = ["6hr", "12hr", "24hr", "48hr"]          # or your existing time_order
condition_order = ["DMSO", "ATMi"]        # or your existing condition_order
class_order = ["mutation", "uncut"]    # NEW

df["time"] = pd.Categorical(df["time"], categories=time_order, ordered=True)
df["condition"] = pd.Categorical(df["condition"], categories=condition_order, ordered=True)
df["class"] = pd.Categorical(df["class"], categories=class_order, ordered=True)

# palette (use yours if you want)
custom_palette = {
    "mutation": "#BBBBBB",   # change if desired
    "uncut": "#EFB54F"
}


# 3) Means for stacked bars

means = (
    df.groupby(["time", "condition", "class"], observed=True)["relative_percentage"]
      .mean()
      .unstack("class")
      .reindex(index=pd.MultiIndex.from_product([time_order, condition_order],
                                                names=["time", "condition"]),
               columns=class_order)
      .fillna(0.0)
)



#%%
rep_palette = {
    1: "#ed3658",
    2: "#36bbac",
    3: "#FFFFFF"
}

fig, axes = plt.subplots(
    nrows=1,
    ncols=4,
    figsize=(7.2, 4.8),
    sharey=False,
    gridspec_kw={"wspace": 0.30}
)

x_pos = {
    "DMSO": 0,
    "ATMi": 1
}

bar_width = 0.65

# Reproducible jitter
rng = np.random.default_rng(42)

for ax, t in zip(axes, time_order):

    for cond in condition_order:

        x = x_pos[cond]
        hrow = means.loc[(t, cond)]
        bottom = 0.0

        # Stacked mutation + uncut bars
        for cls in class_order:

            h = float(hrow.get(cls, 0.0))

            ax.bar(
                x,
                h,
                width=bar_width,
                bottom=bottom,
                color=custom_palette[cls],
                edgecolor="black",
                linewidth=1.5,
                alpha=0.8,
                zorder=2
            )

            bottom += h

        # Replicate points for mutation
        dft = df[
            (df["time"] == t) &
            (df["condition"] == cond) &
            (df["class"] == "mutation")
        ].copy()

        if not dft.empty:

            dft["replicate"] = dft["replicate"].astype(int)

            jitter = rng.uniform(
                -bar_width * 0.12,
                bar_width * 0.12,
                size=len(dft)
            )

            for i, (_, row) in enumerate(dft.iterrows()):

                ax.scatter(
                    x + jitter[i],
                    row["relative_percentage"],
                    s=180,
                    facecolors=rep_palette[row["replicate"]],
                    edgecolors="black",
                    linewidths=1.8,
                    zorder=4,
                    clip_on=False
                )

            vals = dft["relative_percentage"].to_numpy(dtype=float)

            if vals.size > 1:

                mean_val = vals.mean()
                sem = vals.std(ddof=1) / np.sqrt(vals.size)

                ax.errorbar(
                    x,
                    mean_val,
                    yerr=sem,
                    fmt="none",
                    ecolor="black",
                    elinewidth=1.8,
                    capsize=6,
                    capthick=1.8,
                    zorder=5
                )

    # Axis formatting
    ax.set_xticks([x_pos["DMSO"], x_pos["ATMi"]])
    ax.set_xticklabels(["DMSO", "ATMi"], fontsize=11)

    ax.set_title(
        t.replace("hr", " hr"),
        fontsize=14,
        pad=10
    )

    ax.set_xlim(-0.55, 1.55)
    ax.set_ylim(0, 102)
    ax.set_yticks([0, 20, 40, 60, 80, 100])

    ax.axhline(
        50,
        linestyle="--",
        color="grey",
        linewidth=1,
        zorder=1
    )

    ax.spines["left"].set_linewidth(2.0)
    ax.spines["bottom"].set_linewidth(2.0)

    ax.tick_params(
        axis="both",
        width=1.6,
        length=4
    )

    ax.tick_params(
        axis="y",
        labelsize=14
    )

    ax.set_xlabel("")

axes[0].set_ylabel("Relative percentage", fontsize=14)
#axes[1].set_ylabel("Relative percentage", fontsize=14)

sns.despine(fig=fig)

plt.tight_layout()
plt.show()







