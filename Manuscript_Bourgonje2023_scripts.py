#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

#Annotate flagellin peptides from IBD (Bourgonje et al., Immunity 2023) and CFS (Vogl et al., Science Advances 2022) PhIP-Seq datasets
#IBD: supplementary table S1F from Bourgonje et al., Immunity 2023 (age- and sex-matched case control analysis of n=256 patients with CD and n=255 population controls)
ibd = pd.read_excel("Flagellin_paper/immunity_st1f.xlsx")

#Flagellin protein annotation (file "New_table_flagellum.xlsx")
import pandas as pd
flagellins = pd.read_excel("Flagellin_paper/New_table_flagellum.xlsx")
flagellins = flagellins.rename({'peptide_name': 'Antibody_peptide'}, axis=1)
singlecolumn = flagellins["Antibody_peptide"]
df = pd.DataFrame(singlecolumn)
df.to_excel("Flagellin_paper/flagellin_list_gabriel_aug2023.xlsx")


ibd['Flagellin'] = ibd['Antibody_peptide'].isin(df['Antibody_peptide']).map({True: 'Yes', False: 'No'})

#Scatter plot CD vs controls labeled by flagellin peptides
colors_list = ['darkturquoise', 'grey']

sns.set(rc={'figure.figsize':(7,7)}, font_scale = 1.5)
sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
sns.set_style("white")
Plot = sns.scatterplot(x="%LL", y="%CD", hue="Flagellin", hue_order=["Yes", "No"], s=80, palette=colors_list, data=ibd)
plt.xlabel("% of HC-NL individuals in whom\n a peptide is significantly bound", size=18, labelpad=10)
plt.ylabel("% of CD patients in whom\n a peptide is significantly bound", size=18, labelpad=10)
plt.text(x=27, y=90, fontsize=16, s="Flagellin ratio overrepresentation\n$\it{P}$=4.6x$10^-$$^1$$^0$", bbox=dict(facecolor='red', alpha=0))
plt.title("CD vs. HC-NL", size=24, y=1.03, fontweight='bold')
plt.xlim(-1,101)
plt.ylim(-1,101)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
legend_handles, _= Plot.get_legend_handles_labels()
Plot.legend(legend_handles, ['Yes','No'], 
          bbox_to_anchor=(0.25,0.99), 
          title='Flagellin')
plt.setp(Plot.get_legend().get_texts(), fontsize='18')
plt.setp(Plot.get_legend().get_title(), fontsize='18')
Plot.figure.savefig("Flagellin_paper/CDvsHC_remake_flagellin.jpg", dpi=600, bbox_inches='tight')

#Test for differential abundance
df_flag = ibd.loc[ibd["Flagellin"] == 'Yes']
from scipy.stats import mannwhitneyu
stat, p = mannwhitneyu(df_flag['%CD'], df_flag['%LL'])
print('Statistics=%.3f, p=%.13f' % (stat, p))

#Same exercise for CFS dataset (file "exist_df.csv" derived from supplementary info of Vogl et al., Science Advances 2022)
cfs = pd.read_csv("Flagellin_paper/exist_df.csv")
cfs = cfs.set_index("SerumName")

#Filter
def Prevalence(Series):
    Prev = Series.sum()/Series.count()
    return(Prev)
Prev_thresh = 0.05
print("Computing antibody prevalence")
Prevalences = cfs.apply(Prevalence, axis=0)
Vector_keep = Prevalences > Prev_thresh #True/False vector of the features with prevalence > 5%
cfs_f = cfs[Vector_keep.index[Vector_keep]] #Removal of columns that do not pass the prevalence threshold
cfs = cfs.reset_index()

#Label groups (ME/CFS vs. controls)
def check_substring(text):
    if '_26_' in text:
        return 0
    elif '_27_' in text:
        return 1
    else:
        return None

cfs['cohort'] = cfs['SerumName'].apply(check_substring)
cfs = cfs.set_index("SerumName")

#Calculate frequencies and percentages of flagellin peptides
cohort_count = cfs.groupby(['cohort']).count()
cohort_sum = cfs.groupby(['cohort']).sum()
cohort_percentages = (cohort_sum / cohort_count) * 100
cohort_percentages = cohort_percentages.T
cohort_percentages = cohort_percentages.reset_index()
cohort_percentages.to_excel("Flagellin_paper/cohort_percentages_CFSpaper_withoutcutoff.xlsx")
pct = pd.read_excel("Flagellin_paper/cohort_percentages_CFSpaper_withoutcutoff.xlsx")
pct = pct.rename({0: '%ME/CFS'}, axis=1)
pct = pct.rename({1: '%HC-UK'}, axis=1)
del pct['Unnamed: 0']
pct['ratio'] = (pct['%ME/CFS'] / pct['%HC-UK'])
pct = pct.rename({'index': 'Unnamed: 0'}, axis=1)

#Annotate peptides (file "df_info_AT.csv" derived from Vogl et al., Nature Medicine 2021)
information = pd.read_csv("datasetspandas/df_info_AT.csv")
pct_info = pct.merge(information, on=["Unnamed: 0"])
del pct_info['aa_seq']
del pct_info['is_IEDB']
del pct_info['is_pos_cntrl']
del pct_info['is_neg_cntrl']
del pct_info['is_phage']
del pct_info['is_influenza']
del pct_info['is_allergens']
del pct_info['is_genome_editing']
del pct_info['IEDB_DOIDs']
del pct_info['IEDB_comments']
del pct_info['IEDB_organism_name']
del pct_info['IEDB_parent_species']
del pct_info['is_rand_cntrl']
del pct_info['is_VFDB']
del pct_info['is_patho_strain']
del pct_info['is_IgA_coated_strain']
del pct_info['is_probio_strain']
del pct_info['bac_src']
del pct_info['is_gut_microbiome']
pct_info.to_excel("Flagellin_paper/CFSpaper_pct_info_withoutcutoff.xlsx")

#Scatterplot ME/CFS vs. controls labeled by flagellin
figure = pd.read_excel("Flagellin_paper/CFSpaper_pct_info.xlsx")
figure['Flagellin'] = figure['Unnamed: 0'].isin(df['Antibody_peptide']).map({True: 'Yes', False: 'No'})
figure_sorted = figure.sort_values(by='Flagellin', ascending=False)
colors_list = ['darkturquoise', 'grey']
sns.set(rc={'figure.figsize':(7,7)}, font_scale = 1.5)
sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
sns.set_style("white")
plt.scatter(x="%HC-UK", y="%ME/CFS", label="Flagellin", s=50, c='grey', data=figure_sorted[figure_sorted['Flagellin'] == 'No'], zorder=1)
plt.scatter(x="%HC-UK", y="%ME/CFS", label="Flagellin", s=50, c='darkturquoise', data=figure_sorted[figure_sorted['Flagellin'] == 'Yes'], zorder=2)
plt.xlabel("% of HC-UK individuals in whom\n a peptide is significantly bound", size=18, labelpad=10)
plt.ylabel("% of ME/CFS patients in whom\n a peptide is significantly bound", size=18, labelpad=10)
plt.title("ME/CFS vs. HC-UK", size=24, y=1.03, fontweight='bold')
plt.text(x=27, y=90, fontsize=16, s="Flagellin ratio overrepresentation\n$\it{P}$=8.9x$10^-$$^1$$^0$", bbox=dict(facecolor='red', alpha=0))
plt.xlim(-1,101)
plt.ylim(-1,101)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.savefig("Flagellin_paper/MECFSvsHC_remake_flagellin.jpg", dpi=600, bbox_inches='tight')

df_flagcfs = figure_sorted.loc[figure_sorted["Flagellin"] == 'Yes']

#Test for differences across cohorts
from scipy.stats import mannwhitneyu
stat, p = mannwhitneyu(df_flagcfs['%ME/CFS'], df_flagcfs['%HC-UK'])
print('Statistics=%.3f, p=%.13f' % (stat, p))

#Age and sex distributions
figuredata = pd.read_excel("Flagellin_paper/Agesexpanel_cohorts.xlsx") #Aggregated age/sex data from IBD (Bourgonje et al., Immunity 2023; Andreu-Sanchez et al., Immunity 2023) and CFS (Vogl et al., Science Advances 2022) datasets
figuredata = figuredata.set_index('sample_id')

sns.set(rc={'figure.figsize':(12,6)}, font_scale = 1.5)
sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
sns.set_style("white")
custom_palette = ["crimson", "dodgerblue", "purple", "limegreen"]
plot = sns.barplot(data=figuredata, x='Age_group', y="Percentage_age", hue="Cohort", palette=custom_palette)
plot.set_title("Age group distribution", y=1.02, fontsize=24, fontweight='bold')
xticklabels = ["18-25", "26-30", "31-35", "36-40", "41-45", "46-50", "51-55", "56-60", "61-65", "66-70", "71-75", "76-80"]
plot.set_xticklabels(labels=xticklabels, rotation=35, ha='center', fontsize=22)
plot.set_ylabel("% of cohort total", fontsize=22)
plot.set_xlabel("")
plt.ylim(0,32)
plt.legend(bbox_to_anchor=(0.805, 0.35), fontsize=16)
plt.yticks(fontsize=22)
plot.figure.savefig("Flagellin_paper/Agedistributions_panel_part.jpg", dpi=600, bbox_inches='tight')


table = pd.crosstab(figuredata.Cohort, figuredata.Sex)
#Test for differences across cohorts
from scipy.stats import chi2_contingency
from scipy.stats import fisher_exact
chi2, p, dof, expected = chi2_contingency(table)
print(f"chi2 statistic: {chi2:.2g}")
print(f"p-value: {p:.3g}")
print(f"degrees of freedom: {dof}")
print("expected frequencies:")
print(expected)


df_CD = figuredata.loc[figuredata["Cohort"] == 'CD']
df_HCNL = figuredata.loc[figuredata["Cohort"] == 'HC-NL']
df_CFS = figuredata.loc[figuredata["Cohort"] == 'ME/CFS']
df_HCUK = figuredata.loc[figuredata["Cohort"] == 'HC-UK']
stat, p = mannwhitneyu(df_CFS.dropna()['Age'], df_HCUK.dropna()['Age'])
print('Statistics=%.3f, p=%.3f' % (stat, p))

#Sex distribution plot
sns.set(rc={'figure.figsize':(12,6)}, font_scale = 1.5)
sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
sns.set_style("white")
custom_palette = ["crimson", "dodgerblue", "purple", "limegreen"]
plot = sns.barplot(data=figuredata, x='Sex', y="Percentage_sex", hue="Cohort", palette=custom_palette)
plot.set_title("Sex distribution", y=1.02, fontsize=24, fontweight='bold')
xticklabels = ["Females", "Males"]
plot.set_xticklabels(labels=xticklabels, rotation=35, ha='center', fontsize=22)
plot.set_ylabel("% of cohort total", fontsize=22)
plot.set_xlabel("")
plt.legend(loc='upper right', fontsize=16)
plt.ylim(0,102)
plt.yticks(fontsize=22)

x1, x2 = -0.3, -.1
y, h, col = figuredata['Percentage_sex'].max() + 5, 2, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*0.5, y+h, "ns", ha='center', va='bottom', color=col, fontsize=20)

x3, x4 = 0.1, 0.3
y, h, col = figuredata['Percentage_sex'].max() + 5, 2, 'k'
plt.plot([x3, x3, x4, x4], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x3+x4)*0.5, y+h, "ns", ha='center', va='bottom', color=col, fontsize=20)

x5, x6 = 0.7, 0.9
y, h, col = figuredata['Percentage_sex'].max() - 35, 2, 'k'
plt.plot([x5, x5, x6, x6], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x5+x6)*0.5, y+h, "ns", ha='center', va='bottom', color=col, fontsize=20)

x7, x8 = 1.1, 1.3
y, h, col = figuredata['Percentage_sex'].max() - 35, 2, 'k'
plt.plot([x7, x7, x8, x8], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x7+x8)*0.5, y+h, "ns", ha='center', va='bottom', color=col, fontsize=20)

x9, x10 = -0.3, 0.1
y, h, col = figuredata['Percentage_sex'].max() + 12, 2, 'k'
plt.plot([x9, x9, x10, x10], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x9+x10)*0.5, y+h, "ns", ha='center', va='bottom', color=col, fontsize=20)

x11, x12 = 0.7, 1.1
y, h, col = figuredata['Percentage_sex'].max() - 25, 2, 'k'
plt.plot([x11, x11, x12, x12], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x11+x12)*0.5, y+h, "ns", ha='center', va='bottom', color=col, fontsize=20)

x13, x14 = 0.9, 1.3
y, h, col = figuredata['Percentage_sex'].max() - 15, 2, 'k'
plt.plot([x13, x13, x14, x14], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x13+x14)*0.5, y+h, "ns", ha='center', va='bottom', color=col, fontsize=20)

x15, x16 = -0.1, 0.3
y, h, col = figuredata['Percentage_sex'].max() + 19, 2, 'k'
plt.plot([x15, x15, x16, x16], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x15+x16)*0.5, y+h, "ns", ha='center', va='bottom', color=col, fontsize=20)

plot.figure.savefig("Flagellin_paper/Sexdistributions_panel_part.jpg", dpi=600, bbox_inches='tight')


#Import combined dataset (file containing overlapping flagellins between CD and ME/CFS cohorts)
combined = pd.read_excel("Flagellin_paper/CDvsCFScohortcomparisons_NEW.xlsx")

#Import combined dataset (containing all overlapping peptides, so also non-flagellins, between CD and ME/CFS cohorts)
all_peptides = pd.read_excel("Flagellin_paper/CDvsMECFS_all_aug21.xlsx")
colors_list = ['crimson', 'purple']

#CD vs ME/CFS plot only flagellins
sns.set(rc={'figure.figsize':(7,7)}, font_scale = 1.5)
sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
sns.set_style("white")
Plot = sns.regplot(x="%ME/CFS", y="%CD", data=combined, color="#953553")
plt.xlabel("% of ME/CFS patients in whom\n a peptide is significantly bound", size=18, labelpad=10)
plt.ylabel("% of CD patients in whom\n a peptide is significantly bound", size=18, labelpad=10)
plt.title("CD vs. ME/CFS", size=24, y=1.03, fontweight='bold')
plt.xlim(0,100)
plt.ylim(0,100)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.text(x=40, y=80, fontsize=16, s="ρ=0.86\n$\it{P}$<0.001", bbox=dict(facecolor='red', alpha=0))
Plot.figure.savefig("Flagellin_paper/CDvsMECFS_regplot.jpg", dpi=600, bbox_inches='tight')

#CD vs. ME/CFS plot all peptides
sns.set(rc={'figure.figsize':(7,7)}, font_scale = 1.5)
sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
sns.set_style("white")
plt.scatter(x="%ME/CFS", y="%CD", label="Flagellin", s=50, c='grey', data=all_peptides[all_peptides['Flagellin'] == 'No'], zorder=1)
plt.scatter(x="%ME/CFS", y="%CD", label="Flagellin", s=50, c='#953553', data=all_peptides[all_peptides['Flagellin'] == 'Yes'], zorder=2)
plt.xlabel("% of ME/CFS patients in whom\n a peptide is significantly bound", size=18, labelpad=10)
plt.ylabel("% of CD patients in whom\n a peptide is significantly bound", size=18, labelpad=10)
plt.title("CD vs. ME/CFS", size=24, y=1.03, fontweight='bold')
plt.xlim(0,100)
plt.ylim(0,100)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.savefig("Flagellin_paper/CDvsMECFS_regplot_allpeptides.jpg", dpi=600, bbox_inches='tight')

#HCNL vs. HCUK plot all peptides
sns.set(rc={'figure.figsize':(7,7)}, font_scale = 1.5)
sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
sns.set_style("white")
plt.scatter(x="%HC-UK", y="%HC-NL", label="Flagellin", s=50, c='grey', data=all_peptides[all_peptides['Flagellin'] == 'No'], zorder=1)
plt.scatter(x="%HC-UK", y="%HC-NL", label="Flagellin", s=50, c='seagreen', data=all_peptides[all_peptides['Flagellin'] == 'Yes'], zorder=2)
plt.xlabel("% of HC-UK individuals in whom\n a peptide is significantly bound", size=18, labelpad=10)
plt.ylabel("% of HC-NL individuals in whom\n a peptide is significantly bound", size=18, labelpad=10)
plt.title("HC-NL vs. HC-UK", size=24, y=1.03, fontweight='bold')
plt.xlim(0,100)
plt.ylim(0,100)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.savefig("Flagellin_paper/HCNLvsUK_regplot_allpeptides.jpg", dpi=600, bbox_inches='tight')

#HCNL vs. HCUK plot only flagellins
colors_list = ['crimson', 'purple']
sns.set(rc={'figure.figsize':(7,7)}, font_scale = 1.5)
sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
sns.set_style("white")
Plot = sns.regplot(x="%HC-UK", y="%HC-NL", data=combined, color="seagreen")
plt.xlabel("% of HC-UK individuals in whom\n a peptide is significantly bound", size=18, labelpad=10)
plt.ylabel("% of HC-NL individuals in whom\n a peptide is significantly boundv", size=18, labelpad=10)
plt.title("HC-NL vs. HC-UK", size=24, y=1.03, fontweight='bold')
plt.xlim(0,100)
plt.ylim(0,100)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.text(x=25, y=60, fontsize=16, s="ρ=0.77\n$\it{P}$<0.001", bbox=dict(facecolor='red', alpha=0))
Plot.figure.savefig("Flagellin_paper/HCUKvsNL_regplot.jpg", dpi=600, bbox_inches='tight')

#Compute correlation coefficients
from scipy.stats import spearmanr
print("correlation", spearmanr(combined.dropna()["%CD"], combined.dropna()["%ME/CFS"]))
print("correlation", spearmanr(combined.dropna()["%HC-NL"], combined.dropna()["%HC-UK"]))

#Annotate flagellin peptides
combined = combined.rename({'Antibody_peptide': 'Unnamed: 0'}, axis=1)
combined_info = combined.merge(information, on=["Unnamed: 0"])
del combined_info['is_IEDB']
del combined_info['is_pos_cntrl']
del combined_info['is_neg_cntrl']
del combined_info['is_phage']
del combined_info['is_influenza']
del combined_info['is_allergens']
del combined_info['is_genome_editing']
del combined_info['IEDB_DOIDs']
del combined_info['IEDB_comments']
del combined_info['IEDB_organism_name']
del combined_info['IEDB_parent_species']
del combined_info['is_rand_cntrl']
del combined_info['is_VFDB']
del combined_info['is_patho_strain']
del combined_info['is_IgA_coated_strain']
del combined_info['is_probio_strain']
del combined_info['bac_src']
del combined_info['is_gut_microbiome']
combined_info.to_excel("Flagellin_paper/CDvsCFScohortcomparisons_NEW.xlsx")

#Plot relative starting positions vs. flagellin peptide frequencies
combined = pd.read_excel("Flagellin_paper/CDvsCFScohortcomparisons_NEW.xlsx")

#CD vs. starting positions
colors_list = ['crimson', 'purple']
sns.set(rc={'figure.figsize':(7,7)}, font_scale = 1.5)
sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
sns.set_style("white")
Plot = sns.scatterplot(x="relative_pos", y="%CD", data=combined, s=80, color='crimson')
plt.xlabel("Relative starting position to protein length", size=18, labelpad=10)
plt.ylabel("% CD patients in whom\n a peptide is significantly bound", size=18, labelpad=10)
plt.title("CD - domain structure", size=24, y=1.03, fontweight='bold')
plt.xlim(-0.02,1.02)
plt.ylim(0,80)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
Plot.figure.savefig("Flagellin_paper/CD_vs_startingposition.jpg", dpi=600, bbox_inches='tight')

#ME/CFS vs. starting positions
colors_list = ['crimson', 'purple']
sns.set(rc={'figure.figsize':(7,7)}, font_scale = 1.5)
sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
sns.set_style("white")
Plot = sns.scatterplot(x="relative_pos", y="%ME/CFS", data=combined, s=80, color='purple')
plt.xlabel("Relative starting position to protein length", size=18, labelpad=10)
plt.ylabel("% ME/CFS patients in whom\n a peptide is significantly bound", size=18, labelpad=10)
plt.title("ME/CFS - domain structure", size=24, y=1.03, fontweight='bold')
plt.xlim(-0.02,1.02)
plt.ylim(0,80)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
Plot.figure.savefig("Flagellin_paper/MECFS_vs_startingposition.jpg", dpi=600, bbox_inches='tight')

#Ratio CD vs. starting positions
colors_list = ['crimson', 'purple']
sns.set(rc={'figure.figsize':(7,7)}, font_scale = 1.5)
sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
sns.set_style("white")
Plot = sns.regplot(x="relative_pos", y="Ratio_CD", data=combined, scatter_kws={'s':80}, color='crimson')
plt.xlabel("Relative starting position to protein length", size=18, labelpad=10)
plt.ylabel("Ratio CD / HC-NL", size=18, labelpad=10)
plt.title("Ratio CD / HC-NL", size=24, y=1.03, fontweight='bold')
plt.xlim(-0.02, 1.02)
plt.ylim(0,40)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
Plot.figure.savefig("Flagellin_paper/RatioCD_vs_startingposition.jpg", dpi=600, bbox_inches='tight')

#Ratio ME/CFS vs. starting positions
colors_list = ['crimson', 'purple']
sns.set(rc={'figure.figsize':(7,7)}, font_scale = 1.5)
sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
sns.set_style("white")
Plot = sns.regplot(x="relative_pos", y="Ratio_CFS", data=combined, scatter_kws={'s':80}, color='purple')
plt.xlabel("Relative starting position to protein length", size=18, labelpad=10)
plt.ylabel("Ratio ME/CFS / HC-UK", size=18, labelpad=10)
plt.title("Ratio ME/CFS / HC-UK", size=24, y=1.03, fontweight='bold')
plt.xlim(-0.02, 1.02)
plt.ylim(0,40)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
Plot.figure.savefig("Flagellin_paper/RatioCFS_vs_startingposition.jpg", dpi=600, bbox_inches='tight')


## Sequence alignment scripts ## - NVH

#FASTA sequence extraction with InterPro tool
with open('sequences_full_length_flagellins.fa', 'rt') as fasta:
    headers = []
    peptide = []
    seq = []
    for l in fasta.readlines():
        l = l.rstrip()
        if '>' in l:
            species = l.split(', ')[0]
            peptide_num = l.split(', ')[1]
            headers.append(species[1:])
            peptide.append(peptide_num)
        else:
            seq.append(l)
    #seq.append(l)

df = pd.DataFrame({'name id' : headers,
                   'peptide_number' : peptide,
                   'sequence' : seq})
print(df)
files = ['results_domain_range_extraction_interpro.xlsx']
for file in files:
    df_blast = pd.read_excel(file)
    df_blast['peptide_number'] = df_blast['peptide_number'].astype(str)
    print(df_blast)
    merged = pd.merge(df, df_blast, on='peptide_number')
    print(merged)

    with open('domain_sequences_interpro.txt', "w") as text_file:
        for x,y,z,dis,num,org in merged[['sequence', 'start', 'end', 'description', 'peptide_number','species']].values:
            #print(x,y,z,dis,num,org)
            y = y - 1
            z = z - 1
            seq = x[y:z]
            org = '>' + org
            range = '({0},{1})'.format(str(y), str(z))
            start = y
            end = z
            length = z - y
            header = ', '.join([org, num, dis])
            #print(header)
            #print(seq)
            print(org, num, dis, start, end, length, seq, sep='\t', file=text_file)

#FASTA sequence extraction - domain alignments
with open('sequences_full_length_flagellins.fa', 'rt') as fasta:
    headers = []
    peptide = []
    seq = []
    for l in fasta.readlines():
        l = l.rstrip()
        if '>' in l:
            species = l.split(', ')[0]
            peptide_num = l.split(', ')[1]
            headers.append(species[1:])
            peptide.append(peptide_num)
        else:
            seq.append(l)
    #seq.append(l)

df = pd.DataFrame({'name id' : headers,
                   'peptide number' : peptide,
                   'sequence' : seq})
#print(df)

df_blast = pd.read_csv('results_domain_range_extraction_domain_alignments.csv', sep='\t', )

merged = pd.merge(df, df_blast, left_on='peptide number', right_on='num_pep')
#print(merged)

with open('domain_sequences_domain_alignments.txt', "w") as text_file:
    for x,y,z,num,org,job in merged[['sequence', 'start', 'end', 'num_pep', 'name id', 'Job title']].values:
        #print(x,y,z,num,org,job)
        y = y - 1
        z = z - 1
        seq = x[y:z]
        org = '>' + org
        range = '(' + str(y) + ',' + str(z) + ')'
        start = y
        end = z
        header = ', '.join([org, job, num, range])
        #print(header)
        #print(seq)
        print(job, org, num, range, sep='\t', file=text_file)

#Domain range extraction with InterPro
import pandas as pd
import json

# Define a function to extract domain start and end values for a given accession
def extract_domain_info(matches, accession, description, xref_name):
    domain_info = []
    for match in matches:
        if 'signature' in match and 'accession' in match['signature'] and match['signature']['accession'] == accession:
            for location in match['locations']:
                domain_info.append({
                    'start': location['start'],
                    'end': location['end'],
                    'description': description,
                    'xref_name': xref_name
                })
    return domain_info

# Prompt the user for the JSON file name
file_name = input("Enter the JSON file name: ")

try:
    # Read the JSON data from the specified file
    with open(file_name, 'r') as file:
        data = json.load(file)

    # Accessions to extract with descriptions
    accessions_to_extract = {
        'PF00669': 'N-terminus',
        'PF00700': 'C-terminus'
    }

    # Create a list to store domain information
    domain_info_list = []

    # Iterate through all xref entries and extract domain information
    for entry in data['results']:
        matches = entry['matches']
        for xref_entry in entry['xref']:
            xref_name = xref_entry['name']
            for accession, description in accessions_to_extract.items():
                domain_info_list.extend(extract_domain_info(matches, accession, description, xref_name))

    # Create a DataFrame from the extracted domain information
    df = pd.DataFrame(domain_info_list)

    # Sort by the xref name
    df = df.sort_values(by=['xref_name']).reset_index(drop=True)

	#Organize DataFrame
    df = pd.concat([df[['start', 'end', 'description']],df['xref_name'].str.split(', ', expand=True)], axis=1)
    df.columns = ['start', 'end', 'description', 'species', 'peptide_number']

except FileNotFoundError:
    print(f"File not found: {file_name}")
except json.JSONDecodeError:
    print("Invalid JSON format in the file.")
except Exception as e:
    print(f"An error occurred: {str(e)}")


df.to_csv("results_file.csv", sep='\t') 

#Domain range extraction - domain alignments
import os
import sys

from collections import defaultdict

d = defaultdict(list)
job_title_list = []
list_of_dicts = []
c=0
with open('source_file.txt', 'rt') as file:
    for line in file.readlines():
        line = line.rstrip()
        if 'Job Title' in line:
            job_title = line.split(':')[1]
            job_title_list.append(job_title)
            if c==0:
                d = defaultdict(list)
                c+=1
            else:
                list_of_dicts.append(d)
                d = defaultdict(list)
            #d['job_title'].append(job_title)
        if '>' in line:
            species = line.split('>')[1].split(',')[0]
            num_pep = line.split('>')[1].split(', ')[1]
            d['species'].append(species)
            #print(species)
            d['num_pep'].append(num_pep)
        if 'Range' in line:
            start = line.split(': ')[1].split(' to ')[0]
            end = line.split(': ')[1].split(' to ')[1]
            rang = (int(start), int(end))
            #print(rang)
            d['range'].append(rang)
    list_of_dicts.append(d)


import pandas as pd

list_df = []
for i,c in zip(job_title_list, list_of_dicts):
    df = pd.DataFrame({'Job title' : [i for x in c['species']], 'species' : [x for x in c['species']],
               'num_pep' : [x for x in c['num_pep']], 'start' : [x[0] for x in c['range']],
                        'end': [x[1] for x in c['range']]})
    list_df.append(df)

merged = pd.concat(list_df)
print(merged)

merged.to_csv('result_file.tsv', sep='\t')

#Figure 4
#Full-length flagellins (approach A)
full_length = pd.read_excel("Flagellin_paper/Alignment full-length Flagellins.xlsx")

##CFS and IBD on separate tabs
sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
sns.set(rc={'figure.figsize':(5,10)})
sns.set_style("white")
colors_list = ['#870303', '#29BF06', '#AA7A38']
plot = sns.boxplot(x="Category", y="Identity (%)", data=full_length, fliersize=0, palette=colors_list, medianprops={'color':'black'}, boxprops={'edgecolor':'black'}, whiskerprops={'color': 'black'}, capprops={'color':'black'})
plot = sns.stripplot(x="Category", y="Identity (%)", data=full_length, jitter=0.1, color='black')
plt.xlabel('')
plt.ylabel('Identity (%)', fontsize=22)
xticklabels = ["Stimulator", "Silent", "Evader"]
plot.set_xticklabels(labels=xticklabels, rotation=35, ha='center', fontsize=20)
plt.ylim(18,115)
plt.yticks(fontsize=20)

x1, x2 = 0, 1
y, h, col = full_length['Identity (%)'].max() + 7, 0.5, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*0.5, y+h, "ns", ha='center', va='bottom', color=col, fontsize=22)

x1, x2 = 1, 2
y, h, col = full_length['Identity (%)'].max() + 3, 0.5, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*0.5, y+h, "***", ha='center', va='bottom', color=col, fontsize=22)

x1, x2 = 0, 2
y, h, col = full_length['Identity (%)'].max() + 11, 0.5, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*0.5, y+h, "***", ha='center', va='bottom', color=col, fontsize=22)

plot.figure.savefig("Flagellin_paper/IBD_full_length.jpg", dpi=600, bbox_inches='tight')

#Reference motifs alignment (approach B)

motifs = pd.read_excel("Flagellin_paper/Alignment full-length Flagellins to Motifs from Clasen.xlsx")

#IBD and CFS on separate tabs
sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
sns.set(rc={'figure.figsize':(5,10)})
sns.set_style("white")
colors_list = ['#870303', '#29BF06']
plot = sns.boxplot(x="Motif", y="Identity (%)", data=motifs, fliersize=0, palette=colors_list, medianprops={'color':'black'}, boxprops={'edgecolor':'black'}, whiskerprops={'color': 'black'}, capprops={'color':'black'})
plot = sns.stripplot(x="Motif", y="Identity (%)", data=motifs, jitter=0.1, color='black')
plt.xlabel('')
plt.ylabel('Identity (%)', fontsize=22)
xticklabels = ["Stimulator", "Silent"]
plot.set_xticklabels(labels=xticklabels, rotation=35, ha='center', fontsize=20)
plt.ylim(18,115)
plt.yticks(fontsize=20)

x1, x2 = 0, 1
y, h, col = motifs['Identity (%)'].max() + 7, 0.5, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*0.5, y+h, "ns", ha='center', va='bottom', color=col, fontsize=22)

plot.figure.savefig("Flagellin_paper/IBD_referencemotifs.jpg", dpi=600, bbox_inches='tight')

#Domain alignments (approach C)
domains = pd.read_excel("Flagellin_paper/Alignment_flagellin_domains_with_reference_flagellin_domains.xlsx")

##CFS and IBD on separate tabs
sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
sns.set(rc={'figure.figsize':(5,10)})
sns.set_style("white")
colors_list = ['#870303', '#29BF06', '#AA7A38']
plot = sns.boxplot(x="Category", y="Identity (%)", data=domains, fliersize=0, palette=colors_list, medianprops={'color':'black'}, boxprops={'edgecolor':'black'}, whiskerprops={'color': 'black'}, capprops={'color':'black'})
plot = sns.stripplot(x="Category", y="Identity (%)", data=domains, jitter=0.1, color='black')
plt.xlabel('')
plt.ylabel('Identity (%)', fontsize=22)
xticklabels = ["Stimulator", "Silent", "Evader"]
plot.set_xticklabels(labels=xticklabels, rotation=35, ha='center', fontsize=20)
plt.ylim(18,115)
plt.title('N-terminal domain (ME/CFS)', fontsize=24, fontweight='bold', y=1.02)
plt.yticks(fontsize=20)

x1, x2 = 0, 1
y, h, col = domains['Identity (%)'].max() + 7, 0.5, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*0.5, y+h, "ns", ha='center', va='bottom', color=col, fontsize=22)

x1, x2 = 1, 2
y, h, col = domains['Identity (%)'].max() + 3, 0.5, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*0.5, y+h, "***", ha='center', va='bottom', color=col, fontsize=22)

x1, x2 = 0, 2
y, h, col = domains['Identity (%)'].max() + 11, 0.5, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*0.5, y+h, "***", ha='center', va='bottom', color=col, fontsize=22)

plot.figure.savefig("Flagellin_paper/CFS_Nterminus_domains.jpg", dpi=600, bbox_inches='tight')

#Histograms sequence lengths
sequences = pd.read_excel("Flagellin_paper/Sequencelengths.xlsx")

sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
sns.set(rc={'figure.figsize':(8,6)})
sns.set_style("white")
color_list = ["#E7872B", "#825CA6"]
plot = sns.displot(sequences, x="Length", kind="hist", bins=40, hue='Category', palette=color_list, legend=False, fill=True)
plt.legend(title='', labels=['C-terminal length', 'N-terminal length'])
plt.xlabel('Sequence length', fontsize=22)
plt.ylabel('Count', fontsize=22)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.xlim(58,)
plt.legend(["C-terminal length", "N-terminal length"], loc="upper center", bbox_to_anchor=(0.55,1))
plt.savefig("Flagellin_paper/Length_distributions.jpg", dpi=600, bbox_inches='tight')


## Supplementary figure: alignment of antibody-bound flagellin peptides
only_peptides = pd.read_excel("Flagellin_paper/Alignment antibody-bound flagellins peptides.xlsx")

##CFS and IBD on separate tabs
## First identity (%)
sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
sns.set(rc={'figure.figsize':(5,10)})
sns.set_style("white")
colors_list = ['#870303', '#29BF06', '#AA7A38']
plot = sns.boxplot(x="Flagellin type", y="Identity (%)", data=only_peptides, fliersize=0, palette=colors_list, medianprops={'color':'black'}, boxprops={'edgecolor':'black'}, whiskerprops={'color': 'black'}, capprops={'color':'black'})
plot = sns.stripplot(x="Flagellin type", y="Identity (%)", data=only_peptides, jitter=0.1, color='black')
plt.xlabel('')
plt.ylabel('Identity (%)', fontsize=22)
xticklabels = ["Stimulator", "Silent", "Evader"]
plot.set_xticklabels(labels=xticklabels, rotation=35, ha='center', fontsize=20)
plt.ylim(18,118)
plt.yticks(fontsize=20)

x1, x2 = 0, 1
y, h, col = only_peptides['Identity (%)'].max() + 7, 0.5, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*0.5, y+h, "ns", ha='center', va='bottom', color=col, fontsize=22)

x1, x2 = 1, 2
y, h, col = only_peptides['Identity (%)'].max() + 3, 0.5, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*0.5, y+h, "***", ha='center', va='bottom', color=col, fontsize=22)

x1, x2 = 0, 2
y, h, col = only_peptides['Identity (%)'].max() + 11, 0.5, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*0.5, y+h, "***", ha='center', va='bottom', color=col, fontsize=22)

plot.figure.savefig("Flagellin_paper/CFS_onlypeptides_alignment.jpg", dpi=600, bbox_inches='tight')

##CFS e-values
sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
sns.set(rc={'figure.figsize':(5,10)})
sns.set_style("white")
colors_list = ['#870303', '#29BF06', '#AA7A38']
plot = sns.boxplot(x="Flagellin type", y="E-value", data=only_peptides, fliersize=0, palette=colors_list, medianprops={'color':'black'}, boxprops={'edgecolor':'black'}, whiskerprops={'color': 'black'}, capprops={'color':'black'})
plot = sns.stripplot(x="Flagellin type", y="E-value", data=only_peptides, jitter=0.1, color='black')
plt.xlabel('')
plt.ylabel('E-value', fontsize=22)
xticklabels = ["Stimulator", "Silent", "Evader"]
plot.set_xticklabels(labels=xticklabels, rotation=35, ha='center', fontsize=20)
plt.yticks(fontsize=20)
plt.yscale('log')
plt.ylim(10e-45, 10e5)

x1, x2 = 0, 1
y, h, col = only_peptides['E-value'].max() + 0.05, 0.03, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*0.5, y+h, "ns", ha='center', va='bottom', color=col, fontsize=22)

x1, x2 = 1, 2
y, h, col = only_peptides['E-value'].max() + 0.0001, 0.0002, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*0.5, y+h, "ns", ha='center', va='bottom', color=col, fontsize=22)

x1, x2 = 0, 2
y, h, col = only_peptides['E-value'].max() + 80, 50, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*0.5, y+h, "ns", ha='center', va='bottom', color=col, fontsize=22)

plot.figure.savefig("Flagellin_paper/CFS_onlypeptides_alignment_evalues.jpg", dpi=600, bbox_inches='tight')

##IBD evalues
sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
sns.set(rc={'figure.figsize':(5,10)})
sns.set_style("white")
colors_list = ['#870303', '#29BF06', '#AA7A38']
plot = sns.boxplot(x="Flagellin type", y="E-value", data=only_peptides, fliersize=0, palette=colors_list, medianprops={'color':'black'}, boxprops={'edgecolor':'black'}, whiskerprops={'color': 'black'}, capprops={'color':'black'})
plot = sns.stripplot(x="Flagellin type", y="E-value", data=only_peptides, jitter=0.1, color='black')
plt.xlabel('')
plt.ylabel('E-value', fontsize=22)
xticklabels = ["Stimulator", "Silent", "Evader"]
plot.set_xticklabels(labels=xticklabels, rotation=35, ha='center', fontsize=20)
plt.yticks(fontsize=20)
plt.yscale('log')
plt.ylim(10e-45, 10e5)

x1, x2 = 0, 1
y, h, col = only_peptides['E-value'].max() + 30, 20, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*0.5, y+h, "ns", ha='center', va='bottom', color=col, fontsize=22)

x1, x2 = 1, 2
y, h, col = only_peptides['E-value'].max() + 3, 3, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*0.5, y+h, "ns", ha='center', va='bottom', color=col, fontsize=22)

x1, x2 = 0, 2
y, h, col = only_peptides['E-value'].max() + 2800, 2000, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*0.5, y+h, "ns", ha='center', va='bottom', color=col, fontsize=22)

plot.figure.savefig("Flagellin_paper/IBD_onlypeptides_alignment_evalues.jpg", dpi=600, bbox_inches='tight')

## Supplementary figure: highly overrepresented flagellins
abundant = pd.read_excel("Flagellin_paper/Alignment highly overrepresented full-length flagellins.xlsx")

##CFS and IBD on separate tabs
## First identity (%)
sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
sns.set(rc={'figure.figsize':(5,10)})
sns.set_style("white")
colors_list = ['#870303', '#29BF06', '#AA7A38']
plot = sns.boxplot(x="Category", y="Identity (%)", data=abundant, fliersize=0, palette=colors_list, medianprops={'color':'black'}, boxprops={'edgecolor':'black'}, whiskerprops={'color': 'black'}, capprops={'color':'black'})
plot = sns.stripplot(x="Category", y="Identity (%)", data=abundant, jitter=0.1, color='black')
plt.xlabel('')
plt.ylabel('Identity (%)', fontsize=22)
xticklabels = ["Stimulator", "Silent", "Evader"]
plot.set_xticklabels(labels=xticklabels, rotation=35, ha='center', fontsize=20)
plt.ylim(18,100)
plt.yticks(fontsize=20)

x1, x2 = 0, 1
y, h, col = abundant['Identity (%)'].max() + 7, 0.5, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*0.5, y+h, "ns", ha='center', va='bottom', color=col, fontsize=22)

x1, x2 = 1, 2
y, h, col = abundant['Identity (%)'].max() + 3, 0.5, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*0.5, y+h, "***", ha='center', va='bottom', color=col, fontsize=22)

x1, x2 = 0, 2
y, h, col = abundant['Identity (%)'].max() + 11, 0.5, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*0.5, y+h, "***", ha='center', va='bottom', color=col, fontsize=22)

plot.figure.savefig("Flagellin_paper/CFS_highlyoverrepresented_alignment.jpg", dpi=600, bbox_inches='tight')

##CFS evalues
sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
sns.set(rc={'figure.figsize':(5,10)})
sns.set_style("white")
colors_list = ['#870303', '#29BF06', '#AA7A38']
plot = sns.boxplot(x="Category", y="E-value", data=abundant, fliersize=0, palette=colors_list, medianprops={'color':'black'}, boxprops={'edgecolor':'black'}, whiskerprops={'color': 'black'}, capprops={'color':'black'})
plot = sns.stripplot(x="Category", y="E-value", data=abundant, jitter=0.1, color='black')
plt.xlabel('')
plt.ylabel('E-value', fontsize=22)
xticklabels = ["Stimulator", "Silent", "Evader"]
plot.set_xticklabels(labels=xticklabels, rotation=35, ha='center', fontsize=20)
plt.yticks(fontsize=20)
plt.yscale('log')
plt.ylim(10e-201, 10e20)

x1, x2 = 0, 1
y, h, col = abundant['E-value'].max() + 10e4, 30e4, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*0.5, y+h, "ns", ha='center', va='bottom', color=col, fontsize=22)

x1, x2 = 1, 2
y, h, col = abundant['E-value'].max() + 1, 10, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*0.5, y+h, "ns", ha='center', va='bottom', color=col, fontsize=22)

x1, x2 = 0, 2
y, h, col = abundant['E-value'].max() + 10e12, 50e12, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*0.5, y+h, "ns", ha='center', va='bottom', color=col, fontsize=22)

plot.figure.savefig("Flagellin_paper/CFS_highlyoverrepresented_alignment_evalues.jpg", dpi=600, bbox_inches='tight')

#IBD evalues
sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
sns.set(rc={'figure.figsize':(5,10)})
sns.set_style("white")
colors_list = ['#870303', '#29BF06', '#AA7A38']
plot = sns.boxplot(x="Category", y="E-value", data=abundant, fliersize=0, palette=colors_list, medianprops={'color':'black'}, boxprops={'edgecolor':'black'}, whiskerprops={'color': 'black'}, capprops={'color':'black'})
plot = sns.stripplot(x="Category", y="E-value", data=abundant, jitter=0.1, color='black')
plt.xlabel('')
plt.ylabel('E-value', fontsize=22)
xticklabels = ["Stimulator", "Silent", "Evader"]
plot.set_xticklabels(labels=xticklabels, rotation=35, ha='center', fontsize=20)
plt.yticks(fontsize=20)
plt.yscale('log')
plt.ylim(10e-201, 10e20)

x1, x2 = 0, 1
y, h, col = abundant['E-value'].max() + 10e4, 30e4, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*0.5, y+h, "ns", ha='center', va='bottom', color=col, fontsize=22)

x1, x2 = 1, 2
y, h, col = abundant['E-value'].max() + 1, 10, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*0.5, y+h, "ns", ha='center', va='bottom', color=col, fontsize=22)

x1, x2 = 0, 2
y, h, col = abundant['E-value'].max() + 10e12, 50e12, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*0.5, y+h, "ns", ha='center', va='bottom', color=col, fontsize=22)

plot.figure.savefig("Flagellin_paper/IBD_highlyoverrepresented_alignment_evalues.jpg", dpi=600, bbox_inches='tight')

full_length = pd.read_excel("Flagellin_paper/Alignment full-length Flagellins.xlsx")

#IBD evalues full-length flagellins
sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
sns.set(rc={'figure.figsize':(5,10)})
sns.set_style("white")
colors_list = ['#870303', '#29BF06', '#AA7A38']
plot = sns.boxplot(x="Category", y="E value", data=full_length, fliersize=0, palette=colors_list, medianprops={'color':'black'}, boxprops={'edgecolor':'black'}, whiskerprops={'color': 'black'}, capprops={'color':'black'})
plot = sns.stripplot(x="Category", y="E value", data=full_length, jitter=0.1, color='black')
plt.xlabel('')
plt.ylabel('E-value', fontsize=22)
xticklabels = ["Stimulator", "Silent", "Evader"]
plot.set_xticklabels(labels=xticklabels, rotation=35, ha='center', fontsize=20)
plt.yticks(fontsize=20)
plt.yscale('log')
plt.ylim(10e-201, 10e20)

x1, x2 = 0, 1
y, h, col = full_length['E value'].max() + 10e4, 30e4, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*0.5, y+h, "ns", ha='center', va='bottom', color=col, fontsize=22)

x1, x2 = 1, 2
y, h, col = full_length['E value'].max() + 1, 10, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*0.5, y+h, "ns", ha='center', va='bottom', color=col, fontsize=22)

x1, x2 = 0, 2
y, h, col = full_length['E value'].max() + 10e12, 50e12, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*0.5, y+h, "ns", ha='center', va='bottom', color=col, fontsize=22)

plot.figure.savefig("Flagellin_paper/IBD_fulllength_alignment_evalues.jpg", dpi=600, bbox_inches='tight')

#CFS evalues full-length flagellins
sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
sns.set(rc={'figure.figsize':(5,10)})
sns.set_style("white")
colors_list = ['#870303', '#29BF06', '#AA7A38']
plot = sns.boxplot(x="Category", y="E value", data=full_length, fliersize=0, palette=colors_list, medianprops={'color':'black'}, boxprops={'edgecolor':'black'}, whiskerprops={'color': 'black'}, capprops={'color':'black'})
plot = sns.stripplot(x="Category", y="E value", data=full_length, jitter=0.1, color='black')
plt.xlabel('')
plt.ylabel('E-value', fontsize=22)
xticklabels = ["Stimulator", "Silent", "Evader"]
plot.set_xticklabels(labels=xticklabels, rotation=35, ha='center', fontsize=20)
plt.yticks(fontsize=20)
plt.yscale('log')
plt.ylim(10e-201, 10e20)

x1, x2 = 0, 1
y, h, col = full_length['E value'].max() + 10e4, 30e4, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*0.5, y+h, "ns", ha='center', va='bottom', color=col, fontsize=22)

x1, x2 = 1, 2
y, h, col = full_length['E value'].max() + 1, 10, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*0.5, y+h, "ns", ha='center', va='bottom', color=col, fontsize=22)

x1, x2 = 0, 2
y, h, col = full_length['E value'].max() + 10e12, 50e12, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*0.5, y+h, "ns", ha='center', va='bottom', color=col, fontsize=22)

plot.figure.savefig("Flagellin_paper/CFS_fulllength_alignment_evalues.jpg", dpi=600, bbox_inches='tight')

motifs = pd.read_excel("Flagellin_paper/Alignment full-length Flagellins to Motifs from Clasen.xlsx")

#IBD motifs evalues
sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
sns.set(rc={'figure.figsize':(5,10)})
sns.set_style("white")
colors_list = ['#870303', '#29BF06', '#AA7A38']
plot = sns.boxplot(x="Motif", y="E value", data=motifs, fliersize=0, palette=colors_list, medianprops={'color':'black'}, boxprops={'edgecolor':'black'}, whiskerprops={'color': 'black'}, capprops={'color':'black'})
plot = sns.stripplot(x="Motif", y="E value", data=motifs, jitter=0.1, color='black')
plt.xlabel('')
plt.ylabel('E-value', fontsize=22)
xticklabels = ["Stimulator", "Silent"]
plot.set_xticklabels(labels=xticklabels, rotation=35, ha='center', fontsize=20)
plt.yticks(fontsize=20)
plt.yscale('log')
plt.ylim(10e-26, 100)

x1, x2 = 0, 1
y, h, col = motifs['E value'].max() + 0.001, 0.0005, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*0.5, y+h, "ns", ha='center', va='bottom', color=col, fontsize=22)

plot.figure.savefig("Flagellin_paper/IBD_referencemotifs_alignment_evalues.jpg", dpi=600, bbox_inches='tight')

#CFS motifs evalues
sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
sns.set(rc={'figure.figsize':(5,10)})
sns.set_style("white")
colors_list = ['#870303', '#29BF06', '#AA7A38']
plot = sns.boxplot(x="Motif", y="E value", data=motifs, fliersize=0, palette=colors_list, medianprops={'color':'black'}, boxprops={'edgecolor':'black'}, whiskerprops={'color': 'black'}, capprops={'color':'black'})
plot = sns.stripplot(x="Motif", y="E value", data=motifs, jitter=0.1, color='black')
plt.xlabel('')
plt.ylabel('E-value', fontsize=22)
xticklabels = ["Stimulator", "Silent"]
plot.set_xticklabels(labels=xticklabels, rotation=35, ha='center', fontsize=20)
plt.yticks(fontsize=20)
plt.yscale('log')
plt.ylim(10e-26, 100)

x1, x2 = 0, 1
y, h, col = motifs['E value'].max() + 0.001, 0.0005, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*0.5, y+h, "ns", ha='center', va='bottom', color=col, fontsize=22)

plot.figure.savefig("Flagellin_paper/CFS_referencemotifs_alignment_evalues.jpg", dpi=600, bbox_inches='tight')

domains = pd.read_excel("Flagellin_paper/Alignment_flagellin_domains_with_reference_flagellin_domains.xlsx")

#CFS evalues N- and C-terminal domains
sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
sns.set(rc={'figure.figsize':(5,10)})
sns.set_style("white")
colors_list = ['#870303', '#29BF06', '#AA7A38']
plot = sns.boxplot(x="Category", y="E value", data=domains, fliersize=0, palette=colors_list, medianprops={'color':'black'}, boxprops={'edgecolor':'black'}, whiskerprops={'color': 'black'}, capprops={'color':'black'})
plot = sns.stripplot(x="Category", y="E value", data=domains, jitter=0.1, color='black')
plt.xlabel('')
plt.ylabel('E-value', fontsize=22)
xticklabels = ["Stimulator", "Silent", "Evader"]
plot.set_xticklabels(labels=xticklabels, rotation=35, ha='center', fontsize=20)
plt.title('N-terminal domain (CD)', fontsize=24, fontweight='bold', y=1.02)
plt.yticks(fontsize=20)
plt.yscale('log')
plt.ylim(10e-101, 10)

x1, x2 = 0, 1
y, h, col = domains['E value'].max() + 10e-11, 5e-10, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*0.5, y+h, "ns", ha='center', va='bottom', color=col, fontsize=22)

x1, x2 = 1, 2
y, h, col = domains['E value'].max() + 10e-15, 5e-14, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*0.5, y+h, "*", ha='center', va='bottom', color=col, fontsize=22)

x1, x2 = 0, 2
y, h, col = domains['E value'].max() + 10e-7, 5e-6, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*0.5, y+h, "*", ha='center', va='bottom', color=col, fontsize=22)

plot.figure.savefig("Flagellin_paper/IBD_Nterm_domain_alignment_evalues.jpg", dpi=600, bbox_inches='tight')

## Additional suppl figure based on alternative domain sequence alignments (InterPro)
domains_extra = pd.read_excel("Flagellin_paper/Blast_results_InterPro_Domains.xlsx")

#Histograms sequence lengths
sequences_extra = pd.read_excel("Flagellin_paper/Blast_results_InterPro_Domains.xlsx")
sequences_extra

sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
sns.set(rc={'figure.figsize':(8,6)})
sns.set_style("white")
color_list = ["#E7872B", "#825CA6"]
plot = sns.displot(sequences_extra, x="Length", kind="hist", bins=100, hue='Category', palette=color_list, legend=False, fill=True)
plt.legend(title='', labels=['C-terminal length', 'N-terminal length'])
plt.xlabel('Sequence length', fontsize=22)
plt.ylabel('Count', fontsize=22)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.xlim(78,)
plt.legend(["C-terminal length", "N-terminal length"], loc="upper center", bbox_to_anchor=(0.55,1))
plt.savefig("Flagellin_paper/Length_distributions_update.jpg", dpi=600, bbox_inches='tight')


#IBD&CFS evalues N- and C-terminal domains
sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
sns.set(rc={'figure.figsize':(5,10)})
sns.set_style("white")
colors_list = ['#870303', '#29BF06', '#AA7A38']
plot = sns.boxplot(x="Category", y="E value", data=domains_extra, fliersize=0, palette=colors_list, medianprops={'color':'black'}, boxprops={'edgecolor':'black'}, whiskerprops={'color': 'black'}, capprops={'color':'black'})
plot = sns.stripplot(x="Category", y="E value", data=domains_extra, jitter=0.1, color='black')
plt.xlabel('')
plt.ylabel('E-value', fontsize=22)
xticklabels = ["Stimulator", "Silent", "Evader"]
plot.set_xticklabels(labels=xticklabels, rotation=35, ha='center', fontsize=20)
plt.title('C-terminal domain (ME/CFS)', fontsize=24, fontweight='bold', y=1.02)
plt.yticks(fontsize=20)
plt.yscale('log')
plt.ylim(10e-65, 10e5)

x1, x2 = 0, 1
y, h, col = domains_extra['E value'].max() + 10e-2, 5e-1, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*0.5, y+h, "ns", ha='center', va='bottom', color=col, fontsize=22)

x1, x2 = 1, 2
y, h, col = domains_extra['E value'].max() + 10e-4, 5e-3, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*0.5, y+h, "ns", ha='center', va='bottom', color=col, fontsize=22)

x1, x2 = 0, 2
y, h, col = domains_extra['E value'].max() + 10e1, 5e2, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*0.5, y+h, "ns", ha='center', va='bottom', color=col, fontsize=22)

plot.figure.savefig("Flagellin_paper/CFS_Cterm_domain_alignment_evalues_update.jpg", dpi=600, bbox_inches='tight')

## Additional figures based on alternative domain sequence alignments
domains_extra = pd.read_excel("Flagellin_paper/Blast_results_InterPro_Domains.xlsx")

##CFS and IBD on separate tabs
sns.set(rc={"figure.dpi":600, 'savefig.dpi':600})
sns.set(rc={'figure.figsize':(5,10)})
sns.set_style("white")
colors_list = ['#870303', '#29BF06', '#AA7A38']
plot = sns.boxplot(x="Category", y="Identity (%)", data=domains_extra, fliersize=0, palette=colors_list, medianprops={'color':'black'}, boxprops={'edgecolor':'black'}, whiskerprops={'color': 'black'}, capprops={'color':'black'})
plot = sns.stripplot(x="Category", y="Identity (%)", data=domains_extra, jitter=0.1, color='black')
plt.xlabel('')
plt.ylabel('Identity (%)', fontsize=22)
xticklabels = ["Stimulator", "Silent", "Evader"]
plot.set_xticklabels(labels=xticklabels, rotation=35, ha='center', fontsize=20)
plt.ylim(18,115)
plt.title('C-terminal domain (CD)', fontsize=24, fontweight='bold', y=1.02)
plt.yticks(fontsize=20)

x1, x2 = 0, 1
y, h, col = domains_extra['Identity (%)'].max() + 7, 0.5, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*0.5, y+h, "*", ha='center', va='bottom', color=col, fontsize=22)

x1, x2 = 1, 2
y, h, col = domains_extra['Identity (%)'].max() + 3, 0.5, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*0.5, y+h, "***", ha='center', va='bottom', color=col, fontsize=22)

x1, x2 = 0, 2
y, h, col = domains_extra['Identity (%)'].max() + 11, 0.5, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*0.5, y+h, "***", ha='center', va='bottom', color=col, fontsize=22)

plot.figure.savefig("Flagellin_paper/CD_Cterminus_domains_updatenicoOct3.jpg", dpi=600, bbox_inches='tight')
