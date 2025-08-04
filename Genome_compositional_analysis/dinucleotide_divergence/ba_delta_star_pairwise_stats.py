"""
Perform pairwise Wilcoxon rank-sum tests on Karlin δ* values between all
biogenic amine gene types to identify significant differences in
dinucleotide composition deviation.

Input:
- CSV file 'karlin_delta_star_results.csv' with columns including 'Gene_Name' and 'Delta_Star'

Process:
- For each unique pair of gene types, run Wilcoxon rank-sum test on their δ* distributions
- Adjust raw p-values for multiple testing using Benjamini-Hochberg FDR correction
- Save results with test statistics, raw and adjusted p-values, and significance flags

Output:
- CSV file 'ba_delta_star_pairwise_stats.csv' containing pairwise test results

Usage:
- Ensure 'karlin_delta_star_results.csv' is in the working directory
- Run this script to generate pairwise comparison statistics

Requirements:
- pandas
- scipy
- statsmodels
Author:
    Aqib Javaid
"""
import pandas as pd
import itertools
from scipy.stats import ranksums
from statsmodels.stats.multitest import multipletests

# Load the BA gene delta* data
df = pd.read_csv("karlin_delta_star_results.csv")

# Unique gene types
genes = df['Gene_Name'].unique()
pairs = list(itertools.combinations(genes, 2))

# Perform Wilcoxon tests
pvals = []
stats = []
for g1, g2 in pairs:
    group1 = df[df['Gene_Name'] == g1]['Delta_Star']
    group2 = df[df['Gene_Name'] == g2]['Delta_Star']
    stat, p = ranksums(group1, group2)
    stats.append(stat)
    pvals.append(p)

# FDR correction
reject, pvals_corr, _, _ = multipletests(pvals, method='fdr_bh')

# Save results
comparison_df = pd.DataFrame({
    "Gene1": [x[0] for x in pairs],
    "Gene2": [x[1] for x in pairs],
    "Wilcoxon_Stat": stats,
    "p_raw": pvals,
    "p_adj": pvals_corr,
    "Significant": reject
})
comparison_df.to_csv("ba_delta_star_pairwise_stats.csv", index=False)
print("✅ Pairwise comparison saved as ba_delta_star_pairwise_stats.csv")
