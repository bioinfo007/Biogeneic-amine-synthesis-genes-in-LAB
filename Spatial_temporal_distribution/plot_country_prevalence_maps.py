#!/usr/bin/env python3
"""
plot_country_prevalence_maps.py
------------------------------

Description:
    Computes country-wise prevalence (with 95% Wilson confidence intervals) of selected BA biosynthesis genes
    from a gene presence/absence matrix and plots choropleth maps.

Inputs:
    - gene_presence_absence_matrix.csv
        Required columns: at minimum 'Country' and gene presence/absence columns (adc, hdc, agmatinase, odc, tdc).
        Country names should match standard country names recognizable by Plotly.

Configuration (can be edited in script):
    - genes: list of gene column names to plot.
    - min_samples: minimum number of isolates per country to include in the map.

Outputs:
    - One SVG choropleth map per gene, named like: <gene>_prevalence_map_CI.svg

Usage:
    Place this script in the directory containing `gene_presence_absence_matrix.csv` and run:

        ./plot_country_prevalence_maps.py
    or:
        python3 plot_country_prevalence_maps.py

Requirements:
    - Python 3
    - pandas
    - plotly
    - statsmodels

Install dependencies with:
    pip install pandas plotly statsmodels

Notes:
    - Countries with fewer than `min_samples` isolates are excluded for that gene.
    - Prevalence is shown in percentage, with Wilson 95% confidence intervals in hover tooltip.
"""

import pandas as pd
import plotly.express as px
from statsmodels.stats.proportion import proportion_confint

# Load and clean data
df = pd.read_csv("gene_presence_absence_matrix.csv")
df = df[df['Country'].notna()]

genes = ['adc', 'hdc', 'agmatinase', 'odc', 'tdc']
min_samples = 10  # Minimum isolates per country to include

for gene in genes:
    # Group by country and calculate sums and counts
    grouped = df.groupby('Country')[gene].agg(['sum', 'count']).reset_index()
    
    # Filter countries with sufficient sample size
    grouped = grouped[grouped['count'] >= min_samples].copy()
    
    # Calculate prevalence in percentage
    grouped['Prevalence (%)'] = 100 * grouped['sum'] / grouped['count']
    
    # Calculate 95% confidence intervals (Wilson method)
    ci_low, ci_upp = proportion_confint(
        grouped['sum'], grouped['count'], alpha=0.05, method='wilson'
    )
    grouped['CI Lower (%)'] = ci_low * 100
    grouped['CI Upper (%)'] = ci_upp * 100
    
    grouped = grouped.rename(columns={'Country': 'country'})
    
    # Plot choropleth map
    fig = px.choropleth(
        grouped,
        locations='country',
        locationmode='country names',
        color='Prevalence (%)',
        hover_data={
            'sum': True,
            'count': True,
            'Prevalence (%)': ':.2f',
            'CI Lower (%)': ':.2f',
            'CI Upper (%)': ':.2f'
        },
        color_continuous_scale='Reds',
        title=f"Country-wise Prevalence of '{gene}' Gene in LAB Isolates (≥{min_samples} samples)"
    )
    
    fig.update_geos(showcountries=True, projection_type="natural earth")
    fig.update_layout(margin={"r":0,"t":40,"l":0,"b":0})
    
    # Save figure as SVG
    fig.write_image(f"{gene}_prevalence_map_CI.svg", width=1200, height=700)
    
    print(f"Saved SVG map for gene: {gene}")
