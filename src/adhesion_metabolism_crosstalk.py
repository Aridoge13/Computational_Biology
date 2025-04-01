from gseapy import enrichr
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np

# Define gene sets
desmosome_genes = ['DSG1', 'DSP', 'PKP1', 'JUP',
                   'DSC1', 'PKP3']
metabolic_genes = ['PRKAA1', 'SLC2A1', 'HK2',
                   'PFKM', 'PDK1', 'MTOR']

# Run enrichment analyses
print("Running enrichment analyses...")
desmo_enrich = enrichr(gene_list=desmosome_genes,
                       gene_sets=['KEGG_2021_Human',
                                  'GO_Biological_Process_2023'],
                       organism='human')

metab_enrich = enrichr(gene_list=metabolic_genes,
                       gene_sets=['KEGG_2021_Human'],
                       organism='human')

# Process results


def process_enrichment_results(df):
    """Clean and prepare enrichment results"""
    df = df.copy()
    df['-log10(Adj P)'] = -np.log10(df['Adjusted P-value'])
    df['Term'] = df['Term'].str.replace('^GOBP_|^KEGG_', '', regex=True)
    return df


df_desmo = process_enrichment_results(desmo_enrich.results)
df_desmo['source'] = 'Desmosome'

df_metab = process_enrichment_results(metab_enrich.results)
df_metab['source'] = 'Metabolism'

# Find overlapping pathways
common_terms = set(df_desmo['Term']).intersection(set(df_metab['Term']))
print(f"\nShared pathways ({len(common_terms)}):")
for i, term in enumerate(common_terms, 1):
    print(f"{i}. {term}")

# Visualization
# Combine and get top 10 terms from each analysis
top_terms = pd.concat([
    df_desmo.nsmallest(10, 'Adjusted P-value'),
    df_metab.nsmallest(10, 'Adjusted P-value')
])

# Create the plot
plt.figure(figsize=(12, 8))
scatter = sns.scatterplot(
    data=top_terms,
    x='source',
    y='Term',
    size='Odds Ratio',
    hue='-log10(Adj P)',
    palette='viridis',
    sizes=(50, 250),
    alpha=0.8
)

# Customize plot
plt.title('Pathway Enrichment: Desmosome vs. Metabolism Genes', pad=20)
plt.xlabel('Gene Set', labelpad=10)
plt.ylabel('Enriched Pathway', labelpad=10)
plt.xticks(rotation=0)

# Move legend
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

# Highlight shared pathways
for term in common_terms:
    if term in top_terms['Term'].values:
        idx = top_terms[top_terms['Term'] == term].index[0]
        plt.gca().get_yticklabels()[idx].set_color('red')
        plt.gca().get_yticklabels()[idx].set_weight('bold')

plt.tight_layout()

# Save and show
plt.savefig('desmosome_metabolism_crosstalk.png', dpi=300, bbox_inches='tight')
print("\nSaved plot to 'desmosome_metabolism_crosstalk.png'")
plt.show()
