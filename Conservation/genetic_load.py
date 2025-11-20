import pandas as pd
from pathlib import Path
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# --- 1. Genetic Load Calculation Function ---

def calculate_genetic_load(variant_file, output_file="genetic_load_report.tsv"):
    """
    Calculate genetic load metrics from variant data with improved scoring.
    """
    variants = pd.read_csv(variant_file, sep='\t')
    
    print(f"Processing {len(variants)} variants...")
    
    def estimate_impact(consequence, allele_freq, genotype="0/1"):
        """Improved impact scoring with genotype consideration"""
        high_impact = ['stop_gained', 'splice_acceptor', 'splice_donor', 'frameshift']
        moderate_impact = ['missense', 'inframe_indel', 'disruptive_inframe_deletion']
        low_impact = ['synonymous', 'intronic', 'intergenic', 'upstream', 'downstream']
        
        if consequence in high_impact:
            score = 1.0
        elif consequence in moderate_impact:
            score = 0.7
        elif consequence in low_impact:
            score = 0.1
        else:
            score = 0.3
        
        # Better frequency adjustment (rare variants = higher load)
        if allele_freq < 0.001:
            freq_factor = 1.0
        elif allele_freq < 0.01:
            freq_factor = 0.9
        elif allele_freq < 0.05:
            freq_factor = 0.7
        else:
            freq_factor = 0.5
        
        # Genotype adjustment (homozygous = 2x impact)
        if genotype in ["1/1", "1"]:
            genotype_factor = 2.0
        else:
            genotype_factor = 1.0
        
        return score * freq_factor * genotype_factor
    
    # Apply scoring
    variants['IMPACT_SCORE'] = variants.apply(
        lambda x: estimate_impact(
            x['CONSEQUENCE'], 
            x['AF'], 
            x.get('GENOTYPE', '0/1')
        ), axis=1
    )
    
    # Calculate metrics
    total_load = variants['IMPACT_SCORE'].sum()
    mean_load = variants['IMPACT_SCORE'].mean()
    median_load = variants['IMPACT_SCORE'].median()
    std_load = variants['IMPACT_SCORE'].std()
    
    print(f"\n=== Genetic Load Summary ===")
    print(f"Total load: {total_load:.2f}")
    print(f"Mean load per variant: {mean_load:.3f}")
    print(f"Median load per variant: {median_load:.3f}")
    print(f"Std dev: {std_load:.3f}")
    print(f"Number of variants: {len(variants)}")
    
    # Group by sample if available
    sample_load_df = None
    if 'SAMPLE' in variants.columns:
        sample_load_df = variants.groupby('SAMPLE')['IMPACT_SCORE'].agg([
            ('Total_Load', 'sum'),
            ('Mean_Load', 'mean'),
            ('Num_Variants', 'count'),
            ('High_Impact_Variants', lambda x: (variants.loc[x.index, 'IMPACT_SCORE'] > 0.9).sum())
        ]).round(3)
        
        sample_load_df.to_csv(output_file, sep='\t')
        print(f"\n✓ Sample-level report saved to {output_file}")
        print(sample_load_df)
    else:
        variants.to_csv(output_file, sep='\t', index=False)
        print(f"\n✓ Variant-level report saved to {output_file}")
    
    return variants, sample_load_df

# --- 2. Sample Data Creation Function ---

def create_sample_data():
    """
    Create example variant data representing two populations with different genetic loads.
    """
    sample_variants = [
        # Population A - Zoo (Higher inbreeding, more deleterious variants)
        {'CHROM': 'chr1', 'POS': 1000, 'CONSEQUENCE': 'missense', 'AF': 0.01, 'GENOTYPE': '0/1', 'SAMPLE': 'IND001'},
        {'CHROM': 'chr1', 'POS': 2000, 'CONSEQUENCE': 'stop_gained', 'AF': 0.005, 'GENOTYPE': '1/1', 'SAMPLE': 'IND001'},
        {'CHROM': 'chr2', 'POS': 5000, 'CONSEQUENCE': 'stop_gained', 'AF': 0.001, 'GENOTYPE': '0/1', 'SAMPLE': 'IND002'},
        {'CHROM': 'chr2', 'POS': 6000, 'CONSEQUENCE': 'intronic', 'AF': 0.3, 'GENOTYPE': '0/1', 'SAMPLE': 'IND002'},
        
        # Population B - Wild (Lower inbreeding, fewer deleterious variants)
        {'CHROM': 'chr3', 'POS': 7000, 'CONSEQUENCE': 'synonymous', 'AF': 0.45, 'GENOTYPE': '0/1', 'SAMPLE': 'IND003'},
        {'CHROM': 'chr3', 'POS': 8000, 'CONSEQUENCE': 'missense', 'AF': 0.25, 'GENOTYPE': '0/1', 'SAMPLE': 'IND003'},
        {'CHROM': 'chr4', 'POS': 9000, 'CONSEQUENCE': 'intronic', 'AF': 0.1, 'GENOTYPE': '0/1', 'SAMPLE': 'IND004'},
        {'CHROM': 'chr4', 'POS': 10000, 'CONSEQUENCE': 'synonymous', 'AF': 0.02, 'GENOTYPE': '0/1', 'SAMPLE': 'IND004'},
    ]
    
    df = pd.DataFrame(sample_variants)
    Path('data').mkdir(exist_ok=True)
    df.to_csv('data/sample_variants.tsv', sep='\t', index=False)
    print("✓ Sample data created: data/sample_variants.tsv\n")
    return df

# --- 3. Visualization Function ---

def visualize_sample_load(sample_load_df, population_mapping=None):
    """
    Visualize genetic load distribution across populations.
    """
    if sample_load_df is None or 'Total_Load' not in sample_load_df.columns:
        print("Error: Sample load summary not available. Skipping visualization.")
        return
    
    plot_df = sample_load_df.reset_index().rename(columns={'SAMPLE': 'Sample'})
    
    # Add population mapping
    if population_mapping:
        plot_df['Population'] = plot_df['Sample'].map(population_mapping)
    else:
        plot_df['Population'] = 'All Samples'
    
    # Create figure with subplots
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Plot 1: Box plot with swarm overlay
    sns.boxplot(
        x='Population', 
        y='Total_Load', 
        data=plot_df, 
        palette="Set2",
        ax=axes[0]
    )
    sns.swarmplot(
        x='Population', 
        y='Total_Load', 
        data=plot_df, 
        color='black', 
        alpha=0.6,
        size=8,
        ax=axes[0]
    )
    axes[0].set_title('Total Genetic Load by Population', fontsize=12, fontweight='bold')
    axes[0].set_ylabel('Total Weighted Impact Score')
    axes[0].grid(axis='y', linestyle='--', alpha=0.7)
    
    # Plot 2: Mean load per variant
    sns.barplot(
        x='Population',
        y='Mean_Load',
        data=plot_df,
        palette="Set2",
        ax=axes[1]
    )
    axes[1].set_title('Mean Genetic Load per Variant', fontsize=12, fontweight='bold')
    axes[1].set_ylabel('Mean Impact Score')
    axes[1].grid(axis='y', linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig("genetic_load_comparison.png", dpi=300, bbox_inches='tight')
    print("✓ Visualization saved as 'genetic_load_comparison.png'\n")
    plt.show()

# --- 4. Statistical Comparison Function ---

def compare_populations(sample_load_df, population_mapping):
    """
    Perform statistical comparison between populations.
    """
    from scipy import stats
    
    plot_df = sample_load_df.reset_index()
    plot_df['Population'] = plot_df['SAMPLE'].map(population_mapping)
    
    populations = plot_df['Population'].unique()
    
    if len(populations) == 2:
        pop_a = plot_df[plot_df['Population'] == populations[0]]['Total_Load'].values
        pop_b = plot_df[plot_df['Population'] == populations[1]]['Total_Load'].values
        
        t_stat, p_value = stats.ttest_ind(pop_a, pop_b)
        
        print(f"=== Population Comparison ===")
        print(f"{populations[0]}: Mean={pop_a.mean():.2f}, SD={pop_a.std():.2f}")
        print(f"{populations[1]}: Mean={pop_b.mean():.2f}, SD={pop_b.std():.2f}")
        print(f"t-test p-value: {p_value:.4f}")
        if p_value < 0.05:
            print("✓ Significant difference detected (p < 0.05)\n")
        else:
            print("✗ No significant difference (p >= 0.05)\n")

# --- 5. Main Execution Block ---

if __name__ == "__main__":
    variant_file = 'data/sample_variants.tsv'
    
    if not Path(variant_file).exists():
        create_sample_data()
    
    population_map = {
        'IND001': 'Zoo Population',
        'IND002': 'Zoo Population',
        'IND003': 'Wild Population',
        'IND004': 'Wild Population',
    }
    
    # Calculate load
    variants_df, sample_load_df = calculate_genetic_load(variant_file)
    
    # Visualize
    if sample_load_df is not None:
        visualize_sample_load(sample_load_df, population_mapping=population_map)
        compare_populations(sample_load_df, population_map)