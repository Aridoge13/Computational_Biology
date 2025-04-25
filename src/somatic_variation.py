import pandas as pd


class SimpleVCFAnnotator:
    def __init__(self, cosmic_genes_file=None):
        self.cancer_genes = set()
        if cosmic_genes_file:
            with open(cosmic_genes_file, 'r') as f:
                self.cancer_genes = {line.strip()
                                     for line in f if line.strip()}

    def annotate_vcf(self, vcf_path, output_path=None):
        df = pd.read_csv(vcf_path, sep='\t', comment='#',
                         names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'])

        # Extract gene info from INFO (if available)
        df['GENE'] = df['INFO'].str.extract(r'GENE=([^;]+)')
        if self.cancer_genes:
            df['In_Cancer_Gene'] = df['GENE'].isin(self.cancer_genes)

        if output_path:
            df.to_csv(output_path, sep='\t', index=False)
        return df
