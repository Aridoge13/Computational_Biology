import pandas as pd
import argparse
import subprocess

def call_svs_with_sniffles(bam_file, vcf_output="svs.vcf", threads=4):
    """
    Call structural variants using Sniffles from a BAM file
    """
    print(f"Calling SVs from {bam_file} using Sniffles...")
    
    cmd = [
        "sniffles",
        "--input", bam_file,
        "--vcf", vcf_output,
        "--threads", str(threads),
        "--allow-overwrite"
    ]
    
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(f"Sniffles completed! VCF saved to {vcf_output}")
        return vcf_output
    except FileNotFoundError:
        print(" Sniffles not found. Install with: conda install -c bioconda sniffles")
        raise

def parse_sniffles_vcf(vcf_file):
    """Parse Sniffles VCF output into DataFrame"""
    svs = []
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            info_dict = {item.split('=')[0]: item.split('=')[1] for item in fields[7].split(';') if '=' in item}
            
            svs.append({
                'CHROM': fields[0],
                'POS': int(fields[1]),
                'END': int(info_dict.get('END', fields[1])),
                'SVTYPE': info_dict.get('SVTYPE', 'UNKNOWN'),
                'SV_ID': fields[2]
            })
    return pd.DataFrame(svs)

def annotate_svs(sv_df, gff_file=None, output_file="annotated_svs.tsv"):
    """
    Annotate structural variants with basic genomic context
    """
    print(f"Annotating {len(sv_df)} structural variants...")
    
    annotations = []
    for idx, row in sv_df.iterrows():
        sv_size = abs(row['END'] - row['POS'])
        
        if sv_size > 1000:
            impact = "HIGH"
            consequence = "large_deletion" if row['SVTYPE'] == 'DEL' else "large_insertion"
        else:
            impact = "MODERATE" 
            consequence = "small_SV"
        
        annotations.append({
            'SV_ID': row['SV_ID'],
            'CHROM': row['CHROM'],
            'POS': row['POS'],
            'END': row['END'],
            'TYPE': row['SVTYPE'],
            'SIZE': sv_size,
            'CONSEQUENCE': consequence,
            'IMPACT': impact
        })
    
    result_df = pd.DataFrame(annotations)
    result_df.to_csv(output_file, sep='\t', index=False)
    print(f"✓ Results saved to {output_file}")
    return result_df

def main():
    parser = argparse.ArgumentParser(description='SV calling and annotation with Sniffles')
    parser.add_argument('--bam_file', required=True, help='Input BAM file')
    parser.add_argument('--vcf_output', default='svs.vcf', help='Output VCF file')
    parser.add_argument('--annotated_output', default='annotated_svs.tsv', help='Annotated output')
    parser.add_argument('--threads', default=4, type=int, help='Sniffles threads')
    
    args = parser.parse_args()
    
    # Call SVs
    vcf_file = call_svs_with_sniffles(args.bam_file, args.vcf_output, args.threads)
    
    # Parse and annotate
    sv_df = parse_sniffles_vcf(vcf_file)
    print(f"Found {len(sv_df)} SVs")
    annotate_svs(sv_df, output_file=args.annotated_output)

if __name__ == "__main__":
    main()