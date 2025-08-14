import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

def main():
    # Parse arguments
    assoc_results_path = sys.argv[1]
    output_plot = sys.argv[2]
    output_top_hits = sys.argv[3]
    p_threshold = float(sys.argv[4])
    top_n = int(sys.argv[5])

    # Load data
    assoc = pd.read_csv(assoc_results_path, sep='\t')
    
    # Validate columns
    required_cols = ['snp_id', 'chrom', 'pos', 'p_value']
    if not all(col in assoc.columns for col in required_cols):
        raise ValueError(f"Missing required columns: {required_cols}")

    # Clean data
    assoc = assoc.dropna(subset=['chrom', 'pos', 'p_value'])
    assoc['chrom'] = assoc['chrom'].astype(str)
    assoc['-log10p'] = -np.log10(assoc['p_value'])

    # Assign x-axis positions
    chrom_order = sorted(assoc['chrom'].unique(), 
                        key=lambda x: int(x) if x.isdigit() else x)
    chrom_offsets = {chrom: i*1e8 for i, chrom in enumerate(chrom_order)}
    assoc['x_pos'] = assoc.apply(
        lambda row: chrom_offsets[row['chrom']] + row['pos'], axis=1
    )

    # Create plot
    plt.figure(figsize=(14, 6))
    
    # Assign colors by chromosome
    color_map = {chrom: '#1f77b4' if i%2 == 0 else '#ff7f0e' 
                for i, chrom in enumerate(chrom_order)}
    colors = assoc['chrom'].map(color_map)

    plt.scatter(
        x=assoc['x_pos'],
        y=assoc['-log10p'],
        c=colors,
        s=8,
        alpha=0.6
    )

    # Highlight significant hits
    sig_hits = assoc[assoc['p_value'] < p_threshold]
    if not sig_hits.empty:
        plt.scatter(
            x=sig_hits['x_pos'],
            y=sig_hits['-log10p'],
            c='red',
            s=20,
            label=f'p < {p_threshold}'
        )

    # Format plot
    plt.axhline(-np.log10(p_threshold), color='gray', linestyle='--', linewidth=0.8)
    plt.xticks(
        ticks=[chrom_offsets[chrom] + 5e7 for chrom in chrom_order],
        labels=chrom_order,
        rotation=45
    )
    plt.xlabel('Chromosome')
    plt.ylabel('-log10(p-value)')
    plt.title('Manhattan Plot')
    plt.grid(True, alpha=0.2)
    plt.tight_layout()
    plt.savefig(output_plot, dpi=300, bbox_inches='tight')
    plt.close()

    # Save top hits
    top_hits = assoc.sort_values('p_value').head(top_n)
    top_hits.to_csv(output_top_hits, sep='\t', index=False)

if __name__ == "__main__":
    main()