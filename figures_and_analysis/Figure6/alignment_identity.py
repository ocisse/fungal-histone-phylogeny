#!/usr/bin/env python3
"""
Calculate average percent identity for a multiple sequence alignment
"""

from Bio import AlignIO
import numpy as np
import argparse
from itertools import combinations

def calculate_pairwise_identity(seq1, seq2, ignore_gaps=True):
    """
    Calculate percent identity between two aligned sequences
    
    Parameters:
    -----------
    seq1, seq2 : str
        Aligned sequences
    ignore_gaps : bool
        If True, ignore positions where either sequence has a gap
    
    Returns:
    --------
    float : Percent identity (0-100)
    """
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be the same length (aligned)")
    
    matches = 0
    valid_positions = 0
    
    for i in range(len(seq1)):
        aa1, aa2 = seq1[i], seq2[i]
        
        if ignore_gaps and (aa1 == '-' or aa2 == '-'):
            continue
        
        valid_positions += 1
        if aa1 == aa2:
            matches += 1
    
    if valid_positions == 0:
        return 0.0
    
    return (matches / valid_positions) * 100

def calculate_alignment_identity(alignment, ignore_gaps=True, sample_size=None):
    """
    Calculate average percent identity for entire alignment
    
    Parameters:
    -----------
    alignment : Bio.Align.MultipleSeqAlignment
    ignore_gaps : bool
        Whether to ignore gap positions
    sample_size : int or None
        If specified, randomly sample this many sequence pairs
        (useful for very large alignments)
    
    Returns:
    --------
    dict : Statistics including mean, median, min, max identity
    """
    sequences = [str(record.seq) for record in alignment]
    n_seqs = len(sequences)
    
    print(f"Calculating pairwise identities for {n_seqs} sequences...")
    
    # Generate all pairwise combinations
    pairs = list(combinations(range(n_seqs), 2))
    n_pairs = len(pairs)
    
    print(f"Total pairwise comparisons: {n_pairs}")
    
    # Sample if requested (for large alignments)
    if sample_size and sample_size < n_pairs:
        import random
        pairs = random.sample(pairs, sample_size)
        print(f"Sampling {sample_size} pairs...")
    
    identities = []
    
    # Calculate pairwise identities
    for i, (idx1, idx2) in enumerate(pairs):
        if i % 1000 == 0:
            print(f"Processed {i}/{len(pairs)} pairs...")
        
        identity = calculate_pairwise_identity(
            sequences[idx1], 
            sequences[idx2], 
            ignore_gaps=ignore_gaps
        )
        identities.append(identity)
    
    # Calculate statistics
    identities = np.array(identities)
    
    stats = {
        'n_sequences': n_seqs,
        'n_comparisons': len(identities),
        'mean_identity': np.mean(identities),
        'median_identity': np.median(identities),
        'std_identity': np.std(identities),
        'min_identity': np.min(identities),
        'max_identity': np.max(identities),
        'q25_identity': np.percentile(identities, 25),
        'q75_identity': np.percentile(identities, 75)
    }
    
    return stats, identities

def print_statistics(stats):
    """Print alignment identity statistics"""
    print("\n" + "="*50)
    print("ALIGNMENT IDENTITY STATISTICS")
    print("="*50)
    print(f"Number of sequences: {stats['n_sequences']}")
    print(f"Pairwise comparisons: {stats['n_comparisons']}")
    print(f"\nPercent Identity:")
    print(f"  Mean:   {stats['mean_identity']:.2f}%")
    print(f"  Median: {stats['median_identity']:.2f}%")
    print(f"  Std:    {stats['std_identity']:.2f}%")
    print(f"  Min:    {stats['min_identity']:.2f}%")
    print(f"  Max:    {stats['max_identity']:.2f}%")
    print(f"  25th percentile: {stats['q25_identity']:.2f}%")
    print(f"  75th percentile: {stats['q75_identity']:.2f}%")

def plot_identity_distribution(identities, output_file=None):
    """Plot histogram of percent identities"""
    try:
        import matplotlib.pyplot as plt
        
        plt.figure(figsize=(10, 6))
        plt.hist(identities, bins=50, alpha=0.7, color='steelblue', edgecolor='black')
        plt.xlabel('Percent Identity (%)', fontsize=12)
        plt.ylabel('Number of Sequence Pairs', fontsize=12)
        plt.title('Distribution of Pairwise Sequence Identity', fontsize=14)
        plt.grid(alpha=0.3)
        
        # Add statistics text
        mean_id = np.mean(identities)
        plt.axvline(mean_id, color='red', linestyle='--', linewidth=2, 
                   label=f'Mean: {mean_id:.1f}%')
        plt.legend()
        
        if output_file:
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            print(f"Plot saved to: {output_file}")
        else:
            plt.show()
        
    except ImportError:
        print("Matplotlib not available. Skipping plot.")

def main():
    parser = argparse.ArgumentParser(description='Calculate average percent identity for alignment')
    parser.add_argument('alignment', help='Input alignment file (FASTA format)')
    parser.add_argument('--ignore-gaps', action='store_true', default=True,
                       help='Ignore gap positions (default: True)')
    parser.add_argument('--include-gaps', action='store_true',
                       help='Include gap positions in calculation')
    parser.add_argument('--sample-size', type=int,
                       help='Sample this many pairs (for large alignments)')
    parser.add_argument('--output', '-o', 
                       help='Output file for statistics')
    parser.add_argument('--plot', 
                       help='Save histogram plot to this file')
    parser.add_argument('--format', default='fasta',
                       help='Alignment format (default: fasta)')
    
    args = parser.parse_args()
    
    # Handle gap parameter
    ignore_gaps = args.ignore_gaps and not args.include_gaps
    
    # Load alignment
    print(f"Loading alignment from: {args.alignment}")
    try:
        alignment = AlignIO.read(args.alignment, args.format)
        print(f"Loaded {len(alignment)} sequences, length {alignment.get_alignment_length()}")
    except Exception as e:
        print(f"Error loading alignment: {e}")
        return
    
    # Calculate identities
    stats, identities = calculate_alignment_identity(
        alignment, 
        ignore_gaps=ignore_gaps,
        sample_size=args.sample_size
    )
    
    # Print results
    print_statistics(stats)
    
    # Save statistics
    if args.output:
        with open(args.output, 'w') as f:
            f.write("Alignment Identity Statistics\n")
            f.write("="*40 + "\n")
            for key, value in stats.items():
                f.write(f"{key}: {value}\n")
        print(f"\nStatistics saved to: {args.output}")
    
    # Plot distribution
    if args.plot:
        plot_identity_distribution(identities, args.plot)

if __name__ == "__main__":
    main()
