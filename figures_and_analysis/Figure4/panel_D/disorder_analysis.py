#!/usr/bin/env python3
"""
Disorder Prediction for Protein Sequence Alignments
Handles 2K+ sequences efficiently with multiple prediction methods
Enhanced with customizable colors and vector output formats
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from Bio import SeqIO
from collections import Counter
import argparse
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from functools import partial
import warnings
import json
warnings.filterwarnings('ignore')

# Disorder propensity scales (normalized 0-1, higher = more disordered)
DISORDER_SCALES = {
    # Based on various disorder predictors and amino acid properties
    'charge_hydropathy': {
        'A': 0.3, 'C': 0.2, 'D': 0.8, 'E': 0.8, 'F': 0.1,
        'G': 0.6, 'H': 0.5, 'I': 0.1, 'K': 0.9, 'L': 0.1,
        'M': 0.2, 'N': 0.7, 'P': 0.7, 'Q': 0.6, 'R': 0.9,
        'S': 0.5, 'T': 0.4, 'V': 0.1, 'W': 0.1, 'Y': 0.2,
        'X': 0.5, '-': 0.0
    },
    # Flexibility/disorder promoting residues
    'flexibility': {
        'A': 0.4, 'C': 0.3, 'D': 0.7, 'E': 0.7, 'F': 0.2,
        'G': 0.8, 'H': 0.5, 'I': 0.2, 'K': 0.8, 'L': 0.2,
        'M': 0.3, 'N': 0.6, 'P': 0.9, 'Q': 0.6, 'R': 0.8,
        'S': 0.6, 'T': 0.5, 'V': 0.2, 'W': 0.2, 'Y': 0.3,
        'X': 0.5, '-': 0.0
    },
    # Combined empirical scale
    'combined': {
        'A': 0.35, 'C': 0.25, 'D': 0.75, 'E': 0.75, 'F': 0.15,
        'G': 0.70, 'H': 0.50, 'I': 0.15, 'K': 0.85, 'L': 0.15,
        'M': 0.25, 'N': 0.65, 'P': 0.80, 'Q': 0.60, 'R': 0.85,
        'S': 0.55, 'T': 0.45, 'V': 0.15, 'W': 0.15, 'Y': 0.25,
        'X': 0.50, '-': 0.0
    }
}

# Color scheme presets
COLOR_SCHEMES = {
    'default': {
        'disorder': '#DC2626',
        'entropy': '#2563EB',
        'gap': '#4B5563',
        'threshold': '#000000',
        'fill': '#DC2626',
        'fill_alpha': 0.25
    },
    'publication': {
        'disorder': '#000000',
        'entropy': '#666666',
        'gap': '#999999',
        'threshold': '#000000',
        'fill': '#CCCCCC',
        'fill_alpha': 0.4
    },
    'colorblind': {
        'disorder': '#E69F00',
        'entropy': '#0072B2',
        'gap': '#009E73',
        'threshold': '#000000',
        'fill': '#E69F00',
        'fill_alpha': 0.3
    },
    'vibrant': {
        'disorder': '#FF006E',
        'entropy': '#3A86FF',
        'gap': '#8338EC',
        'threshold': '#000000',
        'fill': '#FF006E',
        'fill_alpha': 0.25
    },
    'pastel': {
        'disorder': '#F4A261',
        'entropy': '#2A9D8F',
        'gap': '#E76F51',
        'threshold': '#264653',
        'fill': '#F4A261',
        'fill_alpha': 0.3
    }
}

def parse_alignment(fasta_file):
    """Parse FASTA alignment file efficiently"""
    sequences = []
    names = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(str(record.seq).upper())
        names.append(record.id)
    return sequences, names

def predict_disorder_single_sequence(sequence, scale='combined', window=9):
    """Predict disorder for a single sequence using sliding window"""
    scores = [DISORDER_SCALES[scale].get(aa, 0.5) for aa in sequence]
    
    # Smooth with sliding window
    if window > 1:
        smoothed = []
        half_window = window // 2
        for i in range(len(scores)):
            start = max(0, i - half_window)
            end = min(len(scores), i + half_window + 1)
            smoothed.append(np.mean(scores[start:end]))
        return smoothed
    return scores

def calculate_consensus_disorder(sequences, scale='combined', window=9, method='mean'):
    """
    Calculate consensus disorder across alignment
    
    Methods:
    - mean: Average disorder score at each position
    - median: Median disorder score at each position
    - max: Maximum disorder score (most permissive)
    - voting: Fraction of sequences predicted as disordered (>0.5)
    """
    if not sequences:
        return []
    
    seq_length = len(sequences[0])
    num_seqs = len(sequences)
    
    # Calculate disorder for each sequence
    print(f"Calculating disorder for {num_seqs} sequences...")
    
    # Use parallel processing for speed
    with ThreadPoolExecutor() as executor:
        predict_fn = partial(predict_disorder_single_sequence, scale=scale, window=window)
        disorder_scores = list(executor.map(predict_fn, sequences))
    
    # Combine scores at each position
    consensus = []
    for pos in range(seq_length):
        position_scores = []
        for seq_disorder in disorder_scores:
            if pos < len(seq_disorder) and sequences[disorder_scores.index(seq_disorder)][pos] != '-':
                position_scores.append(seq_disorder[pos])
        
        if position_scores:
            if method == 'mean':
                consensus.append(np.mean(position_scores))
            elif method == 'median':
                consensus.append(np.median(position_scores))
            elif method == 'max':
                consensus.append(np.max(position_scores))
            elif method == 'voting':
                # Fraction predicted as disordered
                consensus.append(sum(1 for s in position_scores if s > 0.5) / len(position_scores))
            else:
                consensus.append(np.mean(position_scores))
        else:
            consensus.append(0.0)
    
    return consensus

def calculate_position_entropy(sequences):
    """Calculate Shannon entropy at each position (measure of variability)"""
    if not sequences:
        return []
    
    seq_length = len(sequences[0])
    entropies = []
    
    for pos in range(seq_length):
        aas = [seq[pos] for seq in sequences if pos < len(seq) and seq[pos] != '-']
        if not aas:
            entropies.append(0.0)
            continue
        
        counts = Counter(aas)
        total = len(aas)
        entropy = 0
        for count in counts.values():
            if count > 0:
                p = count / total
                entropy -= p * np.log2(p)
        
        # Normalize by max entropy (log2(20) for 20 amino acids)
        max_entropy = np.log2(min(20, len(set(aas))))
        entropies.append(entropy / max_entropy if max_entropy > 0 else 0)
    
    return entropies

def calculate_gap_fraction(sequences):
    """Calculate fraction of gaps at each position"""
    if not sequences:
        return []
    
    seq_length = len(sequences[0])
    num_seqs = len(sequences)
    gap_fractions = []
    
    for pos in range(seq_length):
        gaps = sum(1 for seq in sequences if pos < len(seq) and seq[pos] == '-')
        gap_fractions.append(gaps / num_seqs)
    
    return gap_fractions

def load_color_scheme(color_scheme_name=None, custom_colors=None):
    """Load color scheme from preset or custom colors"""
    if custom_colors:
        # Parse custom colors from command line
        try:
            colors = json.loads(custom_colors)
            # Fill in any missing colors with defaults
            default = COLOR_SCHEMES['default'].copy()
            default.update(colors)
            return default
        except json.JSONDecodeError:
            print("Warning: Invalid custom colors JSON, using default scheme")
            return COLOR_SCHEMES['default']
    elif color_scheme_name in COLOR_SCHEMES:
        return COLOR_SCHEMES[color_scheme_name]
    else:
        return COLOR_SCHEMES['default']

def create_comprehensive_disorder_plot(sequences, names, output_prefix='disorder_analysis',
                                      scale='combined', window=9, method='mean',
                                      disorder_threshold=0.5, colors=None, 
                                      format='png', dpi=300, font_family='sans-serif',
                                      line_width=2.0):
    """Create comprehensive disorder analysis plots with custom colors"""
    
    # Load color scheme
    if colors is None:
        colors = COLOR_SCHEMES['default']
    
    # Set matplotlib parameters for vector output
    if format in ['pdf', 'svg']:
        mpl.rcParams['pdf.fonttype'] = 42  # TrueType fonts
        mpl.rcParams['ps.fonttype'] = 42
    
    mpl.rcParams['font.family'] = font_family
    
    print("Calculating consensus disorder...")
    consensus_disorder = calculate_consensus_disorder(sequences, scale, window, method)
    
    print("Calculating position entropy...")
    entropy = calculate_position_entropy(sequences)
    
    print("Calculating gap fractions...")
    gap_fraction = calculate_gap_fraction(sequences)
    
    positions = range(len(consensus_disorder))
    
    # Create main figure with multiple panels
    fig = plt.figure(figsize=(18, 10))
    gs = fig.add_gridspec(4, 2, hspace=0.3, wspace=0.3)
    
    # 1. Consensus disorder plot
    ax1 = fig.add_subplot(gs[0, :])
    ax1.plot(positions, consensus_disorder, linewidth=line_width, 
             color=colors['disorder'], label='Disorder Score')
    ax1.axhline(y=disorder_threshold, color=colors['threshold'], 
                linestyle='--', linewidth=line_width*0.75, 
                label=f'Threshold ({disorder_threshold})')
    ax1.fill_between(positions, consensus_disorder, disorder_threshold,
                     where=np.array(consensus_disorder) >= disorder_threshold,
                     alpha=colors['fill_alpha'], color=colors['fill'], 
                     label='Disordered regions')
    ax1.set_ylabel('Disorder Score', fontsize=11)
    ax1.set_title(f'Consensus Disorder Prediction ({len(sequences)} sequences, {method} method)', 
                  fontsize=13, fontweight='bold')
    ax1.set_ylim(0, 1)
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='upper right')
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    
    # 2. Disorder with entropy overlay
    ax2 = fig.add_subplot(gs[1, :])
    ax2_twin = ax2.twinx()
    
    line1 = ax2.plot(positions, consensus_disorder, linewidth=line_width, 
                     color=colors['disorder'], label='Disorder')
    line2 = ax2_twin.plot(positions, entropy, linewidth=line_width, 
                         color=colors['entropy'], alpha=0.7, label='Entropy')
    
    ax2.set_ylabel('Disorder Score', fontsize=11, color=colors['disorder'])
    ax2_twin.set_ylabel('Sequence Variability (Entropy)', fontsize=11, color=colors['entropy'])
    ax2.set_title('Disorder vs. Sequence Variability', fontsize=13, fontweight='bold')
    ax2.set_ylim(0, 1)
    ax2_twin.set_ylim(0, 1)
    ax2.grid(True, alpha=0.3)
    ax2.spines['top'].set_visible(False)
    ax2_twin.spines['top'].set_visible(False)
    
    lines = line1 + line2
    labels = [l.get_label() for l in lines]
    ax2.legend(lines, labels, loc='upper right')
    
    # 3. Gap analysis
    ax3 = fig.add_subplot(gs[2, :])
    ax3.fill_between(positions, gap_fraction, alpha=0.5, color=colors['gap'])
    ax3.plot(positions, gap_fraction, linewidth=line_width*0.75, color=colors['gap'])
    ax3.set_ylabel('Gap Fraction', fontsize=11)
    ax3.set_xlabel('Alignment Position', fontsize=11)
    ax3.set_title('Gap Distribution Across Alignment', fontsize=13, fontweight='bold')
    ax3.set_ylim(0, 1)
    ax3.grid(True, alpha=0.3)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    
    # 4. Disorder distribution histogram
    ax4 = fig.add_subplot(gs[3, 0])
    ax4.hist(consensus_disorder, bins=50, color=colors['disorder'], alpha=0.7, edgecolor='black')
    ax4.axvline(x=disorder_threshold, color=colors['threshold'], linestyle='--', linewidth=2)
    ax4.set_xlabel('Disorder Score', fontsize=11)
    ax4.set_ylabel('Frequency', fontsize=11)
    ax4.set_title('Distribution of Disorder Scores', fontsize=12, fontweight='bold')
    ax4.grid(True, alpha=0.3, axis='y')
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)
    
    # 5. Summary statistics
    ax5 = fig.add_subplot(gs[3, 1])
    ax5.axis('off')
    
    # Calculate statistics
    disordered_positions = sum(1 for score in consensus_disorder if score >= disorder_threshold)
    disorder_percent = (disordered_positions / len(consensus_disorder)) * 100
    
    high_entropy_positions = sum(1 for e in entropy if e >= 0.7)
    high_gap_positions = sum(1 for g in gap_fraction if g >= 0.3)
    
    stats_text = "SUMMARY STATISTICS\n" + "="*40 + "\n\n"
    stats_text += f"Total positions: {len(consensus_disorder)}\n"
    stats_text += f"Sequences analyzed: {len(sequences)}\n\n"
    stats_text += f"Disorder Analysis:\n"
    stats_text += f"  Mean disorder: {np.mean(consensus_disorder):.3f}\n"
    stats_text += f"  Median disorder: {np.median(consensus_disorder):.3f}\n"
    stats_text += f"  Disordered positions: {disordered_positions} ({disorder_percent:.1f}%)\n\n"
    stats_text += f"Sequence Variability:\n"
    stats_text += f"  Mean entropy: {np.mean(entropy):.3f}\n"
    stats_text += f"  High variability positions: {high_entropy_positions}\n\n"
    stats_text += f"Gap Analysis:\n"
    stats_text += f"  Mean gap fraction: {np.mean(gap_fraction):.3f}\n"
    stats_text += f"  High gap positions (>30%): {high_gap_positions}\n"
    
    ax5.text(0.1, 0.9, stats_text, transform=ax5.transAxes,
            fontsize=10, verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    output_file = f'{output_prefix}_comprehensive.{format}'
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight', format=format)
    print(f"Comprehensive plot saved to {output_file}")
    plt.close()
    
    # Create simplified publication-ready plot
    fig, ax = plt.subplots(figsize=(14, 5))
    
    ax.plot(positions, consensus_disorder, linewidth=line_width*1.25, 
            color=colors['disorder'])
    ax.axhline(y=disorder_threshold, color=colors['threshold'], 
               linestyle='--', linewidth=line_width*0.75, alpha=0.7)
    ax.fill_between(positions, consensus_disorder, disorder_threshold,
                    where=np.array(consensus_disorder) >= disorder_threshold,
                    alpha=colors['fill_alpha'], color=colors['fill'])
    
    ax.set_xlabel('Alignment Position', fontsize=13)
    ax.set_ylabel('Consensus Disorder Score', fontsize=13)
    ax.set_title(f'Intrinsic Disorder Prediction (n={len(sequences)} sequences)', 
                fontsize=14, fontweight='bold')
    ax.set_ylim(0, 1)
    ax.grid(True, alpha=0.3, linestyle=':')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Add text annotation for disorder percentage
    ax.text(0.98, 0.95, f'{disorder_percent:.1f}% disordered', 
           transform=ax.transAxes, fontsize=11,
           verticalalignment='top', horizontalalignment='right',
           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    output_file = f'{output_prefix}_simple.{format}'
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight', format=format)
    print(f"Simple plot saved to {output_file}")
    plt.close()
    
    return consensus_disorder, entropy, gap_fraction

def export_disorder_scores(consensus_disorder, entropy, gap_fraction, output_file='disorder_scores.txt'):
    """Export disorder scores to text file"""
    with open(output_file, 'w') as f:
        f.write("Position\tDisorder_Score\tEntropy\tGap_Fraction\n")
        for i, (d, e, g) in enumerate(zip(consensus_disorder, entropy, gap_fraction), 1):
            f.write(f"{i}\t{d:.4f}\t{e:.4f}\t{g:.4f}\n")
    print(f"Disorder scores exported to {output_file}")

def main():
    parser = argparse.ArgumentParser(
        description='Predict intrinsic disorder for protein sequence alignments',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python disorder_analysis.py alignment.fasta
  
  # Use colorblind-friendly palette with SVG output
  python disorder_analysis.py alignment.fasta --color-scheme colorblind --format svg
  
  # Custom colors (JSON format)
  python disorder_analysis.py alignment.fasta --custom-colors '{"disorder":"#FF0000","entropy":"#0000FF"}'
  
  # Publication-ready grayscale PDF
  python disorder_analysis.py alignment.fasta --color-scheme publication --format pdf --dpi 600
  
  # Adjust styling
  python disorder_analysis.py alignment.fasta --line-width 3.0 --font-family serif

Available color schemes: default, publication, colorblind, vibrant, pastel

Custom colors JSON keys:
  - disorder: Main disorder line color
  - entropy: Entropy line color
  - gap: Gap plot color
  - threshold: Threshold line color
  - fill: Fill area color
  - fill_alpha: Fill transparency (0-1)
        """
    )
    
    parser.add_argument('alignment', type=str, help='Input alignment file (FASTA format)')
    parser.add_argument('--scale', type=str, default='combined',
                       choices=['charge_hydropathy', 'flexibility', 'combined'],
                       help='Disorder propensity scale (default: combined)')
    parser.add_argument('--method', type=str, default='mean',
                       choices=['mean', 'median', 'max', 'voting'],
                       help='Consensus method (default: mean)')
    parser.add_argument('--window', type=int, default=9,
                       help='Smoothing window size (default: 9)')
    parser.add_argument('--threshold', type=float, default=0.5,
                       help='Disorder threshold (default: 0.5)')
    parser.add_argument('--output-prefix', type=str, default='disorder_analysis',
                       help='Output file prefix (default: disorder_analysis)')
    parser.add_argument('--export-scores', action='store_true',
                       help='Export disorder scores to text file')
    
    # Color and formatting options
    parser.add_argument('--color-scheme', type=str, default='default',
                       choices=list(COLOR_SCHEMES.keys()),
                       help='Color scheme preset (default: default)')
    parser.add_argument('--custom-colors', type=str, default=None,
                       help='Custom colors as JSON string')
    parser.add_argument('--format', type=str, default='png',
                       choices=['png', 'pdf', 'svg', 'eps'],
                       help='Output format (default: png). Use pdf/svg for editing')
    parser.add_argument('--dpi', type=int, default=300,
                       help='DPI for raster formats (default: 300)')
    parser.add_argument('--font-family', type=str, default='sans-serif',
                       choices=['sans-serif', 'serif', 'monospace'],
                       help='Font family (default: sans-serif)')
    parser.add_argument('--line-width', type=float, default=2.0,
                       help='Line width for plots (default: 2.0)')
    
    args = parser.parse_args()
    
    # Parse alignment
    print(f"Reading alignment from {args.alignment}...")
    sequences, names = parse_alignment(args.alignment)
    print(f"Loaded {len(sequences)} sequences, alignment length: {len(sequences[0])}")
    
    # Load color scheme
    colors = load_color_scheme(args.color_scheme, args.custom_colors)
    print(f"Using color scheme: {args.color_scheme}")
    
    # Create plots
    consensus_disorder, entropy, gap_fraction = create_comprehensive_disorder_plot(
        sequences, names,
        output_prefix=args.output_prefix,
        scale=args.scale,
        window=args.window,
        method=args.method,
        disorder_threshold=args.threshold,
        colors=colors,
        format=args.format,
        dpi=args.dpi,
        font_family=args.font_family,
        line_width=args.line_width
    )
    
    # Export scores if requested
    if args.export_scores:
        export_disorder_scores(consensus_disorder, entropy, gap_fraction,
                              f'{args.output_prefix}_scores.txt')
    
    print("\nAnalysis complete!")
    print(f"Generated files:")
    print(f"  - {args.output_prefix}_comprehensive.{args.format} (detailed multi-panel plot)")
    print(f"  - {args.output_prefix}_simple.{args.format} (publication-ready plot)")
    if args.export_scores:
        print(f"  - {args.output_prefix}_scores.txt (disorder scores)")
    
    if args.format in ['pdf', 'svg']:
        print(f"\nVector format ({args.format}) generated - fully editable in Illustrator/Inkscape!")

if __name__ == "__main__":
    main()
