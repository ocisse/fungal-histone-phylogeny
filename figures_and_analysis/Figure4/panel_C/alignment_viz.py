#!/usr/bin/env python3
"""
Protein Alignment Divergence Visualizer for Large Datasets
Compatible with CIAlign output and handles 2K+ sequences efficiently
Enhanced with customizable colors and vector output formats
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib as mpl
from matplotlib.colors import ListedColormap, BoundaryNorm
import seaborn as sns
from Bio import SeqIO
from collections import defaultdict
import argparse
from pathlib import Path
import json

# Color scheme presets for alignment visualization
COLOR_SCHEMES = {
    'default': {
        'match': '#10B981',
        'mismatch': '#EF4444',
        'gap': '#E5E7EB',
        'reference_label': '#DC2626',
        'text': '#000000',
        'grid': '#D1D5DB'
    },
    'publication': {
        'match': '#000000',
        'mismatch': '#FFFFFF',
        'gap': '#CCCCCC',
        'reference_label': '#000000',
        'text': '#000000',
        'grid': '#666666'
    },
    'colorblind': {
        'match': '#0072B2',
        'mismatch': '#E69F00',
        'gap': '#F0F0F0',
        'reference_label': '#D55E00',
        'text': '#000000',
        'grid': '#CCCCCC'
    },
    'vibrant': {
        'match': '#06D6A0',
        'mismatch': '#EF476F',
        'gap': '#F8F9FA',
        'reference_label': '#FF006E',
        'text': '#073B4C',
        'grid': '#CED4DA'
    },
    'chemistry': {
        'hydrophobic': '#FCD34D',
        'polar': '#60A5FA',
        'positive': '#F472B6',
        'negative': '#A78BFA',
        'special': '#34D399',
        'gap': '#E5E7EB',
        'reference_label': '#DC2626',
        'text': '#000000',
        'grid': '#D1D5DB'
    }
}

def parse_alignment(fasta_file):
    """Parse FASTA alignment file"""
    sequences = []
    names = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(str(record.seq))
        names.append(record.id)
    return sequences, names

def calculate_divergence_matrix(sequences, ref_idx=0):
    """Calculate divergence of all sequences from reference"""
    ref_seq = sequences[ref_idx]
    divergences = []
    
    for seq in sequences:
        matches = 0
        valid_positions = 0
        
        for i, (ref_aa, aa) in enumerate(zip(ref_seq, seq)):
            if ref_aa != '-' and aa != '-':
                valid_positions += 1
                if ref_aa == aa:
                    matches += 1
        
        div_percent = ((valid_positions - matches) / valid_positions * 100) if valid_positions > 0 else 0
        divergences.append(div_percent)
    
    return divergences

def calculate_position_conservation(sequences, ref_idx=0):
    """Calculate conservation score at each position relative to reference"""
    ref_seq = sequences[ref_idx]
    seq_length = len(ref_seq)
    num_seqs = len(sequences)
    
    conservation_scores = []
    position_divergence = []
    
    for pos in range(seq_length):
        ref_aa = ref_seq[pos]
        if ref_aa == '-':
            conservation_scores.append(np.nan)
            position_divergence.append(np.nan)
            continue
        
        matches = sum(1 for seq in sequences if seq[pos] == ref_aa)
        conservation = (matches / num_seqs) * 100
        conservation_scores.append(conservation)
        position_divergence.append(100 - conservation)
    
    return conservation_scores, position_divergence

def get_aa_color_chemistry(aa, colors):
    """Get color for amino acid based on chemical properties"""
    hydrophobic = 'AILMFWV'
    polar = 'STNQ'
    positive = 'KRH'
    negative = 'DE'
    special = 'CGP'
    
    if aa == '-':
        return colors.get('gap', '#E5E7EB')
    if hydrophobic.find(aa) != -1:
        return colors.get('hydrophobic', '#FCD34D')
    if polar.find(aa) != -1:
        return colors.get('polar', '#60A5FA')
    if positive.find(aa) != -1:
        return colors.get('positive', '#F472B6')
    if negative.find(aa) != -1:
        return colors.get('negative', '#A78BFA')
    if special.find(aa) != -1:
        return colors.get('special', '#34D399')
    return '#9CA3AF'

def load_color_scheme(color_scheme_name=None, custom_colors=None):
    """Load color scheme from preset or custom colors"""
    if custom_colors:
        try:
            colors = json.loads(custom_colors)
            default = COLOR_SCHEMES['default'].copy()
            default.update(colors)
            return colors
        except json.JSONDecodeError:
            print("Warning: Invalid custom colors JSON, using default scheme")
            return COLOR_SCHEMES['default']
    elif color_scheme_name in COLOR_SCHEMES:
        return COLOR_SCHEMES[color_scheme_name]
    else:
        return COLOR_SCHEMES['default']

def create_divergence_heatmap(sequences, names, ref_idx=0, output_file='divergence_heatmap.png', 
                              max_seqs=100, sample_positions=None, colors=None, format='png',
                              dpi=300, show_sequences=True, font_size=8):
    """Create heatmap showing divergence from reference"""
    if colors is None:
        colors = COLOR_SCHEMES['default']
    
    # Set matplotlib parameters for vector output
    if format in ['pdf', 'svg']:
        mpl.rcParams['pdf.fonttype'] = 42
        mpl.rcParams['ps.fonttype'] = 42
    
    ref_seq = sequences[ref_idx]
    
    # Sample sequences if too many
    if len(sequences) > max_seqs:
        print(f"Sampling {max_seqs} sequences from {len(sequences)} total for visualization")
        indices = [ref_idx] + list(np.random.choice(
            [i for i in range(len(sequences)) if i != ref_idx], 
            max_seqs - 1, 
            replace=False
        ))
        indices.sort()
        sequences = [sequences[i] for i in indices]
        names = [names[i] for i in indices]
        ref_idx = 0
    
    # Sample positions if specified
    seq_length = len(ref_seq)
    if sample_positions and sample_positions < seq_length:
        print(f"Sampling {sample_positions} positions from {seq_length} total")
        positions = sorted(np.random.choice(seq_length, sample_positions, replace=False))
    else:
        positions = range(seq_length)
    
    # Create divergence matrix
    matrix = np.zeros((len(sequences), len(positions)))
    for i, seq in enumerate(sequences):
        for j, pos in enumerate(positions):
            if ref_seq[pos] == '-' or seq[pos] == '-':
                matrix[i, j] = -1  # Gap
            elif ref_seq[pos] == seq[pos]:
                matrix[i, j] = 0  # Match
            else:
                matrix[i, j] = 1  # Mismatch
    
    # Create figure
    fig, ax = plt.subplots(figsize=(20, max(8, len(sequences) * 0.15)))
    
    # Custom colormap
    cmap_colors = [colors['gap'], colors['match'], colors['mismatch']]
    cmap = ListedColormap(cmap_colors)
    bounds = [-1.5, -0.5, 0.5, 1.5]
    norm = BoundaryNorm(bounds, cmap.N)
    
    im = ax.imshow(matrix, aspect='auto', cmap=cmap, norm=norm, interpolation='nearest')
    
    # Labels
    ax.set_yticks(range(len(names)))
    ylabels = []
    for i in range(len(names)):
        if i == ref_idx:
            ylabels.append(f"*REF {names[i]}")
        else:
            ylabels.append(names[i])
    ax.set_yticklabels(ylabels, fontsize=font_size)
    
    # Color reference label
    ax.get_yticklabels()[ref_idx].set_color(colors['reference_label'])
    ax.get_yticklabels()[ref_idx].set_weight('bold')
    
    # X-axis position labels
    if len(positions) <= 100:
        step = 5
    elif len(positions) <= 500:
        step = 25
    else:
        step = 50
    
    xtick_positions = [i for i, pos in enumerate(positions) if pos % step == 0]
    xtick_labels = [pos + 1 for pos in positions if pos % step == 0]
    ax.set_xticks(xtick_positions)
    ax.set_xticklabels(xtick_labels, fontsize=font_size)
    
    ax.set_xlabel('Alignment Position', fontsize=12, fontweight='bold')
    ax.set_ylabel('Sequences', fontsize=12, fontweight='bold')
    ax.set_title(f'Sequence Divergence from Reference: {names[ref_idx]}', 
                 fontsize=14, fontweight='bold', pad=20)
    
    # Remove spines for cleaner look
    for spine in ax.spines.values():
        spine.set_visible(False)
    
    # Colorbar with custom labels
    cbar = plt.colorbar(im, ax=ax, ticks=[-1, 0, 1], fraction=0.02, pad=0.02)
    cbar.ax.set_yticklabels(['Gap', 'Match', 'Mismatch'], fontsize=font_size)
    cbar.outline.set_visible(False)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight', format=format)
    print(f"Heatmap saved to {output_file}")
    plt.close()

def create_conservation_plot(sequences, names, ref_idx=0, output_file='conservation_plot.png',
                            colors=None, format='png', dpi=300, line_width=2.0, 
                            show_aa_labels=True, font_size=10):
    """Create position-wise conservation plot"""
    if colors is None:
        colors = COLOR_SCHEMES['default']
    
    if format in ['pdf', 'svg']:
        mpl.rcParams['pdf.fonttype'] = 42
        mpl.rcParams['ps.fonttype'] = 42
    
    conservation_scores, position_divergence = calculate_position_conservation(sequences, ref_idx)
    ref_seq = sequences[ref_idx]
    positions = range(len(ref_seq))
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 8), sharex=True, 
                                    gridspec_kw={'height_ratios': [3, 1]})
    
    # Conservation plot
    ax1.plot(positions, conservation_scores, linewidth=line_width, 
             color=colors.get('match', '#2563EB'))
    ax1.fill_between(positions, conservation_scores, alpha=0.3, 
                     color=colors.get('match', '#3B82F6'))
    ax1.set_ylabel('Conservation (%)', fontsize=12, fontweight='bold')
    ax1.set_title(f'Position-wise Conservation Relative to Reference: {names[ref_idx]}', 
                  fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3, linestyle=':')
    ax1.set_ylim(0, 100)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    
    # Reference sequence track with color coding
    colors_map = []
    for pos in positions:
        if ref_seq[pos] == '-':
            colors_map.append(colors.get('gap', '#E5E7EB'))
        elif conservation_scores[pos] >= 90:
            colors_map.append('#10B981')  # Highly conserved
        elif conservation_scores[pos] >= 70:
            colors_map.append('#FCD34D')  # Moderately conserved
        else:
            colors_map.append('#EF4444')  # Variable
    
    ax2.bar(positions, [1] * len(positions), color=colors_map, width=1.0, edgecolor='none')
    ax2.set_ylabel('Reference\nSequence', fontsize=11, fontweight='bold')
    ax2.set_xlabel('Alignment Position', fontsize=12, fontweight='bold')
    ax2.set_ylim(0, 1)
    ax2.set_yticks([])
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    
    # Add amino acid labels if requested and not too many positions
    if show_aa_labels and len(positions) <= 100:
        for i, aa in enumerate(ref_seq):
            if aa != '-':
                ax2.text(i, 0.5, aa, ha='center', va='center', 
                        fontsize=font_size, fontfamily='monospace', 
                        fontweight='bold', color='white' if conservation_scores[i] >= 90 else 'black')
    
    # Legend
    legend_elements = [
        mpatches.Patch(facecolor='#10B981', label='Highly conserved (≥90%)'),
        mpatches.Patch(facecolor='#FCD34D', label='Moderately conserved (70-90%)'),
        mpatches.Patch(facecolor='#EF4444', label='Variable (<70%)'),
        mpatches.Patch(facecolor=colors.get('gap', '#E5E7EB'), label='Gap in reference')
    ]
    ax2.legend(handles=legend_elements, loc='upper right', fontsize=10, frameon=False)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight', format=format)
    print(f"Conservation plot saved to {output_file}")
    plt.close()

def create_divergence_distribution(sequences, names, ref_idx=0, 
                                   output_file='divergence_distribution.png',
                                   colors=None, format='png', dpi=300, bins=50):
    """Create distribution plot of sequence divergences"""
    if colors is None:
        colors = COLOR_SCHEMES['default']
    
    if format in ['pdf', 'svg']:
        mpl.rcParams['pdf.fonttype'] = 42
        mpl.rcParams['ps.fonttype'] = 42
    
    divergences = calculate_divergence_matrix(sequences, ref_idx)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # Histogram
    ax1.hist(divergences, bins=bins, color=colors.get('mismatch', '#3B82F6'), 
             alpha=0.7, edgecolor='black', linewidth=0.5)
    ax1.axvline(divergences[ref_idx], color=colors.get('reference_label', 'red'), 
                linestyle='--', linewidth=2, 
                label=f'Reference ({names[ref_idx]})')
    ax1.set_xlabel('Divergence from Reference (%)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Number of Sequences', fontsize=12, fontweight='bold')
    ax1.set_title('Distribution of Sequence Divergence', fontsize=14, fontweight='bold')
    ax1.legend(frameon=False)
    ax1.grid(True, alpha=0.3, linestyle=':')
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    
    # Box plot
    bp = ax2.boxplot(divergences, vert=True, patch_artist=True,
                     boxprops=dict(facecolor=colors.get('mismatch', '#3B82F6'), 
                                  alpha=0.7, linewidth=1.5),
                     medianprops=dict(color='red', linewidth=2),
                     whiskerprops=dict(linewidth=1.5),
                     capprops=dict(linewidth=1.5))
    ax2.set_ylabel('Divergence from Reference (%)', fontsize=12, fontweight='bold')
    ax2.set_title('Divergence Summary Statistics', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3, axis='y', linestyle=':')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.set_xticklabels([''])
    
    # Add statistics text
    stats_text = f'Mean: {np.mean(divergences):.2f}%\n'
    stats_text += f'Median: {np.median(divergences):.2f}%\n'
    stats_text += f'Std Dev: {np.std(divergences):.2f}%\n'
    stats_text += f'Min: {np.min(divergences):.2f}%\n'
    stats_text += f'Max: {np.max(divergences):.2f}%'
    
    ax2.text(1.15, 0.5, stats_text, transform=ax2.transAxes,
            fontsize=10, verticalalignment='center', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight', format=format)
    print(f"Divergence distribution saved to {output_file}")
    plt.close()

def export_divergence_data(sequences, names, ref_idx, output_file='divergence_data.txt'):
    """Export divergence data for custom analysis"""
    divergences = calculate_divergence_matrix(sequences, ref_idx)
    conservation_scores, _ = calculate_position_conservation(sequences, ref_idx)
    
    with open(output_file, 'w') as f:
        # Write sequence divergences
        f.write("# Sequence Divergences\n")
        f.write("Sequence_Name\tDivergence_Percent\n")
        for name, div in zip(names, divergences):
            f.write(f"{name}\t{div:.4f}\n")
        
        # Write position conservation
        f.write("\n# Position Conservation\n")
        f.write("Position\tConservation_Percent\n")
        for i, cons in enumerate(conservation_scores, 1):
            if not np.isnan(cons):
                f.write(f"{i}\t{cons:.4f}\n")
    
    print(f"Divergence data exported to {output_file}")

def main():
    parser = argparse.ArgumentParser(
        description='Visualize protein alignment divergence from a reference sequence',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with SVG output
  python alignment_viz.py alignment.fasta --format svg
  
  # Specify reference by name with custom colors
  python alignment_viz.py alignment.fasta --ref "Sequence_1" --color-scheme colorblind
  
  # Publication-ready PDF with custom styling
  python alignment_viz.py alignment.fasta --format pdf --dpi 600 --color-scheme publication
  
  # Custom colors (JSON format)
  python alignment_viz.py alignment.fasta --custom-colors '{"match":"#00FF00","mismatch":"#FF0000"}'
  
  # Limit visualization size
  python alignment_viz.py alignment.fasta --max-seqs 50 --sample-positions 200

Available color schemes: default, publication, colorblind, vibrant, chemistry

Custom colors JSON keys:
  - match: Color for matching residues
  - mismatch: Color for mismatched residues
  - gap: Color for gaps
  - reference_label: Color for reference sequence label
  - text: General text color
  - grid: Grid line color
        """
    )
    
    parser.add_argument('alignment', type=str, help='Input alignment file (FASTA format)')
    parser.add_argument('--ref', type=str, help='Reference sequence name')
    parser.add_argument('--ref-index', type=int, default=0, 
                       help='Reference sequence index (0-based, default: 0)')
    parser.add_argument('--output-prefix', type=str, default='alignment_viz',
                       help='Output file prefix (default: alignment_viz)')
    parser.add_argument('--max-seqs', type=int, default=100,
                       help='Maximum sequences to show in heatmap (default: 100)')
    parser.add_argument('--sample-positions', type=int, default=None,
                       help='Sample N positions for heatmap (default: all)')
    
    # Color and formatting options
    parser.add_argument('--color-scheme', type=str, default='default',
                       choices=list(COLOR_SCHEMES.keys()),
                       help='Color scheme preset (default: default)')
    parser.add_argument('--custom-colors', type=str, default=None,
                       help='Custom colors as JSON string')
    parser.add_argument('--format', type=str, default='png',
                       choices=['png', 'pdf', 'svg', 'eps'],
                       help='Output format (default: png). Use pdf/svg for Illustrator')
    parser.add_argument('--dpi', type=int, default=300,
                       help='DPI for raster formats (default: 300)')
    parser.add_argument('--line-width', type=float, default=2.0,
                       help='Line width for plots (default: 2.0)')
    parser.add_argument('--font-size', type=int, default=8,
                       help='Font size for labels (default: 8)')
    parser.add_argument('--no-aa-labels', action='store_true',
                       help='Hide amino acid labels on conservation plot')
    parser.add_argument('--export-data', action='store_true',
                       help='Export divergence data to text file')
    
    args = parser.parse_args()
    
    # Parse alignment
    print(f"Reading alignment from {args.alignment}...")
    sequences, names = parse_alignment(args.alignment)
    print(f"Loaded {len(sequences)} sequences, alignment length: {len(sequences[0])}")
    
    # Determine reference index
    ref_idx = args.ref_index
    if args.ref:
        try:
            ref_idx = names.index(args.ref)
            print(f"Using reference sequence: {args.ref} (index {ref_idx})")
        except ValueError:
            print(f"Warning: Reference '{args.ref}' not found, using index {ref_idx}")
    else:
        print(f"Using reference sequence: {names[ref_idx]} (index {ref_idx})")
    
    # Load color scheme
    colors = load_color_scheme(args.color_scheme, args.custom_colors)
    print(f"Using color scheme: {args.color_scheme}")
    print(f"Output format: {args.format}")
    
    # Create visualizations
    print("\nGenerating visualizations...")
    
    create_divergence_heatmap(
        sequences, names, ref_idx, 
        f'{args.output_prefix}_heatmap.{args.format}',
        max_seqs=args.max_seqs,
        sample_positions=args.sample_positions,
        colors=colors,
        format=args.format,
        dpi=args.dpi,
        font_size=args.font_size
    )
    
    create_conservation_plot(
        sequences, names, ref_idx,
        f'{args.output_prefix}_conservation.{args.format}',
        colors=colors,
        format=args.format,
        dpi=args.dpi,
        line_width=args.line_width,
        show_aa_labels=not args.no_aa_labels,
        font_size=args.font_size
    )
    
    create_divergence_distribution(
        sequences, names, ref_idx,
        f'{args.output_prefix}_distribution.{args.format}',
        colors=colors,
        format=args.format,
        dpi=args.dpi
    )
    
    # Export data if requested
    if args.export_data:
        export_divergence_data(sequences, names, ref_idx,
                              f'{args.output_prefix}_data.txt')
    
    # Calculate and print summary statistics
    divergences = calculate_divergence_matrix(sequences, ref_idx)
    print("\n" + "="*60)
    print("DIVERGENCE SUMMARY")
    print("="*60)
    print(f"Reference sequence: {names[ref_idx]}")
    print(f"Total sequences: {len(sequences)}")
    print(f"Alignment length: {len(sequences[0])}")
    print(f"\nDivergence statistics:")
    print(f"  Mean:   {np.mean(divergences):.2f}%")
    print(f"  Median: {np.median(divergences):.2f}%")
    print(f"  Std:    {np.std(divergences):.2f}%")
    print(f"  Range:  {np.min(divergences):.2f}% - {np.max(divergences):.2f}%")
    print("="*60)
    print("\nAll visualizations completed!")
    
    if args.format in ['pdf', 'svg']:
        print(f"\nVector format ({args.format}) generated - fully editable in Illustrator/Inkscape!")

if __name__ == "__main__":
    main()
