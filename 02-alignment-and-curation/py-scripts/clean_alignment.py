#!/usr/bin/env python3
"""
Clean Protein Alignment Visualizer
Publication-quality alignment visualization with dot notation for identical positions
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib as mpl
from matplotlib.patches import Rectangle
from Bio import SeqIO
import argparse
from pathlib import Path

# Color schemes - light and clean
COLOR_SCHEMES = {
    'light_pastels': {
        'hydrophobic': '#FFE5B4',  # Peach
        'polar': '#B4E7FF',        # Light sky blue
        'positive': '#FFB4D4',     # Light pink
        'negative': '#E5B4FF',     # Light lavender
        'special': '#C4FFC4',      # Light mint
        'gap': '#F5F5F5',          # Very light gray
        'identical': '#FFFFFF',     # White (shown as dot)
        'text': '#2C3E50',         # Dark blue-gray
        'background': '#FFFFFF'
    },
    'soft_blue': {
        'hydrophobic': '#E8F4F8',  # Ice blue
        'polar': '#A7C7E7',        # Powder blue
        'positive': '#D4E6F1',     # Light periwinkle
        'negative': '#AED6F1',     # Baby blue
        'special': '#85C1E9',      # Sky blue
        'gap': '#F8F9FA',          # Off white
        'identical': '#FFFFFF',
        'text': '#2C3E50',
        'background': '#FFFFFF'
    },
    'minimal_gray': {
        'hydrophobic': '#F0F0F0',  # Light gray
        'polar': '#E0E0E0',        # Medium-light gray
        'positive': '#D0D0D0',     # Medium gray
        'negative': '#C0C0C0',     # Gray
        'special': '#B0B0B0',      # Darker gray
        'gap': '#FAFAFA',          # Almost white
        'identical': '#FFFFFF',
        'text': '#2C3E50',
        'background': '#FFFFFF'
    },
    'highlight_diff': {
        'hydrophobic': '#FFFFFF',  # White for identical
        'polar': '#FFFFFF',
        'positive': '#FFFFFF',
        'negative': '#FFFFFF',
        'special': '#FFFFFF',
        'gap': '#F5F5F5',
        'different': '#FFE5E5',    # Light red for differences
        'identical': '#FFFFFF',
        'text': '#2C3E50',
        'background': '#FFFFFF'
    },
    'conservation': {
        'conserved': '#E8F5E9',    # Light green
        'similar': '#FFF9C4',      # Light yellow
        'different': '#FFEBEE',    # Light red
        'gap': '#F5F5F5',
        'identical': '#FFFFFF',
        'text': '#2C3E50',
        'background': '#FFFFFF'
    }
}

# Amino acid groups
AA_GROUPS = {
    'hydrophobic': 'AILMFWV',
    'polar': 'STNQ',
    'positive': 'KRH',
    'negative': 'DE',
    'special': 'CGP'
}

def parse_alignment(fasta_file):
    """Parse FASTA alignment file"""
    sequences = []
    names = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(str(record.seq).upper())
        names.append(record.id)
    return sequences, names

def get_aa_group(aa):
    """Get amino acid group"""
    for group, aas in AA_GROUPS.items():
        if aa in aas:
            return group
    return 'special'

def get_aa_color(aa, is_identical, ref_aa, color_scheme, highlight_mode='chemistry'):
    """
    Get color for amino acid based on scheme and whether it's identical to reference
    
    highlight_mode:
    - chemistry: Color by chemical properties
    - difference: Only highlight differences
    - conservation: Color by conservation type
    """
    colors = color_scheme
    
    if aa == '-':
        return colors['gap']
    
    if is_identical:
        return colors['identical']
    
    if highlight_mode == 'chemistry':
        group = get_aa_group(aa)
        return colors.get(group, colors.get('special', '#FFFFFF'))
    
    elif highlight_mode == 'difference':
        return colors.get('different', '#FFE5E5')
    
    elif highlight_mode == 'conservation':
        # Check if similar properties
        if ref_aa == '-':
            return colors.get('different', '#FFEBEE')
        ref_group = get_aa_group(ref_aa)
        aa_group = get_aa_group(aa)
        
        if aa_group == ref_group:
            return colors.get('similar', '#FFF9C4')
        else:
            return colors.get('different', '#FFEBEE')
    
    return '#FFFFFF'

def create_clean_alignment(sequences, names, output_file='alignment.png',
                          color_scheme_name='light_pastels',
                          format='png', dpi=300,
                          reference_idx=0,
                          use_dots=True,
                          highlight_mode='chemistry',
                          show_ruler=True,
                          block_size=10,
                          font_size=9,
                          show_conservation=True):
    """
    Create clean alignment visualization
    """
    if format in ['pdf', 'svg']:
        mpl.rcParams['pdf.fonttype'] = 42
        mpl.rcParams['ps.fonttype'] = 42
    
    colors = COLOR_SCHEMES.get(color_scheme_name, COLOR_SCHEMES['light_pastels'])
    
    n_seqs = len(sequences)
    seq_len = len(sequences[0])
    ref_seq = sequences[reference_idx]
    
    # Calculate figure size
    char_width = 0.25
    char_height = 0.35
    label_width = 3.0
    ruler_height = 0.5 if show_ruler else 0
    conservation_height = 0.3 if show_conservation else 0
    
    fig_width = label_width + (seq_len * char_width) + 1
    fig_height = ruler_height + (n_seqs * char_height) + conservation_height + 1
    
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    ax.set_xlim(-0.5, seq_len + 0.5)
    ax.set_ylim(-0.5, n_seqs + 0.5 + (1 if show_ruler else 0))
    
    # Draw sequences
    for seq_idx, (seq, name) in enumerate(zip(sequences, names)):
        y_pos = n_seqs - seq_idx - 1
        
        # Draw sequence name
        ax.text(-1, y_pos, name, fontsize=font_size, 
               ha='right', va='center',
               fontweight='bold' if seq_idx == reference_idx else 'normal',
               color=colors['text'])
        
        # Mark reference
        if seq_idx == reference_idx:
            ax.text(-0.5, y_pos, '*', fontsize=font_size*1.5, 
                   ha='right', va='center', color='#DC2626', fontweight='bold')
        
        # Draw amino acids
        for pos, aa in enumerate(seq):
            ref_aa = ref_seq[pos]
            is_identical = (aa == ref_aa) and (seq_idx != reference_idx)
            
            # Get color
            bg_color = get_aa_color(aa, is_identical, ref_aa, colors, highlight_mode)
            
            # Draw background
            if bg_color != colors['background']:
                rect = Rectangle((pos - 0.4, y_pos - 0.35), 0.8, 0.7,
                               facecolor=bg_color, edgecolor='none')
                ax.add_patch(rect)
            
            # Add subtle border for blocks
            if pos % block_size == 0 and pos > 0:
                ax.axvline(x=pos - 0.5, color='#E0E0E0', linewidth=0.5, 
                          linestyle='-', alpha=0.5)
            
            # Draw amino acid or dot
            if use_dots and is_identical:
                display_char = '·'  # Middle dot
                text_color = '#999999'
            else:
                display_char = aa
                text_color = colors['text'] if aa != '-' else '#CCCCCC'
            
            ax.text(pos, y_pos, display_char, fontsize=font_size,
                   ha='center', va='center', color=text_color,
                   fontfamily='monospace', fontweight='normal')
    
    # Draw ruler
    if show_ruler:
        ruler_y = n_seqs + 0.3
        for pos in range(0, seq_len, 10):
            ax.text(pos, ruler_y, str(pos + 1), fontsize=font_size - 1,
                   ha='left', va='bottom', color='#666666')
            ax.plot([pos, pos], [ruler_y - 0.1, ruler_y - 0.2], 
                   color='#666666', linewidth=0.5)
        
        # Draw ruler line
        ax.plot([0, seq_len], [ruler_y - 0.15, ruler_y - 0.15], 
               color='#CCCCCC', linewidth=0.5)
    
    # Draw conservation bar
    if show_conservation:
        cons_y = -0.8
        for pos in range(seq_len):
            col = [sequences[i][pos] for i in range(n_seqs)]
            unique = len(set([aa for aa in col if aa != '-']))
            
            if unique == 1:
                cons_color = '#10B981'  # Fully conserved - green
                cons_height = 0.4
            elif unique == 2:
                cons_color = '#FCD34D'  # Moderately conserved - yellow
                cons_height = 0.3
            else:
                cons_color = '#EF4444'  # Variable - red
                cons_height = 0.2
            
            rect = Rectangle((pos - 0.4, cons_y), 0.8, cons_height,
                           facecolor=cons_color, edgecolor='none', alpha=0.6)
            ax.add_patch(rect)
        
        ax.text(-1, cons_y + 0.2, 'Conservation', fontsize=font_size - 1,
               ha='right', va='center', color=colors['text'])
    
    # Remove axes
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_facecolor(colors['background'])
    fig.patch.set_facecolor(colors['background'])
    
    # Add legend
    if highlight_mode == 'chemistry':
        legend_elements = [
            mpatches.Patch(facecolor=colors['hydrophobic'], label='Hydrophobic (AILMFWV)'),
            mpatches.Patch(facecolor=colors['polar'], label='Polar (STNQ)'),
            mpatches.Patch(facecolor=colors['positive'], label='Positive (KRH)'),
            mpatches.Patch(facecolor=colors['negative'], label='Negative (DE)'),
            mpatches.Patch(facecolor=colors['special'], label='Special (CGP)'),
            mpatches.Patch(facecolor=colors['gap'], label='Gap'),
        ]
        if use_dots:
            legend_elements.append(
                mpatches.Patch(facecolor='white', label='Identical to ref (·)')
            )
        
        ax.legend(handles=legend_elements, loc='upper right', 
                 bbox_to_anchor=(1.0, 1.0), fontsize=font_size - 1,
                 frameon=True, fancybox=True, shadow=False)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight', format=format,
               facecolor=colors['background'])
    print(f"Alignment saved to {output_file}")
    plt.close()

def calculate_identity_matrix(sequences, names):
    """Calculate pairwise percent identity matrix"""
    n = len(sequences)
    identity_matrix = np.zeros((n, n))
    
    for i in range(n):
        for j in range(n):
            if i == j:
                identity_matrix[i][j] = 100.0
            else:
                seq1, seq2 = sequences[i], sequences[j]
                matches = 0
                valid_positions = 0
                
                for aa1, aa2 in zip(seq1, seq2):
                    if aa1 != '-' and aa2 != '-':
                        valid_positions += 1
                        if aa1 == aa2:
                            matches += 1
                
                if valid_positions > 0:
                    identity = (matches / valid_positions) * 100
                else:
                    identity = 0.0
                
                identity_matrix[i][j] = identity
    
    return identity_matrix

def export_identity_matrix(identity_matrix, names, output_file='identity_matrix.txt'):
    """Export identity matrix to text file"""
    with open(output_file, 'w') as f:
        # Header
        f.write("# Pairwise Percent Identity Matrix\n")
        f.write("# Values represent percent identity between sequences\n\n")
        
        # Tab-delimited matrix
        f.write("Sequence\t" + "\t".join(names) + "\n")
        for i, name in enumerate(names):
            values = "\t".join(f"{identity_matrix[i, j]:.2f}" for j in range(len(names)))
            f.write(f"{name}\t{values}\n")
    
    print(f"Identity matrix saved to {output_file}")

def create_identity_heatmap(identity_matrix, names, output_file='identity_heatmap.png',
                           format='png', dpi=300, color_scheme_name='light_pastels'):
    """Create visual heatmap of identity matrix"""
    if format in ['pdf', 'svg']:
        mpl.rcParams['pdf.fonttype'] = 42
        mpl.rcParams['ps.fonttype'] = 42
    
    colors = COLOR_SCHEMES.get(color_scheme_name, COLOR_SCHEMES['light_pastels'])
    
    fig, ax = plt.subplots(figsize=(8, 7))
    
    # Create heatmap
    im = ax.imshow(identity_matrix, cmap='RdYlGn', aspect='auto', vmin=0, vmax=100)
    
    # Add text annotations
    for i in range(len(names)):
        for j in range(len(names)):
            text_color = 'white' if identity_matrix[i, j] < 50 else 'black'
            text = ax.text(j, i, f'{identity_matrix[i, j]:.1f}',
                          ha="center", va="center", color=text_color,
                          fontsize=10, fontweight='bold')
    
    # Labels
    ax.set_xticks(range(len(names)))
    ax.set_yticks(range(len(names)))
    ax.set_xticklabels(names, rotation=45, ha='right', fontsize=10)
    ax.set_yticklabels(names, fontsize=10)
    
    # Colorbar
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('Percent Identity (%)', fontsize=11, fontweight='bold')
    
    ax.set_title('Pairwise Percent Identity Matrix', fontsize=13, fontweight='bold', pad=15)
    
    # Grid
    ax.set_xticks(np.arange(len(names)) - 0.5, minor=True)
    ax.set_yticks(np.arange(len(names)) - 0.5, minor=True)
    ax.grid(which='minor', color='white', linestyle='-', linewidth=2)
    
    ax.set_facecolor(colors['background'])
    fig.patch.set_facecolor(colors['background'])
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight', format=format,
               facecolor=colors['background'])
    print(f"Identity heatmap saved to {output_file}")
    plt.close()

def create_compact_alignment(sequences, names, output_file='alignment_compact.png',
                            color_scheme_name='light_pastels',
                            format='png', dpi=300,
                            reference_idx=0,
                            positions_per_line=60):
    """
    Create compact multi-line alignment for long sequences
    """
    if format in ['pdf', 'svg']:
        mpl.rcParams['pdf.fonttype'] = 42
        mpl.rcParams['ps.fonttype'] = 42
    
    colors = COLOR_SCHEMES.get(color_scheme_name, COLOR_SCHEMES['light_pastels'])
    
    n_seqs = len(sequences)
    seq_len = len(sequences[0])
    ref_seq = sequences[reference_idx]
    
    # Calculate number of blocks
    n_blocks = (seq_len + positions_per_line - 1) // positions_per_line
    
    # Figure setup
    block_height = n_seqs * 0.3 + 0.8
    fig_height = n_blocks * block_height + 1
    fig_width = 16
    
    fig, axes = plt.subplots(n_blocks, 1, figsize=(fig_width, fig_height))
    if n_blocks == 1:
        axes = [axes]
    
    for block_idx in range(n_blocks):
        ax = axes[block_idx]
        
        start_pos = block_idx * positions_per_line
        end_pos = min(start_pos + positions_per_line, seq_len)
        block_len = end_pos - start_pos
        
        ax.set_xlim(-0.5, block_len + 0.5)
        ax.set_ylim(-0.5, n_seqs + 0.5)
        
        # Draw position numbers
        ax.text(-1, n_seqs - 0.5, f'{start_pos + 1}-{end_pos}', 
               fontsize=9, ha='right', va='top', 
               color=colors['text'], fontweight='bold')
        
        # Draw sequences for this block
        for seq_idx, (seq, name) in enumerate(zip(sequences, names)):
            y_pos = n_seqs - seq_idx - 1
            
            # Sequence name (only for first block)
            if block_idx == 0:
                ax.text(-1, y_pos, name, fontsize=9, 
                       ha='right', va='center',
                       fontweight='bold' if seq_idx == reference_idx else 'normal',
                       color=colors['text'])
            
            # Draw amino acids for this block
            for i, pos in enumerate(range(start_pos, end_pos)):
                aa = seq[pos]
                ref_aa = ref_seq[pos]
                is_identical = (aa == ref_aa) and (seq_idx != reference_idx)
                
                bg_color = get_aa_color(aa, is_identical, ref_aa, colors, 'chemistry')
                
                if bg_color != colors['background']:
                    rect = Rectangle((i - 0.4, y_pos - 0.35), 0.8, 0.7,
                                   facecolor=bg_color, edgecolor='none')
                    ax.add_patch(rect)
                
                display_char = '·' if is_identical else aa
                text_color = '#999999' if is_identical else colors['text']
                
                ax.text(i, y_pos, display_char, fontsize=9,
                       ha='center', va='center', color=text_color,
                       fontfamily='monospace')
        
        # Remove axes
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.set_facecolor(colors['background'])
    
    fig.patch.set_facecolor(colors['background'])
    plt.tight_layout()
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight', format=format,
               facecolor=colors['background'])
    print(f"Compact alignment saved to {output_file}")
    plt.close()

def main():
    parser = argparse.ArgumentParser(
        description='Create clean, publication-quality protein alignment visualizations',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with light pastels
  python clean_alignment.py sequences.fasta
  
  # Soft blue theme with SVG output
  python clean_alignment.py sequences.fasta --color-scheme soft_blue --format svg
  
  # Highlight only differences (minimal coloring)
  python clean_alignment.py sequences.fasta --highlight-mode difference
  
  # With identity matrix (text file)
  python clean_alignment.py sequences.fasta --identity-matrix
  
  # With identity heatmap visualization
  python clean_alignment.py sequences.fasta --identity-heatmap
  
  # Both matrix and heatmap
  python clean_alignment.py sequences.fasta --identity-matrix --identity-heatmap
  
  # Compact multi-line format for long sequences
  python clean_alignment.py sequences.fasta --compact --positions-per-line 60
  
  # Custom reference and no dots
  python clean_alignment.py sequences.fasta --ref-index 2 --no-dots

Color schemes: light_pastels, soft_blue, minimal_gray, highlight_diff, conservation
Highlight modes: chemistry, difference, conservation
        """
    )
    
    parser.add_argument('alignment', type=str, help='FASTA alignment file')
    parser.add_argument('--ref-index', type=int, default=0,
                       help='Reference sequence index (default: 0)')
    parser.add_argument('--color-scheme', type=str, default='light_pastels',
                       choices=list(COLOR_SCHEMES.keys()),
                       help='Color scheme (default: light_pastels)')
    parser.add_argument('--highlight-mode', type=str, default='chemistry',
                       choices=['chemistry', 'difference', 'conservation'],
                       help='Highlighting mode (default: chemistry)')
    parser.add_argument('--no-dots', action='store_true',
                       help='Show amino acids instead of dots for identical positions')
    parser.add_argument('--no-ruler', action='store_true',
                       help='Hide position ruler')
    parser.add_argument('--no-conservation', action='store_true',
                       help='Hide conservation bar')
    parser.add_argument('--compact', action='store_true',
                       help='Use compact multi-line format')
    parser.add_argument('--positions-per-line', type=int, default=60,
                       help='Positions per line in compact mode (default: 60)')
    parser.add_argument('--format', type=str, default='png',
                       choices=['png', 'pdf', 'svg', 'eps'],
                       help='Output format (default: png)')
    parser.add_argument('--dpi', type=int, default=300,
                       help='DPI for raster formats (default: 300)')
    parser.add_argument('--font-size', type=int, default=9,
                       help='Font size (default: 9)')
    parser.add_argument('--output', type=str, default=None,
                       help='Output filename (default: auto-generated)')
    parser.add_argument('--identity-matrix', action='store_true',
                       help='Generate percent identity matrix (text file)')
    parser.add_argument('--identity-heatmap', action='store_true',
                       help='Generate identity matrix heatmap visualization')
    
    args = parser.parse_args()
    
    # Parse alignment
    print(f"Reading alignment from {args.alignment}...")
    sequences, names = parse_alignment(args.alignment)
    print(f"Loaded {len(sequences)} sequences, length: {len(sequences[0])}")
    
    # Determine output filename
    if args.output:
        output_file = args.output
    else:
        input_path = Path(args.alignment)
        suffix = '_compact' if args.compact else ''
        output_file = f"{input_path.stem}_alignment{suffix}.{args.format}"
    
    print(f"Reference sequence: {names[args.ref_index]}")
    print(f"Color scheme: {args.color_scheme}")
    print(f"Highlight mode: {args.highlight_mode}")
    
    # Calculate identity matrix if requested
    if args.identity_matrix or args.identity_heatmap:
        print("\nCalculating pairwise identity matrix...")
        identity_matrix = calculate_identity_matrix(sequences, names)
        
        if args.identity_matrix:
            input_path = Path(args.alignment)
            matrix_file = f"{input_path.stem}_identity_matrix.txt"
            export_identity_matrix(identity_matrix, names, matrix_file)
        
        if args.identity_heatmap:
            input_path = Path(args.alignment)
            heatmap_file = f"{input_path.stem}_identity_heatmap.{args.format}"
            create_identity_heatmap(identity_matrix, names, heatmap_file,
                                   args.format, args.dpi, args.color_scheme)
    
    # Create visualization
    if args.compact:
        create_compact_alignment(
            sequences, names, output_file,
            args.color_scheme, args.format, args.dpi,
            args.ref_index, args.positions_per_line
        )
    else:
        create_clean_alignment(
            sequences, names, output_file,
            args.color_scheme, args.format, args.dpi,
            args.ref_index,
            use_dots=not args.no_dots,
            highlight_mode=args.highlight_mode,
            show_ruler=not args.no_ruler,
            show_conservation=not args.no_conservation,
            font_size=args.font_size
        )
    
    print(f"\n✓ Visualization complete: {output_file}")

if __name__ == "__main__":
    main()
