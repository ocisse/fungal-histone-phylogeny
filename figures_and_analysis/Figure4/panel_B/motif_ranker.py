#!/usr/bin/env python3
"""
Regular Expression Motif Comparator
Compare motifs defined as regular expressions (e.g., from GLAM2)
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib as mpl
import re
import argparse
from collections import Counter
import json

COLOR_SCHEMES = {
    'default': {'primary': '#2563EB', 'secondary': '#DC2626', 'tertiary': '#10B981'},
    'pastel_blue': {'primary': '#A7C7E7', 'secondary': '#7EB3D9', 'tertiary': '#D4E6F1'},
    'publication': {'primary': '#2C3E50', 'secondary': '#7F8C8D', 'tertiary': '#95A5A6'}
}

def parse_regex_motifs(input_file):
    """
    Parse motif regular expressions from file
    Expected format:
    >Motif_1
    G.[AILV]{3}.G[KR][ST]T
    >Motif_2
    G[AT]...[KR]G[KR][ST]
    """
    motifs = []
    names = []
    
    with open(input_file, 'r') as f:
        current_name = None
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                current_name = line[1:]
            elif line and current_name:
                motifs.append(line)
                names.append(current_name)
                current_name = None
    
    return motifs, names

def regex_to_consensus(regex_pattern):
    """
    Convert regex pattern to consensus sequence
    Examples:
    - . or [A-Z] -> X (any amino acid)
    - [AILV] -> Most common in group (or first)
    - {n} -> Repeat n times
    - Literal characters -> Keep as is
    """
    consensus = ""
    i = 0
    
    while i < len(regex_pattern):
        char = regex_pattern[i]
        
        if char == '.':
            consensus += 'X'
            i += 1
        elif char == '[':
            # Find closing bracket
            j = regex_pattern.find(']', i)
            if j == -1:
                consensus += 'X'
                i += 1
                continue
            
            aa_group = regex_pattern[i+1:j]
            # Use first amino acid as representative
            if aa_group and aa_group[0] != '^':
                consensus += aa_group[0]
            else:
                consensus += 'X'
            i = j + 1
        elif char == '{':
            # Handle repetitions
            j = regex_pattern.find('}', i)
            if j == -1:
                i += 1
                continue
            
            try:
                repeat = int(regex_pattern[i+1:j])
                if consensus:
                    last_char = consensus[-1]
                    consensus = consensus[:-1] + (last_char * repeat)
            except ValueError:
                pass
            i = j + 1
        elif char.isalpha():
            consensus += char
            i += 1
        else:
            i += 1
    
    return consensus

def regex_to_profile(regex_pattern, alphabet='ACDEFGHIKLMNPQRSTVWY'):
    """
    Convert regex to position-specific scoring matrix (profile)
    Returns list of probability distributions for each position
    """
    profile = []
    i = 0
    
    while i < len(regex_pattern):
        char = regex_pattern[i]
        
        if char == '.':
            # Equal probability for all amino acids
            probs = {aa: 1.0/len(alphabet) for aa in alphabet}
            profile.append(probs)
            i += 1
        elif char == '[':
            j = regex_pattern.find(']', i)
            if j == -1:
                probs = {aa: 1.0/len(alphabet) for aa in alphabet}
                profile.append(probs)
                i += 1
                continue
            
            aa_group = regex_pattern[i+1:j]
            # Equal probability within group, zero outside
            probs = {aa: 0.0 for aa in alphabet}
            if aa_group:
                group_aas = [c for c in aa_group if c.isalpha()]
                if group_aas:
                    for aa in group_aas:
                        if aa in alphabet:
                            probs[aa] = 1.0 / len(group_aas)
            profile.append(probs)
            i = j + 1
        elif char == '{':
            j = regex_pattern.find('}', i)
            if j == -1:
                i += 1
                continue
            
            try:
                repeat = int(regex_pattern[i+1:j])
                if profile:
                    last_probs = profile[-1]
                    profile = profile[:-1] + [last_probs] * repeat
            except ValueError:
                pass
            i = j + 1
        elif char.isalpha():
            # Specific amino acid - 100% probability
            probs = {aa: 0.0 for aa in alphabet}
            probs[char] = 1.0
            profile.append(probs)
            i += 1
        else:
            i += 1
    
    return profile

def calculate_profile_similarity(profile1, profile2):
    """
    Calculate similarity between two profiles using Jensen-Shannon divergence
    Returns similarity score (0-100)
    """
    from scipy.spatial.distance import jensenshannon
    
    # Pad to same length
    max_len = max(len(profile1), len(profile2))
    alphabet = list(profile1[0].keys()) if profile1 else list(profile2[0].keys())
    
    # Extend shorter profile with uniform distribution
    uniform = {aa: 1.0/len(alphabet) for aa in alphabet}
    while len(profile1) < max_len:
        profile1.append(uniform)
    while len(profile2) < max_len:
        profile2.append(uniform)
    
    # Calculate position-wise divergence
    divergences = []
    for pos1, pos2 in zip(profile1, profile2):
        vec1 = np.array([pos1[aa] for aa in alphabet])
        vec2 = np.array([pos2[aa] for aa in alphabet])
        
        # Jensen-Shannon divergence (0 = identical, 1 = completely different)
        js_div = jensenshannon(vec1, vec2, base=2)
        divergences.append(js_div)
    
    # Average divergence across positions
    avg_divergence = np.mean(divergences)
    
    # Convert to similarity (0-100 scale)
    similarity = (1 - avg_divergence) * 100
    
    return similarity

def calculate_consensus_similarity(consensus1, consensus2):
    """
    Simple similarity based on consensus sequences
    """
    max_len = max(len(consensus1), len(consensus2))
    
    # Pad shorter sequence
    cons1 = consensus1.ljust(max_len, '-')
    cons2 = consensus2.ljust(max_len, '-')
    
    matches = 0
    valid_positions = 0
    
    for aa1, aa2 in zip(cons1, cons2):
        if aa1 != '-' and aa2 != '-':
            valid_positions += 1
            if aa1 == aa2:
                matches += 1
            elif aa1 == 'X' or aa2 == 'X':
                matches += 0.5  # Partial credit for wildcards
    
    if valid_positions > 0:
        similarity = (matches / valid_positions) * 100
    else:
        similarity = 0
    
    return similarity

def compare_regex_motifs(regex_motifs, method='profile'):
    """
    Compare regex motifs using specified method
    
    Methods:
    - profile: Use position-specific profiles (most accurate)
    - consensus: Use consensus sequences (simpler)
    """
    n = len(regex_motifs)
    similarity_matrix = np.zeros((n, n))
    
    if method == 'profile':
        # Convert to profiles
        profiles = [regex_to_profile(regex) for regex in regex_motifs]
        
        for i in range(n):
            similarity_matrix[i][i] = 100
            for j in range(i+1, n):
                try:
                    similarity = calculate_profile_similarity(profiles[i], profiles[j])
                except:
                    # Fallback to consensus method
                    consensus1 = regex_to_consensus(regex_motifs[i])
                    consensus2 = regex_to_consensus(regex_motifs[j])
                    similarity = calculate_consensus_similarity(consensus1, consensus2)
                
                similarity_matrix[i][j] = similarity
                similarity_matrix[j][i] = similarity
    else:  # consensus
        consensuses = [regex_to_consensus(regex) for regex in regex_motifs]
        
        for i in range(n):
            similarity_matrix[i][i] = 100
            for j in range(i+1, n):
                similarity = calculate_consensus_similarity(consensuses[i], consensuses[j])
                similarity_matrix[i][j] = similarity
                similarity_matrix[j][i] = similarity
    
    return similarity_matrix

def create_regex_comparison_plot(regex_motifs, names, similarity_matrix,
                                 output_file='regex_comparison.png',
                                 colors=None, format='png', dpi=300):
    """Create comprehensive comparison visualization"""
    if colors is None:
        colors = COLOR_SCHEMES['default']
    
    if format in ['pdf', 'svg']:
        mpl.rcParams['pdf.fonttype'] = 42
        mpl.rcParams['ps.fonttype'] = 42
    
    n = len(names)
    fig = plt.figure(figsize=(16, max(10, n * 1.5)))
    gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3, height_ratios=[1.5, 1])
    
    # 1. Similarity heatmap
    ax1 = fig.add_subplot(gs[0, :])
    im = ax1.imshow(similarity_matrix, cmap='RdYlGn', aspect='auto', vmin=0, vmax=100)
    
    # Add values
    for i in range(n):
        for j in range(n):
            text = ax1.text(j, i, f'{similarity_matrix[i, j]:.1f}',
                          ha="center", va="center", 
                          color="black" if similarity_matrix[i, j] > 50 else "white",
                          fontsize=9, fontweight='bold')
    
    ax1.set_xticks(range(n))
    ax1.set_yticks(range(n))
    ax1.set_xticklabels(names, rotation=45, ha='right', fontsize=10)
    ax1.set_yticklabels(names, fontsize=10)
    ax1.set_title('Motif Similarity Matrix', fontsize=14, fontweight='bold', pad=20)
    
    cbar = plt.colorbar(im, ax=ax1, fraction=0.046, pad=0.04)
    cbar.set_label('Similarity (%)', fontsize=11, fontweight='bold')
    
    ax1.set_xticks(np.arange(n) - 0.5, minor=True)
    ax1.set_yticks(np.arange(n) - 0.5, minor=True)
    ax1.grid(which='minor', color='white', linestyle='-', linewidth=2)
    
    # 2. Regex patterns
    ax2 = fig.add_subplot(gs[1, 0])
    ax2.axis('off')
    
    text_content = "REGEX PATTERNS\n" + "="*40 + "\n\n"
    for i, (name, regex) in enumerate(zip(names, regex_motifs)):
        consensus = regex_to_consensus(regex)
        text_content += f"{name}:\n"
        text_content += f"  Regex: {regex}\n"
        text_content += f"  Consensus: {consensus}\n\n"
    
    ax2.text(0.05, 0.95, text_content, transform=ax2.transAxes,
            fontsize=9, verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # 3. Statistics
    ax3 = fig.add_subplot(gs[1, 1])
    ax3.axis('off')
    
    # Calculate stats
    avg_similarity = np.mean(similarity_matrix[np.triu_indices(n, k=1)])
    max_sim_idx = np.unravel_index(
        np.argmax(similarity_matrix - np.eye(n) * 100), 
        similarity_matrix.shape
    )
    min_sim_idx = np.unravel_index(
        np.argmin(similarity_matrix + np.eye(n) * 100), 
        similarity_matrix.shape
    )
    
    stats_text = "SUMMARY STATISTICS\n" + "="*40 + "\n\n"
    stats_text += f"Total motifs: {n}\n\n"
    stats_text += f"Average pairwise similarity:\n  {avg_similarity:.2f}%\n\n"
    stats_text += f"Most similar pair:\n"
    stats_text += f"  {names[max_sim_idx[0]]} vs {names[max_sim_idx[1]]}\n"
    stats_text += f"  {similarity_matrix[max_sim_idx]:.2f}%\n\n"
    stats_text += f"Most divergent pair:\n"
    stats_text += f"  {names[min_sim_idx[0]]} vs {names[min_sim_idx[1]]}\n"
    stats_text += f"  {similarity_matrix[min_sim_idx]:.2f}%\n"
    
    ax3.text(0.05, 0.95, stats_text, transform=ax3.transAxes,
            fontsize=10, verticalalignment='top', fontfamily='monospace',
            bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
    
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight', format=format)
    print(f"Comparison plot saved to {output_file}")
    plt.close()

def export_results(regex_motifs, names, similarity_matrix, output_prefix='regex_analysis'):
    """Export results to text files"""
    # Export consensus sequences
    with open(f'{output_prefix}_consensus.fasta', 'w') as f:
        for name, regex in zip(names, regex_motifs):
            consensus = regex_to_consensus(regex)
            f.write(f">{name}\n{consensus}\n")
    print(f"Consensus sequences saved to {output_prefix}_consensus.fasta")
    
    # Export similarity matrix
    with open(f'{output_prefix}_similarity.txt', 'w') as f:
        f.write("Motif\t" + "\t".join(names) + "\n")
        for i, name in enumerate(names):
            f.write(name + "\t" + "\t".join(f"{similarity_matrix[i, j]:.2f}" for j in range(len(names))) + "\n")
    print(f"Similarity matrix saved to {output_prefix}_similarity.txt")

def main():
    parser = argparse.ArgumentParser(
        description='Compare motifs defined as regular expressions',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Input file format (FASTA-style):
>Motif_1
G.[AILV]{3}.G[KR][ST]T
>Motif_2
G[AT]...[KR]G[KR][ST]

Examples:
  # Compare regex motifs
  python regex_motif_comparator.py glam2_regex.txt
  
  # Use profile-based comparison (more accurate)
  python regex_motif_comparator.py glam2_regex.txt --method profile
  
  # Generate SVG output
  python regex_motif_comparator.py glam2_regex.txt --format svg --export-data
        """
    )
    
    parser.add_argument('motifs', type=str, help='Input file with regex motifs')
    parser.add_argument('--method', type=str, default='consensus',
                       choices=['profile', 'consensus'],
                       help='Comparison method (default: consensus)')
    parser.add_argument('--output-prefix', type=str, default='regex_analysis',
                       help='Output file prefix')
    parser.add_argument('--color-scheme', type=str, default='default',
                       choices=list(COLOR_SCHEMES.keys()))
    parser.add_argument('--format', type=str, default='png',
                       choices=['png', 'pdf', 'svg', 'eps'])
    parser.add_argument('--dpi', type=int, default=300)
    parser.add_argument('--export-data', action='store_true',
                       help='Export consensus sequences and similarity matrix')
    
    args = parser.parse_args()
    
    # Parse input
    print(f"Reading regex motifs from {args.motifs}...")
    regex_motifs, names = parse_regex_motifs(args.motifs)
    print(f"Loaded {len(regex_motifs)} motifs:")
    for name, regex in zip(names, regex_motifs):
        consensus = regex_to_consensus(regex)
        print(f"  {name}: {regex} -> {consensus}")
    
    # Compare motifs
    print(f"\nComparing motifs using {args.method} method...")
    similarity_matrix = compare_regex_motifs(regex_motifs, args.method)
    
    # Create visualization
    print("Generating visualization...")
    create_regex_comparison_plot(
        regex_motifs, names, similarity_matrix,
        f'{args.output_prefix}_comparison.{args.format}',
        COLOR_SCHEMES[args.color_scheme], args.format, args.dpi
    )
    
    # Export if requested
    if args.export_data:
        export_results(regex_motifs, names, similarity_matrix, args.output_prefix)
    
    print("\n" + "="*60)
    print("ANALYSIS COMPLETE")
    print("="*60)
    avg_sim = np.mean(similarity_matrix[np.triu_indices(len(names), k=1)])
    print(f"Average pairwise similarity: {avg_sim:.2f}%")
    print("="*60)

if __name__ == "__main__":
    main()
