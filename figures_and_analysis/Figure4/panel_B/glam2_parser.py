#!/usr/bin/env python3
"""
GLAM2 Output Parser
Extract and compare representative sequences from GLAM2 alignment output
"""

import re
import argparse
from collections import defaultdict
from pathlib import Path

def parse_glam2_output(glam2_file):
    """
    Parse GLAM2 output file and extract motifs
    Returns list of motifs, each containing aligned sequences
    """
    motifs = []
    current_motif = None
    in_alignment = False
    
    with open(glam2_file, 'r') as f:
        lines = f.readlines()
    
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        
        # Detect motif header (various GLAM2 formats)
        if re.match(r'(Motif\s+\d+|Score|Aligned)', line, re.IGNORECASE):
            if current_motif and current_motif['sequences']:
                motifs.append(current_motif)
            
            current_motif = {
                'name': f"Motif_{len(motifs) + 1}",
                'sequences': [],
                'consensus': '',
                'score': 0
            }
            in_alignment = True
        
        # Parse aligned sequences
        # Format: SeqName Position AlignedSeq EndPos Strand Score
        # Example: Spec1349|KAJ 20 LATK....AA.RKTA...AT...G.GVKKPH 38 + 66.1
        if in_alignment and line:
            parts = line.split()
            
            # Check if this is a sequence line (has sequence ID, positions, sequence, score)
            if len(parts) >= 4:
                # Try to identify sequence line
                # Usually: ID start_pos sequence end_pos [+/-] [score]
                seq_id = parts[0]
                
                # Find the aligned sequence (typically contains uppercase letters and dots)
                seq_pattern = None
                for part in parts[1:]:
                    if re.match(r'^[A-Z.]+$', part) and len(part) > 5:
                        seq_pattern = part
                        break
                
                if seq_pattern:
                    # Extract score if available
                    score = 0
                    try:
                        if len(parts) >= 5:
                            score = float(parts[-1])
                    except ValueError:
                        pass
                    
                    current_motif['sequences'].append({
                        'id': seq_id,
                        'aligned': seq_pattern,
                        'ungapped': seq_pattern.replace('.', ''),
                        'score': score
                    })
            
            # Check if line is consensus (typically lowercase or mixed with spaces)
            elif re.search(r'[A-Z]{2,}', line) and ' ' in line:
                # Might be consensus line - store it
                consensus = ''.join(c for c in line if c.isalpha())
                if len(consensus) > 5:
                    current_motif['consensus'] = consensus
        
        i += 1
    
    # Add last motif
    if current_motif and current_motif['sequences']:
        motifs.append(current_motif)
    
    return motifs

def build_consensus_from_alignment(sequences):
    """
    Build consensus sequence from aligned sequences
    Using most common amino acid at each position
    """
    if not sequences:
        return ""
    
    # Get alignment length
    max_len = max(len(seq['aligned']) for seq in sequences)
    
    consensus = []
    for pos in range(max_len):
        aas = []
        for seq in sequences:
            if pos < len(seq['aligned']):
                aa = seq['aligned'][pos]
                if aa != '.':
                    aas.append(aa)
        
        if aas:
            # Most common amino acid
            from collections import Counter
            most_common = Counter(aas).most_common(1)[0][0]
            consensus.append(most_common)
    
    return ''.join(consensus)

def get_representative_sequences(motifs, method='consensus'):
    """
    Get representative sequence for each motif
    
    Methods:
    - consensus: Build consensus from alignment
    - top_score: Use highest scoring sequence
    - first: Use first sequence
    """
    representatives = []
    
    for motif in motifs:
        if method == 'consensus':
            if motif['consensus']:
                seq = motif['consensus']
            else:
                seq = build_consensus_from_alignment(motif['sequences'])
            rep_id = f"{motif['name']}_consensus"
        
        elif method == 'top_score':
            # Find highest scoring sequence
            top_seq = max(motif['sequences'], key=lambda x: x['score'])
            seq = top_seq['ungapped']
            rep_id = f"{motif['name']}_{top_seq['id']}"
        
        else:  # first
            seq = motif['sequences'][0]['ungapped']
            rep_id = f"{motif['name']}_{motif['sequences'][0]['id']}"
        
        representatives.append({
            'id': rep_id,
            'sequence': seq,
            'motif_name': motif['name'],
            'num_sequences': len(motif['sequences'])
        })
    
    return representatives

def export_fasta(representatives, output_file):
    """Export representative sequences to FASTA"""
    with open(output_file, 'w') as f:
        for rep in representatives:
            f.write(f">{rep['id']}\n")
            f.write(f"{rep['sequence']}\n")
    print(f"Exported {len(representatives)} representative sequences to {output_file}")

def export_summary(motifs, representatives, output_file):
    """Export summary of motifs"""
    with open(output_file, 'w') as f:
        f.write("GLAM2 Motif Summary\n")
        f.write("="*80 + "\n\n")
        
        for motif, rep in zip(motifs, representatives):
            f.write(f"{motif['name']}\n")
            f.write(f"  Number of sequences: {len(motif['sequences'])}\n")
            f.write(f"  Representative: {rep['sequence']}\n")
            f.write(f"  Length: {len(rep['sequence'])}\n")
            
            if motif['sequences']:
                top_seq = max(motif['sequences'], key=lambda x: x['score'])
                f.write(f"  Top scoring sequence: {top_seq['id']} (score: {top_seq['score']:.2f})\n")
            
            f.write(f"\n  Aligned sequences:\n")
            for seq in motif['sequences'][:5]:  # Show top 5
                f.write(f"    {seq['id']:20s} {seq['aligned']} (score: {seq['score']:.2f})\n")
            
            if len(motif['sequences']) > 5:
                f.write(f"    ... and {len(motif['sequences']) - 5} more\n")
            
            f.write("\n" + "-"*80 + "\n\n")
    
    print(f"Summary saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(
        description='Parse GLAM2 output and extract representative sequences',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Extract consensus sequences (output in same directory as input)
  python glam2_parser.py cdhit_cluster_fastas/cluster_10.fasta.glam2_out/glam2.txt
  
  # Use top-scoring sequence from each motif
  python glam2_parser.py path/to/glam2.txt --method top_score
  
  # Custom output directory
  python glam2_parser.py path/to/glam2.txt --output-dir ./results
  
  # Custom output name
  python glam2_parser.py path/to/glam2.txt --output my_motifs

This will create files in the same directory as the input:
  - cluster_10_motifs_representatives.fasta (for use with motif_ranker.py)
  - cluster_10_motifs_summary.txt (detailed motif information)
        """
    )
    
    parser.add_argument('glam2_file', type=str, help='GLAM2 output file')
    parser.add_argument('--method', type=str, default='consensus',
                       choices=['consensus', 'top_score', 'first'],
                       help='Method to select representative (default: consensus)')
    parser.add_argument('--output', type=str, default=None,
                       help='Output file prefix (default: same dir as input with _motifs suffix)')
    parser.add_argument('--output-dir', type=str, default=None,
                       help='Output directory (default: same directory as input file)')
    
    args = parser.parse_args()
    
    # Determine output directory and prefix
    from pathlib import Path
    input_path = Path(args.glam2_file)
    
    if args.output_dir:
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
    else:
        # Use same directory as input file
        output_dir = input_path.parent
    
    if args.output:
        output_prefix = args.output
    else:
        # Create output prefix from input filename
        # e.g., cluster_10.fasta.glam2_out/glam2.txt -> cluster_10_motifs
        if 'cluster_' in str(input_path):
            cluster_name = re.search(r'cluster_\d+', str(input_path))
            if cluster_name:
                output_prefix = f"{cluster_name.group()}_motifs"
            else:
                output_prefix = input_path.stem + "_motifs"
        else:
            output_prefix = input_path.stem + "_motifs"
    
    # Construct full output paths
    output_fasta = output_dir / f"{output_prefix}_representatives.fasta"
    output_summary = output_dir / f"{output_prefix}_summary.txt"
    
    print(f"Parsing GLAM2 output from {args.glam2_file}...")
    print(f"Output directory: {output_dir}")
    print(f"Output prefix: {output_prefix}")
    motifs = parse_glam2_output(args.glam2_file)
    
    print(f"\nFound {len(motifs)} motifs:")
    for motif in motifs:
        print(f"  {motif['name']}: {len(motif['sequences'])} sequences")
    
    print(f"\nExtracting representatives using '{args.method}' method...")
    representatives = get_representative_sequences(motifs, args.method)
    
    print("\nRepresentative sequences:")
    for rep in representatives:
        print(f"  {rep['id']}: {rep['sequence']}")
    
    # Export files
    export_fasta(representatives, str(output_fasta))
    export_summary(motifs, representatives, str(output_summary))
    
    print("\n" + "="*80)
    print("EXTRACTION COMPLETE")
    print("="*80)
    print(f"Output files:")
    print(f"  {output_fasta}")
    print(f"  {output_summary}")
    print(f"\nNext steps:")
    print(f"  1. Review: {output_summary}")
    print(f"  2. Compare motifs:")
    print(f"     python motif_ranker.py {output_fasta}")
    print("="*80)

if __name__ == "__main__":
    main()
