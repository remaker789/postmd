#!/usr/bin/env python
"""
Extract specific elements from XYZ file while preserving XYZ format.
For Windows systems.

Usage:
    python extract_xyz.py -i input.xyz -e "Na K" [-o output.xyz]
"""

import argparse
import os
import sys

def parse_args():
    parser = argparse.ArgumentParser(description='Extract specific elements from XYZ file')
    parser.add_argument('-i', '--input', required=True, help='Input XYZ file')
    parser.add_argument('-e', '--elements', required=True, help='Space-separated list of elements to extract')
    parser.add_argument('-o', '--output', help='Output XYZ file (optional)')
    return parser.parse_args()

def generate_output_name(input_file, elements):
    """Generate default output filename from input filename and elements"""
    base, ext = os.path.splitext(input_file)
    elements_str = '_'.join(elements.split())
    return f"{base}_{elements_str}{ext}"

def extract_xyz(input_file, elements, output_file=None):
    """Extract specified elements from XYZ file"""
    if output_file is None:
        output_file = generate_output_name(input_file, elements)
    
    elements = set(elements.split())
    
    try:
        with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
            while True:
                # Read number of atoms
                try:
                    natoms = int(fin.readline())
                except (ValueError, EOFError):
                    break
                
                # Read comment line
                comment = fin.readline()
                
                # Read and filter atom coordinates
                selected_atoms = []
                for _ in range(natoms):
                    line = fin.readline()
                    if not line:
                        break
                    atom = line.split()
                    if atom[0] in elements:
                        selected_atoms.append(line)
                
                # Write frame if any atoms were selected
                if selected_atoms:
                    fout.write(f"{len(selected_atoms)}\n")
                    fout.write(comment)
                    for atom in selected_atoms:
                        fout.write(atom)
        
        print(f"Extracted XYZ file has been saved to {output_file}")
        
    except IOError as e:
        print(f"Error processing files: {str(e)}", file=sys.stderr)
        sys.exit(1)

def main():
    args = parse_args()
    
    # Check if input file exists
    if not os.path.isfile(args.input):
        print(f"Error: Input file '{args.input}' does not exist.", file=sys.stderr)
        sys.exit(1)
    
    extract_xyz(args.input, args.elements, args.output)

if __name__ == "__main__":
    main() 