#!/usr/bin/env python
"""Extract specific elements from XYZ file while preserving XYZ format."""

import os
import sys

def generate_output_name(input_file, elements):
    """Generate default output filename from input filename and elements"""
    base, ext = os.path.splitext(input_file)
    elements_str = '_'.join(elements.split())
    return f"{base}_{elements_str}{ext}"

def extract_xyz(input_file, elements, output_file=None):
    """
    Extract specified elements from XYZ file.
    
    Parameters
    ----------
    input_file : str
        Path to input XYZ file
    elements : str
        Space-separated list of elements to extract
    output_file : str, optional
        Path to output XYZ file. If not provided, will generate automatically.
    """
    if output_file is None:
        output_file = generate_output_name(input_file, elements)
    
    if not os.path.isfile(input_file):
        raise FileNotFoundError(f"Input file '{input_file}' does not exist.")
    
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
        raise IOError(f"Error processing files: {str(e)}") 