#!/usr/bin/env python
"""
postmd: A toolkit for post-processing molecular dynamics simulations.
"""

import argparse
import sys
from . import lmp_snippets
from . import plumed_snippets
from . import cp2k_snippets

def main_parser():
    parser = argparse.ArgumentParser(
        description="postmd: A toolkit for post-processing molecular dynamics simulations."
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # LAMMPS snippets
    lmp_parser = subparsers.add_parser('lmp', help='LAMMPS snippets')
    lmp_parser.add_argument('--pre', type=str, default="all", 
                           help='Print pre-processing snippets')
    lmp_parser.add_argument('--post', type=str, default="all", 
                           help='Print post-processing snippets')
    
    # PLUMED snippets
    plumed_parser = subparsers.add_parser('plumed', help='PLUMED snippets')
    plumed_parser.add_argument('--pre', type=str, default="all", 
                              help='Print pre-processing snippets')
    plumed_parser.add_argument('--post', type=str, default="all", 
                              help='Print post-processing snippets')
    
    # CP2K snippets
    cp2k_parser = subparsers.add_parser('cp2k', help='CP2K snippets')
    cp2k_parser.add_argument('--pre', type=str, default="all", 
                            help='Print pre-processing snippets')
    cp2k_parser.add_argument('--post', type=str, default="all", 
                            help='Print post-processing snippets')
    
    return parser

def main():
    parser = main_parser()
    args = parser.parse_args()
    
    if args.command is None:
        parser.print_help()
        sys.exit(0)
        
    if args.command == 'lmp':
        if args.pre:
            lmp_snippets.handle_pre(args.pre)
        if args.post:
            lmp_snippets.handle_post(args.post)
    
    elif args.command == 'plumed':
        if args.pre:
            plumed_snippets.handle_pre(args.pre)
        if args.post:
            plumed_snippets.handle_post(args.post)
            
    elif args.command == 'cp2k':
        if args.pre:
            cp2k_snippets.handle_pre(args.pre)
        if args.post:
            cp2k_snippets.handle_post(args.post)

if __name__ == "__main__":
    main()
