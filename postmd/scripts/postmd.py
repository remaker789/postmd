#!/usr/bin/env python
"""
postmd: A toolkit for post-processing molecular dynamics simulations.
"""

import argparse
import sys
import platform
from . import lmp_snippets
from . import plumed_snippets
from . import cp2k_snippets
from ..dump.extract_xyz import extract_xyz  # 导入 extract_xyz 函数

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
    
    # 添加 extract 子命令
    extract_parser = subparsers.add_parser('extract', help='Extract elements from trajectory files')
    extract_parser.add_argument('-i', '--input', required=True, help='Input XYZ file')
    extract_parser.add_argument('-e', '--elements', required=True, help='Space-separated list of elements to extract')
    extract_parser.add_argument('-o', '--output', help='Output XYZ file (optional)')
    
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

    elif args.command == 'extract':
        extract_xyz(args.input, args.elements, args.output)
        # 根据系统类型选择提取方法
        # if platform.system() == 'Windows':
            # extract_xyz(args.input, args.elements, args.output)
        # else:  # Linux or macOS
        #     import subprocess
        #     import os
        #     script_path = os.path.join(os.path.dirname(__file__), '..', 'dump', 'extract_xyz.sh')
        #     cmd = ['/bin/bash', script_path, '-i', args.input, '-e', args.elements]
        #     if args.output:
        #         cmd.extend(['-o', args.output])
        #     subprocess.run(cmd)

if __name__ == "__main__":
    main()
