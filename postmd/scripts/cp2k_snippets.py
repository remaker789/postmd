import argparse



def main_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="""
    cp2k-snippets is a convenient script that show the input and post-processing 
    snippets of LAMMPS, including MSD, friction coefficient, trajectory, etc. 
    To see the options, type "lmp-snippets -h"."""
    )

    parser.add_argument(
        "--pre", type=str, default="all", help="print the pre-processing snippets of LAMMPS."
    )
    
    parser.add_argument(
        "--post", type=str, default="all", help="print the post-processing snippets of LAMMPS."
    )
    
    return parser
    

def main():
    # print("Description\n------------")
    parser = main_parser()
    # try:
    #     import argcomplete

    #     argcomplete.autocomplete(parser)
    # except ImportError:
    #     # argcomplete not present.
    #     pass

    args = parser.parse_args()

    # try:
    #     getattr(args, "func")
    # except AttributeError:
    #     parser.print_help()
    #     sys.exit(0)
    # args.func(args)

    if args.pre.lower() == "all":
        pass
        # 动态地调用所有以'pre'开头的函数
        # for func_name, func in pre_functions.items():
        #     func()

        


if __name__ == "__main__":
    main()