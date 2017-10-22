import argparse
import sys

import pycortex.command.view as view


def main(argv):
    parser = argparse.ArgumentParser('pycortex')

    # Copied from https://stackoverflow.com/a/26217957/528691
    subparsers = parser.add_subparsers()
    view.add_subparser_to(subparsers)

    args = parser.parse_args(argv)
    if getattr(args, 'func', False):
        args.func(args)
    else:
        parser.print_help()


if __name__ == '__main__':
    main(sys.argv[1:])
