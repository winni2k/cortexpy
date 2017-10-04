import argparse
import sys

from pycortex.commands.print_ import print_contig


def main(argv):
    parser = argparse.ArgumentParser('pycortex')

    # Copied from https://stackoverflow.com/a/26217957/528691
    subparsers = parser.add_subparsers()
    parser_print = subparsers.add_parser('print', help='Print contig')
    parser_print.add_argument('--graph')
    parser_print.add_argument('--record')

    parser_print.set_defaults(func=print_contig)

    args = parser.parse_args(argv)
    if args.func:
        args.func(args)


if __name__ == '__main__':
    main(sys.argv[1:])
