import argparse
import sys

from pycortex.command.view import view


def main(argv):
    parser = argparse.ArgumentParser('pycortex')

    # Copied from https://stackoverflow.com/a/26217957/528691
    subparsers = parser.add_subparsers()
    view_cmd = subparsers.add_parser('view', help='Show contig')
    view_cmd.add_argument('graph')
    view_cmd.add_argument('--record')
    view_cmd.add_argument('--output-type', default='term', choices=['term', 'png'])
    view_cmd.add_argument('--output')

    view_cmd.set_defaults(func=view)

    args = parser.parse_args(argv)
    if args.func:
        args.func(args)


if __name__ == '__main__':
    main(sys.argv[1:])
