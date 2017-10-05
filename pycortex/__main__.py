import argparse
import sys

from pycortex.command.show import show


def main(argv):
    parser = argparse.ArgumentParser('pycortex')

    # Copied from https://stackoverflow.com/a/26217957/528691
    subparsers = parser.add_subparsers()
    show_cmd = subparsers.add_parser('show', help='Show contig')
    show_cmd.add_argument('graph')
    show_cmd.add_argument('--record')
    show_cmd.add_argument('--output-type', default='term', choices=['term', 'png'])
    show_cmd.add_argument('--output')

    show_cmd.set_defaults(func=show)

    args = parser.parse_args(argv)
    if args.func:
        args.func(args)


if __name__ == '__main__':
    main(sys.argv[1:])
