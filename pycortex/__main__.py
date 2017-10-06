import argparse
import sys

import pycortex.command.view as view


def main(argv):
    parser = argparse.ArgumentParser('pycortex')

    # Copied from https://stackoverflow.com/a/26217957/528691
    subparsers = parser.add_subparsers()
    view_cmd = subparsers.add_parser('view', help='Show contig')
    view_cmd.add_argument('graph')
    view_cmd.add_argument('--record')
    view_cmd.add_argument('--output-type', default='term',
                          choices=[v.name for v in view.ViewChoice])
    view_cmd.add_argument('--output')

    view_cmd.set_defaults(func=view.view)

    args = parser.parse_args(argv)
    if args.func:
        args.func(args)


if __name__ == '__main__':
    main(sys.argv[1:])
