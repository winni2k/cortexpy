import logging
import sys

logger = logging.getLogger('cortexpy')


def main(argv=sys.argv):
    from . import __version__
    import importlib
    import argparse
    subcommands = {
        'view': 'cortexpy.command.view.view',
        'assemble': 'cortexpy.command.assemble.assemble',
        'traverse': 'cortexpy.command.traverse.traverse',
        'subgraph': 'cortexpy.command.subgraph.subgraph',
        'prune': 'cortexpy.command.prune.prune',
    }
    parser = argparse.ArgumentParser(prog='cortexpy')
    parser.add_argument('--version', action='version',
                        version='%(prog)s version {}'.format(__version__))
    parser.add_argument('subcommand', choices=sorted(subcommands.keys()),
                        help='cortexpy sub-command')
    parser.add_argument('args', nargs=argparse.REMAINDER, help='sub-command arguments')
    args = parser.parse_args(argv[1:])

    package_string, method_string = subcommands[args.subcommand].rsplit('.', 1)
    module = importlib.import_module(package_string)
    return getattr(module, method_string)(args.args)
