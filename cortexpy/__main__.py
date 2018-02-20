"""
cortexpy

Usage: cortexpy [--help] <command> [<args>...]

Options:
   -h, --help  Display this help message

Allowed cortexpy commands are:
   view       View a cortex graph
   assemble   Assemble from a cortex graph
   traverse   Traverse a cortex graph

See 'cortexpy <command>' for more information on a specific command.

"""
import sys
import importlib
from docopt import docopt
from cortexpy import VERSION_STRING
import logging

logger = logging.getLogger('cortexpy')


def main_without_argv():
    return main(sys.argv[1:])


def main(argv):
    args = docopt(__doc__, argv=argv, version=VERSION_STRING, options_first=True)
    commands = {
        'view': 'cortexpy.command.view.view',
        'assemble': 'cortexpy.command.assemble.assemble',
        'traverse': 'cortexpy.command.traverse.traverse',
        'prune': 'cortexpy.command.prune.prune',
    }
    argv = [args['<command>']] + args['<args>']
    if args['<command>'] in commands.keys():
        package_string, method_string = commands[args['<command>']].rsplit('.', 1)
        module = importlib.import_module(package_string)
        run_command(getattr(module, method_string), argv)
    else:
        logger.error(
            "'{}' is not a cortexpy command. See 'cortexpy --help'.".format(args['<command>'])
            )
        return 1


def run_command(function, argv):
    from schema import SchemaError

    try:
        function(argv)
    except SchemaError as e:
        logger.error('Input argument error for arguments: {}'.format(argv))
        logger.error(e)
        return e
    return 0


if __name__ == '__main__':
    exit(main_without_argv())
