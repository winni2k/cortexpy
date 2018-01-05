"""
Usage: cortexpy [--help] <command> [<args>...]

Options:
   -h, --help  Display this help message

Allowed cortexpy commands are:
   view       View a cortex graph

See 'cortexpy <command>' for more information on a specific command.

"""
import sys
from docopt import docopt
from schema import SchemaError

from cortexpy import VERSION_STRING


def main_without_argv():
    return main(sys.argv[1:])


def main(argv):
    args = docopt(__doc__, argv=argv, version=VERSION_STRING, options_first=True)

    argv = [args['<command>']] + args['<args>']
    if args['<command>'] == 'view':
        import cortexpy.command.view
        try:
            cortexpy.command.view.view(argv)
        except SchemaError as e:
            print(e)
            return e
    else:
        print("'{}' is not a cortexpy command. See 'cortexpy --help'.".format(args['<command>']))
        return 1
    return 0


if __name__ == '__main__':
    exit(main_without_argv())
