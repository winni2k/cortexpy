"""
cortexpy prune

Usage:
  cortexpy prune --out <graph> <graph>

Options:
    -h, --help          Display this help message.
    -o, --out <graph>   Output graph.

Description:
    graph            A cortexpy graph.  '-' redirects to or from stdout.
"""


def prune(argv):
    from docopt import docopt
    from cortexpy import VERSION_STRING

    args = docopt(__doc__, argv=argv, version=VERSION_STRING)

    import networkx as nx
    import sys

    if args['--out'] == '-':
        output = sys.stdout.buffer
    else:
        output = open(args['--out'], 'wb')

    if args['<graph>'] == '-':
        graph = nx.read_gpickle(sys.stdin.buffer)
    else:
        graph = nx.read_gpickle(args['<graph>'])

    nx.write_gpickle(graph, output)
