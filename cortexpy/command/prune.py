"""
cortexpy prune

Usage:
  cortexpy prune --out <graph> <graph> [-t <n>] [-v | -s]

Options:
    -h, --help             Display this help message.
    -v, --verbose          Increase log level to debug
    -s, --silent           Decrease log level to warnings and errors

    -o, --out <graph>      Output graph.
    -t, --remove-tips <n>  Remove tips less than <n>. [default: 0]

Description:
    graph            A cortexpy graph.  '-' redirects to or from stdout.
"""
from cortexpy.graph import Interactor
from cortexpy.utils import get_graph_stream_iterator


def validate_prune(args):
    from schema import Schema, Use

    schema = Schema({
        '--remove-tips': Use(int),
        str: object,
    })
    return schema.validate(args)


def prune(argv):
    from docopt import docopt
    from cortexpy import VERSION_STRING

    args = docopt(__doc__, argv=argv, version=VERSION_STRING)
    args = validate_prune(args)

    from cortexpy.logging_config import configure_logging_from_args
    configure_logging_from_args(args)

    import logging
    logger = logging.getLogger('cortexpy.prune')

    import networkx as nx
    import sys

    if args['--out'] == '-':
        output = sys.stdout.buffer
    else:
        output = open(args['--out'], 'wb')

    if args['<graph>'] == '-':
        graphs = get_graph_stream_iterator(sys.stdin.buffer)
    else:
        graphs = get_graph_stream_iterator(open(args['<graph>'], 'rb'))

    logger.info('Removing tips shorter than {} k-mers'.format(args['--remove-tips']))
    for graph_idx, graph in enumerate(graphs):
        logger.info('Processing graph {}'.format(graph_idx))
        Interactor(graph, colors=None).prune_tips_less_than(args['--remove-tips'])
        nx.write_gpickle(graph, output)
