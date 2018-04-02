from cortexpy.graph import Interactor
from cortexpy.utils import get_graph_stream_iterator


def prune(argv):
    import argparse
    from .shared import get_shared_argsparse
    shared_parser = get_shared_argsparse()
    parser = argparse.ArgumentParser('cortexpy prune', parents=[shared_parser])
    parser.add_argument('-t', '--remove-tips', default=0, type=int,
                        help='Remove tips shorter than this number')
    parser.add_argument('graph', help="Input cortexpy graph.  '-' reads from stdin")
    args = parser.parse_args(argv)

    from cortexpy.logging_config import configure_logging_from_args
    configure_logging_from_args(args)

    import logging
    logger = logging.getLogger('cortexpy.prune')

    import networkx as nx
    import sys

    if args.out == '-':
        output = sys.stdout.buffer
    else:
        output = open(args.out, 'wb')

    if args.graph == '-':
        graphs = get_graph_stream_iterator(sys.stdin.buffer)
    else:
        graphs = get_graph_stream_iterator(open(args.graph, 'rb'))

    logger.info('Removing tips shorter than {} k-mers'.format(args.remove_tips))
    for graph_idx, graph in enumerate(graphs):
        logger.info('Processing graph {}'.format(graph_idx))
        Interactor(graph, colors=None).prune_tips_less_than(args.remove_tips)
        nx.write_gpickle(graph, output)
