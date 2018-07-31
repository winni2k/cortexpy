def prune(argv):
    import argparse
    from .shared import get_shared_argsparse
    shared_parser = get_shared_argsparse()
    parser = argparse.ArgumentParser('cortexpy prune', parents=[shared_parser])
    parser.add_argument('-t', '--remove-tips', default=0, type=int,
                        help='Remove tips shorter than this number')
    parser.add_argument('graph', help="Input cortexpy graph.  '-' reads from stdin")
    args = parser.parse_args(argv)

    from cortexpy.logging_config import configure_logging_from_args_and_get_logger
    logger = configure_logging_from_args_and_get_logger(args, 'cortexpy.prune')

    from cortexpy.graph.interactor import Interactor
    from cortexpy.graph.parser.streaming import load_cortex_graph
    from cortexpy.graph.serializer.kmer import dump_colored_de_bruijn_graph_to_cortex

    import sys
    if args.out == '-':
        output = sys.stdout.buffer
    else:
        output = open(args.out, 'wb')

    logger.info('Loading de Bruijn graph')
    if args.graph == '-':
        graph = load_cortex_graph(sys.stdin.buffer)
    else:
        graph = load_cortex_graph(open(args.graph, 'rb'))

    logger.info('Removing tips shorter than {} k-mers'.format(args.remove_tips))
    graph = Interactor(graph).prune_tips_less_than(args.remove_tips).graph
    dump_colored_de_bruijn_graph_to_cortex(graph, output)
