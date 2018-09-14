def prune(argv):
    import argparse
    from .shared import get_shared_argparse
    shared_parser = get_shared_argparse()
    parser = argparse.ArgumentParser('cortexpy prune', parents=[shared_parser])
    parser.add_argument('-t', '--remove-tips', required=True, type=int,
                        help='Remove tips shorter than this number')
    parser.add_argument('graph', help="Input cortexpy graph.  '-' reads from stdin")
    args = parser.parse_args(argv)

    from cortexpy.logging_config import configure_logging_from_args_and_get_logger
    logger = configure_logging_from_args_and_get_logger(args, 'cortexpy.prune')

    if args.remove_tips < 2:
        logger.error('--remove-tips (%s) needs to be greater than 1', args.remove_tips)
        return 1

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

    logger.info(f'Loaded {len(graph)} kmers')

    graph = Interactor(graph).prune_tips_less_than(args.remove_tips).graph
    dump_colored_de_bruijn_graph_to_cortex(graph, output)
