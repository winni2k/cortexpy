import cortexpy.constants


def traverse(argv):
    import argparse
    from cortexpy.graph import traversal
    from .shared import get_shared_argsparse
    shared_parser = get_shared_argsparse()
    parser = argparse.ArgumentParser(
        'cortexpy traverse', parents=[shared_parser],
        description="""
        Traverse a cortex graph starting from each k-mer in an initial_contig and return the
        subgraph as a Cortex graph."""
    )
    parser.add_argument('--graphs', nargs='+',
                        required=True,
                        help="Input cortexpy graphs."
                             "  Multiple graphs can be specified and are joined on-the-fly.")
    parser.add_argument('initial_contig', help="Initial contig from which to start traversal")
    parser.add_argument('--orientation',
                        type=cortexpy.constants.EngineTraversalOrientation,
                        choices=[o.name for o in cortexpy.constants.EngineTraversalOrientation],
                        default=cortexpy.constants.EngineTraversalOrientation.both,
                        help='Traversal orientation')
    parser.add_argument('-c', '--colors',
                        nargs='+',
                        type=int,
                        help="""Colors to traverse.  May take multiple color numbers separated by
                        a space.  The traverser will follow all colors
                        specified.  Will follow all colors if not specified.
                        """, default=None)
    parser.add_argument('--initial-fasta', action='store_true',
                        help='Treat initial_contig as a file in FASTA format')
    parser.add_argument('--max-nodes', type=int, default=None,
                        help='Maximum number of nodes to traverse (int).'
                             '  Die without output if max nodes is exceeded')
    parser.add_argument('--logging-interval', type=int, default=90,
                        help='Logging interval.  [default: %(default)s]')
    parser.add_argument('--cache-size', type=int, default=0, help='Number of kmers to cache')
    parser.add_argument('--binary-search-cache-size', type=int, default=0,
                        help='Number of kmers to cache for binary search')
    args = parser.parse_args(argv)

    from cortexpy.logging_config import configure_logging_from_args_and_get_logger
    logger = configure_logging_from_args_and_get_logger(args, 'cortexpy.traverse')

    import sys
    from cortexpy.graph.serializer.kmer import dump_colored_de_bruijn_graph_to_cortex
    from cortexpy.graph.parser.random_access_collection import RandomAccessCollection
    from cortexpy.constants import EngineTraversalOrientation
    from cortexpy.graph.traversal.engine import Engine
    from contextlib import ExitStack
    with ExitStack() as stack:
        if args.out == '-':
            output = sys.stdout.buffer
        else:
            output = stack.enter_context(open(args.out, 'wb'))

        from cortexpy.graph.parser.random_access import RandomAccess
        if len(args.graphs) == 1:
            ra_parser = RandomAccess(
                stack.enter_context(open(args.graphs[0], 'rb')),
                kmer_cache_size=args.cache_size
            )
        else:
            ra_parser = RandomAccessCollection(
                [RandomAccess(stack.enter_context(open(graph_path, 'rb')),
                                       kmer_cache_size=args.cache_size)
                 for graph_path in args.graphs])
        engine = Engine(
            ra_parser,
            orientation=EngineTraversalOrientation[args.orientation.name],
            max_nodes=args.max_nodes,
            logging_interval=args.logging_interval
        )

        if args.colors is not None:
            engine.traversal_colors = args.colors
        else:
            engine.traversal_colors = tuple(list(range(engine.ra_parser.num_colors)))
        logger.info('Traversing colors: ' + ','.join([str(c) for c in engine.traversal_colors]))

        if args.initial_fasta:
            engine.traverse_from_each_kmer_in_fasta(args.initial_contig)
        else:
            engine.traverse_from_each_kmer_in(args.initial_contig)

        dump_colored_de_bruijn_graph_to_cortex(engine.graph, output)
