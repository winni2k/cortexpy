def traverse(argv):
    import argparse
    from cortexpy.command.shared import get_shared_argparse
    shared_parser = get_shared_argparse()

    parser = argparse.ArgumentParser(
        prog='cortexpy traverse', parents=[shared_parser],
        description="""
        Traverse all simple paths between all sources and targets of an input graph.

        Input is a cortex graph. Output is a FASTA.

        This tool also allows the creation of a JSON representation of a CORTEX graph that is consistent 
        with seed strings by using the --to-json and --seed-strings arguments.

        If a links file is supplied, then branches consistent with the links will be preferred in
        the traversal. 
        """
    )
    parser.add_argument('graph', help="cortex graph. Slurp graph from stdin is '-'.")
    parser.add_argument('--to-json', action='store_true')
    parser.add_argument('--seed-strings', nargs='*', default=[],
                        help="Strings with seed kmers from which to start contig traversal. "
                             "Multiple strings can be specified.")
    parser.add_argument('--color', type=int, help='Restrict view to single color')
    parser.add_argument('--max-paths', type=int, default=0,
                        help='Return exit status 64 if more than this '
                             'number of paths are encountered. '
                             '0 turns off this check.')
    parser.add_argument('--graph-index', type=int, default=0,
                        help='Graph index to be added to description of all output paths')
    parser.add_argument('--extra-start-kmer',
                        help='Disconnect this k-mer from incoming k-mers before '
                             'candidate transcript creation. '
                             'This argument may fail if not used together with --seed-strings.')
    parser.add_argument('--links-file', help='gzipped Mccortex-style links file for graph')
    args = parser.parse_args(argv)

    from cortexpy.logging_config import configure_logging_from_args_and_get_logger
    logger = configure_logging_from_args_and_get_logger(args, 'cortexpy.traverse')

    import sys
    import gzip
    from cortexpy.graph.interactor import Interactor
    from cortexpy.graph.serializer.serializer import Serializer
    from cortexpy.graph.parser.streaming import load_cortex_graph
    from cortexpy.links import Links
    from . import get_exit_code_yaml_path
    import yaml

    EXIT_CODES = yaml.load(open(get_exit_code_yaml_path(), 'rt'))

    if args.out == '-':
        output = sys.stdout
    else:
        output = open(args.out, 'wt')

    logger.info(f'Loading graph: %s', args.graph)
    if args.graph == '-':
        graph = load_cortex_graph(sys.stdin.buffer)
    else:
        graph = load_cortex_graph(open(args.graph, 'rb'))
    logger.info(f'Loaded {len(graph)} kmers')

    consistent_graph = None
    if args.seed_strings:
        seed_kmer_strings = strings_to_kmer_strings(args.seed_strings, graph.graph['kmer_size'])
        logger.info(
            f'Making graph consistent with {len(seed_kmer_strings)} kmers from --seed-strings')
        consistent_graph = Interactor(graph) \
            .make_graph_nodes_consistent(seed_kmer_strings) \
            .graph

    if args.to_json:
        logger.info('Writing JSON representation of graph to STDOUT')
        if consistent_graph:
            graph = consistent_graph
        print(Serializer(graph).to_json())
        return

    if not consistent_graph:
        logger.info('Making graph consistent')
        consistent_graph = Interactor.from_graph(graph) \
            .make_graph_nodes_consistent() \
            .graph

    if args.extra_start_kmer:
        if args.extra_start_kmer not in graph:
            logger.error(f'Could not find extra start kmer ({args.extra_start_kmer}) in graph')
            return 1

    links = None
    if args.links_file is not None:
        logger.info(f'Loading links file {args.links_file}')
        links = Links.from_binary_stream(gzip.open(args.links_file, 'rb'))
    seq_record_generator = Interactor(consistent_graph) \
        .all_simple_paths(args.extra_start_kmer, links=links)
    seq_record_generator = annotated_seq_records(seq_record_generator, graph_idx=args.graph_index)
    if args.max_paths > 0:
        logger.info('Exiting after element %s', args.max_paths)
        seq_record_generator = raise_after_nth_element(seq_record_generator, args.max_paths)
    logger.info('Writing seq records to %s', args.out)
    try:
        for record in seq_record_generator:
            output.write(record.format('fasta'))
    except IndexError:
        logger.error('Max paths (%s) exceeded', args.max_paths)
        return EXIT_CODES['MAX_PATH_EXCEEDED']


def annotated_seq_records(seq_record_generator, *, graph_idx):
    for rec in seq_record_generator:
        rec.id = 'g{}_p{}'.format(graph_idx, rec.id)
        yield rec


def raise_after_nth_element(iterator, n):
    for idx, val in enumerate(iterator):
        if idx == n:
            raise IndexError
        yield val


def strings_to_kmer_strings(strings, kmer_size):
    kmers = []
    for string in strings:
        assert len(string) >= kmer_size
        for pos in range(0, len(string) - kmer_size + 1):
            kmers.append(string[pos:(pos + kmer_size)])
    return kmers
