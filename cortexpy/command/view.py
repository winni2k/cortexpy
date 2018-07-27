class ArgparseError(ValueError):
    """Is raised if the args object is misconfigured"""
    pass


def view(argv):
    import argparse
    parser = argparse.ArgumentParser(prog='cortexpy view')
    subcommands = {
        'graph': view_graph,
        'contig': view_contig,
        'traversal': view_traversal,
    }
    parser.add_argument('subcommand', choices=sorted(subcommands.keys()),
                        help='cortexpy view sub-command')
    parser.add_argument('args', nargs=argparse.REMAINDER, help='sub-command arguments')
    args = parser.parse_args(argv)
    return subcommands[args.subcommand](args.args)


def view_traversal(argv):
    import argparse
    from .shared import get_shared_argsparse
    shared_parser = get_shared_argsparse()

    parser = argparse.ArgumentParser(prog='cortexpy view traversal', parents=[shared_parser])
    parser.add_argument('traversal', help="cortexpy traversal in Python pickle format."
                                          " Read traversal from stdin traversal is '-'.")
    parser.add_argument('--to-json', action='store_true')
    parser.add_argument('--kmers', action='store_true')
    parser.add_argument('--seed-strings', nargs='*', default=[],
                        help="Strings with seed kmers from which to start contig traversal. "
                             "Multiple strings can be specified.")
    parser.add_argument('--color', type=int, help='Restrict view to single color')
    args = parser.parse_args(argv)

    from cortexpy.logging_config import configure_logging_from_args_and_get_logger
    logger = configure_logging_from_args_and_get_logger(args, 'cortexpy.view')

    assert not args.kmers or not args.to_json
    from Bio import SeqIO
    import sys
    from contextlib import ExitStack
    from cortexpy.graph.interactor import Interactor, Contigs
    from cortexpy.graph.serializer.serializer import Serializer
    from cortexpy.graph.serializer import KmerGraph
    from cortexpy.graph.parser.streaming import load_cortex_graph

    with ExitStack() as stack:
        if args.out == '-':
            output = sys.stdout
        else:
            output = stack.enter_context(open(args.out, 'wt'))

        logger.info(f'Loading graph: {args.traversal}')
        if args.traversal == '-':
            graph = load_cortex_graph(sys.stdin.buffer)
        else:
            with open(args.traversal, 'rb') as fh:
                graph = load_cortex_graph(fh)
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
            logger.info('Writing JSON representation of graph Unitigs to STDOUT')
            if consistent_graph:
                graph = consistent_graph
            print(Serializer(graph).to_json())
            return

        if args.kmers:
            seq_record_generator = KmerGraph(graph).to_seq_records()
        else:
            if not consistent_graph:
                logger.info('Making graph consistent')
                consistent_graph = Interactor.from_graph(graph) \
                    .make_graph_nodes_consistent() \
                    .graph
            seq_record_generator = Contigs(consistent_graph, args.color).all_simple_paths()
        seq_record_generator = annotated_seq_records(seq_record_generator, graph_idx="x")
        logger.info('Writing seq records to {args.out}')
        SeqIO.write(seq_record_generator, output, 'fasta')


def view_graph(argv):
    import argparse
    parser = argparse.ArgumentParser(prog='cortexpy view graph')
    parser.add_argument('graph', help="cortex graph")
    args = parser.parse_args(argv)

    print_cortex_file(open(args.graph, 'rb'))


def view_contig(argv):
    import argparse
    parser = argparse.ArgumentParser(prog='cortexpy view contig')
    parser.add_argument('graph', help="cortex graph")
    parser.add_argument('contig', help='contig to explore inside graph')
    parser.add_argument('--to-json', action='store_true')
    args = parser.parse_args(argv)

    from cortexpy.graph import ContigRetriever
    from cortexpy.graph.serializer.serializer import Serializer

    contig_retriever = ContigRetriever(open(args.graph, 'rb'))
    if args.to_json:
        serializer = Serializer(contig_retriever.get_kmer_graph(args.contig))
        print(serializer.to_json())
    else:
        print_contig(contig_retriever, args.contig)


def print_cortex_file(graph_handle):
    from cortexpy.graph.parser.streaming import kmer_generator_from_stream

    for kmer in kmer_generator_from_stream(graph_handle):
        print(kmer_to_cortex_jdk_print_string(kmer))


def print_contig(contig_retriever, contig):
    contig_kmers = contig_retriever.get_kmers(contig)
    for kmer, kmer_string in contig_kmers:
        print(kmer_to_cortex_jdk_print_string(kmer, alt_kmer_string=kmer_string))


def kmer_to_cortex_jdk_print_string(kmer, alt_kmer_string=None):
    is_revcomp = bool(alt_kmer_string is not None and kmer.kmer != alt_kmer_string)

    edge_set_strings = (edge_set.to_str(as_revcomp=is_revcomp) for edge_set in kmer.edges)
    to_print = [str(kmer.kmer)]
    if alt_kmer_string is not None:
        to_print.append(': ' + alt_kmer_string)
    to_print.append(' ' + ' '.join(map(str, kmer.coverage)))
    to_print.append(' ' + ' '.join(edge_set_strings))
    return ''.join(to_print)


def annotated_seq_records(seq_record_generator, *, graph_idx):
    for rec in seq_record_generator:
        rec.id = 'g{}_p{}'.format(graph_idx, rec.id)
        yield rec


def strings_to_kmer_strings(strings, kmer_size):
    kmers = []
    for string in strings:
        assert len(string) >= kmer_size
        for pos in range(0, len(string) - kmer_size + 1):
            kmers.append(string[pos:(pos + kmer_size)])
    return kmers
