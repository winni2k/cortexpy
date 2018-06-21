import logging

logger = logging.getLogger('cortexpy.view')


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
    parser = argparse.ArgumentParser(prog='cortexpy view traversal')
    parser.add_argument('traversal', help="cortexpy traversal in Python pickle format."
                                          " Read traversal from stdin traversal is '-'.")
    parser.add_argument('--to-json', action='store_true')
    parser.add_argument('--kmers', action='store_true')
    parser.add_argument('--seed-strings', nargs='*', default=[],
                        help="Strings with seed kmers from which to start contig traversal. "
                             "Multiple strings can be specified.")
    parser.add_argument('--color', type=int, help='Restrict view to single color')
    args = parser.parse_args(argv)

    from Bio import SeqIO
    import sys
    from cortexpy.graph import interactor
    from cortexpy.graph import serializer
    from cortexpy.graph.serializer import unitig
    from cortexpy.graph.parser.streaming import load_de_bruijn_graph

    if args.traversal == '-':
        graph = load_de_bruijn_graph(sys.stdin.buffer)
    else:
        with open(args.traversal, 'rb') as fh:
            graph = load_de_bruijn_graph(fh)

    if args.to_json:
        print(unitig.Serializer(graph).to_json())
    else:
        if args.kmers:
            seq_record_generator = serializer.KmerGraph(graph).to_seq_records()
        else:
            seed_kmer_strings = strings_to_kmer_strings(args.seed_strings, graph.graph['kmer_size'])
            consistent_graph = interactor.Interactor(graph,
                                                     colors=None).make_graph_nodes_consistent(
                seed_kmer_strings).graph
            seq_record_generator = interactor.Contigs(consistent_graph,
                                                      args.color).all_simple_paths()
        seq_record_generator = annotated_seq_records(seq_record_generator, graph_idx="x")
        SeqIO.write(seq_record_generator, sys.stdout, 'fasta')


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
    from cortexpy.graph.serializer import unitig

    contig_retriever = ContigRetriever(open(args.graph, 'rb'))
    if args.to_json:
        serializer = unitig.Serializer(contig_retriever.get_kmer_graph(args.contig))
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
