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
    parser.add_argument('--color', type=int, help='Restrict view to single color')
    parser.add_argument('--subgraphs', default=None,
                        help='Do not merge subgraphs and store under this prefix')
    args = parser.parse_args(argv)

    from Bio import SeqIO
    import sys
    from cortexpy.graph import interactor
    from cortexpy.graph import serializer
    import networkx as nx
    from cortexpy.utils import get_graph_stream_iterator

    if args.traversal == '-':
        graphs = get_graph_stream_iterator(sys.stdin.buffer)
    else:
        graphs = get_graph_stream_iterator(open(args.traversal, 'rb'))

    if args.to_json:
        if args.subgraphs:
            for graph_idx, graph in enumerate(graphs):
                with open('{}_{}.json'.format(args.subgraphs, graph_idx), 'w') as fh:
                    fh.write(serializer.Serializer(graph).to_json())
        else:
            graph = nx.compose_all(list(graphs))
            print(serializer.Serializer(graph).to_json())
    else:
        for graph_idx, graph in enumerate(graphs):
            if args.kmers:
                seq_record_generator = serializer.Kmers(graph).to_seq_records()
            else:
                seq_record_generator = interactor.Contigs(graph, args.color).all_simple_paths()
            seq_record_generator = annotated_seq_records(seq_record_generator, graph_idx=graph_idx)
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
    from cortexpy.graph.serializer import Serializer

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
