from enum import Enum

from cortexpy.graph import ContigRetriever
from cortexpy.graph.parser.streaming import kmer_generator_from_stream
from cortexpy.graph.serializer import Serializer


class ViewChoice(Enum):
    term = 0
    json = 1


class ArgparseError(ValueError):
    """Is raised if the args object is misconfigured"""
    pass


def add_subparser_to(subparsers):
    parser = subparsers.add_parser('view', help='Show contig')
    parser.add_argument('graph')
    parser.add_argument('--record', help='Returns record with unitigs collapsed in json format')
    parser.add_argument('--output-type', default='term',
                        choices=[v.name for v in ViewChoice])
    parser.set_defaults(func=view)


def view(args):
    with open(args.graph, 'rb') as graph_handle:
        if args.record is None and args.output_type == ViewChoice.term.name:
            print_cortex_file(graph_handle)
        else:
            contig_retriever = ContigRetriever(graph_handle)
            if args.output_type == ViewChoice.term.name:
                print_contig(contig_retriever, args.record)
            else:
                serializer = Serializer(contig_retriever.get_kmer_graph(args.record),
                                        colors=contig_retriever.colors)
                if args.output_type == ViewChoice.json.name:
                    print(serializer.to_json())
                else:
                    raise ArgparseError


def print_cortex_file(graph_handle):
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
