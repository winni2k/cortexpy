from enum import Enum

import logging

import attr
import sys

from Bio import SeqIO

from cortexpy.graph import ContigRetriever, parser as g_parser, traversal, interactor
from cortexpy.graph.parser.streaming import kmer_generator_from_stream
from cortexpy.graph.serializer import Serializer

logger = logging.getLogger('cortexpy.view')


class ViewContigOutputFormat(Enum):
    term = 0
    json = 1


class ViewTraversalOutputFormat(Enum):
    json = 1
    fasta = 2


class ViewTraversalOutputType(Enum):
    kmers = 0
    contigs = 1


class ArgparseError(ValueError):
    """Is raised if the args object is misconfigured"""
    pass


def add_subparser_to(subparsers):
    parent_parser = subparsers.add_parser('view', help='Display a subset of a Cortex graph')
    subparsers = parent_parser.add_subparsers()
    child_parsers = []

    graph_parser = subparsers.add_parser('graph', help='List graph records')
    contig_parser = subparsers.add_parser('contig', help='Display a contig from a Cortex graph')
    traversal = ViewTraversal(subparsers)

    child_parsers.append(graph_parser)
    child_parsers.append(contig_parser)
    child_parsers.append(traversal.parser)

    for parser in child_parsers:
        parser.add_argument('graph', help='Input Cortex graph')

    graph_parser.set_defaults(func=view_graph)
    contig_parser.set_defaults(func=view_contig)

    contig_parser.add_argument('contig', help='Contig to retrieve records for')
    contig_parser.add_argument('--output-format', default='term',
                               choices=[v.name for v in ViewContigOutputFormat])

    traversal.add_arguments()


@attr.s(slots=True)
class ViewTraversal(object):
    subparsers = attr.ib()
    parser = attr.ib(init=False)

    def __attrs_post_init__(self):
        self.parser = self.subparsers.add_parser('traversal',
                                                 help='Display a Cortex graph traversal')

    def add_arguments(self):
        parser = self.parser
        parser.set_defaults(func=view_traversal)
        parser.add_argument('initial_contig',
                            help='Start traversal from this kmer (contig to be added)')
        parser.add_argument('--orientation', default='both',
                            choices=[o.name for o in traversal.EngineTraversalOrientation])
        parser.add_argument('--color', default=0, type=int)
        parser.add_argument('--max-nodes', default=1000, type=int)
        parser.add_argument('--output-type', default=ViewTraversalOutputType.kmers.name)
        parser.add_argument('--output-format', default=ViewTraversalOutputFormat.fasta.name,
                            choices=[v.name for v in ViewTraversalOutputFormat])


def get_file_handle(graph):
    return open(graph, 'rb')


def view_graph(args):
    print_cortex_file(get_file_handle(args.graph))


def view_contig(args):
    contig_retriever = ContigRetriever(get_file_handle(args.graph))
    if args.output_format == ViewContigOutputFormat.term.name:
        print_contig(contig_retriever, args.contig)
    else:
        serializer = Serializer(contig_retriever.get_kmer_graph(args.contig))
        if args.output_format == ViewContigOutputFormat.json.name:
            print(serializer.to_json())
        else:
            raise ArgparseError


def view_traversal(args):
    traverser = traversal.Engine(
        g_parser.RandomAccess(get_file_handle(args.graph)),
        color=args.color,
        orientation=traversal.EngineTraversalOrientation[args.orientation],
        max_nodes=args.max_nodes,
    )
    graph = traverser.traverse_from_each_kmer_in(args.initial_contig)
    serializer = Serializer(graph)
    if args.output_format == ViewTraversalOutputFormat.json.name:
        serializer.annotate_graph_edges = True
        print(serializer.to_json())
    elif args.output_format == ViewTraversalOutputFormat.fasta.name:
        serializer.annotate_graph_edges = False
        serializer.collapse_unitigs = False
        if args.output_type == ViewTraversalOutputType.contigs.name:
            seq_record_generator = interactor.Contigs(graph, args.color).all_simple_paths()
        else:
            seq_record_generator = serializer.to_seq_records()
        SeqIO.write(seq_record_generator, sys.stdout, 'fasta')


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
