"""
cortexpy view

Usage:
  cortexpy view graph <graph>
  cortexpy view contig [--to-json] <graph> <contig>
  cortexpy view traversal [--to-json --kmers --color <color>] <traversal>

Options:
    -h, --help                        Display this help message
    --kmers                           View graph as k-mers
    --color <color>                   Restrict view to single color

Subcommands:
    graph       Print all kmers in a cortex graph.
    contig      Retrieve all kmers in a contig from a cortex graph.
    traversal   View a cortexpy traversal in Python pickle format. Reads traversal from stdin if
                <traversal> is '-'.

"""
from docopt import docopt
import logging

logger = logging.getLogger('cortexpy.view')


class ArgparseError(ValueError):
    """Is raised if the args object is misconfigured"""
    pass


def view(argv):
    from cortexpy import VERSION_STRING

    args = docopt(__doc__, argv=argv, version=VERSION_STRING)
    if args['graph']:
        view_graph(args)
    elif args['contig']:
        view_contig(args)
    elif args['traversal']:
        view_traversal(args)
    else:
        raise ArgparseError


def validate_view_traversal(args):
    from schema import Schema, Use, Or

    schema = Schema({
        '--color': Or(None, Use(int)),
        str: object,
    })
    return schema.validate(args)


def view_traversal(args):
    from Bio import SeqIO
    import sys
    from cortexpy.graph import interactor
    from cortexpy.graph import serializer

    args = validate_view_traversal(args)

    import networkx as nx

    if args['<traversal>'] == '-':
        graph = nx.read_gpickle(sys.stdin.buffer)
    else:
        graph = nx.read_gpickle(args['<traversal>'])

    if args['--to-json']:
        print(serializer.Serializer(graph).to_json())
    else:
        kmer_serializer = serializer.Kmers(graph)
        if args['--kmers']:
            seq_record_generator = kmer_serializer.to_seq_records()
        else:
            seq_record_generator = interactor.Contigs(
                graph, args['--color']
            ).all_simple_paths()
        SeqIO.write(seq_record_generator, sys.stdout, 'fasta')


def view_graph(args):
    print_cortex_file(open(args['<graph>'], 'rb'))


def view_contig(args):
    from cortexpy.graph import ContigRetriever
    from cortexpy.graph.serializer import Serializer

    contig_retriever = ContigRetriever(open(args['<graph>'], 'rb'))
    if args['--to-json']:
        serializer = Serializer(contig_retriever.get_kmer_graph(args['<contig>']))
        print(serializer.to_json())
    else:
        print_contig(contig_retriever, args['<contig>'])


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
