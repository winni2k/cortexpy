"""
Usage:
  cortexpy view graph <graph> [--help]
  cortexpy view contig <graph> <contig> [--help] [--output-format=<format>]
  cortexpy view traversal <graph> <initial_contig> [options] [--output-format=<format>]

Options:
    -h, --help                    Display this help message.
    --orientation <orientation>   Traversal orientation [default: both].
    -c --colors <colors>          Colors to traverse [default: 0].
           May take multiple color numbers separated by a comma (example: '0,2,3').
           The traverser will follow all colors specified.
    --max-nodes <n>               Maximum number of nodes to traverse [default: 1000].
    --output-type <type>          Output type [default: kmers].
    --output-format <format>      Output format [default: term].

Subcommands:
    graph       Print all kmers in a cortex graph
    contig      Retrieve all kmers in a contig from a cortex graph
    traversal   Traverse a cortex graph starting from each kmer in an initial_contig
                and return the subgraph
"""
from docopt import docopt
from enum import Enum
import logging
import attr

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


def view(argv):
    from cortexpy import VERSION_STRING

    args = docopt(__doc__, argv=argv, version=VERSION_STRING)
    if args['graph']:
        view_graph(args)
    elif args['contig']:
        view_contig(args)
    elif args['traversal']:
        ViewTraversal.run(args)
    else:
        raise ArgparseError


@attr.s(slots=True)
class ViewTraversal(object):
    subparsers = attr.ib()
    parser = attr.ib(init=False)

    def __attrs_post_init__(self):
        self.parser = self.subparsers.add_parser('traversal',
                                                 help='Display a Cortex graph traversal')

    @classmethod
    def validate(cls, args):
        from cortexpy.graph import traversal
        from schema import Schema, Use

        schema = Schema({
            '--orientation': (
                lambda x: x in traversal.constants.EngineTraversalOrientation.__members__.keys()
            ),
            '--colors': Use(lambda colors: [int(color) for color in colors.split(',')]),
            '--max-nodes': Use(int),
            '--output-type': lambda x: x in ViewTraversalOutputType.__members__.keys(),
            '--output-format': lambda x: x in ViewTraversalOutputFormat.__members__.keys(),
            str: object,
        })
        return schema.validate(args)

    @classmethod
    def run(cls, args):
        from Bio import SeqIO
        import sys
        from cortexpy.graph import parser as g_parser, traversal, interactor
        from cortexpy.graph import serializer

        args = cls.validate(args)
        graph = traversal.Engine(
            g_parser.RandomAccess(get_file_handle(args['<graph>'])),
            traversal_colors=args['--colors'],
            orientation=traversal.constants.EngineTraversalOrientation[
                args['--orientation']],
            max_nodes=args['--max-nodes'],
        ).traverse_from_each_kmer_in(args['<initial_contig>']).graph
        if args['--output-format'] == ViewTraversalOutputFormat.json.name:
            print(serializer.Serializer(graph).to_json())
        elif args['--output-format'] == ViewTraversalOutputFormat.fasta.name:
            kmer_serializer = serializer.Kmers(graph)
            if args['--output-type'] == ViewTraversalOutputType.contigs.name:
                seq_record_generator = interactor.Contigs(
                    graph, int(args['--colors'][0])
                ).all_simple_paths()
            else:
                seq_record_generator = kmer_serializer.to_seq_records()
            SeqIO.write(seq_record_generator, sys.stdout, 'fasta')


def get_file_handle(graph):
    return open(graph, 'rb')


def view_graph(args):
    print_cortex_file(get_file_handle(args['<graph>']))


def view_contig(args):
    from cortexpy.graph import ContigRetriever
    from cortexpy.graph.serializer import Serializer

    contig_retriever = ContigRetriever(get_file_handle(args['<graph>']))
    if args['--output-format'] == ViewContigOutputFormat.term.name:
        print_contig(contig_retriever, args['<contig>'])
    elif args['--output-format'] == ViewContigOutputFormat.json.name:
        serializer = Serializer(contig_retriever.get_kmer_graph(args['<contig>']))
        print(serializer.to_json())
    else:
        raise ArgparseError


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
