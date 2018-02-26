"""
cortexpy traverse

Usage:
  cortexpy traverse <graph> <initial-contigs> --out <file> [options]

Options:
    -h, --help                    Display this help message.
    --out <file>                  Output graph. '-' prints to stdout
    --orientation <orientation>   Traversal orientation [default: both].
    -c --colors <colors>          Colors to traverse [default: 0].
           May take multiple color numbers separated by a comma (example: '0,2,3').
           The traverser will follow all colors specified.
    --max-nodes <n>               Maximum number of nodes to traverse [default: 1000].
    --initial-fasta               Treat <initial-contigs> as fasta
    --subgraphs                  Emit traversal as sequence of networkx subgraphs

Description:
    Traverse a cortex graph starting from each k-mer in an initial_contig and return the subgraph as
    a Python pickle object.
"""
import sys
from docopt import docopt
from schema import Schema, Use, Optional
import logging

logger = logging.getLogger('cortexpy.traverse')


def validate(args):
    from cortexpy.graph import traversal
    schema = Schema({
        Optional('--orientation'): (
            lambda x: x in traversal.constants.EngineTraversalOrientation.__members__.keys()
        ),
        Optional('--colors'): Use(lambda colors: [int(color) for color in colors.split(',')]),
        Optional('--max-nodes'): Use(int),
        str: object,
    })
    return schema.validate(args)


def traverse(argv):
    from cortexpy import VERSION_STRING

    args = docopt(__doc__, argv=argv, version=VERSION_STRING)
    args = validate(args)
    if args['--out'] == '-':
        output = sys.stdout.buffer
    else:
        output = open(args['--out'], 'wb')

    import networkx as nx
    from cortexpy.graph import parser as g_parser, traversal
    engine = traversal.Engine(
        g_parser.RandomAccess(open(args['<graph>'], 'rb')),
        traversal_colors=args['--colors'],
        orientation=traversal.constants.EngineTraversalOrientation[
            args['--orientation']],
        max_nodes=args['--max-nodes'],
    )
    if args['--initial-fasta']:
        engine.traverse_from_each_kmer_in_fasta(args['<initial-contigs>'])
    else:
        engine.traverse_from_each_kmer_in(args['<initial-contigs>'])

    output_graph = engine.graph
    if args['--subgraphs'] and len(output_graph) > 0:
        for subgraph in nx.weakly_connected_component_subgraphs(output_graph):
            nx.write_gpickle(subgraph, output)
    else:
        nx.write_gpickle(output_graph, output)
