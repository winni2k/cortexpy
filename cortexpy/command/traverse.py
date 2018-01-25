"""
Usage:
  cortexpy traverse graph <graph> <initial_contig> [options]

Options:
    -h, --help                    Display this help message.
    --orientation <orientation>   Traversal orientation [default: both].
    -c --colors <colors>          Colors to traverse [default: 0].
           May take multiple color numbers separated by a comma (example: '0,2,3').
           The traverser will follow all colors specified.
    --max-nodes <n>               Maximum number of nodes to traverse [default: 1000].

Description:
    Traverse a cortex graph starting from each kmer in an initial_contig and return the subgraph as
    a Python pickle object.
"""
import sys
import networkx as nx
from schema import Schema, Use
from cortexpy.graph import parser as g_parser, traversal, interactor


def validate(args):
    schema = Schema({
        '--orientation': (
            lambda x: x in traversal.constants.EngineTraversalOrientation.__members__.keys()
        ),
        '--colors': Use(lambda colors: [int(color) for color in colors.split(',')]),
        '--max-nodes': Use(int),
        str: object,
    })
    return schema.validate(args)


def traverse(args):
    args = validate(args)
    graph = traversal.Engine(
        g_parser.RandomAccess(open(args['<graph>'], 'rb')),
        traversal_colors=args['--colors'],
        orientation=traversal.constants.EngineTraversalOrientation[
            args['--orientation']],
        max_nodes=args['--max-nodes'],
    ).traverse_from_each_kmer_in(args['<initial_contig>']).graph
    nx.write_gpickle(graph, sys.stdout.buffer)
