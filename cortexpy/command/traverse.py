"""
cortexpy traverse

Usage:
  cortexpy traverse --graph=<graph>... --out=<file> [options] [-v | -s] <initial-contigs>

Options:
    -h, --help                    Display this help message.
    -v, --verbose                 Increase log level to debug
    -s, --silent                  Decrease log level to warnings and errors

    --out <file>                  Output graph. '-' prints to stdout
    --graph <graph>               Input mccortex graph.  Multiple graphs can be specified and are
                                  joined on-the-fly.
    --orientation <orientation>   Traversal orientation [default: both].
    -c --colors <colors>          Colors to traverse.  May take multiple color numbers separated by
                                  a comma (example: '0,2,3').  The traverser will follow all colors
                                  specified.  Will follow all colors if not specified.
    --initial-fasta               Treat <initial-contigs> as fasta
    --subgraphs                   Emit traversal as sequence of networkx subgraphs
    --max-nodes <n>               Maximum number of nodes to traverse (int).
                                  Die without output if max nodes is exceeded.
    --logging-interval <seconds>  Logging interval [default: 90]

Description:
    Traverse a cortex graph starting from each k-mer in an initial_contig and return the subgraph as
    a Python pickle object.
"""


def validate(args):
    """Validate inputs for traverse command"""
    from schema import Schema, Use, Or
    from cortexpy.graph import traversal

    schema = Schema({
        '--orientation': (
            lambda x: x in traversal.constants.EngineTraversalOrientation.__members__.keys()
        ),
        '--colors': Or(None,
                       Use(lambda colors: tuple([int(color) for color in colors.split(',')]))),
        '--max-nodes': Or(None, Use(int)),
        '--logging-interval': Use(int),
        str: object,
    })
    return schema.validate(args)


def traverse(argv):
    from cortexpy import VERSION_STRING
    from docopt import docopt

    args = docopt(__doc__, argv=argv, version=VERSION_STRING)
    args = validate(args)

    from cortexpy.logging_config import configure_logging_from_args
    configure_logging_from_args(args)

    import logging
    logger = logging.getLogger('cortexpy.traverse')

    import sys
    from contextlib import ExitStack
    with ExitStack() as stack:
        if args['--out'] == '-':
            output = sys.stdout.buffer
        else:
            output = stack.enter_context(open(args['--out'], 'wb'))

        import networkx as nx
        from cortexpy.graph import parser as g_parser, traversal
        if len(args['--graph']) == 1:
            ra_parser = g_parser.RandomAccess(stack.enter_context(open(args['--graph'][0], 'rb')))
        else:
            ra_parser = g_parser.RandomAccessCollection(
                [g_parser.RandomAccess(stack.enter_context(open(graph_path, 'rb'))) for graph_path
                 in args['--graph']])
        engine = traversal.Engine(
            ra_parser,
            orientation=traversal.constants.EngineTraversalOrientation[
                args['--orientation']],
            max_nodes=args['--max-nodes'],
            logging_interval=args['--logging-interval']
        )

        if args['--colors'] is not None:
            engine.traversal_colors = args['--colors']
        else:
            engine.traversal_colors = tuple(list(range(engine.ra_parser.num_colors)))
        logger.info('Traversing colors: ' + ','.join([str(c) for c in engine.traversal_colors]))

        if args['--initial-fasta']:
            engine.traverse_from_each_kmer_in_fasta(args['<initial-contigs>'])
        else:
            engine.traverse_from_each_kmer_in(args['<initial-contigs>'])

        output_graph = engine.graph
        if args['--subgraphs'] and len(output_graph) > 0:
            for subgraph in sorted(nx.weakly_connected_component_subgraphs(output_graph),
                                   key=lambda g: len(g), reverse=True):
                nx.write_gpickle(subgraph, output)
        else:
            nx.write_gpickle(output_graph, output)
