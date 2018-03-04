"""
cortexpy traverse

Usage:
  cortexpy traverse <graph> <initial-contigs> --out <file> [options] [-v | -s]

Options:
    -h, --help                    Display this help message.
    -v, --verbose                 Increase log level to debug
    -s, --silent                  Decrease log level to warnings and errors
    --out <file>                  Output graph. '-' prints to stdout
    --orientation <orientation>   Traversal orientation [default: both].
    -c --colors <colors>          Colors to traverse.
           May take multiple color numbers separated by a comma (example: '0,2,3').
           The traverser will follow all colors specified.
           Will follow all colors if not specified.
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


log_level_conversion = {
    0: 'NOTSET',
    10: 'DEBUG',
    20: 'INFO',
    30: 'WARNING',
    40: 'ERROR',
    50: 'CRITICAL',
}


def traverse(argv):
    from cortexpy import VERSION_STRING
    from docopt import docopt

    args = docopt(__doc__, argv=argv, version=VERSION_STRING)
    args = validate(args)

    import logging
    if args['--verbose']:
        log_level = logging.DEBUG
    elif args['--silent']:
        log_level = logging.WARNING
    else:
        log_level = logging.INFO
    logging.basicConfig(level=log_level)
    logger = logging.getLogger('cortexpy.traverse')
    logger.info('Log level is {}'.format(log_level_conversion[logger.getEffectiveLevel()]))

    import sys
    from contextlib import ExitStack
    with ExitStack() as stack:
        if args['--out'] == '-':
            output = sys.stdout.buffer
        else:
            output = stack.enter_context(open(args['--out'], 'wb'))

        import networkx as nx
        from cortexpy.graph import parser as g_parser, traversal
        input_graph = stack.enter_context(open(args['<graph>'], 'rb'))
        engine = traversal.Engine(
            g_parser.RandomAccess(input_graph),
            orientation=traversal.constants.EngineTraversalOrientation[
                args['--orientation']],
            max_nodes=args['--max-nodes'],
            logging_interval=args['--logging-interval']
        )

        if args['--colors'] is not None:
            engine.traversal_colors = args['--colors']
        else:
            engine.traversal_colors = tuple(list(range(engine.ra_parser.header.num_colors)))
        logger.info('Traversing colors: ' + ','.join([str(c) for c in engine.traversal_colors]))

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