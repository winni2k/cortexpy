"""
cortexpy assemble

Usage:
  cortexpy assemble <graph> <start-sequences-fasta> [options]

Options:
    -h, --help                    Display this help message.
    -c --colors <colors>          Colors to traverse [default: 0].
                           May take multiple color numbers separated by a comma (example: '0,2,3').
                           The traverser will follow all colors specified.
    --max-nodes <n>               Maximum number of nodes to traverse [default: 1000].

Description:
    Assemble all possible transcripts in <graph> from all k-mers in <start-sequences> and print the
    resulting transcripts as a FASTA to stdout. All specified colors are traversed and collapsed
    before output.

    graph            A cortex graph
    start-sequences  A FASTA file with sequences from which to start assembly.
"""


def validate(args):
    from schema import Schema, Use

    schema = Schema({
        '--colors': Use(lambda colors: [int(color) for color in colors.split(',')]),
        '--max-nodes': Use(int),
        str: object,
    })
    return schema.validate(args)


def assemble(argv):
    import sys
    from docopt import docopt
    from Bio import SeqIO
    from cortexpy import VERSION_STRING
    from cortexpy.graph import traversal, parser, interactor

    args = docopt(__doc__, argv=argv, version=VERSION_STRING)
    args = validate(args)

    traverser = traversal.Engine(
        parser.RandomAccess(open(args['<graph>'], 'rb')),
        traversal_colors=args['--colors'],
        orientation=traversal.constants.EngineTraversalOrientation.both,
        max_nodes=args['--max-nodes'],
    )
    traverser.traverse_from_each_kmer_in_fasta(args['<start-sequences-fasta>'])
    graph = traverser.graph

    seq_record_generator = interactor.Contigs(graph, int(args['--colors'][0])).all_simple_paths()

    SeqIO.write(seq_record_generator, sys.stdout, 'fasta')
