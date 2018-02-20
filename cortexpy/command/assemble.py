"""
cortexpy assemble

Usage:
  cortexpy assemble <graph> <start-sequences-fasta> [options]

Options:
    -h, --help               Display this help message.
    -c --color <color>       Restrict assembly to single color (int)
    --max-nodes <n>          Maximum number of nodes to traverse [default: 1000].

Description:
    Assemble all possible transcripts in <graph> from all k-mers in <start-sequences> and print the
    resulting transcripts as a FASTA to stdout. All specified colors are traversed and collapsed
    before output.

    graph            A cortex graph
    start-sequences  A FASTA file with sequences from which to start assembly.
"""


def validate(args):
    from schema import Schema, Use, Optional, Or

    schema = Schema({
        Optional('--color'): Or(None, Use(int)),
        Optional('--max-nodes'): Use(int),
        str: object,
    })
    return schema.validate(args)


def assemble(argv):
    from docopt import docopt
    from cortexpy import VERSION_STRING
    args = docopt(__doc__, argv=argv, version=VERSION_STRING)
    args = validate(args)

    import sys
    from Bio import SeqIO
    from cortexpy.graph import traversal, parser, interactor

    random_access = parser.RandomAccess(open(args['<graph>'], 'rb'))
    if args['--color'] is None:
        colors = list(range(random_access.header.num_colors))
    else:
        colors = [args['--color']]
    traverser = traversal.Engine(
        random_access,
        traversal_colors=colors,
        orientation=traversal.constants.EngineTraversalOrientation.both,
        max_nodes=args['--max-nodes'],
    )
    traverser.traverse_from_each_kmer_in_fasta(args['<start-sequences-fasta>'])

    seq_record_generator = interactor.Contigs(traverser.graph).all_simple_paths()

    SeqIO.write(seq_record_generator, sys.stdout, 'fasta')
