def assemble(argv):
    import argparse
    parser = argparse.ArgumentParser(prog='cortexpy assemble', description="""
    Assemble all possible transcripts in <graph> from all k-mers in <start-sequences> and print the
    resulting transcripts as a FASTA to stdout. All specified colors are traversed and collapsed
    before output.
    """)
    parser.add_argument('graph', help='cortex graph')
    parser.add_argument('start_sequences_fasta', help='FASTA file with sequences to start from')
    parser.add_argument('--color', type=int, help='Restrict view to single color')
    parser.add_argument('--max-nodes', type=int, default=1000,
                        help='Maximum number of nodes to traverse [default: %(default)s]')
    args = parser.parse_args(argv)

    import sys
    from Bio import SeqIO
    from cortexpy.utils import kmerize_fasta
    from cortexpy.graph.interactor import Interactor, Contigs
    from cortexpy.graph.parser.random_access import RandomAccess
    from cortexpy.constants import EngineTraversalOrientation
    from cortexpy.graph.traversal.engine import Engine

    random_access = RandomAccess(open(args.graph, 'rb'))
    if args.color is None:
        colors = list(range(random_access.num_colors))
    else:
        colors = [args.color]
    traverser = Engine(
        random_access,
        traversal_colors=colors,
        orientation=EngineTraversalOrientation.both,
        max_nodes=args.max_nodes,
    )
    traverser.traverse_from_each_kmer_in_fasta(args.start_sequences_fasta)
    kmers = kmerize_fasta(args.start_sequences_fasta, traverser.ra_parser.kmer_size)
    interactor = Interactor.from_graph(traverser.graph).make_graph_nodes_consistent(
        seed_kmer_strings=kmers)

    seq_record_generator = Contigs(interactor.graph).all_simple_paths()

    SeqIO.write(seq_record_generator, sys.stdout, 'fasta')
