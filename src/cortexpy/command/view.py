def traverse_deprecated(*args, **kwargs):
    import warnings

    from cortexpy.command.traverse import traverse
    warnings.warn('"cortexpy view traverse" is deprecated. Please use "cortexpy traverse" instead.',
                  DeprecationWarning,
                  stacklevel=2)
    return traverse(*args, **kwargs)


def view(argv):
    import argparse
    parser = argparse.ArgumentParser(prog='cortexpy view')
    subcommands = {
        'graph': view_graph,
        'contig': view_contig,
        'traversal': traverse_deprecated,
    }
    parser.add_argument('subcommand', choices=sorted(subcommands.keys()),
                        help='cortexpy view sub-command')
    parser.add_argument('args', nargs=argparse.REMAINDER, help='sub-command arguments')
    args = parser.parse_args(argv)
    return subcommands[args.subcommand](args.args)


def view_graph(argv):
    import argparse

    parser = argparse.ArgumentParser(prog='cortexpy view graph')
    parser.add_argument('graph', help="cortex graph")
    parser.add_argument('--kmers', action='store_true')
    args = parser.parse_args(argv)

    import sys
    with open(args.graph, 'rb') as fh:
        if args.kmers:
            write_kmer_strings(fh, sys.stdout)
        else:
            print_cortex_file(fh)


def view_contig(argv):
    import argparse
    parser = argparse.ArgumentParser(prog='cortexpy view contig')
    parser.add_argument('graph', help="cortex graph")
    parser.add_argument('contig', help='contig to explore inside graph')
    parser.add_argument('--to-json', action='store_true')
    args = parser.parse_args(argv)

    from cortexpy.graph.contig_retriever import ContigRetriever
    from cortexpy.graph.serializer.serializer import Serializer

    contig_retriever = ContigRetriever(open(args.graph, 'rb'))
    if args.to_json:
        serializer = Serializer(contig_retriever.get_kmer_graph(args.contig))
        print(serializer.to_json())
    else:
        print_contig(contig_retriever, args.contig)


def write_kmer_strings(input, output):
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    from cortexpy.graph.parser.streaming import kmer_string_generator_from_stream

    seqs = (SeqRecord(Seq(s), id=str(i), description="") for i, s in
            enumerate(kmer_string_generator_from_stream(input)))
    SeqIO.write(seqs, output, "fasta")


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
