from enum import Enum

import networkx as nx

from pycortex.graph import ContigRetriever
from pycortex.graph.parser.constants import NUM_TO_LETTER
from pycortex.graph.parser.streaming import kmer_generator_from_stream
from pycortex.utils import revcomp


class ViewChoice(Enum):
    term = 0
    png = 1


def check_view_arguments(args):
    if args.output_type != ViewChoice.term.name:
        if not args.output:
            raise RuntimeError("Need to specify --output if --output-type is not 'term'")


def view(args):
    check_view_arguments(args)
    with open(args.graph, 'rb') as graph_handle:
        if args.output_type == ViewChoice.term.name:
            print_contig(graph_handle, record=args.record)
        elif args.output_type == ViewChoice.png.name:
            plot_contig(graph_handle, args.output, args.record)


def plot_contig(graph_handle, output, record):
    kmer_graph = ContigRetriever(graph_handle).get_kmer_graph(record)
    nx.draw(kmer_graph)


def print_contig(graph_handle, *, record=None):
    if record is None:
        for kmer in kmer_generator_from_stream(graph_handle):
            print(kmer_to_cortex_jdk_print_string(kmer))
    else:
        contig_kmers = ContigRetriever(graph_handle=graph_handle).get_kmers(record)
        for kmer, kmer_string in contig_kmers:
            print(kmer_to_cortex_jdk_print_string(kmer, alt_kmer_string=kmer_string))


def kmer_to_cortex_jdk_print_string(kmer, alt_kmer_string=None):
    if kmer is None:
        revcomp_kmer = revcomp(alt_kmer_string)
        if revcomp_kmer > alt_kmer_string:
            revcomp_kmer = alt_kmer_string
        return '{}: {} missing'.format(revcomp_kmer, alt_kmer_string)
    is_revcomp = bool(alt_kmer_string is not None and kmer.kmer != alt_kmer_string)

    edge_set_strings = [edge_set_as_string(edge_set, is_revcomp=is_revcomp) for edge_set in
                        kmer.edges]
    to_print = [str(kmer.kmer)]
    if alt_kmer_string is not None:
        to_print.append(': ' + alt_kmer_string)
    to_print.append(' ' + ' '.join(map(str, kmer.coverage)))
    to_print.append(' ' + ' '.join(edge_set_strings))
    return ''.join(to_print)


def edge_set_as_string(edge_set, is_revcomp=False):
    letters = []

    if is_revcomp:
        num_to_letter = list(reversed(NUM_TO_LETTER))
    else:
        num_to_letter = NUM_TO_LETTER

    for idx, edge in enumerate(edge_set):
        letter = num_to_letter[idx % 4]
        if idx < 4:
            letter = letter.lower()
        if edge:
            letters.append(letter)
        else:
            letters.append('.')

    if is_revcomp:
        incoming, outgoing = letters[:4], letters[4:]
        incoming, outgoing = list(reversed(incoming)), list(reversed(outgoing))
        letters = outgoing + incoming

    return ''.join(letters)
