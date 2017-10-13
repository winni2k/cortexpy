from enum import Enum

from pycortex.graph import ContigRetriever
from pycortex.graph.parser.constants import NUM_TO_LETTER
from pycortex.graph.parser.streaming import kmer_generator_from_stream
from pycortex.graph.serializer import Serializer
from pycortex.utils import revcomp


class ViewChoice(Enum):
    term = 0
    json = 1


class ArgparseError(ValueError):
    """Is raised if the args object is misconfigured"""
    pass


def add_subparser_to(subparsers):
    parser = subparsers.add_parser('view', help='Show contig')
    parser.add_argument('graph')
    parser.add_argument('--record')
    parser.add_argument('--output-type', default='term',
                        choices=[v.name for v in ViewChoice])
    parser.add_argument('--collapse-kmer-unitigs', action='store_true')
    parser.set_defaults(func=view)


def view(args):
    with open(args.graph, 'rb') as graph_handle:
        if args.record is None and args.output_type == ViewChoice.term.name:
            print_cortex_file(graph_handle)
        else:
            contig_retriever = ContigRetriever(graph_handle)
            if args.output_type == ViewChoice.term.name:
                print_contig(contig_retriever, args.record)
            else:
                serializer = Serializer(
                    contig_retriever.get_kmer_graph(args.record),
                    collapse_kmer_unitigs=args.collapse_kmer_unitigs)
                if args.output_type == ViewChoice.json.name:
                    print(serializer.to_json())
                else:
                    raise ArgparseError


def print_cortex_file(graph_handle):
    for kmer in kmer_generator_from_stream(graph_handle):
        print(kmer_to_cortex_jdk_print_string(kmer))


def print_contig(contig_retriever, contig):
    contig_kmers = contig_retriever.get_kmers(contig)
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
