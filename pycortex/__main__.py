import sys
import argparse
import attr

from pycortex.cortex_graph import CortexGraphRandomAccessParser, \
    cortex_kmer_as_cortex_jdk_print_string


@attr.s(slots=True)
class CortexGraphContigRetriever(object):
    graph_parser = attr.ib(default=None)
    fh = attr.ib(default=None)

    def __attrs_post_init__(self):
        if self.graph_parser is None:
            assert self.fh is not None
            self.graph_parser = CortexGraphRandomAccessParser(self.fh)

    def get_kmers_for_contig(self, contig):
        kmer_size = self.graph_parser.header.kmer_size
        assert len(contig) >= kmer_size
        kmers = []
        for kmer_start in range(len(contig) - kmer_size + 1):
            kmer_string = contig[kmer_start:(kmer_start + kmer_size)]
            kmer = self.graph_parser.get_kmer_for_string(kmer_string)
            kmers.append((kmer, kmer_string))
        return kmers


def print_contig(args):
    if args.record:
        with open(args.graph, 'rb') as fh:
            contig_retriever = CortexGraphContigRetriever(fh=fh)
            contig_kmers = contig_retriever.get_kmers_for_contig(args.record)
            if len(contig_kmers) == 1:
                print(cortex_kmer_as_cortex_jdk_print_string(contig_kmers[0][0]))
            else:
                for kmer, kmer_string in contig_kmers:
                    print(cortex_kmer_as_cortex_jdk_print_string(kmer, alt_kmer_string=kmer_string))


def main(argv):
    parser = argparse.ArgumentParser('pycortex')

    # Copied from https://stackoverflow.com/a/26217957/528691
    subparsers = parser.add_subparsers()
    parser_print = subparsers.add_parser('print', help='Print contig')
    parser_print.add_argument('--graph')
    parser_print.add_argument('--record')

    parser_print.set_defaults(func=print_contig)

    args = parser.parse_args(argv)
    if args.func:
        args.func(args)


if __name__ == '__main__':
    main(sys.argv[1:])
