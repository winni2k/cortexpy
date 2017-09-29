import sys
import argparse
import attr

from pycortex.cortex_graph import CortexGraphRandomAccessParser, \
    cortex_kmer_as_cortex_jdk_print_string


@attr.s(slots=True)
class CortexGraphContigRetriever(object):
    cgp = attr.ib()

    def get_kmers_for_contig(self, contig):
        pass


def print_contig(args):
    if args.record:
        with open(args.graph, 'rb') as fh:
            cg = CortexGraphRandomAccessParser(fh)
            kmer = cg.get_kmer(args.record)

            print(cortex_kmer_as_cortex_jdk_print_string(kmer))


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
