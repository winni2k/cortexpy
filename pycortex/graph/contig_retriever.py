import attr
import networkx as nx

from pycortex.graph import parser as parser


@attr.s(slots=True)
class ContigRetriever(object):
    graph_handle = attr.ib()

    def get_kmers(self, contig):
        graph_parser = parser.RandomAccess(self.graph_handle)
        kmer_size = graph_parser.header.kmer_size
        assert len(contig) >= kmer_size
        kmers = []
        for kmer_start in range(len(contig) - kmer_size + 1):
            kmer_string = contig[kmer_start:(kmer_start + kmer_size)]
            try:
                kmer = graph_parser.get_kmer_for_string(kmer_string)
            except KeyError:
                kmer = None
            kmers.append((kmer, kmer_string))
        return kmers

    def get_kmer_graph(self, contig):
        g = nx.DiGraph()
        for kmer, kmer_string in self.get_kmers(contig):
            g.add_node(kmer_string)
        return g
