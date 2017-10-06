import attr
import networkx as nx

from pycortex.graph import parser as parser
from pycortex.utils import revcomp


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
        kmer_graph = nx.DiGraph()
        for kmer, kmer_string in self.get_kmers(contig):
            kmer_graph.add_node(kmer_string)
            is_revcomp = kmer.kmer != kmer_string
            incoming_kmers = set()
            outgoing_kmers = set()
            for edge_set in kmer.edges:
                for incoming_kmer in edge_set.get_incoming_kmers(kmer.kmer):
                    incoming_kmers.add(incoming_kmer)
                for outgoing_kmer in edge_set.get_outgoing_kmers(kmer.kmer):
                    outgoing_kmers.add(outgoing_kmer)
            if is_revcomp:
                incoming_kmers, outgoing_kmers = \
                    {revcomp(kmer) for kmer in outgoing_kmers}, \
                    {revcomp(kmer) for kmer in incoming_kmers}
            for incoming_kmer in incoming_kmers:
                kmer_graph.add_edge(incoming_kmer, kmer_string)
            for outgoing_kmer in outgoing_kmers:
                kmer_graph.add_edge(kmer_string, outgoing_kmer)
        return kmer_graph
