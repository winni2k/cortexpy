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
        assert contig is not None
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
        previous_kmer_tuple = None
        for kmer, kmer_string in self.get_kmers(contig):
            kmer_graph.add_node(kmer_string, kmer=kmer)
            if previous_kmer_tuple is not None:
                if previous_kmer_tuple[0] is None or kmer is None:
                    kmer_graph.add_edge(previous_kmer_tuple[1], kmer_string, is_missing=True)
            previous_kmer_tuple = (kmer, kmer_string)
            if kmer is None:
                continue
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
                kmer_graph.add_edge(incoming_kmer, kmer_string, is_missing=False)
            for outgoing_kmer in outgoing_kmers:
                kmer_graph.add_edge(kmer_string, outgoing_kmer, is_missing=False)
        for kmer_node in kmer_graph:
            kmer_graph.nodes[kmer_node]['repr'] = kmer_node
        return kmer_graph
