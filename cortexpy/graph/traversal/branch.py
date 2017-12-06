import attr
import networkx as nx

from cortexpy.graph.serializer import EdgeTraversalOrientation


@attr.s(slots=True)
class Branch(object):
    ra_parser = attr.ib()
    config = attr.ib()

    def traverse_from(self, kmer_string, *, direction=EdgeTraversalOrientation.original):
        graph = nx.MultiDiGraph()
        kmer = self.ra_parser.get_kmer_for_string(kmer_string)
        graph.add_node(kmer_string, kmer=kmer)

        if direction == EdgeTraversalOrientation.original:
            traversal_direction_kmer_strings = 'get_outgoing_kmers'
        else:
            traversal_direction_kmer_strings = 'get_incoming_kmers'

        edges = set()
        for edge_set in kmer.edges:
            edges |= set(getattr(edge_set, traversal_direction_kmer_strings)(kmer_string))
        while len(edges) == 1:
            kmer_string = edges.pop()
            kmer = self.ra_parser.get_kmer_for_string(kmer_string)
            graph.add_node(kmer_string, kmer=kmer)
        return graph
