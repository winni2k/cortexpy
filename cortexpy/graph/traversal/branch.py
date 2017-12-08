import attr
import networkx as nx

from cortexpy.graph.serializer import EdgeTraversalOrientation


@attr.s(slots=True)
class Branch(object):
    ra_parser = attr.ib()
    traversal_color = attr.ib(0)

    def traverse_from(self, kmer_string, *, orientation=EdgeTraversalOrientation.original):
        graph = nx.MultiDiGraph()
        kmer = self.ra_parser.get_kmer_for_string(kmer_string)
        graph.add_node(kmer_string, kmer=kmer)
        while True:
            edge_set = kmer.edges[self.traversal_color].oriented(orientation)
            neighbor_edge_set = edge_set
            if kmer.kmer != kmer_string:
                neighbor_edge_set = kmer.edges[self.traversal_color].oriented(
                    EdgeTraversalOrientation.other(orientation))
            if neighbor_edge_set.num_neighbor() != 1:
                break
            prev_kmer_string = kmer_string
            kmer_string = edge_set.neighbor_kmer_strings(prev_kmer_string)[0]
            kmer = self.ra_parser.get_kmer_for_string(kmer_string)
            graph.add_node(kmer_string, kmer=kmer)
            graph.add_edge(prev_kmer_string, kmer_string, key=self.traversal_color)
        return graph
