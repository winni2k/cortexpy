import attr
import networkx as nx

from cortexpy.graph.serializer import EdgeTraversalOrientation


@attr.s(slots=True)
class Branch(object):
    ra_parser = attr.ib()
    config = attr.ib()

    def traverse_from(self, kmer_string, *, orientation=EdgeTraversalOrientation.original):
        traversal_color = self.config.traversal_color
        graph = nx.MultiDiGraph()
        kmer = self.ra_parser.get_kmer_for_string(kmer_string)
        graph.add_node(kmer_string, kmer=kmer)
        while True:
            if kmer.kmer != kmer_string:
                kmer_orientation = EdgeTraversalOrientation.other(orientation)
            else:
                kmer_orientation = orientation
            edge_set = kmer.edges[traversal_color].oriented(kmer_orientation)
            if edge_set.num_neighbor() != 1:
                break
            kmer_string = edge_set.neighbor_kmer_strings(kmer_string)[0]
            kmer = self.ra_parser.get_kmer_for_string(kmer_string)
            graph.add_node(kmer_string, kmer=kmer)
        return graph
