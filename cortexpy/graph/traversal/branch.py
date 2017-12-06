import attr
import networkx as nx

from cortexpy.graph.serializer import EdgeTraversalOrientation


@attr.s(slots=True)
class Branch(object):
    ra_parser = attr.ib()
    config = attr.ib()

    def traverse_from(self, kmer_string, *, direction=EdgeTraversalOrientation.original):
        graph = nx.MultiDiGraph()
        try:
            kmer = self.ra_parser.get_kmer_for_string(kmer_string)
            graph.add_node(kmer_string, kmer=kmer)
        except KeyError:
            pass
        return graph
