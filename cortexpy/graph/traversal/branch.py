import attr
import networkx as nx

from cortexpy.graph.serializer import EdgeTraversalOrientation


@attr.s(slots=True)
class Branch(object):
    ra_parser = attr.ib()
    config = attr.ib()

    def traverse_from(self, kmer_string, *, direction=EdgeTraversalOrientation.original):
        graph = nx.MultiDiGraph()
        return graph
        # while True:
        #     ra_parser[kmer_string]
