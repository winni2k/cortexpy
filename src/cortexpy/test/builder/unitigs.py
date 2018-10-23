import attr
import networkx as nx


@attr.s(slots=True)
class UnitigBuilder:
    graph = attr.ib(attr.Factory(nx.DiGraph))

    def __getattr__(self, item):
        return getattr(self.graph, item)

    def with_kmer_size(self, k):
        self.graph.graph['kmer_size'] = k
        return self

    def add_node(self, node_id, unitig_string):
        self.graph.add_node(node_id, unitig=unitig_string)
        return self

    def build(self):
        return self.graph
