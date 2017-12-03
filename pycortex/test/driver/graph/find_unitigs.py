import attr
from pycortex.graph.serializer import UnitigCollapser
from pycortex.test.builder.graph.networkx import add_kmers_to_graph, NetworkxGraphBuilder
from pycortex.test.expectation.unitig_graph import GraphWithUnitigExpectation


@attr.s(slots=True)
class FindUnitigsTestDriver(object):
    builder = attr.ib(attr.Factory(NetworkxGraphBuilder))
    test_coverage = attr.ib(True)
    graph = attr.ib(init=False)

    def __attrs_post_init__(self):
        self.builder.with_colors(0)
        self.graph = self.builder.graph

    def without_test_coverage(self):
        self.test_coverage = False
        return self

    def build(self):
        self.graph = self.builder.build()
        self.graph = add_kmers_to_graph(self.graph)
        return self

    def run(self):
        self.build()
        collapser = UnitigCollapser(self.graph, colors=list(self.builder.colors))
        collapser.find_unitigs()
        return GraphWithUnitigExpectation(collapser.unitig_graph)
