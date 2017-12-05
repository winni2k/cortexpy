import attr

from cortexpy.graph.serializer import UnitigFinder
from cortexpy.test.builder.graph.networkx import NetworkxGraphBuilder


@attr.s(slots=True)
class UnitigFinderBuilder(object):
    builder = attr.ib(attr.Factory(NetworkxGraphBuilder))
    test_coverage = attr.ib(True)

    def __attrs_post_init__(self):
        self.builder.with_colors(0)

    @property
    def graph(self):
        return self.builder.graph

    def without_test_coverage(self):
        self.test_coverage = False
        return self

    def with_colors(self, *colors):
        self.builder.with_colors(*colors)
        return self

    def build(self):
        graph = self.builder.build()
        return UnitigFinder(graph, colors=list(self.builder.colors),
                            test_coverage=self.test_coverage)
