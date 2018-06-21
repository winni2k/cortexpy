import attr

from cortexpy.graph.serializer.unitig import UnitigFinder
from cortexpy.test.builder.graph import colored_de_bruijn


@attr.s(slots=True)
class UnitigFinderBuilder(object):
    builder = attr.ib(attr.Factory(colored_de_bruijn.ColoredDeBruijnGraphBuilder))
    test_coverage = attr.ib(True)

    def __attrs_post_init__(self):
        self.builder.with_colors(0)

    @property
    def graph(self):
        return self.builder

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
