import attr
from cortexpy.test.builder.unitig_finder import UnitigFinderBuilder
from cortexpy.test.expectation.unitig_graph import GraphWithUnitigExpectation


@attr.s(slots=True)
class FindUnitigsTestDriver(object):
    builder = attr.ib(attr.Factory(UnitigFinderBuilder))
    test_coverage = attr.ib(True)
    finder = attr.ib(init=False)

    @property
    def graph(self):
        return self.builder.graph

    @property
    def graph_builder(self):
        return self.builder.builder

    def with_colors(self, *colors):
        self.builder.with_colors(*colors)
        return self

    def without_test_coverage(self):
        self.builder.without_test_coverage()
        return self

    def build(self):
        self.finder = self.builder.build()
        return self.finder

    def run(self):
        self.build()
        return GraphWithUnitigExpectation(self.finder.find_unitigs())
