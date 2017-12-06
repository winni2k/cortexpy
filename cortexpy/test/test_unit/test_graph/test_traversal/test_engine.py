import attr

import cortexpy.graph.parser
import cortexpy.test.builder as builder
import cortexpy.test.expectation as expectation


@attr.s(slots=True)
class EngineTestDriver(object):
    graph_builder = attr.ib(attr.Factory(builder.Graph))

    def with_kmer(self, *args):
        self.graph_builder.with_kmer(*args)
        return self

    def with_kmer_size(self, n):
        self.graph_builder.with_kmer_size(n)
        return self

    def run(self):
        random_access_parser = cortexpy.graph.parser.RandomAccess(self.graph_builder.build())
        graph = cortexpy.graph.traversal.Engine(random_access_parser).traverse()
        return expectation.graph.KmerGraphExpectation(graph)


class Test(object):
    def test_on_empty_graph_returns_empty_graph(self):
        # given
        driver = EngineTestDriver().with_kmer_size(3)

        # when
        expect = driver.run()

        # then
        (expect
         .has_n_nodes(0)
         .has_n_edges(0))
