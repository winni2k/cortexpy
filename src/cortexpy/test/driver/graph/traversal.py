import attr

from cortexpy.constants import EngineTraversalOrientation
from cortexpy.graph.parser.random_access import RandomAccess
from cortexpy.graph.traversal.engine import Engine
from cortexpy.test import builder as builder, expectation as expectation


@attr.s(slots=True)
class EngineTestDriver(object):
    graph_builder = attr.ib(attr.Factory(builder.Graph))
    start_kmer_string = attr.ib(None)
    start_string = attr.ib(None)
    max_nodes = attr.ib(1000)
    traversal_orientation = attr.ib(EngineTraversalOrientation.original)
    traverser = attr.ib(None)
    traversal_colors = attr.ib((0,))
    ra_constructor = attr.ib(RandomAccess)

    def with_kmer(self, *args):
        self.graph_builder.with_kmer(*args)
        return self

    def with_kmer_size(self, n):
        self.graph_builder.with_kmer_size(n)
        return self

    def with_num_colors(self, n):
        self.graph_builder.with_num_colors(n)
        return self

    def with_start_kmer_string(self, start_kmer_string):
        self.start_kmer_string = start_kmer_string
        return self

    def with_start_string(self, start_string):
        self.start_string = start_string
        return self

    def with_max_nodes(self, max_nodes):
        self.max_nodes = max_nodes
        return self

    def with_traversal_orientation(self, orientation):
        self.traversal_orientation = EngineTraversalOrientation[orientation]
        return self

    def with_traversal_colors(self, *colors):
        self.traversal_colors = colors
        return self

    def with_ra_constructor(self, constructor):
        self.ra_constructor = constructor
        return self

    def run(self):
        random_access_parser = self.ra_constructor(self.graph_builder.build())
        self.traverser = Engine(random_access_parser,
                                traversal_colors=self.traversal_colors,
                                max_nodes=self.max_nodes,
                                orientation=self.traversal_orientation)
        assert (self.start_string is None) != (self.start_kmer_string is None)
        if self.start_string:
            self.traverser.traverse_from_each_kmer_in(self.start_string)
        else:
            self.traverser.traverse_from(self.start_kmer_string)
        return expectation.graph.KmerGraphExpectation(self.traverser.graph)
