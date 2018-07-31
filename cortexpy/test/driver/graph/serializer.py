import attr

from cortexpy.graph.interactor import Interactor
from cortexpy.graph.contig_retriever import ContigRetriever
from cortexpy.graph.parser.random_access import RandomAccess
from cortexpy.graph.serializer import unitig
from cortexpy.graph.traversal.engine import Engine
from cortexpy.test import builder as builder
from cortexpy.test.expectation.kmer import CollapsedKmerUnitgGraphExpectation


@attr.s(slots=True)
class SerializerTestDriver(object):
    graph_builder = attr.ib(attr.Factory(builder.Graph))
    contig_to_retrieve = attr.ib(None)
    retriever = attr.ib(None)
    traverse = attr.ib(False)
    retrieve = attr.ib(False)
    traversal_start_kmer = attr.ib(None)
    traversal_colors = attr.ib((0,))

    def with_kmer_size(self, n):
        self.graph_builder.with_kmer_size(n)
        return self

    def with_kmer(self, *args):
        self.graph_builder.with_kmer(*args)
        return self

    def traverse_with_start_kmer_and_colors(self, start_kmer, *colors):
        self.traverse = True
        self.traversal_start_kmer = start_kmer
        self.traversal_colors = colors
        return self

    def retrieve_contig(self, contig):
        self.retrieve = True
        self.contig_to_retrieve = contig
        return self

    def run(self):
        if self.retrieve:
            self.retriever = ContigRetriever(self.graph_builder.build())
            return self.retriever.get_kmer_graph(self.contig_to_retrieve)
        elif self.traverse:
            traverser = Engine(RandomAccess(self.graph_builder.build()),
                               traversal_colors=self.traversal_colors)
            graph = traverser.traverse_from(self.traversal_start_kmer).graph
            return Interactor(graph) \
                .make_graph_nodes_consistent([self.traversal_start_kmer]) \
                .graph
        else:
            raise Exception("Need to load a command")


@attr.s(slots=True)
class CollapseKmerUnitigsTestDriver(object):
    serializer_driver = attr.ib(attr.Factory(SerializerTestDriver))

    def __getattr__(self, name):
        serializer_method = getattr(self.serializer_driver, name)

        def method(*args):
            serializer_method(*args)
            return self

        return method

    def run(self):
        kmer_graph = self.serializer_driver.run()
        collapser = unitig \
            .UnitigCollapser(kmer_graph) \
            .collapse_kmer_unitigs()
        return CollapsedKmerUnitgGraphExpectation(collapser.unitig_graph)
