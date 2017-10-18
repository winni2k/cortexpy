import attr

import pycortex.test.builder as builder
import pycortex.graph as graph
from pycortex.test.expectation.kmer import KmerNodeExpectation


@attr.s(slots=True)
class KmerGraphExpectation(object):
    graph = attr.ib()

    def with_node(self, node):
        assert node in self.graph
        return KmerNodeExpectation(self.graph.node[node])

    def has_n_nodes(self, n):
        assert len(self.graph) == n
        return self

    def has_n_edges(self, n):
        assert len(self.graph.edges) == n
        return self

    def has_nodes(self, *nodes):
        assert set(self.graph.nodes) == set(nodes)
        return self

    def has_edges(self, *edges):
        assert set(self.graph.edges) == set(edges)
        return self

    def has_edge(self, source, target):
        assert (source, target) in self.graph.edges
        return KmerGraphEdgeExpectation(self.graph[source][target])


@attr.s(slots=True)
class KmerGraphEdgeExpectation(object):
    edge = attr.ib()

    def is_missing(self):
        assert self.edge['is_missing']
        return self

    def is_not_missing(self):
        assert not self.edge['is_missing']
        return self


class TestGetKmerGraph(object):
    def test_with_no_kmer_returns_missing_kmer(self):
        # given
        graph_builder = (builder.Graph()
                         .with_kmer_size(3))
        retriever = graph.ContigRetriever(graph_builder.build())

        # when
        expect = KmerGraphExpectation(retriever.get_kmer_graph('AAA'))

        # then
        (expect.has_n_nodes(1)
         .has_n_edges(0)
         .with_node('AAA').is_missing())

    def test_with_one_kmer_returns_one_kmer(self):
        # given
        graph_builder = (builder.Graph()
                         .with_kmer_size(3))
        graph_builder.with_kmer('AAA', 1, '........')
        retriever = graph.ContigRetriever(graph_builder.build())

        # when
        kmer_graph = retriever.get_kmer_graph('AAA')

        # then
        assert len(kmer_graph.edges) == 0
        assert list(kmer_graph) == ['AAA']

    def test_with_two_linked_kmers_returns_two_kmers(self):
        # given
        graph_builder = (builder.Graph()
                         .with_kmer_size(3))
        graph_builder.with_kmer('AAA', 1, '.....C..')
        graph_builder.with_kmer('AAC', 1, 'a.......')
        retriever = graph.ContigRetriever(graph_builder.build())

        # when
        expect = KmerGraphExpectation(retriever.get_kmer_graph('AAA'))

        # then
        expect.has_nodes('AAA', 'AAC').has_n_edges(1)
        expect.has_edge('AAA', 'AAC').is_not_missing()

    def test_with_three_linked_kmers_and_two_colors_returns_three_kmers(self):
        # given
        graph_builder = (builder.Graph()
                         .with_kmer_size(3)
                         .with_num_colors(2))
        graph_builder.with_kmer('AAA', [1, 1], ['.....C..', '.......T'])
        graph_builder.with_kmer('AAC', [1, 0], ['a.......', '........'])
        graph_builder.with_kmer('AAT', [0, 1], ['........', 'a.......'])
        retriever = graph.ContigRetriever(graph_builder.build())

        # when
        kmer_graph = retriever.get_kmer_graph('AAAC')

        # then
        assert set(kmer_graph) == {'AAA', 'AAC', 'AAT'}
        assert set(kmer_graph.edges) == {('AAA', 'AAC'), ('AAA', 'AAT')}

    def test_with_two_unlinked_kmers_request_for_three_returns_three_linked_kmers(self):
        # given
        graph_builder = (builder.Graph()
                         .with_kmer_size(3))
        graph_builder.with_kmer('AAA', 1, '........')
        graph_builder.with_kmer('ACC', 1, '........')
        retriever = graph.ContigRetriever(graph_builder.build())

        # when
        expect = KmerGraphExpectation(retriever.get_kmer_graph('AAACC'))

        # then
        expect.has_nodes('AAA', 'AAC', 'ACC')
        expect.has_n_edges(2)
        expect.has_edge('AAA', 'AAC').is_missing()
        expect.has_edge('AAC', 'ACC').is_missing()

    def test_with_two_neighboring_unlinked_kmers_returns_two_kmers_linked_by_missing_edge(self):
        # given
        graph_builder = (builder.Graph()
                         .with_kmer_size(3))
        graph_builder.with_kmer('AAA', 1, '........')
        graph_builder.with_kmer('AAC', 1, '........')
        retriever = graph.ContigRetriever(graph_builder.build())

        # when
        expect = KmerGraphExpectation(retriever.get_kmer_graph('AAAC'))

        # then
        expect.has_nodes('AAA', 'AAC')
        expect.has_n_edges(1)
        expect.has_edge('AAA', 'AAC').is_missing()


class TestGetKmerGraphRevcomp(object):
    def test_with_one_kmer_returns_kmer(self):
        # given
        graph_builder = (builder.Graph()
                         .with_kmer_size(3))
        graph_builder.with_kmer('AAA', 1, '........')
        retriever = graph.ContigRetriever(graph_builder.build())

        # when
        kmer_graph = retriever.get_kmer_graph('TTT')

        # then
        assert len(kmer_graph.edges) == 0
        assert list(kmer_graph.nodes) == ['TTT']

    def test_with_two_linked_kmers_returns_two_kmers(self):
        # given
        graph_builder = (builder.Graph()
                         .with_kmer_size(3))
        graph_builder.with_kmer('AAA', 1, '.....C..')
        graph_builder.with_kmer('AAC', 1, 'a.......')
        retriever = graph.ContigRetriever(graph_builder.build())

        # when
        kmer_graph = retriever.get_kmer_graph('GTTT')

        # then
        assert set(kmer_graph.nodes) == {'GTT', 'TTT'}
        assert list(kmer_graph.edges) == [('GTT', 'TTT')]

    def test_with_three_linked_kmers_and_two_colors_returns_three_kmers(self):
        # given
        graph_builder = (builder.Graph()
                         .with_kmer_size(3)
                         .with_num_colors(2))
        graph_builder.with_kmer('AAA', [1, 1], ['.....C..', '.......T'])
        graph_builder.with_kmer('AAC', [1, 0], ['a.......', '........'])
        graph_builder.with_kmer('AAT', [0, 1], ['........', 'a.......'])
        retriever = graph.ContigRetriever(graph_builder.build())

        # when
        kmer_graph = retriever.get_kmer_graph('GTTT')

        # then
        assert set(kmer_graph.nodes) == {'TTT', 'GTT', 'ATT'}
        assert set(kmer_graph.edges) == {('GTT', 'TTT'), ('ATT', 'TTT')}

        assert list(kmer_graph.nodes['TTT']['kmer'].coverage) == [1, 1]
        assert list(kmer_graph.nodes['GTT']['kmer'].coverage) == [1, 0]
