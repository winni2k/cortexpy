import networkx as nx

from cortexpy.graph import interactor
from cortexpy.test.builder.graph.colored_de_bruijn import ColoredDeBruijnGraphBuilder


class TestMultiDiGraph(object):
    def test_prunes_three_tips_of_length_1(self):
        # given
        graph = nx.MultiDiGraph()
        graph.add_path([0, 1, 2])
        graph.add_path([0, 1, 3])

        # when
        graph = interactor.Interactor(graph, colors=None).prune_tips_less_than(2).graph

        # then
        assert 1 == len(graph)
        assert {1} == set(graph.nodes)

    def test_prunes_one_tip_of_length_1(self):
        # given
        graph = nx.MultiDiGraph()
        graph.add_path([0, 1, 2, 3])
        graph.add_path([0, 1, 4, 5])

        # when
        graph = interactor.Interactor(graph, colors=None).prune_tips_less_than(2).graph

        # then
        assert set(range(1, 6)) == set(graph.nodes)
        assert 5 == len(graph)


class TestColoredDeBruijnGraph(object):
    def test_prunes_three_tips_of_length_1(self):
        # given
        b = ColoredDeBruijnGraphBuilder()
        b.with_colors(0, 1)
        kmers = ['AAA', 'AAC', 'ACC', 'ACG']
        for kmer in kmers:
            b.add_node(kmer)
        b.add_edge(kmers[0], kmers[1], 0)
        b.add_edge(kmers[0], kmers[1], 1)
        b.add_edge(kmers[1], kmers[2], 0)
        b.add_edge(kmers[1], kmers[3], 1)
        graph = b.build()

        # when
        graph = interactor.Interactor(graph, colors=None).prune_tips_less_than(2).graph

        # then
        assert 1 == len(graph)
        assert {'AAC'} == set(graph.nodes)

    def test_prunes_one_tip_of_length_1(self):
        # given
        b = ColoredDeBruijnGraphBuilder()
        b.with_colors(0, 1)
        kmers = ['AAA', 'AAC', 'ACC', 'CCC', 'ACG', 'CGC']
        for kmer in kmers:
            b.add_node(kmer)
        b.add_edge(kmers[0], kmers[1], 0)
        b.add_edge(kmers[0], kmers[1], 1)
        b.add_edge(kmers[1], kmers[2], 0)
        b.add_edge(kmers[2], kmers[3], 0)
        b.add_edge(kmers[1], kmers[4], 1)
        b.add_edge(kmers[4], kmers[5], 1)
        graph = b.build()

        # when
        graph = interactor.Interactor(graph, colors=None).prune_tips_less_than(2).graph

        # then
        assert set(kmers[1:]) == set(graph.nodes)
        assert 5 == len(graph)
