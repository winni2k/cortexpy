import networkx as nx

from cortexpy.graph import interactor, Interactor
from cortexpy.test.builder.graph.colored_de_bruijn import (
    ColoredDeBruijnGraphBuilder,
    get_cdb_builder,
)


class TestMultiDigraph(object):
    def test_only_follows_one_color_with_color_specified(self):
        # given
        graph = nx.MultiDiGraph()
        graph.add_edge('AAA', 'AAT', key=0)
        graph.add_edge('AAT', 'ATA', key=1)

        # when
        paths = list(interactor.Contigs(graph, color=0).all_simple_paths())

        # then
        assert ['AAAT'] == [str(p.seq) for p in paths]

    def test_follows_two_colors_with_no_color_specified(self):
        # given
        graph = nx.MultiDiGraph()
        graph.add_edge('AAA', 'AAT', key=0)
        graph.add_edge('AAT', 'ATA', key=1)

        # when
        paths = list(interactor.Contigs(graph).all_simple_paths())

        # then
        assert {'AAATA'} == set([str(p.seq) for p in paths])

    def test_follows_three_colors_with_no_color_specified(self):
        # given
        graph = nx.MultiDiGraph()
        graph.add_edge('AAA', 'AAT', key=0)
        graph.add_edge('AAT', 'ATA', key=1)
        graph.add_edge('AAT', 'ATC', key=2)

        # when
        paths = list(interactor.Contigs(graph).all_simple_paths())

        # then
        assert {'AAATA', 'AAATC'} == set([str(p.seq) for p in paths])


class TestDBG(object):

    def test_in_y_graph_finds_two_paths(self):
        # given
        b = ColoredDeBruijnGraphBuilder()
        b.add_path('CAA', 'AAA')
        b.add_path('TAA', 'AAA')
        cdb = b.build()

        # when
        paths = list(interactor.Contigs(cdb).all_simple_paths())

        # then
        assert {'CAAA', 'TAAA'} == set([str(p.seq) for p in paths])

    def test_in_y_graph_finds_two_paths_of_revcomp(self):
        # given
        b = get_cdb_builder()
        b.with_kmer('CGC 0 .......T')
        b.with_kmer('AGC 0 a....CG.')
        b.with_kmer('AAG 0 .....C..')
        b.with_kmer('GCC 0 a.......')
        cdb = b.build()

        cdb = Interactor(cdb, colors=None).make_graph_nodes_consistent(['AAG']).graph

        # when
        paths = list(interactor.Contigs(cdb).all_simple_paths())

        # then
        assert {'AAGCG', 'AAGCC'} == set([str(p.seq) for p in paths])
