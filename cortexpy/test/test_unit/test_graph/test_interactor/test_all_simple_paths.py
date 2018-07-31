from cortexpy.graph.interactor import Interactor, Contigs
from cortexpy.test.builder.graph.cortex import (
    CortexGraphBuilder,
    get_cortex_builder,
)


class TestCortex(object):
    def test_only_follows_one_color_with_color_specified(self):
        # given
        b = CortexGraphBuilder()
        b.with_colors(0, 1)
        b.add_edge('AAA', 'AAT', color=0)
        b.add_edge('AAT', 'ATA', color=1)
        b.make_consistent('AAA')
        graph = b.build()

        # when
        paths = list(Contigs(graph, color=0).all_simple_paths())

        # then
        assert ['AAAT'] == [str(p.seq) for p in paths]

    def test_follows_two_colors_with_no_color_specified(self):
        # given
        b = CortexGraphBuilder()
        b.with_colors(0, 1)
        b.add_edge('AAA', 'AAT', color=0)
        b.add_edge('AAT', 'ATA', color=1)
        b.make_consistent('AAA')
        graph = b.build()

        # when
        paths = list(Contigs(graph).all_simple_paths())

        # then
        assert {'AAATA'} == set([str(p.seq) for p in paths])

    def test_follows_three_colors_with_no_color_specified(self):
        # given
        b = CortexGraphBuilder()
        b.with_colors(0, 1, 2)
        b.add_edge('AAA', 'AAT', color=0)
        b.add_edge('AAT', 'ATA', color=1)
        b.add_edge('AAT', 'ATC', color=2)
        b.make_consistent('AAA')
        graph = b.build()

        # when
        paths = list(Contigs(graph).all_simple_paths())

        # then
        assert {'AAATA', 'AAATC'} == set([str(p.seq) for p in paths])

    def test_in_y_graph_finds_two_paths(self):
        # given
        b = CortexGraphBuilder()
        b.add_path('CAA', 'AAA')
        b.add_path('TAA', 'AAA')
        b.make_consistent('AAA')
        cdb = b.build()

        # when
        paths = list(Contigs(cdb).all_simple_paths())

        # then
        assert {'CAAA', 'TAAA'} == set([str(p.seq) for p in paths])

    def test_in_y_graph_finds_two_paths_of_revcomp(self):
        # given
        b = get_cortex_builder()
        b.with_kmer('CGC 1 .......T')
        b.with_kmer('AGC 1 a....CG.')
        b.with_kmer('AAG 1 .....C..')
        b.with_kmer('GCC 1 a.......')
        cdb = b.build()

        cdb = Interactor(cdb).make_graph_nodes_consistent(['AAG']).graph

        # when
        paths = list(Contigs(cdb).all_simple_paths())

        # then
        assert {'AAGCG', 'AAGCC'} == set([str(p.seq) for p in paths])
