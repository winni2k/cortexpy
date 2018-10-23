from cortexpy.graph.interactor import Interactor
from cortexpy.test.builder.graph.cortex import (
    CortexGraphBuilder,
    get_cortex_builder,
    LinksBuilder,
)


class TestCortex(object):
    def test_emits_one_single_color_unitig(self):
        # given
        b = CortexGraphBuilder()
        b.with_colors(0)
        b.add_edge('AAA', 'AAT', color=0)
        b.make_consistent('AAA')
        graph = b.build()

        # when
        paths = list(Interactor(graph).all_simple_paths())

        # then
        assert ['AAAT'] == [str(p.seq) for p in paths]

    def test_only_follows_one_color_with_color_specified(self):
        # given
        b = CortexGraphBuilder()
        b.with_colors(0, 1)
        b.add_edge('AAA', 'AAT', color=0)
        b.add_edge('AAT', 'ATA', color=1)
        b.make_consistent('AAA')
        graph = b.build()

        # when
        paths = list(Interactor(graph).keep_color(0).all_simple_paths())

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
        paths = list(Interactor(graph).all_simple_paths())

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
        paths = list(Interactor(graph).all_simple_paths())

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
        paths = list(Interactor(cdb).all_simple_paths())

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
        paths = list(Interactor(cdb).all_simple_paths())

        # then
        assert ['AAGCC', 'AAGCG'] == sorted([str(p.seq) for p in paths])


class TestLinks:
    def test_with_link_for_y_graph_emits_one_path(self):
        # given
        b = CortexGraphBuilder()
        b.with_kmer_size(3)
        b.add_path('AAA', 'AAC')
        b.add_path('AAA', 'AAT')
        b.make_consistent('AAA')
        cdb = b.build()

        links = LinksBuilder() \
            .with_link_for_kmer('F 1 1 C', 'AAA') \
            .build()

        # when
        paths = list(Interactor(cdb).all_simple_paths(links=links))

        # then
        assert ['AAAC'] == [str(p.seq) for p in paths]

    def test_bubble_and_y_with_two_links_returns_two_transcripts(self):
        # given
        links = LinksBuilder() \
            .with_link_for_kmer('F 2 1 CT', 'AAA') \
            .with_link_for_kmer('F 1 1 A', 'CCC') \
            .build()

        b = CortexGraphBuilder()
        b.with_kmer_size(3)
        b.add_path('AAA', 'AAC', 'ACC', 'CCC', 'CCA')
        b.add_path('AAA', 'AAG', 'AGC', 'GCC', 'CCC', 'CCT')
        b.make_consistent('AAA')
        cdb = b.build()

        # when
        paths = list(Interactor(cdb).all_simple_paths(links=links))

        # then
        assert ['AAACCCA', 'AAACCCT'] == sorted([str(p.seq) for p in paths])

    def test_bubble_and_y_with_two_links_returns_three_transcripts(self):
        # given
        links = LinksBuilder() \
            .with_link_for_kmer('F 2 1 CT', 'AAA') \
            .with_link_for_kmer('F 1 1 G', 'AAA') \
            .build()

        b = CortexGraphBuilder()
        b.with_kmer_size(3)
        b.add_path('AAA', 'AAC', 'ACC', 'CCC', 'CCA')
        b.add_path('AAA', 'AAG', 'AGC', 'GCC', 'CCC', 'CCT')
        b.make_consistent('AAA')
        cdb = b.build()

        # when
        paths = list(Interactor(cdb).all_simple_paths(links=links))

        # then
        assert ['AAACCCT', 'AAAGCCCA', 'AAAGCCCT'] == sorted([str(p.seq) for p in paths])
