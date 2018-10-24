import copy

import pytest

from cortexpy.links import LinkWalker, UnitigLinkWalker, LinkedGraphTraverser
from cortexpy.test.builder.graph.cortex import LinksBuilder
from cortexpy.test.builder.unitigs import UnitigBuilder


class TestLinkGroup_GetLinkJunctions:
    @pytest.mark.parametrize('is_lexlo', [True, False])
    def test_with_single_link_returns_one_base(self, is_lexlo):
        # given
        b = LinksBuilder()
        b.with_link_for_kmer('F 3 1 ACC', 'AAA')

        # when
        links = b.build()

        # then
        junctions = list(links.body['AAA'].get_link_junctions_in_kmer_orientation(is_lexlo))
        if is_lexlo:
            assert ['ACC'] == junctions
        else:
            assert [] == junctions


class TestWalker:
    @pytest.mark.parametrize('kmer', ('AAA', 'TTT'))
    def test_with_single_link_returns_one_base(self, kmer):
        # given
        b = LinksBuilder()
        b.with_link_for_kmer('F 3 1 ACC', 'AAA')
        links = b.build()

        # when
        walker = LinkWalker.from_links(links).load_kmer(kmer)

        # then
        if kmer == 'AAA':
            assert 1 == walker.n_junctions
            assert ['A'] == list(walker.next_junction_bases())
            walker.choose_branch('A')
            assert 1 == walker.n_junctions
            walker.choose_branch('C')
            assert 1 == walker.n_junctions
            walker.choose_branch('C')
            assert 0 == walker.n_junctions
        else:
            assert 0 == walker.n_junctions
            assert [] == list(walker.next_junction_bases())

    @pytest.mark.parametrize('kmer', ('AAA', 'TTT'))
    def test_with_single_reverse_link_returns_one_base(self, kmer):
        # given
        b = LinksBuilder()
        b.with_link_for_kmer('R 3 1 GTT', 'AAA')
        links = b.build()

        # when
        walker = LinkWalker.from_links(links).load_kmer(kmer)

        # then
        if kmer == 'TTT':
            assert 1 == walker.n_junctions
            assert ['G'] == list(walker.next_junction_bases())
            walker.choose_branch('G')
            assert 1 == walker.n_junctions
            walker.choose_branch('T')
            assert 1 == walker.n_junctions
            walker.choose_branch('T')
            assert 0 == walker.n_junctions
        else:
            assert 0 == walker.n_junctions
            assert [] == list(walker.next_junction_bases())

    @pytest.mark.parametrize('kmer', ('AAA', 'TTT'))
    def test_with_reverse_and_forward_link_returns_one_base_each(self, kmer):
        # given
        b = LinksBuilder()
        b.with_link_for_kmer('F 1 1 A', 'AAA')
        b.with_link_for_kmer('R 1 1 A', 'AAA')
        links = b.build()

        # when
        walker = LinkWalker.from_links(links).load_kmer(kmer)

        # then
        assert ['A'] == list(walker.next_junction_bases())

    def test_with_two_reverse_and_forward_links_returns_two_bases_each(self):
        # given
        b = LinksBuilder()
        b.with_link_for_kmer('F 1 1 A', 'AAA')
        b.with_link_for_kmer('F 1 1 C', 'AAA')
        b.with_link_for_kmer('R 1 1 G', 'AAA')
        b.with_link_for_kmer('R 1 1 T', 'AAA')

        links = b.build()

        # when
        walker = LinkWalker.from_links(links)

        # then
        assert ['A', 'C'] == list(walker.load_kmer('AAA').next_junction_bases())
        assert ['G', 'T'] == list(walker.clear().load_kmer('TTT').next_junction_bases())

    def test_with_single_link_returns_each_junctions_base(self):
        # given
        b = LinksBuilder()
        b.with_link_for_kmer('F 3 1 ACT', 'AAA')
        links = b.build()

        # when
        walker = LinkWalker.from_links(links).load_kmer('AAA')

        # then
        assert ['A'] == list(walker.next_junction_bases())
        assert ['C'] == list(walker.choose_branch('A').next_junction_bases())
        assert ['T'] == list(walker.choose_branch('C').next_junction_bases())
        assert [] == list(walker.choose_branch('T').next_junction_bases())

    def test_raises_when_base_that_does_not_exist_is_chosen(self):
        # given
        b = LinksBuilder()
        b.with_link_for_kmer('F 1 1 A', 'AAA')
        links = b.build()

        # when
        walker = LinkWalker.from_links(links).load_kmer('AAA')

        # then
        with pytest.raises(KeyError):
            walker.choose_branch('C')


class TestUnitigLinkWalker:
    def test_y_graph_with_one_link_returns_one_node(self):
        # given
        links = LinksBuilder() \
            .with_link_for_kmer('F 1 1 C', 'AAA') \
            .build()

        b = UnitigBuilder()
        b.add_node(0, 'AAA')
        b.add_node(1, 'AAC')
        b.add_node(2, 'AAG')

        b.add_edge(0, 1)
        b.add_edge(0, 2)
        unitigs = b.build()

        # when
        walker = UnitigLinkWalker.from_links_unitigs_kmer_size_unitig(links, unitigs, 3, 0)

        # then
        assert [1] == list(walker.link_successors())
        assert [1] == list(walker.successors())
        with pytest.raises(KeyError):
            walker.choose(2)
        walker.choose(1)
        print(list(walker.unitigs.successors(1)))
        with pytest.raises(ValueError):
            list(walker.link_successors())
        assert [] == list(walker.successors())

    def test_two_ys_with_two_links_returns_one_and_two_nodes(self):
        # given
        links = LinksBuilder() \
            .with_link_for_kmer('F 2 1 CT', 'AAA') \
            .with_link_for_kmer('F 1 1 A', 'CCC') \
            .build()

        b = UnitigBuilder()
        b.add_node(0, 'AAA')
        b.add_node(1, 'AACCC')
        b.add_node(2, 'AAGCC')
        b.add_node(3, 'CCC')
        b.add_node(4, 'CCA')
        b.add_node(5, 'CCT')

        b.add_edge(0, 1)
        b.add_edge(0, 2)
        b.add_edge(1, 3)
        b.add_edge(2, 3)
        b.add_edge(3, 4)
        b.add_edge(3, 5)
        unitigs = b.build()

        # when
        walker = UnitigLinkWalker.from_links_unitigs_kmer_size_unitig(links, unitigs, 3, 0)

        # then
        assert [1] == list(walker.link_successors())
        walker.choose(1)
        with pytest.raises(ValueError):
            list(walker.link_successors())
        walker.choose(3)
        assert [4, 5] == sorted(walker.link_successors())
        walker.choose(5)
        with pytest.raises(ValueError):
            list(walker.link_successors())

    def test_two_ys_with_one_link_returns_both_nodes_the_second_time(self):
        # given
        links = LinksBuilder() \
            .with_link_for_kmer('F 1 1 C', 'AAA') \
            .build()

        b = UnitigBuilder()
        b.add_node(0, 'AAA')
        b.add_node(1, 'AACCC')
        b.add_node(2, 'AAGCC')
        b.add_node(3, 'CCC')
        b.add_node(4, 'CCA')
        b.add_node(5, 'CCT')

        b.add_edge(0, 1)
        b.add_edge(0, 2)
        b.add_edge(1, 3)
        b.add_edge(2, 3)
        b.add_edge(3, 4)
        b.add_edge(3, 5)
        unitigs = b.build()

        # when
        walker = UnitigLinkWalker.from_links_unitigs_kmer_size_unitig(links, unitigs, 3, 0)

        # then
        assert [1] == list(walker.successors())
        walker.choose(1)
        assert [3] == list(walker.successors())
        walker.choose(3)
        assert [4, 5] == sorted(walker.successors())
        walker.choose(5)
        assert [] == list(walker.successors())

    def test_raises_when_link_does_match_graph(self):
        # given
        links = LinksBuilder() \
            .with_link_for_kmer('F 1 1 T', 'AAA') \
            .build()

        b = UnitigBuilder()
        b.add_node(0, 'AAA')
        b.add_node(1, 'AAC')
        b.add_node(2, 'AAG')

        b.add_edge(0, 1)
        b.add_edge(0, 2)
        unitigs = b.build()

        # when
        walker = UnitigLinkWalker.from_links_unitigs_kmer_size_unitig(links, unitigs, 3, 0)

        # then
        with pytest.raises(ValueError, match='Links do not appear to match unitigs'):
            list(walker.link_successors())


class TestUnitigLinkWalker_copy:
    def test_y_graph_with_one_link_returns_copy_of_walker(self):
        # given
        links = LinksBuilder() \
            .with_link_for_kmer('F 1 1 C', 'AAA') \
            .build()

        b = UnitigBuilder()
        b.add_node(0, 'AAA')
        b.add_node(1, 'AAC')
        b.add_node(2, 'AAG')

        b.add_edge(0, 1)
        b.add_edge(0, 2)
        unitigs = b.build()
        walker = UnitigLinkWalker.from_links_unitigs_kmer_size_unitig(links, unitigs, 3, 0)

        # when
        new_walker = copy.copy(walker)
        new_walker.choose(1)

        # then
        assert [1] == list(walker.successors())
        assert [] == list(new_walker.successors())


class TestLinkedGraphTraverser:
    def test_y_graph_with_one_link_returns_one_node(self):
        # given
        links = LinksBuilder() \
            .with_link_for_kmer('F 1 1 C', 'AAA') \
            .build()

        b = UnitigBuilder()
        b.add_node(0, 'AAA')
        b.add_node(1, 'AAC')
        b.add_node(2, 'AAG')

        b.add_edge(0, 1)
        b.add_edge(0, 2)
        unitigs = b.build()

        # when
        traverser = LinkedGraphTraverser.from_graph_and_link_walker(
            unitigs,
            UnitigLinkWalker.from_links_unitigs_kmer_size_unitig(links, unitigs, 3, 0)
        )

        # then
        for node in range(3):
            assert node in traverser
        assert [1] == list(traverser[0])
        assert [] == list(traverser[1])
