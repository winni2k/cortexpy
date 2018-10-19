import pytest

from cortexpy.links import LinkWalker
from cortexpy.test.builder.graph.cortex import LinksBuilder


class TestLinkGroup_GetLinkJunctions:
    @pytest.mark.parametrize('is_lexlo,is_forward',
                             [
                                 (True, True),
                                 (True, False),
                                 (False, True),
                                 (False, False),
                             ])
    def test_with_single_link_returns_one_base(self, is_lexlo, is_forward):
        # given
        b = LinksBuilder()
        b.with_link_for_kmer('F 3 1 ACC', 'AAA')

        # when
        links = b.build()

        # then
        junctions = list(
            links.body['AAA'].get_link_junctions(is_lexlo, in_kmer_orientation=is_forward))
        if is_lexlo == is_forward:
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
        walker = LinkWalker(links).load_kmer(kmer)

        # then
        if kmer == 'AAA':
            assert 1 == walker.n_junctions
            assert ['A'] == list(walker.next_junction_bases())
            walker.choose_junction('A')
            assert 1 == walker.n_junctions
            walker.choose_junction('C')
            assert 1 == walker.n_junctions
            walker.choose_junction('C')
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
        walker = LinkWalker(links).load_kmer(kmer)

        # then
        assert [kmer[0]] == list(walker.next_junction_bases())

    def test_with_two_reverse_and_forward_links_returns_two_bases_each(self):
        # given
        b = LinksBuilder()
        b.with_link_for_kmer('F 1 1 A', 'AAA')
        b.with_link_for_kmer('F 1 1 C', 'AAA')
        b.with_link_for_kmer('R 1 1 G', 'AAA')
        b.with_link_for_kmer('R 1 1 T', 'AAA')

        links = b.build()

        # when
        walker = LinkWalker(links)

        # then
        assert ['A', 'C'] == list(walker.load_kmer('AAA').next_junction_bases())
        assert ['C', 'A'] == list(walker.clear().load_kmer('TTT').next_junction_bases())

    def test_with_single_link_returns_each_junctions_base(self):
        # given
        b = LinksBuilder()
        b.with_link_for_kmer('F 3 1 ACT', 'AAA')
        links = b.build()

        # when
        walker = LinkWalker(links).load_kmer('AAA')

        # then
        assert ['A'] == list(walker.next_junction_bases())
        assert ['C'] == list(walker.choose_junction('A').next_junction_bases())
        assert ['T'] == list(walker.choose_junction('C').next_junction_bases())
        assert [] == list(walker.choose_junction('T').next_junction_bases())

    def test_raises_when_base_that_does_not_exist_is_chosen(self):
        # given
        b = LinksBuilder()
        b.with_link_for_kmer('F 1 1 A', 'AAA')
        links = b.build()

        # when
        walker = LinkWalker(links).load_kmer('AAA')

        # then
        with pytest.raises(KeyError):
            walker.choose_junction('C')
