import pytest
from Bio.Seq import reverse_complement
from hypothesis import given, assume, settings
from hypothesis import strategies as s

from pycortex.graph.parser.kmer import EmptyKmerBuilder, connect_kmers
from pycortex.test.builder.graph.kmer import dna_sequences, kmer_strings


class TestBuildKmer(object):
    @given(s.data())
    @settings(max_examples=10)
    def test_requires_odd_length(self, data):
        # given
        kmer_string = data.draw(dna_sequences())
        assume(len(kmer_string) % 2 == 0)

        # when/then
        with pytest.raises(ValueError):
            EmptyKmerBuilder().build_or_get(kmer_string)

    @given(s.data())
    @settings(max_examples=10)
    def test_requires_length_3_or_more(self, data):
        # given
        kmer_string = data.draw(dna_sequences(max_size=2))

        # when/then
        with pytest.raises(ValueError):
            EmptyKmerBuilder().build_or_get(kmer_string)

    @given(s.data())
    @settings(max_examples=10)
    def test_creates_revcomp(self, data):
        # given
        kmer_string = data.draw(kmer_strings())

        # when
        kmer = EmptyKmerBuilder().build_or_get(kmer_string)

        # then
        if reverse_complement(kmer_string) < kmer_string:
            assert kmer.kmer < kmer_string
        else:
            assert kmer.kmer == kmer_string

    def test_raises_on_negative_num_colors(self):
        with pytest.raises(ValueError):
            EmptyKmerBuilder(-1)

    @given(s.data())
    @settings(max_examples=10)
    def test_returns_same_kmer_if_built_twice(self, data):
        # given
        kmer_string = data.draw(kmer_strings())
        builder = EmptyKmerBuilder()

        # when
        kmer1 = builder.build_or_get(kmer_string)
        kmer2 = builder.build_or_get(kmer_string)

        # then
        assert kmer1 is kmer2


class TestAddColor(object):
    @given(s.data(), s.integers(min_value=0, max_value=7))
    @settings(max_examples=10)
    def test_increases_color_count(self, data, num_colors):
        # given
        kmer = EmptyKmerBuilder(num_colors).build_or_get(data.draw(kmer_strings()))

        # when
        kmer.append_color()

        # then
        assert kmer.num_colors == num_colors + 1


class TestConnectKmers(object):
    def test_raises_if_equal_kmers_are_not_the_same_object(self):
        # given
        kmer1 = EmptyKmerBuilder(1).build_or_get('AAA')
        kmer2 = EmptyKmerBuilder(1).build_or_get('AAA')
        assert kmer1 is not kmer2

        # when/then
        with pytest.raises(ValueError):
            connect_kmers(kmer1, kmer2, 0)

    def test_connects_revcomp_kmer_to_itself(self):
        # given
        kmer = EmptyKmerBuilder(1).build_or_get('TAA')

        # when
        connect_kmers(kmer, kmer, 0)

        # then
        assert kmer.edges[0].is_edge('t')
        for letter in 'acgACGT':
            assert not kmer.edges[0].is_edge(letter)

    def test_connects_kmer_to_itself(self):
        # given
        kmer = EmptyKmerBuilder(1).build_or_get('AAA')

        # when
        connect_kmers(kmer, kmer, 0)

        # then
        assert kmer.edges[0].is_edge('a')
        assert kmer.edges[0].is_edge('A')

    def test_connects_kmer_to_itself_for_one_sided_connection(self):
        # given
        kmer = EmptyKmerBuilder(1).build_or_get('TAA')

        # when
        connect_kmers(kmer, kmer, 0)

        # then
        assert kmer.edges[0].is_edge('t')
        for letter in 'acgACGT':
            assert not kmer.edges[0].is_edge(letter)

    def test_raises_if_first_kmer_cannot_be_connected_to_second_kmer(self):
        # given
        builder = EmptyKmerBuilder(1)
        kmer1 = builder.build_or_get('AAA')
        kmer2 = builder.build_or_get('ACC')

        # when/then
        with pytest.raises(ValueError):
            connect_kmers(kmer1, kmer2, 0)

    def test_aag_act_connection(self):
        # given
        builder = EmptyKmerBuilder(1)
        kmer1 = builder.build_or_get('AAG')
        kmer2 = builder.build_or_get('ACT')

        # when
        connect_kmers(kmer1, kmer2, 0)

        # then
        assert kmer1.edges[0].is_edge('T')
        assert kmer2.edges[0].is_edge('T')

    @given(s.data(), s.text('ACGT', min_size=1, max_size=1), s.integers(min_value=0, max_value=7),
           s.booleans(), s.booleans())
    def test_works_on_edit_distance_one_forward_in_any_order(self, data, letter, num_colors,
                                                             add_letter_to_end,
                                                             add_letter_to_revcomp):
        # given
        builder = EmptyKmerBuilder(num_colors)
        first_kmer = builder.build_or_get(data.draw(kmer_strings()))
        assert first_kmer.kmer < reverse_complement(first_kmer.kmer)

        if add_letter_to_revcomp:
            first_kmer_string = reverse_complement(first_kmer.kmer)
        else:
            first_kmer_string = first_kmer.kmer

        if add_letter_to_end:
            second_kmer_string = first_kmer_string[1:] + letter
        else:
            second_kmer_string = letter + first_kmer_string[:-1]

        second_kmer = builder.build_or_get(second_kmer_string)

        for color in range(num_colors):
            # when
            connect_kmers(first_kmer, second_kmer, color)

            # then
            if add_letter_to_revcomp:
                first_kmer_edge_letter = reverse_complement(letter)
            else:
                first_kmer_edge_letter = letter.lower()

            if add_letter_to_end:
                first_kmer_edge_letter = first_kmer_edge_letter.swapcase()
                second_kmer_edge_letter = first_kmer_string[0].lower()
            else:
                second_kmer_edge_letter = first_kmer_string[-1]

            if reverse_complement(second_kmer_string) < second_kmer_string:
                second_kmer_edge_letter = reverse_complement(second_kmer_edge_letter).swapcase()

            assert first_kmer.edges[color].is_edge(first_kmer_edge_letter)
            assert second_kmer.edges[color].is_edge(second_kmer_edge_letter)
