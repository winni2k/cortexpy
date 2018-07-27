import pytest
import math
from Bio.Seq import reverse_complement, Seq
from hypothesis import given, assume, settings, strategies as s

from cortexpy import edge_set
from cortexpy.graph.parser.kmer import (
    EmptyKmerBuilder, connect_kmers, StringKmerConverter, Kmer,
    KmerData,
    disconnect_kmers,
)
from cortexpy.test.builder.graph.body import KmerRecord
from cortexpy.test.builder.graph.kmer import dna_sequences, kmer_strings, kmer_records


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

    def test_raises_on_non_lexlo_kmer(self):
        # when
        rec = KmerRecord('AAA', [1], [edge_set.empty()])
        kmer = Kmer.from_kmer_data(KmerData(rec.to_bytestring(), kmer_size=3, num_colors=1))
        with pytest.raises(AttributeError):
            kmer.kmer = reverse_complement(kmer.kmer)


class TestAddColor(object):
    @given(s.data(), s.integers(min_value=1, max_value=7), s.integers(min_value=0, max_value=10))
    @settings(max_examples=10)
    def test_increases_color_count(self, data, num_colors, increment_color):
        # given
        kmer = EmptyKmerBuilder(num_colors).build_or_get(data.draw(kmer_strings()))

        # when
        kmer.increment_color_coverage(increment_color)

        # then
        assert kmer.num_colors == max(num_colors, increment_color + 1)
        assert len(kmer.coverage) == kmer.num_colors
        assert isinstance(kmer.coverage[increment_color], int)
        assert len(kmer.edges) == kmer.num_colors
        assert list(kmer.edges) == [edge_set.empty() for _ in range(kmer.num_colors)]
        if increment_color >= num_colors:
            assert 1 == kmer.coverage[increment_color]
        else:
            assert 2 == kmer.coverage[increment_color]


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

    @given(s.data(), s.text('ACGT', min_size=1, max_size=1), s.integers(min_value=1, max_value=7),
           s.booleans(), s.booleans())
    def test_works_on_edit_distance_one_forward_in_any_order(self, data, letter, num_colors,
                                                             add_letter_to_end,
                                                             add_letter_to_revcomp):
        # given
        builder = EmptyKmerBuilder(num_colors)
        kmer1 = builder.build_or_get(data.draw(kmer_strings()))
        assert kmer1.kmer < reverse_complement(kmer1.kmer)

        if add_letter_to_revcomp:
            kmer1_string = reverse_complement(kmer1.kmer)
        else:
            kmer1_string = kmer1.kmer

        if add_letter_to_end:
            kmer2_string = kmer1_string[1:] + letter
        else:
            kmer2_string = letter + kmer1_string[:-1]

        kmer2 = builder.build_or_get(kmer2_string)

        if add_letter_to_revcomp:
            kmer1_edge_letter = reverse_complement(letter)
        else:
            kmer1_edge_letter = letter.lower()

        if add_letter_to_end:
            kmer1_edge_letter = kmer1_edge_letter.swapcase()
            kmer2_edge_letter = kmer1_string[0].lower()
        else:
            kmer2_edge_letter = kmer1_string[-1]

        if reverse_complement(kmer2_string) < kmer2_string:
            kmer2_edge_letter = reverse_complement(kmer2_edge_letter).swapcase()

        for color in range(num_colors):
            # when
            connect_kmers(kmer1, kmer2, color)

            # then
            assert kmer1.edges[color].is_edge(kmer1_edge_letter)
            assert kmer2.edges[color].is_edge(kmer2_edge_letter)

            # let's also test connected check
            if not add_letter_to_revcomp == add_letter_to_end:
                assert kmer1.has_outgoing_edge_to_kmer_in_color(kmer2, color)
            else:
                assert kmer1.has_incoming_edge_from_kmer_in_color(kmer2, color)

            # and let's check disconnect kmers
            disconnect_kmers(kmer1, kmer2, [color])
            assert not kmer1.edges[color].is_edge(kmer1_edge_letter)
            assert not kmer2.edges[color].is_edge(kmer2_edge_letter)


class TestHasEdgeTo(object):
    @given(s.integers(min_value=0, max_value=2))
    def test_is_true_for_two_connected_kmers(self, color):
        # given
        kmer1 = EmptyKmerBuilder(3).build_or_get('AAA')
        kmer2 = EmptyKmerBuilder(3).build_or_get('AAT')
        kmer1.edges[color].add_edge('T')
        kmer2.edges[color].add_edge('a')

        # when/then
        assert kmer1.has_outgoing_edge_to_kmer_in_color(kmer2, color)
        assert kmer2.has_incoming_edge_from_kmer_in_color(kmer1, color)

    @given(s.booleans())
    def test_raises_if_connection_is_unidirectional(self, connect_first):
        # given
        kmer1 = EmptyKmerBuilder(1).build_or_get('AAA')
        kmer2 = EmptyKmerBuilder(1).build_or_get('AAT')
        if connect_first:
            kmer1.edges[0].add_edge('T')
        else:
            kmer2.edges[0].add_edge('a')

        # when/then
        with pytest.raises(ValueError):
            kmer1.has_outgoing_edge_to_kmer_in_color(kmer2, 0)
        with pytest.raises(ValueError):
            kmer2.has_incoming_edge_from_kmer_in_color(kmer1, 0)

    def test_is_false_for_unconnected_neighbor_kmers(self):
        # given
        kmer1 = EmptyKmerBuilder(1).build_or_get('AAA')
        kmer2 = EmptyKmerBuilder(1).build_or_get('AAT')

        # when/then
        assert not kmer1.has_outgoing_edge_to_kmer_in_color(kmer2, 0)
        assert not kmer2.has_incoming_edge_from_kmer_in_color(kmer1, 0)

    def test_raises_for_for_non_neighbor_kmers(self):
        # given
        kmer1 = EmptyKmerBuilder(1).build_or_get('AAA')
        kmer2 = EmptyKmerBuilder(1).build_or_get('ACC')

        # when/then
        with pytest.raises(ValueError):
            kmer1.has_outgoing_edge_to_kmer_in_color(kmer2, 0)
        with pytest.raises(ValueError):
            kmer2.has_incoming_edge_from_kmer_in_color(kmer1, 0)

    @given(s.integers(min_value=0, max_value=2), s.integers(min_value=0, max_value=2))
    def test_is_false_for_non_connection_color(self, color_to_link, color_to_check):
        # given
        assume(color_to_check != color_to_link)
        kmer1 = EmptyKmerBuilder(3).build_or_get('AAA')
        kmer2 = EmptyKmerBuilder(3).build_or_get('AAT')
        kmer1.edges[color_to_link].add_edge('T')
        kmer2.edges[color_to_link].add_edge('a')

        # when/then
        assert not kmer1.has_outgoing_edge_to_kmer_in_color(kmer2, color_to_check)
        assert not kmer2.has_incoming_edge_from_kmer_in_color(kmer1, color_to_check)


class TestStringKmerConverter(object):
    @given(s.lists(s.sampled_from(list('ACGT')), min_size=1, max_size=65))
    def test_kmer_sizes(self, letters):
        # given
        string_kmer = ''.join(letters)
        converter = StringKmerConverter(kmer_size=len(letters))

        # when
        uints = converter.to_uints(string_kmer)

        # then
        assert math.ceil(len(letters) / 32) == len(uints)

    @given(s.integers(min_value=0, max_value=3))
    def test_converts_least_significant_letter(self, letter_idx):
        # given
        kmer_string = 'A' * 31 + list('ACGT')[letter_idx]
        converter = StringKmerConverter(kmer_size=len(kmer_string))

        # when
        uints = converter.to_uints(kmer_string)

        # then
        assert letter_idx == uints[0]
        assert 1 == len(uints)

    @given(s.integers(min_value=0, max_value=3))
    def test_converts_most_significant_letter(self, letter_idx):
        # given
        kmer_string = list('ACGT')[letter_idx] + 'A' * 31
        converter = StringKmerConverter(kmer_size=len(kmer_string))

        # when
        uints = converter.to_uints(kmer_string)

        # then
        assert int(letter_idx << 62) == uints[0]
        assert 1 == len(uints)

    @given(s.integers(min_value=0, max_value=3),
           s.integers(min_value=1, max_value=129))
    def test_with_partially_filled_container_converts_least_significant_letter(self, letter_idx,
                                                                               kmer_size):
        # given
        kmer_string = 'A' * (kmer_size - 1) + list('ACGT')[letter_idx]
        converter = StringKmerConverter(kmer_size=kmer_size)

        # when
        uints = converter.to_uints(kmer_string)

        # then
        assert letter_idx == uints[-1]
        assert math.ceil(kmer_size / 32) == len(uints)

    @given(s.data(),
           s.integers(min_value=0, max_value=129),
           s.integers(min_value=1, max_value=3),
           s.booleans()
           )
    def test_converts_to_raw(self, data, kmer_size, num_colors, from_bioseq):
        # given
        assume(kmer_size % 2 == 1)
        kmer = Kmer.from_kmer_data(
            KmerData(data.draw(kmer_records(kmer_size, num_colors)).to_bytestring(),
                     kmer_size=kmer_size,
                     num_colors=num_colors))
        converter = StringKmerConverter(kmer.kmer_size)

        # when
        if from_bioseq:
            raw_kmer = converter.to_raw(Seq(kmer.kmer))
        else:
            raw_kmer = converter.to_raw(kmer.kmer)

        # then
        assert kmer.get_raw_kmer() == raw_kmer
