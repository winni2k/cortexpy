import attr
import numpy as np
import pytest
from hypothesis import given
from hypothesis import strategies as s

import cortexpy.graph.parser.random_access as parser
import cortexpy.test.builder as builder
from cortexpy.graph.parser.random_access_collection import RandomAccessCollection
from cortexpy.test.builder.graph.body import KmerRecord, as_edge_set
from cortexpy.test.builder.graph.kmer import kmer_records


@attr.s(slots=True)
class GraphCollection(object):
    n_colors_per_graph = attr.ib()
    kmer_size = attr.ib()
    graph_builders = attr.ib(init=False)

    def __attrs_post_init__(self):
        self.graph_builders = [
            builder.Graph().with_kmer_size(self.kmer_size).with_num_colors(n_graph_colors)
            for n_graph_colors in self.n_colors_per_graph]

    def with_kmer_record(self, kmer):
        colors_so_far = 0
        for graph_idx, n_colors in enumerate(self.n_colors_per_graph):
            last_color_idx = n_colors + colors_so_far
            graph_kmer = KmerRecord(kmer.kmer,
                                    kmer.coverage[colors_so_far:last_color_idx],
                                    kmer.edges[colors_so_far:last_color_idx])
            self.graph_builders[graph_idx].with_kmer_record(graph_kmer)
            colors_so_far += n_colors
        return self

    def with_kmer_for_graph(self, graph_idx, kmer_string, color_coverage=1, edges='........'):
        self.graph_builders[graph_idx].with_kmer(kmer_string, color_coverage=color_coverage,
                                                 edges=edges)
        return self

    def build(self):
        return RandomAccessCollection(
            ra_parsers=[parser.RandomAccess(builder.build()) for builder in self.graph_builders])


class TestDunderGetitemDunder(object):
    @given(s.data(),
           s.integers(min_value=1, max_value=3),
           s.lists(s.integers(min_value=1, max_value=3), min_size=1, max_size=3),
           s.integers(min_value=0, max_value=3))
    def test_record_retrieval(self, data, base_kmer_size, num_colors_per_graph, n_kmers):
        # given
        kmer_size = base_kmer_size * 2 + 1
        num_colors = sum(num_colors_per_graph)
        collection_builder = GraphCollection(n_colors_per_graph=num_colors_per_graph,
                                             kmer_size=kmer_size)

        expected_kmers = []
        seen_kmers = set()
        for _ in range(n_kmers):
            kmer = data.draw(kmer_records(kmer_size, num_colors))
            while kmer.kmer in seen_kmers:
                kmer = data.draw(kmer_records(kmer_size, num_colors))
            seen_kmers.add(kmer.kmer)
            collection_builder.with_kmer_record(kmer)
            expected_kmers.append(kmer)

        collection = collection_builder.build()

        # when
        for expected_kmer in expected_kmers:
            kmer = collection[expected_kmer.kmer]

            # then
            assert expected_kmer.kmer == kmer.kmer
            assert expected_kmer.coverage == kmer.coverage
            assert expected_kmer.edges == kmer.edges

    def test_raises_on_missing_kmer(self):
        # given
        collection_builder = GraphCollection(n_colors_per_graph=[1, 1],
                                             kmer_size=3)

        collection = collection_builder.build()

        # when
        with pytest.raises(KeyError):
            collection['AAA']

    def test_does_not_raise_on_partially_missing_kmer(self):
        # given
        collection_builder = GraphCollection(n_colors_per_graph=[1, 2],
                                             kmer_size=3)
        collection_builder.with_kmer_for_graph(0, 'AAA', color_coverage=1, edges='....A...')
        collection_builder.with_kmer_for_graph(1, 'CCC', color_coverage=(2, 3),
                                               edges=('........', '...t....'))
        collection = collection_builder.build()

        # when
        kmer1 = collection['AAA']
        kmer2 = collection['CCC']

        # then
        assert 'AAA' == kmer1.kmer
        assert (1, 0, 0) == tuple(kmer1.coverage)
        assert ('....A...', '........', '........') == tuple(str(e) for e in kmer1.edges)

        # then
        assert 'CCC' == kmer2.kmer
        assert np.all((0, 2, 3) == kmer2.coverage)
        assert ('........', '........', '...t....') == tuple(str(e) for e in kmer2.edges)


class TestGetKmerForString(object):
    def test_gets_aaa_for_ttt_query(self):
        # given
        graph_builder = builder.Graph()
        graph_builder.with_kmer_size(3)
        graph_builder.with_num_colors(1)

        expected_kmer = KmerRecord('AAA', [1], [as_edge_set('........')])
        graph_builder.with_kmer_record(expected_kmer)

        cg = RandomAccessCollection([parser.RandomAccess(graph_builder.build())])

        # when
        assert expected_kmer.kmer == cg.get_kmer_for_string('AAA').kmer
        assert expected_kmer.kmer == cg.get_kmer_for_string('TTT').kmer
