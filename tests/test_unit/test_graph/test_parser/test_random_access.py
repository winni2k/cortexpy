import io
import random
from unittest import mock

import numpy as np
import pytest
from hypothesis import given, assume, HealthCheck, settings
from hypothesis import strategies as s

import cortexpy.graph.parser.random_access as parser
import cortexpy.graph.serializer.kmer as kmer_serializer
import cortexpy.test.builder as builder
from cortexpy.graph.parser.header import Header
from cortexpy.graph.parser.random_access import KmerUintSequence
from cortexpy.test.builder.graph.body import KmerRecord, as_edge_set
from cortexpy.test.builder.graph.kmer import kmer_records
from cortexpy.utils import lexlo


class TestDunderGetitemDunder:
    RAClass = parser.RandomAccess

    @settings(suppress_health_check=[HealthCheck.too_slow])
    @given(s.data(),
           s.integers(min_value=1, max_value=4).map(lambda i: i * 2 - 1),
           s.integers(min_value=1, max_value=3),
           s.integers(min_value=0, max_value=3))
    def test_record_retrieval(self, data, kmer_size, num_colors, n_kmers):
        # given
        assume(n_kmers <= 4 ** (kmer_size - 1))
        graph_builder = (builder.Graph()
                         .with_kmer_size(kmer_size)
                         .with_num_colors(num_colors))

        expected_kmers = []
        seen_kmers = set()
        for _ in range(n_kmers):
            kmer = data.draw(kmer_records(kmer_size, num_colors))
            while kmer.kmer in seen_kmers:
                kmer = data.draw(kmer_records(kmer_size, num_colors))
            seen_kmers.add(kmer.kmer)
            graph_builder.with_kmer_record(kmer)
            expected_kmers.append(kmer)

        cg = self.RAClass(graph_builder.build())

        # when
        for expected_kmer in expected_kmers:
            kmer = cg[expected_kmer.kmer]

            # then
            assert expected_kmer.kmer == kmer.kmer
            assert np.all(expected_kmer.coverage == kmer.coverage)
            for expected, actual in zip(expected_kmer.edges, kmer.edges):
                assert expected == actual

    def test_raises_on_missing_kmer(self):
        # given
        graph_builder = builder.Graph()
        graph_builder.with_kmer_size(3)

        cg = self.RAClass(graph_builder.build())

        # when
        with pytest.raises(KeyError):
            cg['AAA']

    def test_raises_on_even_kmer_size(self):
        # given
        graph_builder = builder.Graph()
        graph_builder.with_kmer_size(2)

        # when
        with pytest.raises(ValueError):
            self.RAClass(graph_builder.build())


class TestSlurpedDunderGetitemDunder(TestDunderGetitemDunder):
    RAClass = parser.SlurpedRandomAccess.from_handle


class TestGetKmerForString(object):
    @pytest.mark.parametrize('RAClass',
                             (parser.RandomAccess, parser.SlurpedRandomAccess.from_handle))
    def test_gets_aaa_for_ttt_query(self, RAClass):
        # given
        graph_builder = builder.Graph()
        graph_builder.with_kmer_size(3)
        graph_builder.with_num_colors(1)

        expected_kmer = KmerRecord('AAA', [1], [as_edge_set('........')])
        graph_builder.with_kmer_record(expected_kmer)

        cg = RAClass(graph_builder.build())

        # when
        assert expected_kmer.kmer == cg.get_kmer_for_string('AAA').kmer
        assert expected_kmer.kmer == cg.get_kmer_for_string('TTT').kmer


class TestDunderIterDunder:
    RAClass = parser.RandomAccess

    @settings(suppress_health_check=[HealthCheck.too_slow])
    @given(s.data(),
           s.integers(min_value=1, max_value=7).map(lambda i: i * 2 - 1),
           s.integers(min_value=1, max_value=7),
           s.integers(min_value=0, max_value=4),
           s.booleans())
    def test_parses_records(self, data, kmer_size, num_colors, n_kmers, test_serializer):
        # given
        assume(n_kmers <= (4 ** kmer_size) / 4)

        graph_builder = (builder.Graph()
                         .with_kmer_size(kmer_size)
                         .with_num_colors(num_colors))

        expected_kmers = []
        seen = set()
        for _ in range(n_kmers):
            kmer = data.draw(kmer_records(kmer_size, num_colors))
            while kmer.kmer in seen:
                kmer = data.draw(kmer_records(kmer_size, num_colors))
            seen.add(kmer.kmer)
            graph_builder.with_kmer_record(kmer)
            expected_kmers.append(kmer)
        ra_parser = self.RAClass(graph_builder.build())

        if test_serializer:
            for real_kmer in ra_parser.values():
                buffer = io.BytesIO()
                real_kmer.dump(buffer)
                assert real_kmer._kmer_data._data == buffer.getvalue()

            sample_names = ra_parser.sample_names
            buffer = io.BytesIO()
            key_list = list(ra_parser.keys())
            random.shuffle(key_list)
            kmer_serializer \
                .Kmers(keys=key_list,
                       val_callable=lambda k: ra_parser[k],
                       kmer_size=kmer_size,
                       num_colors=num_colors,
                       sample_names=sample_names) \
                .dump(buffer)
            buffer.seek(0)
            ra_parser = self.RAClass(buffer)

        # when
        for expected_kmer in expected_kmers:
            kmer = ra_parser[expected_kmer.kmer]

            # then
            assert expected_kmer.kmer == kmer.kmer
            assert np.all(expected_kmer.coverage == kmer.coverage)
            assert expected_kmer.edges == kmer.edges

    def test_gets_aaa(self):
        # given
        graph_builder = (builder.Graph()
                         .with_kmer_size(3)
                         .with_num_colors(1))

        expected_kmer = KmerRecord('AAA', (1,), [as_edge_set('........')])
        graph_builder.with_kmer_record(expected_kmer)

        cg = self.RAClass(graph_builder.build())

        # when
        for kmer in cg.values():
            assert expected_kmer.kmer == kmer.kmer
            assert np.all(expected_kmer.coverage == kmer.coverage)
            assert expected_kmer.edges == kmer.edges

    def test_works_with_no_kmers(self):
        # given
        graph_builder = (builder.Graph()
                         .with_kmer_size(3)
                         .with_num_colors(1))

        cg = self.RAClass(graph_builder.build())

        # when/then
        assert list(cg) == []


class TestSlurpedDunderIterDunder(TestDunderIterDunder):
    RAClass = parser.SlurpedRandomAccess.from_handle


class TestKmerUintSequence(object):
    @given(s.data(),
           s.integers(min_value=1, max_value=129),
           s.integers(min_value=1, max_value=5))
    def test_index(self, data, kmer_size, n_kmers):
        # given
        assume(kmer_size % 2 == 1)
        num_colors = 1
        graph_builder = (builder.Graph()
                         .with_kmer_size(kmer_size)
                         .with_num_colors(num_colors))

        expected_kmers = []
        seen_kmers = set()
        for _ in range(n_kmers):
            kmer = data.draw(kmer_records(kmer_size, num_colors))
            while kmer.kmer in seen_kmers:
                kmer = data.draw(kmer_records(kmer_size, num_colors))
            seen_kmers.add(kmer.kmer)
            graph_builder.with_kmer_record(kmer)
            expected_kmers.append(kmer)
        expected_kmers = sorted(expected_kmers)

        graph_stream = graph_builder.build()
        header_stream = graph_builder.header.build()
        header = Header.from_stream(header_stream)

        # when
        sequence = KmerUintSequence(graph_handle=graph_stream,
                                    body_start=len(header_stream.getvalue()),
                                    header=header,
                                    n_records=len(expected_kmers))
        # then
        for idx, expected_kmer in enumerate(expected_kmers):
            # then
            assert idx == sequence.index_kmer_string(expected_kmer.kmer)


class TestCacheComplexity:
    @given(s.integers(min_value=0, max_value=16))
    def test_regular_ra_parser(self, num_kmers):
        # given
        kmer_size = 11
        b = builder.Graph() \
            .with_kmer_size(kmer_size) \
            .with_num_colors(1)
        seen_kmers = set()
        for _ in range(num_kmers):
            kmer_string = lexlo(''.join([random.choice('ACGT') for _ in range(kmer_size)]))
            while kmer_string in seen_kmers:
                kmer_string = lexlo(''.join([random.choice('ACGT') for _ in range(kmer_size)]))
            seen_kmers.add(kmer_string)
            b.with_kmer(kmer_string)
        fh = b.build()

        ra = parser.RandomAccess(fh, kmer_cache_size=None)
        for k_string in list(ra):
            ra[k_string]
        with mock.patch.object(fh, 'read', wraps=fh.read) as mocked_read:
            # when
            for seen_kmer in sorted(seen_kmers):
                ra[seen_kmer]

            # then
            assert 2 * num_kmers == mocked_read.call_count

    @given(s.integers(min_value=0, max_value=16))
    def test_slurping_ra_parser(self, num_kmers):
        # given
        kmer_size = 11
        b = builder.Graph() \
            .with_kmer_size(kmer_size) \
            .with_num_colors(1)
        seen_kmers = set()
        for _ in range(num_kmers):
            kmer_string = lexlo(''.join([random.choice('ACGT') for _ in range(kmer_size)]))
            while kmer_string in seen_kmers:
                kmer_string = lexlo(''.join([random.choice('ACGT') for _ in range(kmer_size)]))
            seen_kmers.add(kmer_string)
            b.with_kmer(kmer_string)
        fh = b.build()

        fh.seek(0)
        with mock.patch.object(fh, 'seek', wraps=fh.seek) as mocked_seek:
            # when
            ra = parser.SlurpedRandomAccess.from_handle(fh)

            # then
            assert 0 == mocked_seek.call_count

            for seen_kmer in sorted(seen_kmers):
                ra[seen_kmer]
            assert 0 == mocked_seek.call_count

        fh.seek(0)
        num_header_reads = 10
        num_eof_reads = 1
        with mock.patch.object(fh, 'read', wraps=fh.read) as mocked_read:
            # when
            ra = parser.SlurpedRandomAccess.from_handle(fh)

            # then
            assert num_kmers + num_eof_reads + num_header_reads == mocked_read.call_count

            for seen_kmer in sorted(seen_kmers):
                ra[seen_kmer]
            assert num_kmers + num_eof_reads + num_header_reads == mocked_read.call_count
