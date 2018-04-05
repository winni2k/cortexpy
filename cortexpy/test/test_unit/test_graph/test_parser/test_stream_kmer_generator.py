import random
from unittest import mock

import numpy as np
from hypothesis import given
from hypothesis import strategies as s

from cortexpy.graph.parser.streaming import kmer_generator_from_stream_and_header
from cortexpy.test.builder.graph.body import Body, KmerRecord
from cortexpy.test.builder.graph.kmer import kmers
from cortexpy.test.mock.graph import Header
from cortexpy.utils import lexlo


class TestStreamKmerGenerator(object):
    @given(s.data(),
           s.integers(min_value=1, max_value=260),
           s.integers(min_value=0, max_value=10),
           s.integers(min_value=0, max_value=4))
    def test_parses_records(self, data, kmer_size, num_colors, n_kmers):
        # given
        builder = Body(kmer_size=kmer_size)

        expected_kmers = []
        for _ in range(n_kmers):
            kmer = data.draw(kmers(kmer_size, num_colors))
            builder.with_kmer_record(kmer)
            expected_kmers.append(kmer)

        header = Header(kmer_size, builder.kmer_container_size, num_colors)

        # when
        for kmer, expected_kmer in zip(
            kmer_generator_from_stream_and_header(builder.build(), header),
            expected_kmers
        ):
            # then
            assert expected_kmer.kmer == kmer.kmer
            assert np.all(expected_kmer.coverage == kmer.coverage)
            assert expected_kmer.edges == kmer.edges

    @given(s.integers(min_value=0, max_value=16))
    def test_complexity(self, n_kmers):
        # given
        num_colors = 0
        kmer_size = 11
        expected_num_calls = n_kmers + 1
        builder = Body(kmer_size=kmer_size)

        for _ in range(n_kmers):
            kmer_string = lexlo(''.join([random.choice('ACGT') for _ in range(kmer_size)]))
            builder.with_kmer_record(KmerRecord(kmer_string, [], []))

        header = Header(kmer_size, builder.kmer_container_size, num_colors)

        # when
        fh = builder.build()
        with mock.patch.object(fh, 'read', wraps=fh.read) as mocked_seek:
            # when
            for _ in kmer_generator_from_stream_and_header(fh, header):
                pass

            # then
            assert expected_num_calls == mocked_seek.call_count

    def test_parses_aac_kmer(self):
        kmer_container_size = 1
        kmer_size = 3
        num_colors = 0
        header = Header(kmer_size, kmer_container_size, num_colors)
        builder = Body(kmer_container_size)

        expected_kmer = KmerRecord('AAC', [], [])
        builder.with_kmer_record(expected_kmer)

        kmer = next(kmer_generator_from_stream_and_header(builder.build(), header))

        assert expected_kmer.kmer == kmer.kmer
        assert np.all(expected_kmer.coverage == kmer.coverage)
        assert expected_kmer.edges == kmer.edges
