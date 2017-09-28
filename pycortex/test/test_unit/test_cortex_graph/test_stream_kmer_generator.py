import attr

from hypothesis import given
from hypothesis import strategies as s

from pycortex.cortex_graph import kmer_generator_from_stream
from pycortex.test.builders.graph_body_builder import CortexGraphBodyBuilder, KmerRecord, kmers


@attr.s(slots=True)
class CortexGraphHeaderStub(object):
    kmer_size = attr.ib()
    kmer_container_size = attr.ib()
    num_colors = attr.ib()


class TestStreamKmerGenerator(object):
    @given(s.data(),
           s.integers(min_value=1),
           s.integers(min_value=0),
           s.integers(min_value=0, max_value=4))
    def test_parses_records(self, data, kmer_size, num_colors, n_kmers):
        # given
        builder = CortexGraphBodyBuilder(kmer_size=kmer_size)

        expected_kmers = []
        for _ in range(n_kmers):
            kmer = data.draw(kmers(kmer_size, num_colors))
            builder.with_kmer_record(kmer)
            expected_kmers.append(kmer)

        header = CortexGraphHeaderStub(kmer_size, builder.kmer_container_size, num_colors)

        # when
        for kmer, expected_kmer in zip(
                kmer_generator_from_stream(builder.build(), header), expected_kmers):
            # then
            assert expected_kmer.kmer == kmer.kmer
            assert expected_kmer.coverage == kmer.coverage
            assert expected_kmer.edges == kmer.edges

    def test_parses_AAC_kmer(self):
        kmer_container_size = 1
        kmer_size = 3
        num_colors = 0
        header = CortexGraphHeaderStub(kmer_size, kmer_container_size, num_colors)
        builder = CortexGraphBodyBuilder(kmer_container_size)

        expected_kmer = KmerRecord('AAC', tuple(), tuple())
        builder.with_kmer_record(expected_kmer)

        kmer = next(kmer_generator_from_stream(builder.build(), header))

        assert expected_kmer.kmer == kmer.kmer
        assert expected_kmer.coverage == kmer.coverage
        assert expected_kmer.edges == kmer.edges
