from itertools import repeat
import cortexpy.test.builder as builder
from cortexpy.graph.parser.header import Header
from cortexpy.graph.parser.random_access import RandomAccess
from cortexpy.graph.parser.streaming import kmer_generator_from_stream
from cortexpy.test.builder.graph.body import KmerRecord, as_edge_set
import logging

logger = logging.getLogger(__name__)


class TestHeaderFromStream(object):
    def test_parses_a_graph_header(self, tmpdir):
        # given
        sample_name = 'sample_0'
        dna_sequence = 'ACGTT'
        kmer_size = 3

        mc_builder = (builder.Mccortex()
                      .with_dna_sequence(dna_sequence, name=sample_name)
                      .with_kmer_size(kmer_size))

        output_graph = mc_builder.build(tmpdir)

        expected_header_attributes = {'version': 6,
                                      'kmer_size': kmer_size,
                                      'record_size': 13,
                                      'kmer_container_size': 1,
                                      'num_colors': 1,
                                      'mean_read_lengths': (len(dna_sequence),),
                                      'total_sequences': (len(dna_sequence),),
                                      'sample_names': (sample_name.encode(),)}

        # when
        header = Header.from_stream(open(output_graph, 'rb'))

        # then
        for key, value in expected_header_attributes.items():
            assert getattr(header, key) == value


class TestKmerGeneratorFromStream(object):
    def test_parses_a_graph(self, tmpdir):
        # given
        kmer_size = 3
        mc_builder = (builder.Mccortex()
                      .with_dna_sequence('ACGTT')
                      .with_kmer_size(kmer_size))

        expected_kmers = [
            KmerRecord('AAC', (1,), [as_edge_set('......G.')]),
            KmerRecord('ACG', (2,), [as_edge_set('a......T')]),
        ]

        # when
        output_graph = mc_builder.build(tmpdir)

        kmer_generator = kmer_generator_from_stream(open(output_graph, 'rb'))

        # then
        actual_kmers = list(kmer_generator)
        for kmer in actual_kmers:
            logger.info(kmer)
        for expected_kmer, kmer in zip(expected_kmers, actual_kmers):
            assert kmer.kmer == expected_kmer.kmer
            assert kmer.coverage == expected_kmer.coverage
            assert kmer.edges == expected_kmer.edges
        assert len(actual_kmers) == len(expected_kmers)

    def test_parses_a_graph_with_kmer_size_9(self, tmpdir):
        # given
        kmer_size = 9
        mc_builder = (builder.Mccortex()
                      .with_dna_sequence('ACGTTCCCC')
                      .with_kmer_size(kmer_size))

        expected_kmers = [
            KmerRecord('ACGTTCCCC', (1,), [as_edge_set('........')]),
        ]

        # when
        output_graph = mc_builder.build(tmpdir)

        kmer_generator = kmer_generator_from_stream(open(output_graph, 'rb'))

        # then
        actual_kmers = list(kmer_generator)
        for kmer in actual_kmers:
            logger.info(kmer)
        for expected_kmer, kmer in zip(expected_kmers, actual_kmers):
            assert kmer.kmer == expected_kmer.kmer
            assert kmer.coverage == expected_kmer.coverage
            assert kmer.edges == expected_kmer.edges

    def test_parses_a_graph_with_kmer_size_32(self, tmpdir):
        # given
        kmer_size = 33
        contig = ''.join(list(repeat('A', kmer_size)))
        mc_builder = (builder.Mccortex()
                      .with_dna_sequence(contig)
                      .with_kmer_size(kmer_size))

        expected_kmers = [
            KmerRecord(contig, (1,), [as_edge_set('........')]),
        ]

        # when
        output_graph = mc_builder.build(tmpdir)

        kmer_generator = kmer_generator_from_stream(open(output_graph, 'rb'))

        # then
        actual_kmers = list(kmer_generator)
        for kmer in actual_kmers:
            logger.info(kmer)
        for expected_kmer, kmer in zip(expected_kmers, actual_kmers):
            assert kmer.kmer == expected_kmer.kmer
            assert kmer.coverage == expected_kmer.coverage
            assert kmer.edges == expected_kmer.edges


class TestRandomAccess(object):
    def test_retrieves_kmer_by_random_access(self, tmpdir):
        # given
        kmer_size = 3
        output_graph = (builder.Mccortex()
                        .with_dna_sequence('ACGTTT')
                        .with_kmer_size(kmer_size)
                        .build(tmpdir))

        expected = KmerRecord('AAC', (1,), [as_edge_set('A.....G.')])
        cg = RandomAccess(open(output_graph, 'rb'))

        # when
        actual = cg['AAC']

        # then
        logger.info(actual)

        assert actual.kmer == expected.kmer
        assert actual.coverage == expected.coverage
        assert actual.edges == expected.edges
