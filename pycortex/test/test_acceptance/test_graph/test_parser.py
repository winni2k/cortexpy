from itertools import repeat

import pycortex.graph.parser as parser
import pycortex.test.builder as builder
from pycortex.test import runner
from pycortex.test.builder.graph.body import KmerRecord, as_edge_set, print_kmer


class TestHeaderFromStream(object):
    def test_parses_a_graph_header(self, tmpdir):
        # given
        sample_name = b'sample_0'
        dna_sequence = 'ACGTT'
        kmer_size = 3

        mc_builder = (builder.Mccortex()
                      .with_dna_sequence(sample_name, dna_sequence)
                      .with_kmer_size(kmer_size))

        expected_header = parser.Header(version=6,
                                        kmer_size=kmer_size,
                                        record_size=13,
                                        kmer_container_size=1,
                                        num_colors=1,
                                        mean_read_lengths=(len(dna_sequence),),
                                        mean_total_sequence=(len(dna_sequence),),
                                        sample_names=(sample_name,))

        # when
        output_graph = mc_builder.build(tmpdir)

        runner.Mccortex(kmer_size).view(output_graph)
        header = parser.header.from_stream(open(output_graph, 'rb'))

        # then
        assert header == expected_header


class TestKmerGeneratorFromStream(object):
    def test_parses_a_graph(self, tmpdir):
        # given
        kmer_size = 3
        mc_builder = (builder.Mccortex()
                      .with_dna_sequence(b'sample_0', 'ACGTT')
                      .with_kmer_size(kmer_size))

        expected_kmers = [
            KmerRecord('AAC', (1,), (as_edge_set('......G.'),)),
            KmerRecord('ACG', (2,), (as_edge_set('A......T'),)),
        ]

        # when
        output_graph = mc_builder.build(tmpdir)

        runner.Mccortex(kmer_size).view(output_graph)
        kmer_generator = parser.kmer_generator_from_stream(open(output_graph, 'rb'))

        # then
        actual_kmers = list(kmer_generator)
        for kmer in actual_kmers:
            print_kmer(kmer)
        for expected_kmer, kmer in zip(expected_kmers, actual_kmers):
            assert kmer.kmer == expected_kmer.kmer
            assert kmer.coverage == expected_kmer.coverage
            assert kmer.edges == expected_kmer.edges

    def test_parses_a_graph_with_kmer_size_9(self, tmpdir):
        # given
        kmer_size = 9
        mc_builder = (builder.Mccortex()
                      .with_dna_sequence(b'sample_0', 'ACGTTCCCC')
                      .with_kmer_size(kmer_size))

        expected_kmers = [
            KmerRecord('ACGTTCCCC', (1,), (as_edge_set('........'),)),
        ]

        # when
        output_graph = mc_builder.build(tmpdir)

        runner.Mccortex(kmer_size).view(output_graph)
        kmer_generator = parser.kmer_generator_from_stream(open(output_graph, 'rb'))

        # then
        actual_kmers = list(kmer_generator)
        for kmer in actual_kmers:
            print_kmer(kmer)
        for expected_kmer, kmer in zip(expected_kmers, actual_kmers):
            assert kmer.kmer == expected_kmer.kmer
            assert kmer.coverage == expected_kmer.coverage
            assert kmer.edges == expected_kmer.edges

    def test_parses_a_graph_with_kmer_size_32(self, tmpdir):
        # given
        kmer_size = 33
        contig = ''.join(list(repeat('A', kmer_size)))
        mc_builder = (builder.Mccortex()
                      .with_dna_sequence(b'sample_0', contig)
                      .with_kmer_size(kmer_size))

        expected_kmers = [
            KmerRecord(contig, (1,), (as_edge_set('........'),)),
        ]

        # when
        output_graph = mc_builder.build(tmpdir)

        runner.Mccortex(kmer_size).view(output_graph)
        kmer_generator = parser.kmer_generator_from_stream(open(output_graph, 'rb'))

        # then
        actual_kmers = list(kmer_generator)
        for kmer in actual_kmers:
            print_kmer(kmer)
        for expected_kmer, kmer in zip(expected_kmers, actual_kmers):
            assert kmer.kmer == expected_kmer.kmer
            assert kmer.coverage == expected_kmer.coverage
            assert kmer.edges == expected_kmer.edges


class TestRandomAccess(object):
    def test_retrieves_kmer_by_random_access(self, tmpdir):
        # given
        kmer_size = 3
        mc_builder = (builder.Mccortex()
                      .with_dna_sequence(b'sample_0', 'ACGTTT')
                      .with_kmer_size(kmer_size))

        expected = KmerRecord('AAC', (1,), (as_edge_set('A.....G.'),))
        output_graph = mc_builder.build(tmpdir)
        runner.Mccortex(kmer_size).view(output_graph)
        cg = parser.RandomAccess(open(output_graph, 'rb'))

        # when
        actual = cg['AAC']

        # then
        print_kmer(actual)

        assert actual.kmer == expected.kmer
        assert actual.coverage == expected.coverage
        assert actual.edges == expected.edges
