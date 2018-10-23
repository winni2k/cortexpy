import gzip
import logging
from itertools import repeat

import pytest

import cortexpy.test.builder.mccortex as builder
from cortexpy.graph.parser.header import Header
from cortexpy.links import Links
from cortexpy.graph.parser.random_access import RandomAccess
from cortexpy.graph.parser.streaming import kmer_generator_from_stream
from cortexpy.test.builder.graph.body import KmerRecord, as_edge_set
from cortexpy.test.expectation.links import LinksExpectation

logger = logging.getLogger(__name__)


class TestHeaderFromStream:
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


class TestKmerGeneratorFromStream:
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


class TestRandomAccess:
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


class TestLinkParsing:
    def test_turner_et_al_example(self, tmpdir):
        # given
        b = builder.MccortexGraphLinks()
        b.with_dna_sequence('ACTGATTTCGATGCGATGCGATGCCACGGTGG')
        b.with_kmer_size(5)
        b.with_link_dna_sequence('TTTCGATGCGATGCGATGCCACG')

        # when
        expect = LinksExpectation(Links.from_binary_stream(gzip.open(b.build(tmpdir)[1], 'rb')))

        # then
        expect.has_link_group_for_kmer('ATCGA').has_links('R 3 1 GGC')
        expect.has_link_group_for_kmer('ATCGC').has_links('R 2 1 GC', 'R 1 1 C')
        expect.has_link_group_for_kmer('ATGCG').has_links('R 2 1 CA', 'R 1 1 A')
        expect.has_link_group_for_kmer('ATGCC').has_links('R 3 1 CCA')
        expect.has_n_link_groups(4)

    def test_simple_tangle_has_four_links(self, tmpdir):
        # given
        b = builder.MccortexGraphLinks()
        b.with_kmer_size(5)
        b.with_dna_sequence('CAAAACCCCC')
        b.with_dna_sequence('TAAAACCCCT')
        b.with_link_dna_sequence('CAAAACCCCT')
        b.with_link_dna_sequence('TAAAACCCCC')

        # when
        expect = LinksExpectation(Links.from_binary_stream(gzip.open(b.build(tmpdir)[1], 'rb')))

        # then
        expect.has_link_group_for_kmer('CAAAA').has_links('F 1 1 T')
        expect.has_link_group_for_kmer('TAAAA').has_links('F 1 1 C')
        expect.has_link_group_for_kmer('CCCCC').has_links('R 1 1 A')
        expect.has_link_group_for_kmer('AGGGG').has_links('F 1 1 G')
        expect.has_n_link_groups(4)

    def test_tangle_with_starting_non_lexlo_kmer_has_four_links(self, tmpdir):
        # given
        b = builder.MccortexGraphLinks()
        b.with_kmer_size(5)
        b.with_dna_sequence('CAAAACCCCC')
        b.with_dna_sequence('TTTAACCCCT')
        b.with_link_dna_sequence('CAAAACCCCT')
        b.with_link_dna_sequence('TTTAACCCCC')

        # when
        expect = LinksExpectation(Links.from_binary_stream(gzip.open(b.build(tmpdir)[1], 'rb')))

        # then
        expect.has_link_group_for_kmer('AAACC').has_links('F 1 1 T')
        expect.has_link_group_for_kmer('GGTTA').has_links('R 1 1 C')
        expect.has_link_group_for_kmer('CCCCC').has_links('R 1 1 A')
        expect.has_link_group_for_kmer('AGGGG').has_links('F 1 1 T')
        expect.has_n_link_groups(4)

    @pytest.mark.xfail(reason='Not implemented in mccortex')
    def test_bubble_has_four_links(self, tmpdir):
        # given
        b = builder.MccortexGraphLinks()
        b.with_kmer_size(5)
        b.with_dna_sequence('AAAAACAACCC')
        b.with_dna_sequence('AAAAATAACCC')
        b.with_link_dna_sequence('AAAAACAACCC')
        b.with_link_dna_sequence('AAAAATAACCC')

        # when
        expect = LinksExpectation(Links.from_binary_stream(gzip.open(b.build(tmpdir)[1], 'rb')))

        # then
        expect.has_link_group_for_kmer('AAAAA').has_links('F 1 1 C', 'F 1 1 T')
        expect.has_link_group_for_kmer('AACCC').has_links('R 1 1 G', 'R 1 1 A')

    @pytest.mark.xfail(reason='Not implemented in mccortex')
    def test_two_bubbles_have_two_links_traversing_them(self, tmpdir):
        # given
        b = builder.MccortexGraphLinks()
        b.with_kmer_size(5)
        b.with_dna_sequence('AAAAACAACCCACCCCC')
        b.with_dna_sequence('AAAAAGAACCCTCCCCC')
        b.with_link_dna_sequence('AAAAACAACCCTCCCCC')
        b.with_link_dna_sequence('AAAAAGAACCCACCCCC')

        # when
        expect = LinksExpectation(Links.from_binary_stream(gzip.open(b.build(tmpdir)[1], 'rb')))

        # then
        expect.has_link_group_for_kmer('AAAAA').has_links('F 2 1 CT', 'F 2 1 GA')
