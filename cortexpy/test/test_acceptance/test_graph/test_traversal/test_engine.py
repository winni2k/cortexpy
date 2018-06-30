from cortexpy.graph import traversal
from cortexpy.graph.parser import RandomAccess
from cortexpy.graph.traversal import EngineTraversalOrientation
from cortexpy.test import builder
from cortexpy.test.expectation import KmerGraphExpectation


class TestTraverseFrom(object):
    def test_with_bubble_and_two_colors_returns_all_kmers(self, tmpdir):
        # given
        kmer_size = 3
        output_graph = (builder.Mccortex(kmer_size)
                        .with_dna_sequence('AAACAAG')
                        .with_dna_sequence('AAATAAG')
                        .with_dna_sequence('AAATAAG', name='sample_1')
                        .build(tmpdir))

        traverser = traversal.Engine(
            RandomAccess(open(output_graph, 'rb')),
            traversal_colors=(0,),
            orientation=EngineTraversalOrientation.both
        )

        # when
        expect = KmerGraphExpectation(traverser.traverse_from('ACA').graph, sort_edges=True)

        # then
        expect.has_node('AAA').has_coverages(2, 1)
        expect.has_node('AAC').has_coverages(1, 0)
        expect.has_node('ACA').has_coverages(1, 0)
        expect.has_node('CAA').has_coverages(1, 0)
        expect.has_node('AAG').has_coverages(2, 1)
        expect.has_node('AAT').has_coverages(1, 1)
        expect.has_node('ATA').has_coverages(1, 1)
        expect.has_node('TAA').has_coverages(1, 1)

        expect.has_edges(
            'AAA AAC 0',
            'AAC ACA 0',
            'ACA CAA 0',
            'CAA AAG 0',
            'AAA AAT 0',
            'AAT ATA 0',
            'ATA TAA 0',
            'TAA AAG 0',
            'AAA AAT 1',
            'AAT ATA 1',
            'ATA TAA 1',
            'AAG TAA 1')


class TestTraverseFromEachKmerInIterable(object):
    def test_with_two_subgraphs_returns_all_kmers(self, tmpdir):
        # given
        kmer_size = 3
        output_graph = (builder.Mccortex(kmer_size)
                        .with_dna_sequence('AAAT')
                        .with_dna_sequence('GGGC')
                        .build(tmpdir))

        traverser = traversal.Engine(
            RandomAccess(open(output_graph, 'rb')),
            traversal_colors=(0,),
            orientation=EngineTraversalOrientation.both
        )

        # when
        expect = KmerGraphExpectation(
            traverser.traverse_from_each_kmer_in_iterable(['AAA', 'GGG']).graph
        )

        # then
        expect.has_edges('AAA AAT 0', 'CCC GCC 0')
