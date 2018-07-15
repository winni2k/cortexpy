from cortexpy.graph.cortex import CortexGraphMapping
from cortexpy.graph.parser.kmer import EmptyKmerBuilder
from cortexpy.test.builder.graph.cortex import get_cortex_graph_mapping_builder


class Test(object):
    def test_loads_kmer(self):
        # given
        b = get_cortex_graph_mapping_builder()
        b.with_kmer('AAA 0 ........')

        # when
        cgm = b.build()

        # then
        assert type(cgm) == CortexGraphMapping
        assert 'AAA' in cgm
        assert [0] == list(cgm['AAA'].coverage)

    def test_deletes_kmer(self):
        # given
        b = get_cortex_graph_mapping_builder()
        b.with_kmer('AAA 0 ........')

        # when
        cgm = b.build()
        del cgm['AAA']

        # then
        assert 'AAA' not in cgm
        assert 0 == len(cgm)

    def test_sets_kmer(self):
        # given
        b = get_cortex_graph_mapping_builder()
        b.with_kmer_size(3)
        kmer_builder = EmptyKmerBuilder()

        # when
        cgm = b.build()
        cgm['AAA'] = kmer_builder.build('AAA')

        # then
        assert 'AAA' in cgm
        assert 1 == len(cgm)

    def test_overwrites_preexisting_kmer(self):
        # given
        kmer_builder = EmptyKmerBuilder()
        b = get_cortex_graph_mapping_builder()
        b.with_kmer('AAA 0 ........')

        # when
        cgm = b.build()
        cgm['AAA'] = kmer_builder.build('AAA')

        # then
        assert 1 == len(cgm)
        assert [] == list(cgm['AAA'].coverage)
