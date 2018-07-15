from hypothesis import given, strategies

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

    @given(strategies.booleans())
    def test_disconnects_two_kmers(self, kmer_cache_off):
        # given
        b = get_cortex_graph_mapping_builder()
        if kmer_cache_off:
            b.with_kmer_cache_size(0)
        b.with_kmer('AAA 0 .....C..')
        b.with_kmer('AAC 0 a.......')

        # when
        cgm = b.build()
        cgm.disconnect_kmers(cgm['AAA'], cgm['AAC'], [0])

        # then
        assert 2 == len(cgm)
        assert '........' == cgm['AAA'].edges[0].to_str()
        assert '........' == cgm['AAC'].edges[0].to_str()

    @given(strategies.booleans())
    def test_connects_two_kmers(self, kmer_cache_off):
        # given
        b = get_cortex_graph_mapping_builder()
        if kmer_cache_off:
            b.with_kmer_cache_size(0)
        b.with_kmer('AAA 0 ........')
        b.with_kmer('AAC 0 ........')

        # when
        cgm = b.build()
        cgm.connect_kmers(cgm['AAA'], cgm['AAC'], 0)

        # then
        assert 2 == len(cgm)
        assert '.....C..' == cgm['AAA'].edges[0].to_str()
        assert 'a.......' == cgm['AAC'].edges[0].to_str()
