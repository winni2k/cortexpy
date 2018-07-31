from hypothesis import given, strategies as strat

from cortexpy.graph.interactor import Interactor
from cortexpy.test.builder.graph.cortex import get_cortex_builder
from cortexpy.test.builder.graph.kmer import kmer_strings
from cortexpy.test.expectation import KmerGraphExpectation
from cortexpy.utils import lexlo


def test_revcomps_a_kmer():
    # given
    b = get_cortex_builder()
    b.with_kmer('AAA 1 ........')
    cdb = b.build()

    # when
    expect = KmerGraphExpectation(
        Interactor(cdb).make_graph_nodes_consistent({'TTT'}).graph)

    # then
    expect.has_node('TTT')
    expect.has_n_nodes(1)


@given(strat.data(),
       strat.integers(min_value=1, max_value=5),
       strat.integers(min_value=3, max_value=7))
def test_revcomps_many_kmers(data, num_kmers, kmer_size):
    # given
    kmers = {}
    for _ in range(num_kmers):
        kmer_string = data.draw(kmer_strings(min_size=kmer_size, max_size=kmer_size))
        kmers[lexlo(kmer_string)] = kmer_string

    b = get_cortex_builder()
    for kmer in kmers.keys():
        b.with_kmer('{} 1 ........'.format(kmer))
    cdb = b.build()

    # when
    expect = KmerGraphExpectation(
        Interactor(cdb).make_graph_nodes_consistent(set(kmers.values())).graph)

    # then
    for kmer_string in kmers.values():
        expect.has_node(kmer_string)
    expect.has_n_nodes(len(kmers))


def test_revcomps_path():
    # given
    b = get_cortex_builder()
    b.with_kmer('CGC 1 .......T')
    b.with_kmer('AGC 1 ......G.')

    cdb = b.build()

    for seed, expected_nodes in [('CGC', ['CGC', 'GCT']),
                                 ('GCT', ['CGC', 'GCT']),
                                 ('AGC', ['AGC', 'GCG']),
                                 ('GCG', ['AGC', 'GCG'])]:
        # when
        expect = KmerGraphExpectation(
            Interactor(cdb).make_graph_nodes_consistent([seed]).graph)

        # then
        expect.has_nodes(*expected_nodes)
        expect.has_n_nodes(2)


def test_keys_y_graph():
    # given
    b = get_cortex_builder()
    b.with_kmer('CGC 1 .......T')
    b.with_kmer('AGC 1 a....CG.')
    b.with_kmer('AAG 1 .....C..')
    b.with_kmer('GCC 1 a.......')

    expected_nodes1 = ['CGC', 'GCT', 'CTT', 'GGC']
    expected_nodes2 = ['AAG', 'AGC', 'GCG', 'GCC']
    for expected_nodes in [expected_nodes1, expected_nodes2]:
        for seed in expected_nodes:
            cdb = b.build()

            # when
            expect = KmerGraphExpectation(
                Interactor(cdb).make_graph_nodes_consistent([seed]).graph)

            # then
            expect.has_nodes(*expected_nodes)
