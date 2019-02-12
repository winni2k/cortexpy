import numpy as np
from hypothesis import strategies as s, assume

from cortexpy.edge_set import EdgeSet
from cortexpy.test.builder.graph.body import KmerRecord
from cortexpy.test.constants import MAX_UINT


@s.composite
def dna_sequences(draw, min_size=0, max_size=7):
    return draw(s.text(alphabet='ACGT', min_size=min_size, max_size=max_size))


@s.composite
def kmer_strings(draw, min_size=3, max_size=7):
    assert min_size >= 3
    dna_seq = draw(dna_sequences(min_size=min_size, max_size=max_size))
    assume(len(dna_seq) % 2 == 1)
    return dna_seq


@s.composite
def kmer_records(draw, kmer_size, num_colors, kmer_strings=dna_sequences):
    kmer = draw(kmer_strings(min_size=kmer_size, max_size=kmer_size))
    coverage = tuple(
        draw(s.lists(s.integers(min_value=1, max_value=MAX_UINT), min_size=num_colors, max_size=num_colors)))
    edges = np.array(
        draw(
            s.lists(
                s.lists(s.integers(min_value=0, max_value=1), min_size=8, max_size=8),
                min_size=num_colors,
                max_size=num_colors)
        ),
        dtype=np.uint8
    )
    edges = [EdgeSet(np.concatenate((e[:4], e[::-1][:4]))) for e in edges]
    return KmerRecord(kmer, coverage, edges)
