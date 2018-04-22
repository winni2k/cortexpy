import attr
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import struct

from cortexpy.graph.parser import Header
from cortexpy.graph.parser.kmer import calc_kmer_container_size


@attr.s(slots=True)
class KmerGraph(object):
    """Serializes kmer graphs."""
    graph = attr.ib()

    def to_seq_records(self):
        return (
            SeqRecord(Seq(str(node)), id=str(node_idx)) for node_idx, node in
            enumerate(self.graph.nodes())
        )


NOT_SET = 0


@attr.s(slots=True)
class ColorInformationBlock(object):
    is_clipped = attr.ib(False)
    are_low_cov_unitigs_removed = attr.ib(False)
    are_low_cov_kmers_removed = attr.ib(False)
    is_cleaned_against_another_graph = attr.ib(False)
    cov_threshold_on_unitigs = attr.ib(0)
    cov_threshold_on_kmers = attr.ib(0)
    name_of_graph_cleaned_against = attr.ib(b'')

    def dump(self, buffer):
        assert isinstance(self.name_of_graph_cleaned_against, bytes)
        string_length = len(self.name_of_graph_cleaned_against)
        binary = struct.pack('4?3I',
                                 self.is_clipped,
                                 self.are_low_cov_unitigs_removed,
                                 self.are_low_cov_kmers_removed,
                                 self.is_cleaned_against_another_graph,
                                 self.cov_threshold_on_unitigs,
                                 self.cov_threshold_on_kmers,
                                 string_length)
        buffer.write(binary)
        buffer.write(self.name_of_graph_cleaned_against)


@attr.s(slots=True)
class Kmers(object):
    """Serializes kmers to cortex binary format"""
    kmers = attr.ib()
    kmer_size = attr.ib(None)
    sample_names = attr.ib(None)
    kmer_container_size = attr.ib(None)
    num_colors = attr.ib(None)
    n_kmers = attr.ib(init=False)

    def __attrs_post_init__(self):
        self.kmers = list(self.kmers)
        self.n_kmers = len(self.kmers)
        if self.kmer_size is None:
            self.kmer_size = self.kmers[0].kmer_size
        if self.kmer_container_size is None:
            self.kmer_container_size = calc_kmer_container_size(self.kmer_size)
        if self.num_colors is None:
            self.num_colors = self.kmers[0].num_colors
        if self.sample_names is None:
            assert self.num_colors > 0
            self.sample_names = [str(i).encode() for i in range(self.num_colors)]

    @property
    def header(self):
        color_info_blocks = [ColorInformationBlock() for _ in range(self.num_colors)]
        return Header(kmer_size=self.kmer_size,
                      kmer_container_size=self.kmer_container_size,
                      num_colors=self.num_colors,
                      sample_names=self.sample_names,
                      color_info_blocks=color_info_blocks)

    def dump(self, buffer):
        """to a filehandle"""
        self.header.dump(buffer)
        for kmer in self.kmers:
            kmer.dump(buffer)
