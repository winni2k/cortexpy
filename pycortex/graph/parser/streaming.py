import attr

from pycortex.graph.parser.constants import UINT64_T
from pycortex.graph.parser.header import Header
from pycortex.kmer import Kmer


@attr.s(slots=True)
class Streaming(object):
    fh = attr.ib()
    header = attr.ib(init=False)

    def __attrs_post_init__(self):
        self.header = Header.from_stream(self.fh)

    def kmers(self):
        return kmer_generator_from_stream(self.fh, self.header)


def kmer_generator_from_stream(stream, cortex_header):
    record_size = cortex_header.kmer_container_size * UINT64_T + 5 * cortex_header.num_colors

    raw_record = stream.read(record_size)
    while raw_record != b'':
        yield Kmer(raw_record,
                   cortex_header.kmer_size,
                   cortex_header.num_colors,
                   cortex_header.kmer_container_size)
        raw_record = stream.read(record_size)
