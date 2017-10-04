from pycortex.graph.parser.constants import UINT64_T
from pycortex.graph.parser.header import header_from_stream
from pycortex.kmer import Kmer


def kmer_generator_from_stream(stream):
    header = header_from_stream(stream)
    return kmer_generator_from_stream_and_header(stream, header)


def kmer_generator_from_stream_and_header(stream, header):
    record_size = header.kmer_container_size * UINT64_T + 5 * header.num_colors

    raw_record = stream.read(record_size)
    while raw_record != b'':
        yield Kmer(raw_record,
                   header.kmer_size,
                   header.num_colors,
                   header.kmer_container_size)
        raw_record = stream.read(record_size)
