import struct
from struct import unpack
from collections import namedtuple
from pycortex.graph.parser.constants import CORTEX_MAGIC_WORD, CORTEX_VERSION, UINT8_T, UINT32_T, \
    UINT64_T


class HeaderParserError(Exception):
    """This error is thrown if the header is not parseable."""


Header = namedtuple('Header', [
    'version',
    'kmer_size',
    'kmer_container_size',
    'num_colors',
    'mean_read_lengths',
    'mean_total_sequence',
    'sample_names',
    'record_size',
])


def header_from_stream(stream):
    magic_word = unpack('cccccc', stream.read(6))

    if magic_word != CORTEX_MAGIC_WORD:
        raise HeaderParserError(
            "Saw magic word {} but was expecting {}".format(magic_word, CORTEX_MAGIC_WORD))

    header = unpack('IIII', stream.read(16))

    version = header[0]
    if version != CORTEX_VERSION:
        raise HeaderParserError(
            "Saw version {} but was expecting {}".format(version, CORTEX_VERSION)
        )

    kmer_size = header[1]
    if kmer_size <= 0:
        raise HeaderParserError(
            "Saw kmer size {} but was expecting value > 0".format(kmer_size)
        )

    kmer_container_size = header[2]
    if kmer_container_size <= 0:
        raise HeaderParserError(
            "Saw kmer bits {} but was expecting value > 0".format(kmer_size)
        )

    num_colors = header[3]
    if num_colors <= 0:
        raise HeaderParserError(
            "Saw number of colors {} but was expecting value > 0".format(kmer_size)
        )

    mean_read_lengths = unpack(
        '{}I'.format(num_colors), stream.read(struct.calcsize('I') * num_colors)
    )

    mean_total_sequence = unpack(
        '{}L'.format(num_colors), stream.read(struct.calcsize('L') * num_colors)
    )

    sample_names = []
    for _ in range(num_colors):
        sample_name_length_string = stream.read(struct.calcsize('I'))
        snlength = unpack('I', sample_name_length_string)[0]
        sample_name = unpack('{}c'.format(snlength), stream.read(snlength))
        sample_names.append(b''.join(sample_name))
    sample_names = tuple(sample_names)

    _ = unpack('16c', stream.read(16))  # error_rate

    for _ in range(num_colors):
        color_info_block_string = stream.read(4 + 3 * struct.calcsize('I'))
        color_info_block = unpack('ccccIII', color_info_block_string)
        stream.read(color_info_block[6])

    concluding_magic_word = unpack('cccccc', stream.read(6))

    if concluding_magic_word != magic_word:
        raise HeaderParserError(
            "Concluding magic word {} != starting magic word {}".format(concluding_magic_word,
                                                                        magic_word))

    record_size = UINT64_T * kmer_container_size + (UINT32_T + UINT8_T) * num_colors
    return Header(version=version,
                  kmer_size=kmer_size,
                  kmer_container_size=kmer_container_size,
                  num_colors=num_colors,
                  mean_read_lengths=mean_read_lengths,
                  mean_total_sequence=mean_total_sequence,
                  sample_names=sample_names,
                  record_size=record_size)
