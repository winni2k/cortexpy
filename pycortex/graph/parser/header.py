import struct
import attr
from struct import unpack

from pycortex.graph.parser.constants import CORTEX_MAGIC_WORD, CORTEX_VERSION, UINT8_T, UINT32_T, \
    UINT64_T


@attr.s(slots=True)
class Header(object):
    _version = attr.ib(CORTEX_VERSION)
    _kmer_size = attr.ib(None)
    _kmer_container_size = attr.ib(None)
    _num_colors = attr.ib(None)
    mean_read_lengths = attr.ib(None)
    mean_total_sequence = attr.ib(None)
    sample_names = attr.ib(None)

    @property
    def record_size(self):
        return UINT64_T * self.kmer_container_size + (UINT32_T + UINT8_T) * self.num_colors

    @property
    def version(self):
        return self._version

    @version.setter
    def version(self, value):
        if value != CORTEX_VERSION:
            raise ValueError('Version is not 6')
        self._version = value

    @property
    def kmer_size(self):
        return self._kmer_size

    @kmer_size.setter
    def kmer_size(self, value):
        if value <= 0:
            raise ValueError('Kmer size < 1')
        self._kmer_size = value

    @property
    def kmer_container_size(self):
        return self._kmer_container_size

    @kmer_container_size.setter
    def kmer_container_size(self, value):
        if value <= 0:
            raise ValueError('Kmer container size < 1')
        self._kmer_container_size = value

    @property
    def num_colors(self):
        return self._num_colors

    @num_colors.setter
    def num_colors(self, value):
        if value <= 0:
            raise ValueError("Number of colors < 1")
        self._num_colors = value


class HeaderParserError(Exception):
    """This error is thrown if the header is not parseable.

    Returns a partially filled header object for inspection.
    """


@attr.s(slots=True)
class HeaderFromStreamBuilder(object):
    header = attr.ib(attr.Factory(Header))

    def fill_first_four_params(self, binary_string):
        params = unpack('4I', binary_string)
        self.header.version = params[0]
        self.header.kmer_size = params[1]
        self.header.kmer_container_size = params[2]
        self.header.num_colors = params[3]

    def fill_sample_names_from_stream(self, stream):
        sample_names = []
        for _ in range(self.header.num_colors):
            sample_name_length_string = stream.read(struct.calcsize('I'))
            sample_name_length = unpack('I', sample_name_length_string)[0]
            sample_name = unpack('{}c'.format(sample_name_length), stream.read(sample_name_length))
            sample_names.append(b''.join(sample_name))
        self.header.sample_names = tuple(sample_names)


def from_stream(stream):
    magic_word = unpack('cccccc', stream.read(6))
    if magic_word != CORTEX_MAGIC_WORD:
        raise HeaderParserError(
            "Saw magic word {} but was expecting {}".format(magic_word, CORTEX_MAGIC_WORD)
        )

    builder = HeaderFromStreamBuilder()
    builder.fill_first_four_params(stream.read(16))
    header = builder.header

    header.mean_read_lengths = unpack(
        '{}I'.format(header.num_colors), stream.read(struct.calcsize('I') * header.num_colors)
    )

    header.mean_total_sequence = unpack(
        '{}L'.format(header.num_colors), stream.read(struct.calcsize('L') * header.num_colors)
    )

    builder.fill_sample_names_from_stream(stream)
    header = builder.header

    _ = unpack('16c', stream.read(16))  # error_rate

    for _ in range(header.num_colors):
        color_info_block_string = stream.read(4 + 3 * struct.calcsize('I'))
        color_info_block = unpack('4c3I', color_info_block_string)
        stream.read(color_info_block[6])

    concluding_magic_word = unpack('6c', stream.read(6))

    if concluding_magic_word != magic_word:
        raise HeaderParserError(
            'Concluding magic word {} != starting magic word {}\nUnparsed: {}\nParsed header: {}'.format(
                concluding_magic_word,
                magic_word,
                stream.read(1000),
                header
            )
        )

    return header
