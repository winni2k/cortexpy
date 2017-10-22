import struct
from struct import unpack
import attr

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


@attr.s(slots=True)
class HeaderFromStreamBuilder(object):
    stream = attr.ib()
    header = attr.ib(attr.Factory(Header))

    def _extract_magic_word(self):
        return unpack('6c', self.stream.read(6))

    def extract_initial_magic_word(self):
        magic_word = self._extract_magic_word()
        if magic_word != CORTEX_MAGIC_WORD:
            raise ValueError(
                "Saw initial magic word {} but was expecting {}".format(magic_word,
                                                                        CORTEX_MAGIC_WORD)
            )
        return self

    def extract_concluding_magic_word(self):
        concluding_magic_word = self._extract_magic_word()
        if concluding_magic_word != CORTEX_MAGIC_WORD:
            raise ValueError(
                (
                    'Concluding magic word {} != starting magic word {}\n'
                    'Unparsed: {}\n'
                    'Parsed header: {}'
                ).format(
                    concluding_magic_word,
                    CORTEX_MAGIC_WORD,
                    self.stream.read(1000),
                    self.header
                )
            )
        return self.header

    def fill_first_four_params(self):
        params = unpack('4I', self.stream.read(16))
        self.header.version = params[0]
        self.header.kmer_size = params[1]
        self.header.kmer_container_size = params[2]
        self.header.num_colors = params[3]
        return self

    def extract_mean_read_lengths(self):
        self.header.mean_read_lengths = unpack(
            '{}I'.format(self.header.num_colors),
            self.stream.read(struct.calcsize('I') * self.header.num_colors)
        )
        return self

    def extract_mean_total_sequence(self):
        self.header.mean_total_sequence = unpack(
            '{}L'.format(self.header.num_colors),
            self.stream.read(struct.calcsize('L') * self.header.num_colors)
        )
        return self

    def extract_sample_names_from_stream(self):
        sample_names = []
        for _ in range(self.header.num_colors):
            sample_name_length_string = self.stream.read(struct.calcsize('I'))
            sample_name_length = unpack('I', sample_name_length_string)[0]
            sample_name = unpack('{}c'.format(sample_name_length),
                                 self.stream.read(sample_name_length))
            sample_names.append(b''.join(sample_name))
        self.header.sample_names = tuple(sample_names)
        return self

    def extract_error_rate(self):
        unpack('16c', self.stream.read(16))  # error_rate
        return self

    def extract_color_info_blocks(self):
        for _ in range(self.header.num_colors):
            color_info_block_string = self.stream.read(4 + 3 * struct.calcsize('I'))
            color_info_block = unpack('4c3I', color_info_block_string)
            self.stream.read(color_info_block[6])
        return self


def from_stream(stream):
    return (HeaderFromStreamBuilder(stream)
            .extract_initial_magic_word()
            .fill_first_four_params()
            .extract_mean_read_lengths()
            .extract_mean_total_sequence()
            .extract_sample_names_from_stream()
            .extract_error_rate()
            .extract_color_info_blocks()
            .extract_concluding_magic_word())
