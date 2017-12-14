import struct
from struct import unpack
import attr

from cortexpy.graph.parser.constants import CORTEX_MAGIC_WORD, CORTEX_VERSION, UINT8_T, UINT32_T, \
    UINT64_T


def none_or_greater_than_zero(_, attribute, value):
    if value is not None and value <= 0:
        raise ValueError("'{}' has to be greater than 0!".format(attribute.name))


@attr.s(slots=True, frozen=True)
class Header(object):
    version = attr.ib(CORTEX_VERSION, validator=[attr.validators.in_([CORTEX_VERSION])])
    kmer_size = attr.ib(None, validator=[none_or_greater_than_zero])
    kmer_container_size = attr.ib(None, validator=[none_or_greater_than_zero])
    num_colors = attr.ib(None, validator=[none_or_greater_than_zero])
    mean_read_lengths = attr.ib(None)
    mean_total_sequence = attr.ib(None)
    sample_names = attr.ib(None)
    error_rate = attr.ib(None)
    color_info_blocks = attr.ib(attr.Factory(list))

    @property
    def record_size(self):
        return UINT64_T * self.kmer_container_size + (UINT32_T + UINT8_T) * self.num_colors

    @property
    def colors(self):
        return tuple(range(self.num_colors))


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
        self.header = attr.evolve(self.header,
                                  version=params[0],
                                  kmer_size=params[1],
                                  kmer_container_size=params[2],
                                  num_colors=params[3])
        return self

    def extract_mean_read_lengths(self):
        self.header = attr.evolve(
            self.header,
            mean_read_lengths=unpack(
                '{}I'.format(self.header.num_colors),
                self.stream.read(struct.calcsize('I') * self.header.num_colors)
            )
        )
        return self

    def extract_mean_total_sequence(self):
        self.header = attr.evolve(
            self.header,
            mean_total_sequence=unpack(
                '{}L'.format(self.header.num_colors),
                self.stream.read(struct.calcsize('L') * self.header.num_colors)
            )
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
        self.header = attr.evolve(self.header, sample_names=tuple(sample_names))
        return self

    def extract_error_rate(self):
        error_rates = []
        for _ in range(self.header.num_colors):
            error_rates.append(unpack('16c', self.stream.read(16)))
        self.header = attr.evolve(self.header, error_rate=tuple(error_rates))
        return self

    def extract_color_info_blocks(self):
        for _ in range(self.header.num_colors):
            color_info_block_string = self.stream.read(4 + 3 * struct.calcsize('I'))
            color_info_block = unpack('4c3I', color_info_block_string)
            cleaned_graph_name = self.stream.read(color_info_block[6])
            self.header.color_info_blocks.append((color_info_block, cleaned_graph_name))
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
