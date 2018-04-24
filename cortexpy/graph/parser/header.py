import struct
from struct import unpack
import attr

from cortexpy.graph.parser.constants import (
    CORTEX_MAGIC_WORD, CORTEX_VERSION, UINT8_T, UINT32_T,
    UINT64_T,
    ERROR_RATE_SIZE,
)


def none_or_greater_than_zero(_, attribute, value):
    if value is not None and value <= 0:
        raise ValueError("'{}' has to be greater than 0!".format(attribute.name))


def greater_than_zero(_, attribute, value):
    if value < 1:
        raise ValueError("'{}' has to be greater than 0!".format(attribute.name))


def odd(_, attribute, value):
    if value % 2 == 0:
        raise ValueError("'{}' has to be odd!".format(attribute.name))


@attr.s(slots=True, frozen=True)
class Header(object):
    version = attr.ib(CORTEX_VERSION, validator=[attr.validators.in_([CORTEX_VERSION])])
    kmer_size = attr.ib(1, validator=[greater_than_zero, odd])
    kmer_container_size = attr.ib(None, validator=[none_or_greater_than_zero])
    num_colors = attr.ib(None, validator=[none_or_greater_than_zero])
    _mean_read_lengths = attr.ib(None)
    _total_sequences = attr.ib(None)
    sample_names = attr.ib(None)
    _error_rates = attr.ib(None)
    color_info_blocks = attr.ib(attr.Factory(list))

    @property
    def record_size(self):
        return UINT64_T * self.kmer_container_size + (UINT32_T + UINT8_T) * self.num_colors

    @property
    def colors(self):
        return tuple(range(self.num_colors))

    @property
    def mean_read_lengths(self):
        if self._mean_read_lengths is None:
            return [0 for _ in range(self.num_colors)]
        return self._mean_read_lengths

    @property
    def total_sequences(self):
        if self._total_sequences is None:
            return [0 for _ in range(self.num_colors)]
        return self._total_sequences

    @property
    def error_rates(self):
        if self._error_rates is None:
            empty_bytes = bytes([0 for _ in range(ERROR_RATE_SIZE)])
            return [empty_bytes for _ in range(self.num_colors)]
        return self._error_rates

    def dump(self, buffer):
        magic_word = b''.join(CORTEX_MAGIC_WORD)
        assert struct.calcsize('L') == UINT64_T
        assert struct.calcsize('I') == UINT32_T
        buffer.write(magic_word)
        buffer.write(struct.pack('4I',
                                 self.version,
                                 self.kmer_size,
                                 self.kmer_container_size,
                                 self.num_colors))
        for mean_read_length in self.mean_read_lengths:
            buffer.write(struct.pack('I', mean_read_length))
        for total_sequence in self.total_sequences:
            buffer.write(struct.pack('L', total_sequence))
        for sample_name in self.sample_names:
            if isinstance(sample_name, str):
                sample_name = sample_name.encode()
            buffer.write(struct.pack('I', len(sample_name)))
            buffer.write(sample_name)
        for error_rate in self.error_rates:
            buffer.write(error_rate)
        for info_block in self.color_info_blocks:
            info_block.dump(buffer)
        buffer.write(magic_word)


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

    def extract_total_sequences(self):
        assert struct.calcsize('L') == UINT64_T
        self.header = attr.evolve(
            self.header,
            total_sequences=unpack(
                '{}L'.format(self.header.num_colors),
                self.stream.read(UINT64_T * self.header.num_colors)
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

    def extract_error_rates(self):
        error_rates = []
        for _ in range(self.header.num_colors):
            error_rates.append(self.stream.read(16))
        self.header = attr.evolve(self.header, error_rates=tuple(error_rates))
        return self

    def extract_color_info_blocks(self):
        for _ in range(self.header.num_colors):
            color_info_block_string = self.stream.read(4 + 3 * struct.calcsize('I'))
            color_info_block = unpack('4c3I', color_info_block_string)
            cleaned_graph_name = self.stream.read(color_info_block[6])
            self.header.color_info_blocks.append((color_info_block, cleaned_graph_name))
        return self


def from_stream(stream):
    header = (HeaderFromStreamBuilder(stream)
              .extract_initial_magic_word()
              .fill_first_four_params()
              .extract_mean_read_lengths()
              .extract_total_sequences()
              .extract_sample_names_from_stream()
              .extract_error_rates()
              .extract_color_info_blocks()
              .extract_concluding_magic_word())
    return header
