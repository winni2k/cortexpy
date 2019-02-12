import pytest
from hypothesis import assume, given
from hypothesis import strategies as s

import cortexpy.graph.parser.header as header
from cortexpy.test.builder.graph.header import Header, ColorInformationBlock
from cortexpy.test.constants import (
    MAX_UINT, MAX_ULONG, UINT8_T, UINT32_T,
    UINT64_T,
)


@s.composite
def color_information_blocks(draw):
    bools = [draw(s.binary(min_size=1, max_size=1)) for _ in range(4)]
    uint32_ts = [draw(s.integers(min_value=0, max_value=MAX_UINT)) for _ in range(2)]
    name_size = draw(s.integers(min_value=0, max_value=3))
    name = draw(s.binary(min_size=name_size, max_size=name_size))

    return ColorInformationBlock(*bools, *uint32_ts, name_size, name)


class TestHeaderParser(object):
    @given(s.binary())
    def test_raises_on_incorrect_magic_word(self, magic_word):
        assume(magic_word != b'CORTEX')

        fh = Header().with_magic_word(magic_word).build()

        with pytest.raises(ValueError) as excinfo:
            header.Header.from_stream(fh)

        assert 'Saw initial magic word' in str(excinfo.value)

    @given(s.integers(min_value=0, max_value=MAX_UINT))
    def test_raises_on_incorrect_version(self, version):
        assume(version != 6)

        fh = Header().with_version(version).build()

        with pytest.raises(ValueError) as excinfo:
            header.Header.from_stream(fh)

        assert "'version' must be in [6]" in str(excinfo.value)

    def test_raises_on_invalid_kmer_size(self):
        fh = Header().with_kmer_size(0).build()

        with pytest.raises(ValueError) as excinfo:
            header.Header.from_stream(fh)

        assert "'kmer_size' has to be greater than 0" in str(excinfo.value)

    def test_raises_on_invalid_kmer_container_size(self):
        fh = Header().with_kmer_size(3).with_kmer_container_size(0).build()

        with pytest.raises(ValueError) as excinfo:
            header.Header.from_stream(fh)

        assert "'kmer_container_size' has to be greater than 0" in str(excinfo.value)

    def test_raises_on_invalid_num_colors(self):
        fh = (Header()
              .with_kmer_size(3)
              .with_kmer_container_size(1)
              .with_num_colors(0)
              .build())

        with pytest.raises(ValueError) as excinfo:
            header.Header.from_stream(fh)

        assert "'num_colors' has to be greater than 0" in str(excinfo.value)

    @given(s.integers(min_value=1, max_value=10))
    def test_raises_when_concluding_magic_word_is_wrong(self, num_colors):
        fh = (Header()
              .with_num_colors(num_colors)
              .with_mean_read_lengths([0 for _ in range(num_colors + 1)])
              .build())

        with pytest.raises(ValueError) as excinfo:
            header.Header.from_stream(fh)

        assert 'Concluding magic word' in str(excinfo.value)

    @given(s.data(),
           s.integers(min_value=1, max_value=3),
           s.integers(min_value=1, max_value=100))
    def test_loads_entire_header_successfully(self, data, num_colors, kmer_size):
        # given
        assume(kmer_size % 2 == 1)
        kmer_container_size = data.draw(s.integers(min_value=1, max_value=5))

        read_lengths = data.draw(
            s.lists(elements=s.integers(min_value=0, max_value=MAX_UINT),
                    min_size=num_colors, max_size=num_colors))

        total_sequence = data.draw(
            s.lists(elements=s.integers(min_value=0, max_value=MAX_ULONG),
                    min_size=num_colors, max_size=num_colors))

        color_names = data.draw(s.lists(
            elements=s.binary(min_size=1, max_size=3), min_size=num_colors, max_size=num_colors)
        )

        cgb = (Header()
               .with_kmer_size(kmer_size)
               .with_kmer_container_size(kmer_container_size)
               .with_num_colors(num_colors)
               .with_mean_read_lengths(read_lengths)
               .with_total_sequence(total_sequence)
               .with_color_names(*color_names)
               .with_error_rates(*[data.draw(s.binary(min_size=16, max_size=16))
                                   for _ in range(num_colors)]))

        for _ in color_names:
            cgb.with_color_information_block(data.draw(color_information_blocks()))

        fh = cgb.build()

        # then
        cgh = header.Header.from_stream(fh)

        assert cgh.version == 6
        assert cgh.kmer_size == kmer_size
        assert cgh.kmer_container_size == kmer_container_size
        assert cgh.num_colors == num_colors
        assert cgh.mean_read_lengths == tuple(read_lengths)
        assert cgh.total_sequences == tuple(total_sequence)
        assert cgh.sample_names == tuple(color_names)
        assert cgh.record_size == UINT64_T * kmer_container_size + (UINT32_T + UINT8_T) * num_colors
