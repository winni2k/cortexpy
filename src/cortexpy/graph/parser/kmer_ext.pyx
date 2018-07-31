# distutils: language=c++

from libcpp.string cimport string
from libcpp.vector cimport vector

cdef size_t SIZE_OF_INT64 = 8
cdef string NUM_TO_LETTER_LIST = b'ACGT'

def raw_kmer_to_bytes(unsigned kmer_size, const unsigned char[:] kmer_bytes not None):
    assert kmer_size > 0
    assert kmer_size <= kmer_bytes.shape[0] * 4
    cdef vector[char] letters
    cdef vector[char] four_letters
    four_letters.resize(4,0)
    cdef char kmer_byte
    cdef size_t ulong_idx, ulong_offset, byte_offset, pair_offset, pair_idx
    for ulong_idx in range(kmer_bytes.shape[0]//SIZE_OF_INT64):
        ulong_offset = ulong_idx*SIZE_OF_INT64
        for byte_offset in reversed(range(SIZE_OF_INT64)):
            kmer_byte = kmer_bytes[ulong_offset + byte_offset]
            for pair_idx in range(4):
                four_letters[3-pair_idx] = NUM_TO_LETTER_LIST[kmer_byte & 0x3]
                kmer_byte >>= 2
            for kmer_byte in four_letters:
                letters.push_back(kmer_byte)
    return bytes(letters[(letters.size() - kmer_size):])

    # the above code was reimplemented from this
    # kmer_as_uint64ts = np.frombuffer(raw_kmer, dtype='<u8')
    # big_endian_kmer = kmer_as_uint64ts.astype('>u8')
    # kmer_as_bits = np.unpackbits(np.frombuffer(big_endian_kmer.tobytes(), dtype=np.uint8))
    # kmer = (kmer_as_bits.reshape(-1, 2) * np.array([2, 1])).sum(1)
    # return NUM_TO_LETTER[kmer[(len(kmer) - self.kmer_size):]]

def raw_kmer_to_string(int kmer_size, kmer_bytes):
    return raw_kmer_to_bytes(kmer_size, kmer_bytes).decode('utf8')

def raw_kmer_to_list(int kmer_size, kmer_bytes):
    return list(raw_kmer_to_string(kmer_size, kmer_bytes))

def raw_edges_to_list(const unsigned char[:] edge_bytes not None):
    cdef char e_byte
    cdef int i, e_byte_idx
    cdef vector[int] edge_set
    edge_set.resize(8, 0)
    tuples = []
    for e_byte_idx in range(edge_bytes.shape[0]):
        e_byte = edge_bytes[e_byte_idx]
        for i in range(8):
            edge_set[7-i] = e_byte & 0x1
            e_byte >>= 1
        tuples.append(tuple(edge_set))
    return tuples

    # edge_bytes = np.frombuffer(self._data[start:], dtype=np.uint8)
    # edge_sets = np.unpackbits(edge_bytes)
    # edge_sets = edge_sets.reshape(-1, 8)
    # edge_sets = [EdgeSet(tuple(edge_set.tolist())) for edge_set in edge_sets]
    # self._edges = edge_sets

def raw_to_coverage(const unsigned char[:] buffer not None, size_t offset, size_t num_colors):
    cdef unsigned coverage
    cdef vector[unsigned] coverages
    coverages.reserve(num_colors)
    for color in range(num_colors):
        coverage = (buffer[offset]<<0) | (buffer[offset+1]<<8) | (buffer[offset+2]<<16) | (buffer[offset+3]<<24)
        coverages.push_back(coverage)
        offset += 4
    return tuple(coverages)

    # originally:
    # start = self.kmer_container_size_in_uint64ts * UINT64_T
    # coverage_raw = self._data[start:(start + self.num_colors * UINT32_T)]
    # fmt_string = ''.join(['I' for _ in range(self.num_colors)])
    # self._coverage = unpack(fmt_string, coverage_raw)
