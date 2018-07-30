# distutils: language=c++

from struct import unpack

from libcpp.string cimport string
from libcpp.vector cimport vector

def raw_kmer_to_bytes(int kmer_size, kmer_bytes):
    assert kmer_size > 0
    assert kmer_size <= len(kmer_bytes) * 4
    cdef string NUM_TO_LETTER_LIST = b'ACGT'
    format_string = ''.join(['<'] + ['Q' for _ in range(len(kmer_bytes) / 8)])
    container = unpack(format_string, kmer_bytes)
    kmer_as_bytes = b''.join([uint.to_bytes(8, 'big') for uint in container])
    cdef vector[char] letters
    for kmer_byte in kmer_as_bytes:
        letters.push_back(NUM_TO_LETTER_LIST[(kmer_byte >> 6) & 0x3])
        letters.push_back(NUM_TO_LETTER_LIST[(kmer_byte >> 4) & 0x3])
        letters.push_back(NUM_TO_LETTER_LIST[(kmer_byte >> 2) & 0x3])
        letters.push_back(NUM_TO_LETTER_LIST[kmer_byte & 0x3])
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
