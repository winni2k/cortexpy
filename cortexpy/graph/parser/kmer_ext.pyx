from struct import unpack

from cortexpy.graph.parser.constants import NUM_TO_LETTER_LIST

def raw_kmer_to_letters(kmer_size, kmer_bytes):
    format_string = ''.join(['<'] + ['Q' for _ in range(len(kmer_bytes) / 8)])
    container = unpack(format_string, kmer_bytes)
    kmer_as_bytes = b''.join([uint.to_bytes(8,'big') for uint in container])
    letters=[]
    for kmer_byte in kmer_as_bytes:
        letters.append(NUM_TO_LETTER_LIST[(kmer_byte >> 6)  & 0x3])
        letters.append(NUM_TO_LETTER_LIST[(kmer_byte >> 4)  & 0x3])
        letters.append(NUM_TO_LETTER_LIST[(kmer_byte >> 2)  & 0x3])
        letters.append(NUM_TO_LETTER_LIST[kmer_byte & 0x3])
    return letters[(len(letters) - kmer_size):]

    # kmer_as_uint64ts = np.frombuffer(raw_kmer, dtype='<u8')
    # big_endian_kmer = kmer_as_uint64ts.astype('>u8')
    # kmer_as_bits = np.unpackbits(np.frombuffer(big_endian_kmer.tobytes(), dtype=np.uint8))
    # kmer = (kmer_as_bits.reshape(-1, 2) * np.array([2, 1])).sum(1)
    # return NUM_TO_LETTER[kmer[(len(kmer) - self.kmer_size):]]
