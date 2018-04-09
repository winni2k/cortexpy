import numpy as np

CORTEX_MAGIC_WORD = (b'C', b'O', b'R', b'T', b'E', b'X')
CORTEX_VERSION = 6
UINT8_T = np.uint8().itemsize
UINT32_T = np.uint32().itemsize
UINT64_T = np.uint64().itemsize
BITS_IN_BYTE = 8
LETTERS_PER_BYTE = BITS_IN_BYTE // 2
NUM_TO_LETTER = np.array(['A', 'C', 'G', 'T'])
LETTER_TO_NUM = bytes.maketrans(b'ACGT', b'0123')
ASCII_OFFSET_OF_ZERO = 48
