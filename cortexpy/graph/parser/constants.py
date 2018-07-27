import numpy as np

CORTEX_MAGIC_WORD = (b'C', b'O', b'R', b'T', b'E', b'X')
CORTEX_VERSION = 6
ERROR_RATE_SIZE = 16
UINT8_T = np.uint8().itemsize
UINT32_T = np.uint32().itemsize
UINT64_T = np.uint64().itemsize
BITS_IN_BYTE = 8
LETTERS_PER_BYTE = BITS_IN_BYTE // 2
NUM_TO_LETTER = np.array(list('ACGT'))
LETTER_TO_NUM = bytes.maketrans(b'ACGT', bytes([0x00, 0x01, 0x02, 0x03]))
ASCII_OFFSET_OF_ZERO = 48
NUM_LETTERS_PER_UINT = UINT64_T * LETTERS_PER_BYTE
NUM_TO_BITS = np.array([[0, 0], [0, 1], [1, 0], [1, 1]])
