import struct

MAX_UINT = 2 ** (struct.calcsize('I') * 8) - 1
MAX_ULONG = 2 ** (struct.calcsize('L') * 8) - 1
UINT8_T = 1
UINT32_T = 4
UINT64_T = 8
