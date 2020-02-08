from enum import Enum


class EdgeTraversalOrientation(Enum):
    original = 0
    reverse = 1

    @classmethod
    def other(cls, orientation):
        if orientation == cls.original:
            return cls.reverse
        return cls.original


class EngineTraversalOrientation(Enum):
    original = 0
    reverse = 1
    both = 2

class EdgeDFSTraversalDirection(Enum):
    forward = 0
    reverse = 1
