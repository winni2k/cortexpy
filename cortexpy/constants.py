from enum import Enum


class EdgeTraversalOrientation(Enum):
    original = 0
    reverse = 1

    @classmethod
    def other(cls, orientation):
        if orientation == cls.original:
            return cls.reverse
        return cls.original


class NodeEdgeDirection(Enum):
    """Refers to direction of edges relative to lexlo kmer"""
    incoming = 0
    outgoing = 1

    @classmethod
    def other(cls, orientation):
        if orientation == cls.incoming:
            return cls.outgoing
        return cls.incoming
