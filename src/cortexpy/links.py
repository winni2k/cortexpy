import json
from collections import defaultdict
from enum import Enum

import attr

from cortexpy.utils import lexlo, comp


@attr.s(slots=True)
class LinkWalker:
    links = attr.ib()
    junctions = attr.ib(init=False)

    def __attrs_post_init__(self):
        self.junctions = defaultdict(list)

    @property
    def n_junctions(self):
        return sum(len(juncs) for juncs in self.junctions.values())

    def load_kmer(self, kmer):
        lexlo_kmer = lexlo(kmer)
        is_lexlo = lexlo_kmer == kmer
        try:
            link_group = self.links.body[lexlo_kmer]
        except KeyError:
            pass
        else:
            for junc in link_group.get_link_junctions(is_lexlo, in_kmer_orientation=True):
                self.junctions[junc[0]].append(junc)
        return self

    def choose_junction(self, base):
        if len(self.junctions.keys()) == 0:
            return self
        if base in self.junctions:
            junctions_to_traverse = self.junctions[base]
            self.clear()
            for junc in junctions_to_traverse:
                if len(junc) > 1:
                    self.junctions[junc[1]].append(junc[1:])
            return self
        raise KeyError('Invalid junction choice. Valid junction choices are: %s',
                       self.junctions.keys())

    def next_junction_bases(self):
        return self.junctions.keys()

    def clear(self):
        self.junctions.clear()
        return self


class LinkOrientation(Enum):
    F = 0
    R = 1

    @classmethod
    def other(cls, orientation):
        if orientation == cls.F:
            return cls.R
        else:
            return cls.F


@attr.s(slots=True)
class Links:
    header = attr.ib()
    body = attr.ib()

    @classmethod
    def from_binary_stream(cls, stream):
        header = LinksHeader.from_binary_stream(stream)
        if header.json['graph']['num_colours'] != 1:
            raise NotImplementedError
        body = LinksBody.from_binary_stream(stream)

        return cls(header, body)


@attr.s(slots=True)
class LinksHeader:
    json = attr.ib()

    @classmethod
    def from_binary_stream(cls, stream):
        lines = []
        bases = list(b'ACTG')
        print(bases)
        last_line = False
        while not last_line:
            line = stream.readline()
            peek = stream.peek(1)
            # this is a really nasty way of determining when to stop reading but I've got nothing better -,-
            if peek[0] in bases:
                last_line = True
            if line.startswith((b'#', b'\n')):
                continue
            lines.append(line.decode().rstrip())
        return cls(json.loads(''.join(lines)))


@attr.s(slots=True)
class LinksBody:

    @classmethod
    def from_binary_stream(cls, stream):
        body_dict = {}
        for group in link_groups(useful_lines(stream)):
            body_dict[group.kmer] = group
        return body_dict


@attr.s(slots=True)
class LinkGroup:
    kmer = attr.ib()
    coverage = attr.ib()
    link_lines = attr.ib(attr.Factory(list))

    # def junction_bases(self, is_lexlo, orientation=LinkOrientation.F):
    #     if not is_lexlo:
    #         orientation = LinkOrientation.other(orientation)
    #     bases = []
    #     for line in self.link_lines:
    #         if line.orientation != orientation:
    #             continue
    #         bases.append(line.juncs[0])
    #     return bases

    def get_link_junctions(self, is_lexlo, in_kmer_orientation=True):
        """kmer orientation is from the perspective of the potentially non-lexlo kmer"""
        if is_lexlo == in_kmer_orientation:
            orientation = LinkOrientation.F
        else:
            orientation = LinkOrientation.R
        for line in self.link_lines:
            if line.orientation != orientation:
                continue
            junctions = line.juncs
            if orientation == LinkOrientation.R:
                junctions = comp(junctions)
            yield junctions


@attr.s(slots=True)
class LinkGroupTraverser:
    link_group = attr.ib()
    age = attr.ib(0)

    def __getattr__(self, item):
        return getattr(self.link_group, item)

    def age_links(self, n_junctions=1):
        pass


@attr.s(slots=True)
class LinkLine:
    orientation = attr.ib()
    num_juncs = attr.ib()
    juncs = attr.ib()
    counts = attr.ib(attr.Factory(list))

    @classmethod
    def from_list(cls, fields):
        num_juncs = int(fields[1])
        juncs = fields[3]
        assert len(juncs) == num_juncs
        return cls(orientation=LinkOrientation[fields[0]],
                   num_juncs=num_juncs,
                   counts=[int(fields[2])],
                   juncs=juncs)


def link_groups(lines):
    lg = None
    for line in lines:
        line = line.decode()
        fields = line.rstrip().split()
        if 2 == len(fields):
            if lg is not None:
                yield lg
            lg = LinkGroup(fields[0], int(fields[1]))
            continue
        lg.link_lines.append(LinkLine.from_list(fields))
    if lg is not None:
        yield lg


def useful_lines(lines):
    for line in lines:
        if line == b'\n' or line.startswith(b'#'):
            continue
        yield line
