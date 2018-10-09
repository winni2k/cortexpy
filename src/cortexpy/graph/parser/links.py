import json
from enum import Enum

import attr


class LinkOrientation(Enum):
    F = 0
    R = 1


@attr.s(slots=True)
class Links:
    header = attr.ib()
    body = attr.ib()

    @classmethod
    def from_stream(cls, stream):
        header = LinksHeader.from_stream(stream)
        if header.json['graph']['num_colours'] != 1:
            raise NotImplementedError
        body = LinksBody.from_stream(stream)

        return cls(header, body)


@attr.s(slots=True)
class LinksHeader:
    json = attr.ib()

    @classmethod
    def from_stream(cls, stream):
        lines = []
        for line in stream:
            lines.append(line)
            if line.startswith('}'):
                break
        return cls(json.loads(''.join(lines)))


@attr.s(slots=True)
class LinksBody:

    @classmethod
    def from_stream(cls, stream):
        body_dict = {}

        for group in link_groups(useful_lines(stream)):
            body_dict[group.kmer] = group
        return body_dict


@attr.s(slots=True)
class LinkGroup:
    kmer = attr.ib()
    coverage = attr.ib()
    link_lines = attr.ib(attr.Factory(list))


@attr.s(slots=True)
class LinkLine:
    orientation = attr.ib()
    num_juncs = attr.ib()
    juncs = attr.ib()
    counts = attr.ib(attr.Factory(list))

    @classmethod
    def from_list(cls, fields):
        return cls(orientation=LinkOrientation[fields[0]],
                   num_juncs=int(fields[1]),
                   counts=[int(fields[2])],
                   juncs=fields[3])


def link_groups(lines):
    lg = None
    for line in lines:
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
        if line.startswith('#'):
            continue
        if line == '\n':
            continue
        yield line
