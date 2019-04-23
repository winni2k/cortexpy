import copy
import json
from collections import defaultdict
from collections.abc import Sequence
from enum import Enum
from logging import getLogger

import attr

from cortexpy.utils import lexlo

logger = getLogger('cortexpy.links')


@attr.s(slots=True)
class LinkedGraphTraverser(Sequence):
    """Adapter for linked walkers to be able to work with nx.all_simple_paths"""
    graph = attr.ib()
    walkers = attr.ib(attr.Factory(dict))

    @classmethod
    def from_graph_and_link_walker(cls, graph, link_walker):
        return cls(graph, {link_walker.current_unitig: link_walker})

    def __contains__(self, item):
        return item in self.graph

    def __iter__(self):
        pass

    def __len__(self):
        pass

    def __getitem__(self, item):
        parent_walker = self.walkers[item]
        successors = list(parent_walker.successors())
        if len(successors) == 0:
            return []

        child_walkers = [copy.copy(parent_walker) for _ in range(len(successors))]
        children = []
        for succ, walker in zip(successors, child_walkers):
            walker.choose(succ)
            children.append(walker.current_unitig)
            self.walkers[children[-1]] = walker
        return children


@attr.s(slots=True)
class UnitigLinkWalker:
    """Traverses a unitig graph with links."""
    link_walker = attr.ib()
    unitigs = attr.ib()
    kmer_size = attr.ib()
    current_unitig = attr.ib()

    @classmethod
    def from_links_unitigs_kmer_size_unitig(cls, links, unitigs, kmer_size, unitig):
        obj = cls(LinkWalker.from_links(links), unitigs, kmer_size, unitig)
        logger.debug('Creating UnitigWalker with unitig: %s', unitigs.nodes[obj.current_unitig])
        obj.link_walker.load_kmer(obj._current_unitig_right_kmer())
        return obj

    def successors(self):
        """Returns nodes from links or all available junctions if no link info exists"""
        successors = list(self.unitigs.successors(self.current_unitig))
        if len(successors) < 2:
            return successors
        j_unitigs = list(self.link_successors())
        if len(j_unitigs) != 0:
            return j_unitigs
        return successors

    def link_successors(self):
        """Only returns unitigs based on link information"""
        successors = list(self.unitigs.successors(self.current_unitig))
        if len(successors) < 2:
            raise ValueError(
                'Tried to call link_successors for unitig that has only one successor %s',
                [self.current_unitig, self.unitigs.nodes[self.current_unitig]]
            )
        available_bases = set(self.link_walker.next_junction_bases())
        successors = {s: self._unitig_choice_base(s) for s in successors}
        if not available_bases <= set(successors.values()):
            raise ValueError(
                f"""Links do not appear to match unitigs. Have these links been constructed on a different cortex graph?
                current unitig: {self.unitigs.nodes[self.current_unitig]}
                successors: {[self.unitigs.nodes[s] for s in successors]}
                available bases: {available_bases}
                junctions: {self.link_walker.junctions}"""

            )
        for succ in successors.keys():
            if successors[succ] in available_bases:
                yield succ

    def choose(self, successor):
        """Register the choice of a successor and advance"""
        logger.debug('Choosing next unitig: %s', self.unitigs.nodes[successor])
        next_unitigs = list(self.unitigs.successors(self.current_unitig))
        assert successor in next_unitigs
        if len(next_unitigs) == 0:
            raise ValueError('Cannot choose a successor for node that has no successors')
        if len(next_unitigs) > 1:
            if next(self.link_successors(), None) is not None:
                self.link_walker.choose_branch(self._unitig_choice_base(successor))
        self._advance_to_successor(successor)
        return self

    def __copy__(self):
        return UnitigLinkWalker(copy.copy(self.link_walker),
                                self.unitigs,
                                self.kmer_size,
                                self.current_unitig)

    def _advance_to_successor(self, successor):
        self.current_unitig = successor
        self.link_walker.load_kmer(self._current_unitig_right_kmer())

    def _current_unitig_string(self):
        return self.unitigs.nodes[self.current_unitig]['unitig']

    def _current_unitig_right_kmer(self):
        unitig_string = self._current_unitig_string()
        return unitig_string[(len(unitig_string) - self.kmer_size):]

    def _unitig_choice_base(self, unitig_id):
        return self.unitigs.nodes[unitig_id]['unitig'][self.kmer_size - 1]


@attr.s(slots=True)
class LinkWalker:
    """Manages the loading and walking of links for kmers"""
    links = attr.ib()
    junctions = attr.ib()

    @classmethod
    def from_links(cls, links):
        junctions = defaultdict(list)
        return cls(links, junctions)

    @property
    def n_junctions(self):
        return sum(len(juncs) for juncs in self.junctions.values())

    def load_kmer(self, kmer):
        """Load the link group for a kmer in the orientation of the kmer."""
        lexlo_kmer = lexlo(kmer)
        is_lexlo = lexlo_kmer == kmer
        try:
            link_group = self.links.body[lexlo_kmer]
            logger.debug('Loaded link group for kmer %s: %s', kmer, link_group)
        except KeyError:
            pass
        else:
            for junc in link_group.get_link_junctions_in_kmer_orientation(is_lexlo):
                self.junctions[junc[0]].append(junc)
        return self

    def choose_branch(self, base):
        """Choose a branch and advance all links. Keep only links consistent with branch."""
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
        """Returns the the bases of the branches that can be chosen."""
        return self.junctions.keys()

    def clear(self):
        self.junctions.clear()
        return self

    def __copy__(self):
        return LinkWalker(self.links, copy.copy(self.junctions))


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
        last_line = False
        while not last_line:
            line = stream.readline()
            peek = stream.peek(1)
            # this is a really nasty way of determining when to stop reading but I've got nothing better -,-
            if len(peek) == 0 or peek[0] in bases:
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

    def get_link_junctions_in_kmer_orientation(self, is_lexlo):
        """kmer orientation is from the perspective of the potentially non-lexlo kmer"""
        if is_lexlo:
            orientation = LinkOrientation.F
        else:
            orientation = LinkOrientation.R
        for line in self.link_lines:
            if line.orientation != orientation:
                continue
            yield line.juncs

    def __str__(self):
        elements = [f'<{self.kmer} {self.coverage}: ']
        for line in self.link_lines:
            elements.append('|'.join([str(l) for l in self.link_lines]))
        elements.append('>')
        return ''.join(elements)


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

    def __str__(self):
        return f'{self.orientation.name} {self.num_juncs} {",".join([str(c) for c in self.counts])} {self.juncs}'


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
