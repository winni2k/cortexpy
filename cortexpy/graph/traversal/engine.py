from enum import Enum

import attr
import collections
import networkx as nx

import cortexpy.graph.traversal
from cortexpy.graph.serializer import SERIALIZER_GRAPH, EdgeTraversalOrientation


@attr.s(slots=True)
class BranchTraversalSetup(object):
    start_string = attr.ib()
    orientation = attr.ib()
    connecting_node = attr.ib(None)


class EngineTraversalOrientation(Enum):
    original = 0
    reverse = 1
    both = 2


@attr.s(slots=True)
class Engine(object):
    ra_parser = attr.ib()
    traversal_color = attr.ib(0)
    orientation = attr.ib(EngineTraversalOrientation.original)
    graph = attr.ib(attr.Factory(SERIALIZER_GRAPH))
    max_nodes = attr.ib(1000)
    branch_queue = attr.ib(attr.Factory(collections.deque))
    _orientations = attr.ib(init=False)

    def __attrs_post_init__(self):
        if self.orientation == EngineTraversalOrientation.both:
            self._orientations = set(EdgeTraversalOrientation)
        else:
            self._orientations = {EdgeTraversalOrientation[self.orientation.name]}

    def traverse_from(self, start_string):
        self.graph = SERIALIZER_GRAPH()
        branch_traverser = cortexpy.graph.traversal.Branch(self.ra_parser, self.traversal_color)
        for orientation in self._orientations:
            self.branch_queue.append(BranchTraversalSetup(start_string, orientation=orientation))
        while 0 < len(self.branch_queue) and len(self.graph) < self.max_nodes:
            setup = self.branch_queue.popleft()
            branch = branch_traverser.traverse_from(setup.start_string,
                                                    orientation=setup.orientation,
                                                    parent_graph=self.graph)
            self.graph = nx.union(self.graph, branch.graph)
            self._connect_branch_to_parent_graph(branch, setup)
            self._link_branch_and_enqueue_neighbor_traversals(branch)
        return self.graph

    def _connect_branch_to_parent_graph(self, branch, setup):
        if setup.connecting_node is not None and branch.first_kmer_string is not None:
            self._add_edge_in_orientation(setup.connecting_node, branch.first_kmer_string,
                                          setup.orientation)

    def _link_branch_and_enqueue_neighbor_traversals(self, branch):
        for neighbor_string in branch.neighbor_kmer_strings:
            if neighbor_string in self.graph:
                self._add_edge_in_orientation(branch.last_kmer_string,
                                              neighbor_string,
                                              branch.orientation)
            else:
                for orientation in self._orientations:
                    self.branch_queue.append(
                        BranchTraversalSetup(neighbor_string,
                                             orientation=orientation,
                                             connecting_node=branch.last_kmer_string)
                    )

    def _add_edge_in_orientation(self, first, second, orientation):
        if orientation == EdgeTraversalOrientation.reverse:
            first, second = second, first
        self.graph.add_edge(first, second, key=self.traversal_color)
