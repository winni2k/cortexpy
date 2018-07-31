import attr
from cortexpy.graph.cortex import (
    build_empty_cortex_graph_from_ra_parser, ConsistentCortexDiGraph,
    CortexDiGraph,
)
from cortexpy.constants import EdgeTraversalOrientation, EngineTraversalOrientation

SERIALIZER_GRAPH = CortexDiGraph


class KmerStringAlreadySeen(Exception):
    pass


@attr.s(slots=True)
class Traverser(object):
    ra_parser = attr.ib()
    traversal_color = attr.ib(0)
    graph = attr.ib(attr.Factory(SERIALIZER_GRAPH))
    other_stopping_colors = attr.ib(attr.Factory(set))
    kmer = attr.ib(init=False, default=None)
    kmer_string = attr.ib(init=False)
    prev_kmer = attr.ib(init=False)
    prev_kmer_string = attr.ib(init=False)
    orientation = attr.ib(init=False)
    parent_graph = attr.ib(init=False)

    def __attrs_post_init__(self):
        assert self.traversal_color not in self.other_stopping_colors

    def traverse_from(self, kmer_string, *,
                      orientation=EdgeTraversalOrientation.original,
                      parent_graph=None):
        if parent_graph is None:
            parent_graph = set()
        self.parent_graph = parent_graph
        self.graph = ConsistentCortexDiGraph(
            graph=build_empty_cortex_graph_from_ra_parser(self.ra_parser).graph)
        self.kmer_string = first_kmer_string = kmer_string
        self.orientation = orientation
        self.prev_kmer_string = None

        try:
            self._get_kmer_and_add_kmer_string_to_graph()
        except KmerStringAlreadySeen:
            return Traversed(self.graph, orientation=self.orientation)

        last_oriented_edge_set = self._traverse()

        reverse_neighbor_kmer_strings = set(
            self._get_neighbors(last_oriented_edge_set.other_orientation()))
        if self.prev_kmer_string is not None:
            reverse_neighbor_kmer_strings.remove(self.prev_kmer_string)
        return Traversed(self.graph,
                         orientation=self.orientation,
                         first_kmer_string=first_kmer_string,
                         last_kmer_string=self.kmer_string,
                         neighbor_kmer_strings=self._get_neighbors(last_oriented_edge_set),
                         reverse_neighbor_kmer_strings=list(reverse_neighbor_kmer_strings))

    def _traverse(self):
        while True:
            traversal_edge_set = self.kmer.edges[self.traversal_color].oriented(
                self.orientation)
            if (
                self._get_num_neighbors(traversal_edge_set) != 1
                or self._get_num_neighbors(traversal_edge_set.other_orientation()) > 1
            ):
                return traversal_edge_set

            for stop_color in self.other_stopping_colors:
                stop_color_edge_set = self.kmer.edges[stop_color].oriented(self.orientation)
                if (
                    self._get_num_neighbors(stop_color_edge_set) != 0
                    or self._get_num_neighbors(stop_color_edge_set.other_orientation()) != 0
                ):
                    return traversal_edge_set

            try:
                self._add_next_kmer_string_to_graph_and_get_next_kmer(traversal_edge_set)
            except KmerStringAlreadySeen:
                return traversal_edge_set

    def _get_num_neighbors(self, oriented_edge_set):
        return oriented_edge_set.num_neighbor(self.kmer_string)

    def _get_neighbors(self, oriented_edge_set):
        return oriented_edge_set.neighbor_kmer_strings(self.kmer_string)

    def _add_next_kmer_string_to_graph_and_get_next_kmer(self, oriented_edge_set):
        next_kmer_string = next(oriented_edge_set.neighbor_kmer_strings(self.kmer_string))
        prev_kmer_string = self.kmer_string
        try:
            self.kmer_string = next_kmer_string
            self._get_kmer_and_add_kmer_string_to_graph()
        except KmerStringAlreadySeen:
            self.kmer_string = prev_kmer_string
            raise
        self.prev_kmer_string = prev_kmer_string

    def _get_kmer(self):
        if self.kmer_string in self.graph or self.kmer_string in self.parent_graph:
            raise KmerStringAlreadySeen
        prev_kmer = self.kmer
        self.kmer = self.ra_parser.get_kmer_for_string(self.kmer_string)
        self.prev_kmer = prev_kmer

    def _get_kmer_and_add_kmer_string_to_graph(self):
        self._get_kmer()
        self.graph.add_node(self.kmer_string, kmer=self.kmer)


@attr.s(slots=True)
class Traversed(object):
    """A branch, the result of a branch traversal"""
    graph = attr.ib()
    orientation = attr.ib()
    first_kmer_string = attr.ib(None)
    last_kmer_string = attr.ib(None)
    neighbor_kmer_strings = attr.ib(attr.Factory(list))
    reverse_neighbor_kmer_strings = attr.ib(attr.Factory(list))

    def is_empty(self):
        return len(self.graph) == 0


@attr.s(slots=True, frozen=True, hash=True)
class TraversalSetup(object):
    start_string = attr.ib()
    orientation = attr.ib()
    traversal_color = attr.ib()
    connecting_node = attr.ib(None)


@attr.s(slots=True)
class Queuer(object):
    queue = attr.ib()
    traversal_colors = attr.ib()
    engine_orientation = attr.ib(EngineTraversalOrientation.original)
    _orientations = attr.ib(init=False)
    _seen_traversal_setups = attr.ib(attr.Factory(set))

    def __attrs_post_init__(self):
        if self.engine_orientation == EngineTraversalOrientation.both:
            self._orientations = list(EdgeTraversalOrientation)
        else:
            self._orientations = [EdgeTraversalOrientation[self.engine_orientation.name]]

    def add_from(self, start_string, orientation, connecting_node, traversal_color):
        traversal_setup = TraversalSetup(start_string=start_string,
                                         orientation=orientation,
                                         connecting_node=connecting_node,
                                         traversal_color=traversal_color)
        if traversal_setup not in self._seen_traversal_setups:
            self._seen_traversal_setups.add(traversal_setup)
            self.queue.append(traversal_setup)

    def add_from_branch(self, branch):
        orientation_neighbor_pairs = [
            (branch.orientation, branch.neighbor_kmer_strings)]
        if self.engine_orientation == EngineTraversalOrientation.both:
            orientation_neighbor_pairs.append(
                (EdgeTraversalOrientation.other(branch.orientation),
                 branch.reverse_neighbor_kmer_strings))
        else:
            assert EdgeTraversalOrientation[
                       self.engine_orientation.name] == branch.orientation
        for orientation, neighbor_strings in orientation_neighbor_pairs:
            for neighbor_string in neighbor_strings:
                for color in self.traversal_colors:
                    self.add_from(start_string=neighbor_string,
                                  orientation=orientation,
                                  connecting_node=branch.last_kmer_string,
                                  traversal_color=color)
