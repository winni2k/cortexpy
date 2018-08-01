import logging
from itertools import chain

import attr
import networkx as nx

from cortexpy.constants import EdgeTraversalOrientation
from cortexpy.graph.parser.kmer import revcomp_target_to_match_ref
from cortexpy.graph.traversal.utils import OrientedGraphFuncs
from cortexpy.utils import lexlo

logger = logging.getLogger(__name__)

UNITIG_GRAPH = nx.MultiDiGraph


@attr.s(slots=True, hash=True)
class Unitig(object):
    left_node = attr.ib()
    right_node = attr.ib()
    coverage = attr.ib()
    contig = attr.ib()
    repr = attr.ib()
    graph = attr.ib()
    is_cycle = attr.ib(False)

    @classmethod
    def from_search(cls, search):
        left_node = search.end_nodes[EdgeTraversalOrientation.reverse]
        right_node = search.end_nodes[EdgeTraversalOrientation.original]
        coverage = [tuple(search.graph.node[left_node]['kmer'].coverage)]
        suffices = [left_node[-1]]
        for next_node in search.next_unitig_node_iter(EdgeTraversalOrientation.original,
                                                      start_node=left_node):
            coverage.append(tuple(search.graph.node[next_node]['kmer'].coverage))
            suffices.append(next_node[-1])
        contig = ''.join([left_node[:-1]] + suffices)
        if search.graph.in_degree(left_node) == 0:
            repr = contig
        else:
            repr = ''.join(suffices)
        return cls(left_node=left_node, right_node=right_node,
                   coverage=tuple(coverage), contig=contig, repr=repr, graph=search.unitig_graph)

    @property
    def unitig_edge_colors(self):
        try:
            an_edge = next(iter(self.graph.edges))
        except StopIteration:
            return []
        return list(self.graph.succ[an_edge[0]][an_edge[1]])


def non_missing_colors_of_kmer(kmer):
    return {c[0] for c in filter(lambda c: c[1] > 0, enumerate(kmer.coverage))}


@attr.s(slots=True)
class UnitigSearch(object):
    """Builds a single unitig"""
    graph = attr.ib()
    start_node = attr.ib()
    unitig_graph = attr.ib(attr.Factory(UNITIG_GRAPH))
    unitig_colors = attr.ib(attr.Factory(set))
    end_nodes = attr.ib(dict)

    @classmethod
    def from_node_and_graph(cls, start_node, graph):
        kmer = graph.node[start_node]['kmer']
        unitig_colors = non_missing_colors_of_kmer(kmer)
        end_nodes = {EdgeTraversalOrientation.original: start_node,
                     EdgeTraversalOrientation.reverse: start_node}
        self = cls(graph, unitig_colors=unitig_colors, end_nodes=end_nodes, start_node=start_node)
        self.unitig_graph.add_node(start_node)
        return self

    @unitig_colors.validator
    def check(self, attribute, value):
        assert len(value) > 0

    def next_unitig_node_iter(self, orientation, start_node=None):
        if start_node is None:
            start_node = self.start_node
        next_node = self.get_next_unitig_node(start_node, orientation)
        while next_node is not None and next_node != start_node:
            yield next_node
            next_node = self.get_next_unitig_node(next_node, orientation)

    def get_next_unitig_node(self, node, orientation):
        next_node = self.find_next_node_in_unitig(node, orientation)
        if next_node is None:
            return None
        if self.unitig_colors != non_missing_colors_of_kmer(self.graph.node[next_node]['kmer']):
            return None
        next_node_coming_back = self.find_next_node_in_unitig(
            next_node, EdgeTraversalOrientation.other(orientation)
        )
        if next_node_coming_back is None:
            return None
        return next_node

    def is_unitig_end(self, node, orientation):
        """Returns true if the unitig ends on this node in the specified orientation"""
        if self.get_next_unitig_node(node, orientation) is None:
            return True
        return False

    def find_next_node_in_unitig(self, node, orientation):
        """Try to get what appears from this node to be the next node in the unitig
        in orientation direction. Return None otherwise.
        """
        funcs = OrientedGraphFuncs(self.graph, orientation)
        edges = list(funcs.edges(node, keys=True))
        if len(self.unitig_colors) != len(edges):
            return None
        if self.unitig_colors != {color for _, _, color in edges}:
            return None
        return funcs.other_edge_node(edges[0])


@attr.s(slots=True)
class UnitigFinder(object):
    graph = attr.ib()
    colors = attr.ib(None)
    test_coverage = attr.ib(True)

    @graph.validator
    def check(self, attribute, value):
        if getattr(value, 'is_consistent', False):
            assert value.is_consistent()

    @classmethod
    def from_graph(cls, graph, **kwargs):
        if 'colors' not in kwargs:
            kwargs['colors'] = graph.graph['colors']
        return cls(graph, **kwargs)

    def find_unitigs(self):
        """Finds unitigs of one or more nodes and adds a Unitig object to the edge between the
        left and the right most node of the unitig.

        Unitigs are described in :py:func:`find_unitig_from`.

        returns a graph containing unitig subgraphs
        """

        unitig_graph = UNITIG_GRAPH()
        visited_nodes = set()
        for start_node in self.graph:
            if start_node in visited_nodes:
                continue
            unitig = self.find_unitig_from(start_node)
            unitig_graph_set = set(unitig.graph)
            assert visited_nodes & unitig_graph_set == set(), "{} is not disjoint from {}".format(
                visited_nodes, unitig_graph_set)
            visited_nodes |= unitig_graph_set
            unitig_graph.add_edge(unitig.left_node, unitig.right_node, unitig=unitig)
            for edge in chain(self.graph.in_edges(unitig.left_node, keys=True),
                              self.graph.out_edges(unitig.right_node, keys=True)):
                unitig_graph.add_edge(*edge)
        return unitig_graph

    def find_unitig_from(self, start_node):
        """Find a unitig that contains :param start_node:

        # from https://github.com/mcveanlab/mccortex/wiki/unitig
        Unitig definition: A maximal run of de Bruijn nodes (kmers, each with an orientation)
        where all nodes must have in-degree and out-degree exactly 1, with the exception of
        in-degree of the first node and out-degree of the last node (they may have any number).
        Additionally, a unitig must not visit a kmer more than once.

        In addition, the colors of the in-edges and out-edges must match.
        """
        search = UnitigSearch.from_node_and_graph(start_node, self.graph)
        for orientation in EdgeTraversalOrientation:
            previous_node = start_node
            for next_node in search.next_unitig_node_iter(orientation):
                if next_node in search.unitig_graph:
                    break
                search.end_nodes[orientation] = next_node
                if orientation == EdgeTraversalOrientation.original:
                    source, target = previous_node, next_node
                else:
                    target, source = previous_node, next_node
                for color in search.unitig_colors:
                    search.unitig_graph.add_edge(source, target, color)
                previous_node = next_node
        for node in search.unitig_graph:
            search.unitig_graph.add_node(node, kmer=self.graph.node[node]['kmer'])
        ret_unitig = Unitig.from_search(search)
        self.set_unitig_cycle(ret_unitig)
        return ret_unitig

    def set_unitig_cycle(self, unitig):
        unitig.is_cycle = False
        try:
            flipped_string, is_flipped = revcomp_target_to_match_ref(
                unitig.left_node,
                unitig.right_node,
                rc_is_after_reference_kmer=True
            )
        except ValueError:
            return
        if lexlo(unitig.right_node) == unitig.right_node:
            edge_letter = flipped_string[-1].upper()
        else:
            edge_letter = lexlo(flipped_string[-1]).lower()
        colors = unitig.unitig_edge_colors
        if len(colors) == 0:
            return
        for color in colors:
            if not self.graph.node[unitig.right_node]['kmer'].edges[color].is_edge(edge_letter):
                return
        unitig.is_cycle = True


@attr.s(slots=True)
class UnitigCollapser(object):
    graph = attr.ib()
    test_coverage = attr.ib(True)
    unitig_graph = attr.ib(None)
    unitig_finder = attr.ib(init=False)

    def __attrs_post_init__(self):
        self.unitig_finder = UnitigFinder.from_graph(self.graph, test_coverage=self.test_coverage)

    @graph.validator
    def is_a_multi_di_graph(self, attribute, value):
        assert value.is_multigraph()
        assert value.is_directed()

    def collapse_kmer_unitigs(self):
        """Collapses unitig kmers into a single graph node.

        All nodes have an attribute `repr` added which is the node label with the left k-1 letters
        removed unless the node has no incoming edges. In that case, `repr` is the full kmer.

        All nodes have an attribute `unitig` added, which is the string representation of the unitg.

        Unitigs are described in :py:func:`find_unitig_from`.

        :return graph:
        """

        logger.info(f'Collapsing graph with {len(self.graph)} kmers')
        unitig_graph = self.unitig_finder.find_unitigs()
        out = UNITIG_GRAPH()
        out.graph = self.graph.graph
        for s, t, color, unitig in unitig_graph.edges(data='unitig', keys=True):
            if unitig is None:
                out.add_edge(s, t, color)
                continue
            # todo: use unitig for unitig_string
            unitig_repr, unitig_string = self._summarize_unitig(unitig)
            out.add_node(unitig, repr=unitig_repr, unitig=unitig_string, coverage=unitig.coverage)
            out.add_edge(unitig.left_node, unitig)
            out.add_edge(unitig, unitig.right_node)
            for pred in self.graph.pred[unitig.left_node]:
                for color in self.graph.pred[unitig.left_node][pred]:
                    out.add_edge(pred, unitig.left_node, color)
            for succ in self.graph.succ[unitig.right_node]:
                for color in self.graph.succ[unitig.right_node][succ]:
                    out.add_edge(unitig.right_node, succ, color)
        new_edges = []
        for unitig in out.nodes:
            if isinstance(unitig, Unitig):
                for other_right_node in out.pred[unitig.left_node]:
                    for other_unitig in out.pred[other_right_node]:
                        if isinstance(other_unitig, Unitig) and other_unitig is not unitig:
                            for color in out.pred[unitig.left_node][other_right_node]:
                                new_edges.append((other_unitig, unitig, color))
        for source, target, color in new_edges:
            out.add_edge(source, target, color)
        for node in list(out.nodes):
            if not isinstance(node, Unitig):
                out.remove_node(node)
            elif node.is_cycle:
                for color in node.unitig_edge_colors:
                    out.add_edge(node, node, color)

        self.unitig_graph = out
        logger.info(
            f'Collapsing complete. Collapsed graph contains {len(self.unitig_graph)} unitigs')
        return self

    def _summarize_unitig(self, unitig):
        return unitig.repr, unitig.contig
