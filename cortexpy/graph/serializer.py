import json

import attr
import networkx as nx

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from networkx.readwrite import json_graph

from cortexpy.constants import EdgeTraversalOrientation

SERIALIZER_GRAPH = nx.MultiDiGraph


@attr.s(slots=True)
class Kmers(object):
    """Serializes kmer graphs."""
    graph = attr.ib()

    def to_seq_records(self):
        return (
            SeqRecord(Seq(str(node)), id=str(node_idx)) for node_idx, node in
            enumerate(self.graph.nodes())
        )


@attr.s(slots=True)
class Serializer(object):
    """Converts kmer graphs to unitig graphs."""
    graph = attr.ib()
    collapse_unitigs = attr.ib(True)
    unitig_graph = attr.ib(init=False)
    colors = attr.ib(init=False)

    def __attrs_post_init__(self):
        self.colors = self.graph.graph['colors']

    def to_json(self):
        self._collapse_graph_to_unitigs()
        serializable = json_graph.node_link_data(self.unitig_graph, attrs={'link': 'edges'})
        return json.dumps(serializable)

    def _collapse_graph_to_unitigs(self):
        self._collapse_kmer_graph()
        self._make_unitig_graph_json_representable()

    def _collapse_kmer_graph(self):
        collapser = UnitigCollapser(self.graph).collapse_kmer_unitigs()
        self.unitig_graph = collapser.unitig_graph

    def _make_unitig_graph_json_representable(self):
        """Makes the unitig graph json representable"""
        self.unitig_graph = nx.convert_node_labels_to_integers(self.unitig_graph,
                                                               label_attribute='node_key')
        self.unitig_graph.graph['colors'] = self.colors
        for _, node_data in self.unitig_graph.nodes.items():
            if isinstance(node_data.get('node_key'), nx.Graph):
                node_data['node_key'] = repr(node_data['node_key'])

            coverage = []
            for coverage_list in node_data['coverage']:
                try:
                    coverage_list = coverage_list.tolist()
                except AttributeError:
                    pass
                finally:
                    coverage.append(coverage_list)
            node_data['coverage'] = coverage


@attr.s(slots=True)
class Unitig(object):
    graph = attr.ib()
    left_node = attr.ib()
    right_node = attr.ib()
    _coverage = attr.ib(None, init=False)

    @property
    def coverage(self):
        if self._coverage is None:
            coverage = [self.graph.node[self.left_node]['kmer'].coverage]
            for _, target, _ in nx.edge_dfs(
                self.graph,
                self.left_node,
                orientation=EdgeTraversalOrientation.original.name
            ):
                coverage.append(self.graph.node[target]['kmer'].coverage)
            self._coverage = coverage
        return self._coverage


@attr.s(slots=True)
class OrientedGraphFuncs(object):
    graph = attr.ib()
    orientation = attr.ib()
    color = attr.ib(0)
    edges = attr.ib(init=False)
    other_edge_node = attr.ib(init=False)

    def __attrs_post_init__(self):
        assert self.orientation in EdgeTraversalOrientation
        if self.orientation == EdgeTraversalOrientation.original:
            self.edges = self.graph.out_edges
            self.other_edge_node = lambda edge: edge[1]
        else:
            self.edges = self.graph.in_edges
            self.other_edge_node = lambda edge: edge[0]

    def degree(self, node):
        degree = 0
        for _, _, color in self.edges(node, keys=True, default=0):
            if color == self.color:
                degree += 1
        return degree


def is_unitig_end_of_color(node, graph, orientation, *, color=0):
    """Returns true if node is a unitig end in orientation direction"""
    assert orientation in EdgeTraversalOrientation
    original = OrientedGraphFuncs(graph, orientation, color=color)
    reverse = OrientedGraphFuncs(graph, EdgeTraversalOrientation.other(orientation), color=color)
    if original.degree(node) != 1:
        return True
    edges = list(original.edges(node))
    next_node = original.other_edge_node(edges[0])
    if reverse.degree(next_node) != 1:
        return True
    return False


def replace_unitig_nodes_with_unitig(out_graph, unitig):
    out_graph.add_node(unitig.graph, is_unitig=True,
                       left_node=unitig.left_node,
                       right_node=unitig.right_node,
                       coverage=unitig.coverage)
    for source, _, color in out_graph.in_edges(unitig.left_node, keys=True):
        out_graph.add_edge(source, unitig.graph, color)
        if source == unitig.right_node:
            out_graph.add_edge(unitig.graph, unitig.graph, color)
    for _, target, color in out_graph.out_edges(unitig.right_node, keys=True):
        out_graph.add_edge(unitig.graph, target, color)
        if target == unitig.left_node:
            out_graph.add_edge(unitig.graph, unitig.graph, color)
    out_graph.remove_nodes_from(unitig.graph)
    return out_graph


@attr.s(slots=True)
class UnitigFinder(object):
    graph = attr.ib()
    colors = attr.ib(None)
    test_coverage = attr.ib(True)

    def __attrs_post_init__(self):
        if 'colors' in self.graph.graph:
            self.colors = self.graph.graph['colors']

    def find_unitigs(self):
        """Finds unitigs of one or more nodes and replaces them by a subgraph containing the
        unitig nodes.

        Unitigs are described in :py:func:`find_unitig_from`.

        returns a graph containing unitig subgraphs
        """

        unitig_graph = self.graph.copy()
        visited_nodes = set()
        for start_node in self.graph:
            if start_node in visited_nodes:
                continue
            unitig = self.find_unitig_from(start_node)
            unitig_graph_set = set(unitig.graph)
            assert visited_nodes & unitig_graph_set == set(), "{} is not disjoint from {}".format(
                visited_nodes, unitig_graph_set)
            visited_nodes |= unitig_graph_set

            unitig_graph = replace_unitig_nodes_with_unitig(unitig_graph, unitig)
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
        unitig_graph = SERIALIZER_GRAPH()
        unitig_graph.add_node(start_node)
        end_nodes = {EdgeTraversalOrientation.original: start_node,
                     EdgeTraversalOrientation.reverse: start_node}
        for orientation in EdgeTraversalOrientation:
            if not self.is_unitig_end(start_node, orientation):
                for edge_info in nx.edge_dfs(self.graph, start_node, orientation=orientation.name):
                    source, target, _ = edge_info[:3]
                    if orientation == EdgeTraversalOrientation.reverse:
                        source, target = target, source
                    if target in unitig_graph or self.is_unitig_end(source, orientation):
                        break
                    end_nodes[orientation] = target
                    unitig_graph.add_edge(*edge_info[:3])
        for node in unitig_graph:
            unitig_graph.add_node(node, **self.graph.node[node])
        return Unitig(unitig_graph,
                      end_nodes[EdgeTraversalOrientation.reverse],
                      end_nodes[EdgeTraversalOrientation.original])

    def is_unitig_end(self, node, orientation):
        """Returns true if the unitig ends in any color

        if test_coverage is True, then look at the coverage of nodes either side of a single link.
        If both nodes have zero coverage for that color, then don't consider that color to
        have ended the unitig.
        """
        ends = [is_unitig_end_of_color(node, self.graph, orientation, color=color)
                for color in self.colors]
        if not self.test_coverage:
            return any(ends)
        if all(ends):
            return True
        end_idx, _ = next(filter(lambda end_tup: not end_tup[1], enumerate(ends)))
        original = OrientedGraphFuncs(self.graph, orientation, color=self.colors[end_idx])

        edge = next(filter(lambda inputs: inputs[2] == self.colors[end_idx],
                           ((u, v, c) for u, v, c in original.edges(node, keys=True))))
        other_node = original.other_edge_node(edge)
        for color_idx, color in enumerate(self.colors):
            node_coverage = self.graph.node[node]['kmer'].coverage[color]
            other_node_coverage = self.graph.node[other_node]['kmer'].coverage[color]
            if ends[color_idx]:
                if node_coverage == other_node_coverage == 0:
                    ends[color_idx] = False
        return any(ends)


@attr.s(slots=True)
class UnitigCollapser(object):
    graph = attr.ib()
    test_coverage = attr.ib(True)
    unitig_graph = attr.ib(None)
    unitig_finder = attr.ib(init=False)

    def __attrs_post_init__(self):
        self.unitig_finder = UnitigFinder(self.graph, colors=self.graph.graph['colors'],
                                          test_coverage=self.test_coverage)

    @graph.validator
    def is_a_multi_di_graph(self, attribute, value):  # noqa
        assert isinstance(value, SERIALIZER_GRAPH)

    def collapse_kmer_unitigs(self):
        """Collapses unitig kmers into a single graph node.

        All nodes have an attribute `repr` added which is the node label with the left k-1 letters
        removed unless the node has no incoming edges. In that case, `repr` is the full kmer.

        All nodes have an attribute `unitig` added, which is the string representation of the unitg.

        Unitigs are described in :py:func:`find_unitig_from`.

        :return graph:
        """
        self.unitig_graph = self.unitig_finder.find_unitigs()
        out_graph = self.unitig_graph.copy()
        for unitig_graph, data in self.unitig_graph.nodes.data():
            assert data['is_unitig']
            left_node = data['left_node']
            kmer_suffices = [left_node[-1]]
            seen_nodes = set()
            for _, target, _ in nx.edge_dfs(unitig_graph, source=left_node):
                if target not in seen_nodes:
                    seen_nodes.add(target)
                    kmer_suffices.append(target[-1])
            unitig_repr = ''.join(kmer_suffices)
            unitig_string = left_node[:-1] + unitig_repr
            if self.graph.in_degree(left_node) == 0:
                unitig_repr = unitig_string
            out_graph.nodes[unitig_graph]['repr'] = unitig_repr
            out_graph.nodes[unitig_graph]['unitig'] = ''.join(unitig_string)
        self.unitig_graph = out_graph
        return self
