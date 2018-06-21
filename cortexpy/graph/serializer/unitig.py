import json

import attr
from itertools import chain
import networkx as nx

from networkx.readwrite import json_graph

from cortexpy.constants import EdgeTraversalOrientation
from cortexpy.graph.parser.kmer import revcomp_kmer_string_to_match
from cortexpy.graph.traversal.utils import OrientedGraphFuncs
from cortexpy.utils import lexlo

UNITIG_GRAPH = nx.MultiDiGraph


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
        if 'sample_names' in self.unitig_graph.graph:
            new_names = []
            for name in self.unitig_graph.graph['sample_names']:
                if isinstance(name, bytes):
                    name = name.decode()
                new_names.append(name)
            self.unitig_graph.graph['sample_names'] = new_names

        for _, node_data in self.unitig_graph.nodes.items():
            del node_data['node_key']

            coverage = []
            for coverage_list in node_data['coverage']:
                try:
                    coverage_list = coverage_list.tolist()
                except AttributeError:
                    pass
                finally:
                    coverage.append([int(i) for i in coverage_list])
            node_data['coverage'] = coverage


@attr.s(slots=True, hash=True)
class Unitig(object):
    graph = attr.ib()
    left_node = attr.ib()
    right_node = attr.ib()
    is_cycle = attr.ib(False)
    coverage = attr.ib(init=False)

    def __attrs_post_init__(self):
        coverage = [tuple(self.graph.node[self.left_node]['kmer'].coverage)]
        for _, target, _ in traverse_unitig_edges(
            self.graph,
            self.left_node,
            orientation=EdgeTraversalOrientation.original
        ):
            coverage.append(tuple(self.graph.node[target]['kmer'].coverage))
        self.coverage = tuple(coverage)

    @property
    def unitig_edge_colors(self):
        try:
            an_edge = next(iter(self.graph.edges))
        except StopIteration:
            return []
        return list(self.graph.succ[an_edge[0]][an_edge[1]])


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


def traverse_unitig_edges(graph, start_node, orientation=EdgeTraversalOrientation.original):
    node = start_node
    assert orientation in EdgeTraversalOrientation
    if orientation == EdgeTraversalOrientation.original:
        edges_func = graph.out_edges
    else:
        edges_func = graph.in_edges
    while node:
        edges = set(edges_func(node))
        if len(edges) != 1:
            break
        edge = next(iter(edges))
        yield edge[0], edge[1], None
        if orientation == EdgeTraversalOrientation.original:
            node = edge[1]
        else:
            node = edge[0]


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
        unitig_graph = UNITIG_GRAPH()
        unitig_graph.add_node(start_node)
        end_nodes = {EdgeTraversalOrientation.original: start_node,
                     EdgeTraversalOrientation.reverse: start_node}
        for orientation in EdgeTraversalOrientation:
            if not self.is_unitig_end(start_node, orientation):
                for edge_info in traverse_unitig_edges(self.graph, start_node,
                                                       orientation=orientation):
                    source, target, _ = edge_info[:3]
                    if orientation == EdgeTraversalOrientation.reverse:
                        source, target = target, source
                    if target in unitig_graph or self.is_unitig_end(source, orientation):
                        break
                    end_nodes[orientation] = target
                    unitig_graph.add_edge(*edge_info[:3])
        for node in unitig_graph:
            unitig_graph.add_node(node, **self.graph.node[node])
        ret_unitig = Unitig(unitig_graph,
                            end_nodes[EdgeTraversalOrientation.reverse],
                            end_nodes[EdgeTraversalOrientation.original])
        self.set_unitig_cycle(ret_unitig)
        return ret_unitig

    def set_unitig_cycle(self, unitig):
        unitig.is_cycle = False
        try:
            flipped_string, is_flipped = revcomp_kmer_string_to_match(
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
        unitig_graph = self.unitig_finder.find_unitigs()
        out = UNITIG_GRAPH()
        out.graph = self.graph.graph
        for _, _, unitig in unitig_graph.edges(data='unitig'):
            if unitig is None:
                continue
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
        return self

    def _summarize_unitig(self, unitig):
        left_node = unitig.left_node
        kmer_suffices = [left_node[-1]]
        seen_nodes = set()
        for _, target, _ in nx.edge_dfs(unitig.graph, source=left_node):
            assert target not in seen_nodes
            seen_nodes.add(target)
            kmer_suffices.append(target[-1])
        unitig_repr = ''.join(kmer_suffices)
        unitig_string = left_node[:-1] + unitig_repr
        if self.graph.in_degree(left_node) == 0:
            unitig_repr = unitig_string
        return unitig_repr, unitig_string
