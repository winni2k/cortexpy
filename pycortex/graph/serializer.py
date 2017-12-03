import json
from enum import Enum
import attr
import networkx as nx
from networkx.readwrite import json_graph


@attr.s(slots=True)
class Serializer(object):
    """A class for serializing kmer graphs."""
    graph = attr.ib()
    colors = attr.ib()
    unitig_graph = attr.ib(None)

    def to_json_serializable(self):
        collapser = UnitigCollapser(self.graph, colors=self.colors).collapse_kmer_unitigs()
        self.unitig_graph = collapser.unitig_graph
        self.unitig_graph = nx.convert_node_labels_to_integers(self.unitig_graph,
                                                               label_attribute='node_key')
        self._make_graph_json_representable()

    def to_json(self):
        if self.unitig_graph is None:
            self.to_json_serializable()
        serializable = json_graph.node_link_data(self.unitig_graph, attrs={'link': 'edges'})
        return json.dumps(serializable)

    def _make_graph_json_representable(self):
        """Makes the unitig graph json representable"""
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


class EdgeTraversalOrientation(Enum):
    original = 0
    reverse = 1

    @classmethod
    def other(cls, orientation):
        if orientation == cls.original:
            return cls.reverse
        return cls.original


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


def is_unitig_end(node, graph, orientation, *, colors=(0,), test_coverage=False):
    """Returns true if the unitig ends in any color

    if test_coverage is True, then look at the coverage of nodes either side of a single link.
    If both nodes have zero coverage for that color, then don't consider that color to have ended
    the unitig.
    """
    ends = [is_unitig_end_of_color(node, graph, orientation, color=color) for color in colors]
    if not test_coverage:
        return any(ends)
    if all(ends):
        return True
    end_idx, _ = next(filter(lambda end_tup: not end_tup[1], enumerate(ends)))
    original = OrientedGraphFuncs(graph, orientation, color=colors[end_idx])

    edge = next(filter(lambda inputs: inputs[2] == colors[end_idx],
                       ((u, v, c) for u, v, c in original.edges(node, keys=True))))
    other_node = original.other_edge_node(edge)
    for color_idx, color in enumerate(colors):
        node_coverage = graph.node[node]['kmer'].coverage[color]
        other_node_coverage = graph.node[other_node]['kmer'].coverage[color]
        if ends[color_idx]:
            if node_coverage == other_node_coverage == 0:
                ends[color_idx] = False
    return any(ends)


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


def find_unitig_from(start_node, graph, *, colors=(0,), test_coverage=False):
    """Find a unitig that contains :param start_node:

    # from https://github.com/mcveanlab/mccortex/wiki/unitig
    Unitig definition: A maximal run of de Bruijn nodes (kmers, each with an orientation)
    where all nodes must have in-degree and out-degree exactly 1, with the exception of in-degree
    of the first node and out-degree of the last node (they may have any number).
    Additionally, a unitig must not visit a kmer more than once.
    """
    unitig_graph = nx.MultiDiGraph()
    unitig_graph.add_node(start_node)
    end_nodes = {EdgeTraversalOrientation.original: start_node,
                 EdgeTraversalOrientation.reverse: start_node}
    for orientation in EdgeTraversalOrientation:
        if not is_unitig_end(start_node, graph, orientation, colors=colors,
                             test_coverage=test_coverage):
            for edge_info in nx.edge_dfs(graph, start_node, orientation=orientation.name):
                source, target, _ = edge_info[:3]
                if orientation == EdgeTraversalOrientation.reverse:
                    source, target = target, source
                if target in unitig_graph or is_unitig_end(source, graph, orientation,
                                                           colors=colors,
                                                           test_coverage=test_coverage):
                    break
                end_nodes[orientation] = target
                unitig_graph.add_edge(*edge_info[:3])
    for node in unitig_graph:
        unitig_graph.add_node(node, **graph.node[node])
    return Unitig(unitig_graph,
                  end_nodes[EdgeTraversalOrientation.reverse],
                  end_nodes[EdgeTraversalOrientation.original])


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
class UnitigCollapser(object):
    graph = attr.ib()
    colors = attr.ib((0,))
    test_coverage = attr.ib(True)
    unitig_graph = attr.ib(None)

    @graph.validator
    def is_a_multi_di_graph(self, attribute, value):  # noqa
        assert isinstance(value, nx.MultiDiGraph)

    def collapse_kmer_unitigs(self):
        """Collapses unitig kmers into a single graph node.

        All nodes have an attribute `repr` added which is the node label with the left k-1 letters
        removed unless the node has no incoming edges. In that case, `repr` is the full kmer.

        Unitigs are described in :py:func:`find_unitig_from`.

        :param graph:
        :return graph:
        """
        self.find_unitigs()
        out_graph = self.unitig_graph.copy()
        for unitig_graph, data in self.unitig_graph.nodes.data():
            assert data['is_unitig']
            left_node = data['left_node']
            if self.graph.in_degree(left_node) > 0:
                short_kmer_name = [left_node[-1]]
            else:
                short_kmer_name = [left_node]
            seen_nodes = set()
            for _, target, _ in nx.edge_dfs(unitig_graph, source=left_node):
                if target not in seen_nodes:
                    seen_nodes.add(target)
                    short_kmer_name.append(target[-1])
            short_kmer_name = ''.join(short_kmer_name)
            out_graph.nodes[unitig_graph]['repr'] = short_kmer_name
        self.unitig_graph = out_graph
        return self

    def find_unitigs(self):
        """Finds unitigs of one or more nodes and replaces them by a subgraph containing the
        unitig nodes.

        Unitigs are described in :py:func:`find_unitig_from`.

        returns a graph containing unitig subgraphs
        """

        self.unitig_graph = self.graph.copy()
        visited_nodes = set()
        for start_node in self.graph:
            if start_node in visited_nodes:
                continue
            unitig = find_unitig_from(start_node, self.graph, colors=self.colors,
                                      test_coverage=self.test_coverage)
            unitig_graph_set = set(unitig.graph)
            assert visited_nodes & unitig_graph_set == set(), "{} is not disjoint from {}".format(
                visited_nodes, unitig_graph_set)
            visited_nodes |= unitig_graph_set

            self.unitig_graph = replace_unitig_nodes_with_unitig(self.unitig_graph, unitig)
        return self
