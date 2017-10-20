import json
from enum import Enum
import attr
import networkx as nx
from networkx.readwrite import json_graph


@attr.s(slots=True)
class Serializer(object):
    """A class for serializing kmer graphs."""
    graph = attr.ib()

    def to_json_serializable(self):
        graph = self.graph
        graph = collapse_kmer_unitigs(graph)
        graph = nx.convert_node_labels_to_integers(graph, label_attribute='node_key')
        return _make_graph_json_representable(graph)

    def to_json(self):
        graph = self.to_json_serializable()
        serializable = json_graph.node_link_data(graph, attrs={'link': 'edges'})
        return json.dumps(serializable)


def _make_graph_json_representable(graph):
    """Makes a unitig graph json representable"""
    graph = graph.copy()
    for node, node_data in graph.nodes.items():
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
    for source, target, edge_data in graph.edges.data():
        is_missing = bool(graph.node[source]['is_missing'] or graph.node[target]['is_missing'])
        graph[source][target]['is_missing'] = is_missing
    return graph


def missing_kmers_change_in_edge(graph, edge):
    source_coverage = graph.node[edge[0]]['kmer'].coverage
    target_coverage = graph.node[edge[1]]['kmer'].coverage
    if source_coverage is None and target_coverage is None:
        return False

    if source_coverage is None or target_coverage is None:
        return True

    if len(source_coverage) != len(target_coverage):
        return True

    for source_color_coverage, target_color_coverage in zip(source_coverage, target_coverage):
        if (source_color_coverage == 0) != (target_color_coverage == 0):
            return True
    return False


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
            for source, target in nx.edge_dfs(self.graph, self.left_node,
                                              orientation=EdgeTraversalOrientation.original.name):
                coverage.append(self.graph.node[target]['kmer'].coverage)
            self._coverage = coverage
        return self._coverage

    @property
    def is_missing(self):
        for color_coverage in self.coverage:
            if color_coverage is not None and any(c != 0 for c in color_coverage):
                return False
        return True


class EdgeTraversalOrientation(Enum):
    original = 0
    reverse = 1

    @classmethod
    def other(cls, orientation):
        if orientation == cls.original:
            return cls.reverse
        return cls.original


@attr.s
class OrientedGraphFuncs(object):
    graph = attr.ib()
    orientation = attr.ib()
    degree = attr.ib(init=False)
    edges = attr.ib(init=False)

    def __attrs_post_init__(self):
        assert self.orientation in EdgeTraversalOrientation
        if self.orientation == EdgeTraversalOrientation.original:
            self.degree = self.graph.out_degree
            self.edges = self.graph.out_edges
        else:
            self.degree = self.graph.in_degree
            self.edges = self.graph.in_edges


def is_unitig_end(node, graph, orientation):
    """Returns true if node is a unitig end in orientation direction"""
    assert orientation in EdgeTraversalOrientation
    original = OrientedGraphFuncs(graph, orientation)
    reverse = OrientedGraphFuncs(graph, EdgeTraversalOrientation.other(orientation))
    if original.degree(node) != 1:
        return True
    edge = next(e for e in original.edges(node))
    if orientation != EdgeTraversalOrientation.original:
        other_edge_node = edge[0]
    else:
        other_edge_node = edge[1]
    if missing_kmers_change_in_edge(graph, edge) or reverse.degree(other_edge_node) > 1:
        return True
    return False


def find_unitig_from(start_node, graph):
    """Find a unitig that contains :param start_node:

    # from https://github.com/mcveanlab/mccortex/wiki/unitig
    Unitig definition: A maximal run of de Bruijn nodes (kmers, each with an orientation)
    where all nodes must have in-degree and out-degree exactly 1, with the exception of in-degree
    of the first node and out-degree of the last node (they may have any number).
    Additionally, a unitig must not visit a kmer more than once.
    """
    unitig_graph = nx.DiGraph()
    unitig_graph.add_node(start_node)
    end_nodes = {EdgeTraversalOrientation.original: start_node,
                 EdgeTraversalOrientation.reverse: start_node}
    for orientation in EdgeTraversalOrientation:
        if not is_unitig_end(start_node, graph, orientation):
            for edge_info in nx.edge_dfs(graph, start_node, orientation=orientation.name):
                source, target = edge_info[0], edge_info[1]
                if orientation == EdgeTraversalOrientation.reverse:
                    source, target = target, source
                if target in unitig_graph or is_unitig_end(source, graph, orientation):
                    break
                end_nodes[orientation] = target
                unitig_graph.add_edge(edge_info[0], edge_info[1])
    for node in unitig_graph:
        unitig_graph.add_node(node, **graph.node[node])
    return Unitig(unitig_graph,
                  end_nodes[EdgeTraversalOrientation.reverse],
                  end_nodes[EdgeTraversalOrientation.original])


def replace_unitig_nodes_with_unitig(out_graph, unitig):
    out_graph.add_node(unitig.graph, is_unitig=True,
                       left_node=unitig.left_node,
                       is_missing=unitig.is_missing,
                       right_node=unitig.right_node,
                       coverage=unitig.coverage)
    for source, _ in out_graph.in_edges(unitig.left_node):
        out_graph.add_edge(source, unitig.graph)
        if source == unitig.right_node:
            out_graph.add_edge(unitig.graph, unitig.graph)
    for _, target in out_graph.out_edges(unitig.right_node):
        out_graph.add_edge(unitig.graph, target)
        if target == unitig.left_node:
            out_graph.add_edge(unitig.graph, unitig.graph)
    out_graph.remove_nodes_from(unitig.graph)
    return out_graph


def find_unitigs(graph):
    """Finds unitigs of 2 or more nodes and replaces them by a subgraph containing the
    unitig nodes.

    Unitigs are described in :py:func:`find_unitig_from`.

    returns a graph containing unitig subgraphs
    """

    out_graph = graph.copy()
    visited_nodes = set()
    for start_node in graph:
        if start_node in visited_nodes:
            continue
        unitig = find_unitig_from(start_node, graph)
        unitig_graph_set = set(unitig.graph)
        assert visited_nodes & unitig_graph_set == set()
        visited_nodes |= unitig_graph_set

        out_graph = replace_unitig_nodes_with_unitig(out_graph, unitig)
    return out_graph


def collapse_kmer_unitigs(graph):
    """Collapses unitig kmers into a single graph node.

    All nodes have an attribute `repr` added which is the node label with the left k-1 letters
    removed unless the node has no incoming edges. In that case, `repr` is the full kmer.

    Unitigs are described in :py:func:`find_unitig_from`.

    :param graph:
    :return graph:
    """
    unitig_graph = find_unitigs(graph)
    out_graph = unitig_graph.copy()
    for node, data in unitig_graph.nodes.data():
        assert data['is_unitig']
        left_node = data['left_node']
        if graph.in_degree(left_node) > 0:
            short_kmer_name = [left_node[-1]]
        else:
            short_kmer_name = [left_node]
        for source, target in nx.edge_dfs(node, left_node):
            short_kmer_name.append(target[-1])
        short_kmer_name = ''.join(short_kmer_name)
        out_graph.nodes[node]['repr'] = short_kmer_name
    return out_graph
