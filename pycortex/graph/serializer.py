import json
import attr
import networkx as nx
from networkx.readwrite import json_graph


@attr.s(slots=True)
class Serializer(object):
    """A class for serializing kmer graphs."""
    graph = attr.ib()
    collapse_kmer_unitigs = attr.ib(False)

    def to_json_serializable(self):
        graph = self.graph
        if self.collapse_kmer_unitigs:
            graph = collapse_kmer_unitigs(graph)
        graph = nx.convert_node_labels_to_integers(graph, label_attribute='node_key')
        return make_graph_json_representable(graph)

    def to_json(self):
        graph = self.to_json_serializable()
        serializable = json_graph.node_link_data(graph, attrs={'link': 'edges'})
        return json.dumps(serializable)


def make_graph_json_representable(graph):
    """
    1. Removes the kmer object
    2. Converts graph node keys to their reprs
    3. Sets every node to missing or not missing
    """
    graph = graph.copy()
    for node, node_data in graph.nodes.items():
        kmer = node_data.pop('kmer', None)
        if kmer is not None:
            node_data['coverage'] = list(kmer.coverage)
        if isinstance(node_data.get('node_key'), nx.Graph):
            node_data['node_key'] = repr(node_data['node_key'])
        if 'is_missing' not in node_data:
            if kmer is None:
                node_data['is_missing'] = True
            else:
                raise ValueError
    for source, target, edge_data in graph.edges.data():
        is_missing = bool(graph.node[source]['is_missing'] or graph.node[target]['is_missing'])
        graph[source][target]['is_missing'] = is_missing
    return graph


def find_unitig_from(start_node, graph):
    """Find a unitig that contains :param start_node:

    # from https://github.com/mcveanlab/mccortex/wiki/unitig
    Unitig definition: A maximal run of de Bruijn nodes (kmers, each with an orientation)
    where all nodes must have in-degree and out-degree exactly 1, with the exception of in-degree
    of the first node and out-degree of the last node (they may have any number).
    Additionally, a unitig must not visit a kmer more than once.
    """
    unitig_graph = nx.DiGraph()
    left_node = right_node = start_node
    seen = {start_node}
    for source, target in nx.edge_dfs(graph, start_node, orientation='forward'):
        if graph.in_degree(target) != 1 or target in seen:
            break
        right_node = target
        unitig_graph.add_edge(source, target)
        seen.add(target)
        if graph.out_degree(target) != 1:
            break
    for source, target, orientation in nx.edge_dfs(graph, start_node, orientation='reverse'):
        if graph.out_degree(source) != 1 or source in seen:
            break
        left_node = source
        unitig_graph.add_edge(source, target)
        seen.add(source)
        if graph.in_degree(source) != 1:
            break
    return unitig_graph, left_node, right_node


def find_unitigs(graph):
    """Finds unitigs of 2 or more nodes and replaces them by a subgraph containing the
    unitig nodes.

    Unitigs are described in :py:func:`find_unitig_from`.

    returns a graph containing unitig subgraphs
    """

    unitigable_nodes = set()
    graph = graph.copy()
    graph_no_miss_edges = graph.copy()
    for source, target, data in graph.edges.data():
        if data.get('is_missing'):
            source_missing = bool(graph.node[source].get('is_missing'))
            target_missing = bool(graph.node[target].get('is_missing'))
            if not (source_missing and target_missing):
                graph_no_miss_edges.remove_edge(source, target)
    for node in graph:
        if graph_no_miss_edges.in_degree(node) < 2 and graph_no_miss_edges.out_degree(node) < 2:
            unitigable_nodes.add(node)

    visited_nodes = set()
    for start_node in unitigable_nodes:
        if start_node in visited_nodes:
            continue
        unitig_graph, left_node, right_node = find_unitig_from(start_node, graph_no_miss_edges)
        unitig_graph_set = set(unitig_graph)
        assert visited_nodes & unitig_graph_set == set()
        visited_nodes |= unitig_graph_set
        if len(unitig_graph) > 1:
            graph.add_node(unitig_graph, is_unitig=True, left_node=left_node, right_node=right_node,
                           is_missing=graph_no_miss_edges.node[left_node].get('is_missing', False))
            for source, _ in graph.in_edges(left_node):
                graph.add_edge(source, unitig_graph)
                if source == right_node:
                    graph.add_edge(unitig_graph, unitig_graph)
            for _, target in graph.out_edges(right_node):
                graph.add_edge(unitig_graph, target)
                if target == left_node:
                    graph.add_edge(unitig_graph, unitig_graph)
            graph.remove_nodes_from(unitig_graph)
    return graph


def non_missing_in_degree(graph, node):
    in_degree = graph.in_degree(node)
    for edge in graph.in_edges(node):
        if graph.get_edge_data(*edge).get('is_missing'):
            in_degree -= 1
    return in_degree


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

        if data.get('is_unitig'):
            left_node = data['left_node']
            if non_missing_in_degree(graph, left_node) > 0:
                short_kmer_name = [left_node[-1]]
            else:
                short_kmer_name = [left_node]
            for source, target in nx.edge_dfs(node, left_node):
                short_kmer_name.append(target[-1])
            short_kmer_name = ''.join(short_kmer_name)
        else:
            if non_missing_in_degree(graph, node) > 0:
                short_kmer_name = node[-1]
            else:
                short_kmer_name = node

        out_graph.nodes[node]['repr'] = short_kmer_name
    return out_graph
