import json

import attr
import networkx as nx
from networkx.readwrite import json_graph

from cortexpy.graph.interactor import CortexDiGraph, Interactor
from .unitig import UnitigCollapser


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
        self._make_kmer_graph_consistent()
        self._collapse_kmer_graph()
        self._make_unitig_graph_json_representable()

    def _make_kmer_graph_consistent(self):
        if isinstance(self.graph, CortexDiGraph):
            self.graph = Interactor(self.graph).make_graph_nodes_consistent().graph

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
            coverage = [list(c) for c in zip(*coverage)]
            node_data['coverage'] = coverage
