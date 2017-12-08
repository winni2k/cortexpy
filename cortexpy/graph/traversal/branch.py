import attr

from cortexpy.graph.serializer import EdgeTraversalOrientation, SERIALIZER_GRAPH


@attr.s(slots=True)
class Branch(object):
    ra_parser = attr.ib()
    traversal_color = attr.ib(0)
    graph = attr.ib(attr.Factory(SERIALIZER_GRAPH))
    kmer = attr.ib(init=False)
    kmer_string = attr.ib(init=False)
    orientation = attr.ib(init=False)

    def traverse_from(self, kmer_string, *, orientation=EdgeTraversalOrientation.original):
        self.graph = SERIALIZER_GRAPH()
        self.kmer_string = kmer_string
        self.orientation = orientation
        self._add_kmer_string_to_graph_and_get_kmer()
        while True:
            oriented_edge_set = self.kmer.edges[self.traversal_color].oriented(self.orientation)
            if self._get_num_neighbors(oriented_edge_set) != 1:
                break
            self._add_next_kmer_string_to_graph_and_get_next_kmer(oriented_edge_set)
        return self.graph

    def _get_num_neighbors(self, oriented_edge_set):
        if self.kmer.kmer != self.kmer_string:
            oriented_edge_set = oriented_edge_set.other_orientation()
        return oriented_edge_set.num_neighbor()

    def _add_next_kmer_string_to_graph_and_get_next_kmer(self, oriented_edge_set):
        prev_kmer_string = self.kmer_string
        self.kmer_string = oriented_edge_set.neighbor_kmer_strings(prev_kmer_string)[0]
        self.graph.add_edge(prev_kmer_string, self.kmer_string, key=self.traversal_color)
        self._add_kmer_string_to_graph_and_get_kmer()

    def _add_kmer_string_to_graph_and_get_kmer(self):
        self.kmer = self.ra_parser.get_kmer_for_string(self.kmer_string)
        self.graph.add_node(self.kmer_string, kmer=self.kmer)
