import attr

from cortexpy.constants import EdgeTraversalOrientation


@attr.s(slots=True)
class OrientedGraphFuncs(object):
    graph = attr.ib()
    orientation = attr.ib()
    color = attr.ib(None)
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
