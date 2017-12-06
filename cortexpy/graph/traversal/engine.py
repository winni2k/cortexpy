import attr
import networkx as nx


@attr.s(slots=True)
class Engine(object):
    place_holder = attr.ib()

    def traverse(self):
        return nx.MultiDiGraph()
