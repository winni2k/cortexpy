import attr
import collections
import networkx as nx

import cortexpy.graph.traversal
from cortexpy.graph.serializer import SERIALIZER_GRAPH


@attr.s(slots=True)
class BranchTraversalSetup(object):
    start_string = attr.ib()
    connecting_node = attr.ib(None)


@attr.s(slots=True)
class Engine(object):
    ra_parser = attr.ib()
    traversal_color = attr.ib(0)
    graph = attr.ib(attr.Factory(SERIALIZER_GRAPH))
    max_nodes = attr.ib(1000)
    branch_queue = attr.ib(attr.Factory(collections.deque))

    def traverse_from(self, start_string):
        self.graph = SERIALIZER_GRAPH()
        branch_traverser = cortexpy.graph.traversal.Branch(self.ra_parser, self.traversal_color)
        self.branch_queue.append(BranchTraversalSetup(start_string))
        while 0 < len(self.branch_queue):
            setup = self.branch_queue.popleft()
            new_branch = branch_traverser.traverse_from(setup.start_string, parent_graph=self.graph)
            if self.max_nodes < len(self.graph) + len(new_branch.graph):
                raise OverflowError
            self.graph = nx.union(self.graph, new_branch.graph)
            if setup.connecting_node is not None and new_branch.first_kmer_string is not None:
                self.graph.add_edge(setup.connecting_node,
                                    new_branch.first_kmer_string,
                                    key=self.traversal_color)
            self._enqueue_branch_traversals(new_branch)
        return self.graph

    def _enqueue_branch_traversals(self, new_branch):
        for kmer_string in new_branch.neighbor_kmer_strings:
            if kmer_string in self.graph:
                self.graph.add_edge(new_branch.last_kmer_string,
                                    kmer_string,
                                    key=self.traversal_color)
            else:
                self.branch_queue.append(
                    BranchTraversalSetup(kmer_string,
                                         connecting_node=new_branch.last_kmer_string)
                )
