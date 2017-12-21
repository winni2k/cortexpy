import attr
import networkx as nx
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import logging

logger = logging.getLogger(__name__)


@attr.s(slots=True)
class Interactor(object):
    graph = attr.ib()
    colors = attr.ib()

    def add_edge_to_graph_for_kmer_pair(self, kmer1, kmer2, kmer1_string, kmer2_string):
        first_is_revcomp = bool(kmer1.kmer != kmer1_string)
        for color in self.colors:
            if first_is_revcomp:
                add_edge = kmer1.has_incoming_edge_from_kmer_in_color(kmer2, color)
            else:
                add_edge = kmer1.has_outgoing_edge_to_kmer_in_color(kmer2, color)
            if add_edge:
                self.graph.add_edge(kmer1_string, kmer2_string, key=color)


def make_copy_of_color(graph, color, include_self_refs=False):
    """Makes a copy of graph, but only copies over links with key=color.
    Only copies over nodes that are linked by a link with key=color.
    """
    out_graph = graph.fresh_copy()
    for u, v, key, data in graph.edges(keys=True, data=True):
        if key == color:
            if u == v and not include_self_refs:
                continue
            out_graph.add_edge(u, v, key=key, **data)
    for node in list(out_graph):
        out_graph.add_node(node, **graph.node[node])
    return out_graph


def convert_kmer_path_to_contig(path):
    if len(path) == 0:
        return ''
    logger.info(path)
    contig = [path[0]]
    for kmer in path[1:]:
        contig.append(kmer[-1])
    return ''.join(contig)


@attr.s(slots=True)
class Contigs(object):
    graph = attr.ib()
    color = attr.ib()

    def all_simple_paths(self):
        graph = make_copy_of_color(self.graph, self.color, include_self_refs=False)
        logger.info(graph.nodes)
        idx = 0
        for source in graph.nodes():
            if graph.in_degree(source) > 0:
                continue
            for target in graph.nodes():
                if source == target or graph.out_degree(target) != 0:
                    continue
                for path in nx.all_simple_paths(graph, source=source, target=target):
                    yield SeqRecord(Seq(convert_kmer_path_to_contig(path)), id=str(idx),
                                    description='')
                    idx += 1
