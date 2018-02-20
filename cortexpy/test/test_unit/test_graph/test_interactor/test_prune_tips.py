import networkx as nx

from cortexpy.graph import interactor


def test_prunes_three_tips_of_length_1():
    # given
    graph = nx.MultiDiGraph()
    graph.add_path([0, 1, 2])
    graph.add_path([0, 1, 3])

    # when
    graph = interactor.Interactor(graph, colors=None).prune_tips_less_than(2).graph

    # then
    assert 1 == len(graph)
    assert {1} == set(graph.nodes)


def test_prunes_one_tip_of_length_1():
    # given
    graph = nx.MultiDiGraph()
    graph.add_path([0, 1, 2, 3])
    graph.add_path([0, 1, 4, 5])

    # when
    graph = interactor.Interactor(graph, colors=None).prune_tips_less_than(2).graph

    # then
    assert 5 == len(graph)
    assert set(range(1, 6)) == set(graph.nodes)
