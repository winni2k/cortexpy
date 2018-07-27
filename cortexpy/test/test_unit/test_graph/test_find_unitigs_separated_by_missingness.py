from cortexpy.test.driver.graph.find_unitigs import FindUnitigsTestDriver


class TestWithMissingEdge(object):
    def test_two_edges_and_two_edges_results_in_two_unitigs(self):
        # given
        driver = FindUnitigsTestDriver()
        driver.without_test_coverage()
        builder = driver.builder
        builder.graph.add_edge('AAA', 'AAC')
        builder.graph.add_edge('ACC', 'CCC')

        # when
        expect = driver.run()

        # then
        expect.has_n_nodes(2)
        expect.has_n_unitigs(2)
        expect.has_unitig_with_edges(('AAA', 'AAC')).with_left_node('AAA').with_right_node('AAC')
        expect.has_unitig_with_edges(('ACC', 'CCC')).with_left_node('ACC').with_right_node('CCC')

    def test_two_pairs_of_edges_separated_by_single_edge_results_in_two_unitigs(self):
        # given
        driver = FindUnitigsTestDriver()
        driver.without_test_coverage()
        builder = driver.graph_builder
        builder.with_colors(0, 1)
        driver.graph.add_path(['AAA', 'AAC', 'ACC', 'CCC'], color=0, coverage=1)
        builder.add_edge_with_color('AAA', 'AAC', 1)
        builder.add_edge_with_color('ACC', 'CCC', 1)

        # when
        expect = driver.run()

        # then
        expect.has_n_unitigs(2)
        expect.has_unitig_with_edges(('AAA', 'AAC')).with_left_node('AAA').with_right_node('AAC')
        expect.has_unitig_with_edges(('ACC', 'CCC')).with_left_node('ACC').with_right_node('CCC')

    def test_two_pairs_of_nodes_separated_by_single_edge_results_in_two_unitigs(self):
        # given
        d = FindUnitigsTestDriver()
        d.without_test_coverage()
        d.with_colors(0, 1)
        d.add_path(['AAA', 'AAC'], color=0, coverage=1)

        # when
        expect = d.run()

        # then
        expect.has_n_unitigs(2)
        expect.has_n_edges(3)
        expect.has_unitig_with_one_node('AAA')
        expect.has_unitig_with_one_node('AAC')
