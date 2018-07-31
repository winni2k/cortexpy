import attr
from hypothesis import given, strategies as s

from cortexpy.constants import EdgeTraversalOrientation, EngineTraversalOrientation
from cortexpy.graph.traversal.branch import Traversed, Queuer, SERIALIZER_GRAPH


@attr.s(slots=True)
class TraversedBranchBuilder(object):
    traversed_branch = attr.ib(init=False)

    def __attrs_post_init__(self):
        self.traversed_branch = Traversed(SERIALIZER_GRAPH(),
                                          EdgeTraversalOrientation.original)

    def with_first_kmer(self, kmer):
        self.traversed_branch.first_kmer_string = kmer
        return self

    def with_last_kmer(self, kmer):
        self.traversed_branch.last_kmer_string = kmer
        return self

    def with_neighbors(self, *neighbors):
        self.traversed_branch.neighbor_kmer_strings = neighbors
        return self

    def with_reverse_neighbors(self, *neighbors):
        self.traversed_branch.reverse_neighbor_kmer_strings = neighbors
        return self

    def with_traversal_orientation(self, orientation):
        self.traversed_branch.orientation = EdgeTraversalOrientation[orientation]
        return self

    def build(self):
        return self.traversed_branch


@attr.s(slots=True)
class BranchQueuerDriver(object):
    branches_in_queue = attr.ib(attr.Factory(list))
    branch_builders_to_queue = attr.ib(attr.Factory(list))
    engine_traversal_orientation = attr.ib(EngineTraversalOrientation.original)
    traversal_colors = attr.ib((0,))

    def queue_from_traversed_branch(self):
        traversed_branch_builder_to_queue = TraversedBranchBuilder()
        self.branch_builders_to_queue.append(traversed_branch_builder_to_queue)
        return traversed_branch_builder_to_queue

    def with_engine_traversal_orientation(self, orientation):
        self.engine_traversal_orientation = EngineTraversalOrientation[orientation]
        return self

    def run(self):
        queuer = Queuer(self.branches_in_queue,
                        traversal_colors=self.traversal_colors,
                        engine_orientation=self.engine_traversal_orientation)
        for traversed_branch_builder in self.branch_builders_to_queue:
            traversed_branch = traversed_branch_builder.build()
            queuer.add_from_branch(traversed_branch)
        return BranchQueueExpectation(queuer)


@attr.s(slots=True)
class BranchQueueExpectation(object):
    queuer = attr.ib()
    queue = attr.ib(init=False)

    def __attrs_post_init__(self):
        self.queue = self.queuer.queue

    def has_n_setups(self, n):
        assert n == len(self.queue)
        return self

    def has_nth_setup(self, n, start_string, orientation, connecting_node):
        assert len(self.queue) > n
        setup = self.queue[n]
        assert setup.start_string == start_string
        assert setup.orientation == orientation
        assert setup.connecting_node == connecting_node
        return self


class Test(object):
    @given(s.sampled_from(('original', 'reverse')), s.booleans())
    def test_enqueues_from_a_traversed_branch(self, traversed_branch_orientation_string,
                                              traverse_both):
        # given
        engine_traversal_orientation_string = traversed_branch_orientation_string
        if traverse_both:
            engine_traversal_orientation_string = 'both'
        traversed_branch_orientation = EdgeTraversalOrientation[traversed_branch_orientation_string]

        driver = BranchQueuerDriver()
        (driver
         .with_engine_traversal_orientation(engine_traversal_orientation_string)
         .queue_from_traversed_branch()
         .with_first_kmer(0)
         .with_last_kmer(1)
         .with_neighbors(2)
         .with_reverse_neighbors(3)
         .with_traversal_orientation(traversed_branch_orientation_string))

        expected_setups = [(2, traversed_branch_orientation, 1)]
        if engine_traversal_orientation_string == 'both':
            expected_setups.append(
                (3, EdgeTraversalOrientation.other(traversed_branch_orientation), 1))

        # when
        expect = driver.run()

        # then

        for n, setup in enumerate(expected_setups):
            expect.has_nth_setup(n,
                                 start_string=setup[0],
                                 orientation=setup[1],
                                 connecting_node=setup[2])

    @given(s.sampled_from(('original', 'reverse')), s.booleans())
    def test_does_not_enqueue_from_a_seen_traversed_branch(self,
                                                           traversed_branch_orientation_string,
                                                           traverse_both):
        # given
        engine_traversal_orientation_string = traversed_branch_orientation_string
        if traverse_both:
            engine_traversal_orientation_string = 'both'
        traversed_branch_orientation = EdgeTraversalOrientation[traversed_branch_orientation_string]

        driver = BranchQueuerDriver()
        driver.with_engine_traversal_orientation(engine_traversal_orientation_string)
        expected_setups = []
        for _ in range(2):
            driver.queue_from_traversed_branch() \
                .with_first_kmer(0) \
                .with_last_kmer(1) \
                .with_neighbors(2) \
                .with_reverse_neighbors(3) \
                .with_traversal_orientation(traversed_branch_orientation_string)

        expected_setups.append((2, traversed_branch_orientation, 1))
        if engine_traversal_orientation_string == 'both':
            expected_setups.append(
                (3, EdgeTraversalOrientation.other(traversed_branch_orientation), 1))

        # when
        expect = driver.run()

        # then
        for n, setup in enumerate(expected_setups):
            expect.has_nth_setup(n,
                                 start_string=setup[0],
                                 orientation=setup[1],
                                 connecting_node=setup[2])
        expect.has_n_setups(len(expected_setups))
