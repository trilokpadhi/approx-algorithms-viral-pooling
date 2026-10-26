import itertools
import random
import unittest
from viral_pooling.design import (compatibility, cut, objective, gain, exact_clique_cut,
                                 gpda, validate_design, tabu_clique_cut, is_clique, maxcut_half)


class DesignTests(unittest.TestCase):
    def test_titer_and_epidemiology(self):
        self.assertEqual(compatibility([1, 3, 8], 2), {frozenset((0, 1))})
        self.assertEqual(compatibility([1, 3], 2, [(0, 1)]), set())

    def test_equations_exhaustive(self):
        n = 4
        pairs = list(itertools.combinations(range(n), 2))
        h = {frozenset(e) for e in pairs[::2]}
        for mask in range(1 << len(pairs)):
            g = {frozenset(e) for i, e in enumerate(pairs) if mask >> i & 1}
            for bits in range(1 << n):
                a = frozenset(v for v in range(n) if bits >> v & 1)
                self.assertEqual(objective(a, g, h) >= 0, is_clique(a, g))
                for v in range(n):
                    self.assertEqual(gain(v, a, n, g, h), objective(a ^ {v}, g, h) - objective(a, g, h))

    def test_oracle_known_cut(self):
        n = 6
        g = {frozenset(e) for e in itertools.combinations(range(n), 2)}
        self.assertEqual(cut(exact_clique_cut(n, g, g, 2), g), 8)

    def test_small_graphs_all_feasible(self):
        n = 4
        pairs = list(itertools.combinations(range(n), 2))
        for bits in range(1 << len(pairs)):
            g = {frozenset(e) for i, e in enumerate(pairs) if bits >> i & 1}
            for limit in (1, 2, 3, 4):
                for solver in ('exact', 'tabu'):
                    pools, trace = gpda(n, g, limit, solver)
                    validate_design(n, g, pools, limit)
                    self.assertTrue(all(t['separated'] > 0 for t in trace))
                    self.assertLessEqual(len(pools), n * (n - 1) // 2 + 1)

    def test_tabu_vs_oracle_random(self):
        rng = random.Random(8)
        for _ in range(40):
            n = 7
            pairs = list(itertools.combinations(range(n), 2))
            g = {frozenset(e) for e in pairs if rng.random() < .6}
            h = {frozenset(e) for e in pairs if rng.random() < .6}
            a, _ = tabu_clique_cut(n, g, h, 3)
            oracle = exact_clique_cut(n, g, h, 3)
            self.assertTrue(is_clique(a, g))
            self.assertLessEqual(len(a), 3)
            self.assertLessEqual(cut(a, h), cut(oracle, h))
            self.assertGreater(cut(a, h), 0)
            self.assertGreaterEqual(cut(maxcut_half(n, h), h), len(h) / 2)

    def test_search_bound_fallback(self):
        h = {frozenset(e) for e in itertools.combinations(range(6), 2)}
        a, info = tabu_clique_cut(6, set(), h, 2, max_steps=1, max_retries=1)
        self.assertTrue(info['fallback'])
        self.assertEqual(len(a), 1)
        self.assertGreater(cut(a, h), 0)

    def test_boundary(self):
        self.assertEqual(gpda(0, set(), 1)[0], [])
        self.assertEqual(gpda(1, set(), 1)[0], [frozenset((0,))])
        with self.assertRaises(ValueError):
            gpda(2, set(), 0)
