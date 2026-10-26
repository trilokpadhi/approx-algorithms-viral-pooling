import tempfile
from pathlib import Path
import unittest
from viral_pooling.sequence import Pool, kgem, generalized_intersection, deconvolve, read_fasta, write_fasta


class SequenceTests(unittest.TestCase):
    def test_fasta(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / 'reads.fa'
            write_fasta(path, {'a': 'AC-G', 'b': 'ACTG'})
            self.assertEqual(read_fasta(path), {'a': 'AC-G', 'b': 'ACTG'})
            path.write_text('>a\nAC\n>b\nA\n')
            with self.assertRaises(ValueError):
                read_fasta(path)

    def test_kgem_known_mixture(self):
        reads = ['A' * 306] * 30 + ['C' * 306] * 10
        fit = kgem(reads, 2)
        self.assertTrue(fit['converged'])
        self.assertEqual(len(set(fit['labels'][:30])), 1)
        self.assertNotEqual(fit['labels'][0], fit['labels'][-1])
        self.assertAlmostEqual(max(fit['weights']), .75)
        self.assertTrue(__import__('math').isfinite(fit['log_likelihood']))

    def test_sampling_variants_generalized_intersection(self):
        p = Pool(frozenset((0, 1)), {'a': 'AAAAAA', 'b': 'CCCCCC'})
        q = Pool(frozenset((1, 2)), {'c': 'CCCCCT', 'd': 'GGGGGG'})
        gi = generalized_intersection(p, q, w=1)
        self.assertEqual(set(gi.reads), {'b', 'c'})
        self.assertEqual(gi.samples, frozenset((1,)))

    def test_chapter_example_one(self):
        p = Pool(frozenset((0, 1)), {'a': 'AAAAAA', 'b': 'CCCCCC'})
        q = Pool(frozenset((1, 2)), {'c': 'CCCCCC', 'd': 'GGGGGG'})
        recovered, info = deconvolve([p, q], w=1)
        self.assertEqual(set(recovered[0]), {'a'})
        self.assertEqual(set(recovered[1]), {'b', 'c'})
        self.assertEqual(set(recovered[2]), {'d'})
        self.assertEqual(info['empty_samples'], [])

    def test_empty_retry_and_unseparated(self):
        recovered, info = deconvolve([Pool(frozenset((0,)), {})], max_w=3)
        self.assertEqual(info['w'], 3)
        self.assertEqual(info['empty_samples'], [0])
        with self.assertRaises(ValueError):
            deconvolve([Pool(frozenset((0, 1)), {'x': 'AAAA'})])

    def test_k_cap_and_input_validation(self):
        self.assertEqual(kgem(['AAAA', 'AAAA'], 9)['effective_k'], 1)
        for reads in ([], ['AA', 'AAA'], ['NNNN']):
            with self.assertRaises(ValueError):
                kgem(reads, 2)
