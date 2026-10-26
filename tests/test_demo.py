import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from viral_pooling.__main__ import generate_data, demo
from viral_pooling.sequence import read_fasta


class DemoTests(unittest.TestCase):
    def test_end_to_end_both_solvers(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            generate_data(root / 'data')
            self.assertEqual(len(read_fasta(root / 'data/synthetic_haplotypes.fasta')), 32)
            for solver in ('tabu', 'exact'):
                with contextlib.redirect_stdout(io.StringIO()):
                    demo(root / 'data', root / solver, solver=solver)
                result = json.loads((root / solver / 'result.json').read_text())
                self.assertEqual(result['samples'], 8)
                self.assertEqual(result['accuracy'], 1)
                self.assertEqual(result['recall'], 1)
                self.assertEqual(result['ambiguous_observations'], 0)
                self.assertEqual(result['deconvolution']['empty_samples'], [])
                self.assertEqual(len(list((root / solver).glob('recovered_*.fasta'))), 8)
