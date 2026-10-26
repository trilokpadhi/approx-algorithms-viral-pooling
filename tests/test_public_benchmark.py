import hashlib
import json
import tempfile
from pathlib import Path
import unittest
from unittest.mock import patch
from viral_pooling.benchmark import load_dataset, run_trial, dataset_statistics
from viral_pooling.public_data import parse_frequency_fasta
from viral_pooling.sequence import read_fasta

ROOT = Path(__file__).resolve().parents[1]
PUBLIC = ROOT / 'data' / 'public_hcv'


class PublicDataTests(unittest.TestCase):
    def test_source_checksums_and_alignment_preserve_bases(self):
        provenance = json.loads((PUBLIC / 'provenance.json').read_text())
        metadata = json.loads((PUBLIC / 'samples.json').read_text())
        aligned = read_fasta(PUBLIC / 'aligned.fasta')
        original = {}
        for item in provenance['files']:
            body = (PUBLIC / item['local_filename']).read_bytes()
            self.assertEqual(hashlib.sha256(body).hexdigest(), item['sha256'])
            self.assertEqual('md5:' + hashlib.md5(body).hexdigest(), item['md5'])
            for name, seq, _ in parse_frequency_fasta(PUBLIC / item['local_filename']):
                original[(item['sample_id'], name)] = seq.replace('-', '')
        for item in metadata['haplotype_mapping']:
            self.assertEqual(aligned[item['haplotype']].replace('-', ''),
                             original[(item['source_sample'], item['original_header_id'])])
        self.assertEqual(hashlib.sha256((PUBLIC / 'aligned.fasta').read_bytes()).hexdigest(),
                         metadata['alignment']['sha256'])

    def test_membership_and_frequency_metadata(self):
        sequences, samples = load_dataset(PUBLIC, public=True)
        self.assertEqual([s['source_sample'] for s in samples], ['53', '80', '142', '153', '157', '165', '200', '201'])
        for s in samples:
            self.assertAlmostEqual(sum(s['frequencies']), 1)
        stats = dataset_statistics(sequences, samples)
        self.assertEqual(stats['haplotype_records'], 338)
        self.assertEqual(stats['aligned_length'], 180)
        self.assertEqual(stats['distinct_strings_shared_between_samples'], 1)

    def test_malformed_frequency_metadata(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'source.fasta'
            for header in ('>r', '>r freq=-1', '>r freq=0'):
                path.write_text(header + '\nACGT\n')
                with self.assertRaises(ValueError):
                    parse_frequency_fasta(path)

    def test_metrics_detect_wrong_and_unassigned_reads(self):
        # Two singleton pools, one observation each. Deliberately return the
        # second source's read under the first label and omit the first read.
        sequences = {'a': 'AAAA', 'b': 'CCCC'}
        samples = [{'id': i, 'haplotypes': [key], 'frequencies': [1]}
                   for i, key in enumerate(('a', 'b'))]
        def incorrect(pools, **kwargs):
            wrong = next(iter(pools[1].reads))
            return {0: {wrong: pools[1].reads[wrong]}, 1: {}}, {}
        with patch('viral_pooling.benchmark.deconvolve', side_effect=incorrect):
            result = run_trial(sequences, samples, depth=1)
        self.assertEqual(result['correct_assignments'], 0)
        self.assertEqual(result['unassigned'], 1)
        self.assertEqual(result['accuracy'], 0)
        self.assertEqual(result['recall'], 0)

    def test_public_end_to_end_conservation_and_diagnostics(self):
        sequences, samples = load_dataset(PUBLIC, public=True)
        result = run_trial(sequences, samples, seed=91)
        self.assertEqual(result['observations'], 1040)
        self.assertEqual(result['assignment_entries'] + result['unassigned'], 1040)
        self.assertEqual(result['ambiguous'], 0)
        self.assertEqual(len(result['per_sample']), 8)
        self.assertEqual(sum(s['correct'] for s in result['per_sample']), result['correct_assignments'])
        self.assertTrue(result['deconvolution']['clustering_trace'])
        self.assertTrue(all('effective_k' in t for t in result['deconvolution']['clustering_trace']))
