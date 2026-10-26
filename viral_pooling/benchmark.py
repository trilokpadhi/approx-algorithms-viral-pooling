"""Controlled synthetic/public-HCV comparison; pooling remains simulated."""
import argparse
from itertools import combinations
import json
import math
from pathlib import Path
import random
import statistics
import time
from .design import gpda, validate_design
from .sequence import Pool, deconvolve, hamming, read_fasta


def load_dataset(folder, public=False):
    metadata = json.loads((folder / 'samples.json').read_text())
    sequences = read_fasta(folder / ('aligned.fasta' if public else 'synthetic_haplotypes.fasta'))
    samples = metadata['samples']
    if [s['id'] for s in samples] != list(range(len(samples))):
        raise ValueError('sample IDs must be consecutive integers')
    for sample in samples:
        keys, weights = sample['haplotypes'], sample['frequencies']
        if not keys or len(keys) != len(weights) or any(k not in sequences for k in keys):
            raise ValueError('invalid haplotype membership')
        if any(not math.isfinite(w) or w < 0 for w in weights) or sum(weights) <= 0:
            raise ValueError('invalid sampling weights')
    return sequences, samples


def dataset_statistics(sequences, samples):
    sets = [set(sequences[k] for k in s['haplotypes']) for s in samples]
    distances = [hamming(a, b) for seqs in sets for a, b in combinations(sorted(seqs), 2)]
    shared = sorted(seq for seq in set.union(*sets) if sum(seq in s for s in sets) > 1)
    return {'samples': len(samples), 'haplotype_records': len(sequences),
            'distinct_strings': len(set(sequences.values())),
            'aligned_length': len(next(iter(sequences.values()))),
            'haplotypes_per_sample': [len(s['haplotypes']) for s in samples],
            'distinct_strings_shared_between_samples': len(shared),
            'mean_within_sample_hamming': statistics.mean(distances) if distances else 0,
            'max_within_sample_hamming': max(distances, default=0)}


def run_trial(sequences, samples, seed=91, depth=80, w=2, max_w=4):
    if depth < 1:
        raise ValueError('positive depth required')
    n = len(samples)
    # Controlled design assumption: titers/epidemiological exclusions unavailable.
    g = {frozenset(pair) for pair in combinations(range(n), 2)}
    design, design_trace = gpda(n, g, 4)
    validate_design(n, g, design, 4)
    rng = random.Random(seed)
    pools, truth, reads = [], {}, {}
    for pool_id, members in enumerate(design):
        observations = {}
        for sample in sorted(members):
            source = samples[sample]
            selected = rng.choices(source['haplotypes'], weights=source['frequencies'], k=depth)
            for key in selected:
                obs = f'p{pool_id}_r{len(observations)}'
                observations[obs] = sequences[key]
                truth[obs] = sample
                reads[obs] = sequences[key]
        pools.append(Pool(members, observations))
    started = time.perf_counter()
    recovered, info = deconvolve(pools, w=w, max_w=max_w)
    elapsed = time.perf_counter() - started
    predicted = {key: set() for key in truth}
    confusion = [[0] * n for _ in range(n)]
    for sample, assigned in recovered.items():
        for key in assigned:
            if key not in truth:
                raise ValueError('deconvolution invented an observation')
            predicted[key].add(sample)
            confusion[truth[key]][sample] += 1
    singleton_samples = {next(iter(p)) for p in design if len(p) == 1}
    pooled_keys = {key for key, sample in truth.items() if sample not in singleton_samples}
    correct = sum(truth[key] in labels for key, labels in predicted.items())
    assigned = sum(bool(labels) for labels in predicted.values())
    entries = sum(map(len, predicted.values()))
    per_sample = []
    for sample in range(n):
        expected = {key for key, s in truth.items() if s == sample}
        output = set(recovered.get(sample, {}))
        right = expected & output
        observed_haplotypes = {reads[k] for k in expected}
        recovered_haplotypes = {reads[k] for k in right}
        per_sample.append({'sample': sample, 'source_sample': samples[sample].get('source_sample', str(sample)),
                           'true_observations': len(expected), 'assigned_observations': len(output),
                           'correct': len(right),
                           'precision': len(right) / len(output) if output else 0,
                           'recall': len(right) / len(expected) if expected else 0,
                           'observed_true_haplotypes': len(observed_haplotypes),
                           'correctly_recovered_observed_haplotypes': len(recovered_haplotypes)})
    return {'seed': seed, 'depth_per_sample_per_pool': depth, 'requested_w': w,
            'pool_memberships': [sorted(p) for p in design], 'pool_to_sample_ratio': len(design) / n,
            'observations': len(truth), 'assigned_observations': assigned,
            'assignment_entries': entries, 'correct_assignments': correct,
            'accuracy': correct / entries if entries else 0,
            'recall': correct / len(truth),
            'singleton_samples': sorted(singleton_samples),
            'non_singleton_sample_recall': sum(truth[k] in predicted[k] for k in pooled_keys) / len(pooled_keys) if pooled_keys else None,
            'unassigned': len(truth) - assigned,
            'ambiguous': sum(len(labels) > 1 for labels in predicted.values()),
            'macro_sample_recall': statistics.mean(s['recall'] for s in per_sample),
            'deconvolution_seconds': elapsed,
            'deconvolution': info, 'design_trace': design_trace,
            'per_sample': per_sample, 'confusion_matrix_true_rows_predicted_columns': confusion}


def aggregate(trials):
    return {metric: {'mean': statistics.mean(t[metric] for t in trials),
                     'min': min(t[metric] for t in trials),
                     'max': max(t[metric] for t in trials)}
            for metric in ('accuracy', 'recall', 'macro_sample_recall', 'non_singleton_sample_recall', 'unassigned', 'deconvolution_seconds')}


def write_report(report, path):
    synthetic, public = (report['datasets'][key] for key in ('synthetic', 'public_hcv'))
    lines = [
        '# Public HCV versus synthetic pooling benchmark', '',
        'Real source sequences; simulated pools and source-labelled read observations.', '',
        'Source: [Janiak et al. (2019)](https://doi.org/10.1038/s41598-019-49258-y); '
        '[Zenodo dataset, DOI 10.5281/zenodo.3266795](https://doi.org/10.5281/zenodo.3266795). '
        'Dataset creator: Kamila Caraballo Cortes; CC BY 4.0.', '',
        '## Protocol', '',
        'Selected the first eight chronic genotype-1 files in numeric specimen-ID order '
        'before measuring accuracy: 53, 80, 142, 153, 157, 165, 200, 201. '
        'All source haplotypes were retained, with header frequencies normalized within specimen. '
        'No additional frequency cutoff was applied. Source gaps were stripped; '
        'all sequences were realigned together with MAFFT 7.526, --auto --thread 1.', '',
        'Both datasets: eight samples, T=4, tabu GPDA, all pairs assumed compatible, '
        '80 independently drawn observations per included specimen per pool, '
        'seeds 91–95, initial W=2, retry through W=4 only on an empty recovered sample. '
        'The complete-graph assumption controls the design comparison; measured titers '
        'and epidemiological exclusions were not supplied. Both designs have four pools '
        'and 1,040 observations per run. Ground-truth source IDs are used only for evaluation.', '',
        '## Results', '',
        '| Measure | Synthetic | Public HCV |', '|---|---:|---:|',
    ]
    for label, key in [('Haplotype records', 'haplotype_records'), ('Distinct strings', 'distinct_strings'),
                       ('Alignment columns', 'aligned_length')]:
        lines.append(f'| {label} | {synthetic["statistics"][key]} | {public["statistics"][key]} |')
    for label, key in [('Mean read-assignment accuracy', 'accuracy'), ('Mean read recall', 'recall'),
                       ('Mean per-specimen recall (macro)', 'macro_sample_recall'),
                       ('Mean recall excluding singleton specimen', 'non_singleton_sample_recall')]:
        lines.append(f'| {label} | {100*synthetic["aggregate"][key]["mean"]:.2f}% | {100*public["aggregate"][key]["mean"]:.2f}% |')
    lines += ['', '| Seed | Synthetic correct / total | HCV correct / total | Final HCV W |',
              '|---|---:|---:|---:|']
    for left, right in zip(synthetic['trials'], public['trials']):
        lines.append(f'| {right["seed"]} | {left["correct_assignments"]}/{left["observations"]} | '
                     f'{right["correct_assignments"]}/{right["observations"]} | {right["deconvolution"]["w"]} |')
    fits = [fit for trial in public['trials'] for fit in trial['deconvolution']['clustering_trace']]
    lines += ['', 'Every observation was assigned exactly once in these runs. '
              'Accuracy is correct assignment entries / all assignment entries; recall is '
              'correctly recovered observation IDs / all input IDs. They coincide here. '
              'Repeated seeds are simulation replicates on the same specimens, not independent cohorts.', '',
              f'All {len(fits)} public-data clustering fits across the main runs and retries '
              f'reported outer convergence: {all(f["converged"] for f in fits)}; '
              f'inner convergence: {all(f["inner_converged"] for f in fits)}. '
              'Convergence does not imply biologically correct clustering.', '',
              '### Per-specimen recall, seed 91', '', '| Public specimen | Correct / true observations | Recall |',
              '|---|---:|---:|']
    for sample in public['trials'][0]['per_sample']:
        lines.append(f'| {sample["source_sample"]} | {sample["correct"]}/{sample["true_observations"]} | {100*sample["recall"]:.2f}% |')
    if 'public_w_sensitivity' in report:
        lines += ['', '### Fixed-W sensitivity, seed 91', '',
                  'These runs disable retries to isolate W. The main protocol above retains retries.', '',
                  '| Fixed W | Correct / total | Accuracy | Empty recovered specimen IDs |', '|---|---:|---:|---|']
        for trial in report['public_w_sensitivity']:
            empty = [public['trials'][0]['per_sample'][i]['source_sample'] for i in trial['deconvolution']['empty_samples']]
            lines.append(f'| {trial["requested_w"]} | {trial["correct_assignments"]}/{trial["observations"]} | '
                         f'{100*trial["accuracy"]:.2f}% | {", ".join(empty) or "none"} |')
    lines += ['', '## Limitations and interpretation', '',
              '- These are reconstructed, error-corrected HVR1 haplotypes, not raw FASTQ reads or experimentally pooled specimens. '
              'Public specimen membership is a simulation label, not independently validated laboratory truth.',
              '- Frequency estimates come from source FASTA headers. Independent resampling does not reproduce PCR amplification, '
              'titer imbalance, contamination, dropout, read errors, or pooling-dependent bias.',
              '- Only eight chronic genotype-1 specimens from one study are included. Numeric-ID selection is deterministic '
              'but is not random or representative. No parameters or specimens were selected to maximize the measured score.',
              '- The public alignment has 180 columns versus 306 in the synthetic demo. Diversity, length, '
              'haplotype counts and frequency profiles all differ; the accuracy gap cannot be attributed to one factor.',
              '- One distinct aligned string occurs in multiple public specimens. Identical strings alone cannot uniquely '
              'establish an originating specimen; the simulated observation labels still retain their original source.',
              '- Specimen 201 occupies a singleton pool and is recovered trivially. The table therefore also reports recall '
              'excluding that specimen (synthetic sample 7 is excluded in the matching metric). '
              'The 4/8 pool ratio is a count reduction, not a demonstrated 50% monetary saving: total sequencing depth is not fixed.',
              '- Empty-sample retries do not detect nonempty but incorrect assignments. In seed 91, specimen 80 has very poor '
              'recall despite a nonempty result. Fixed-W sensitivity also shows parameter dependence.',
              '- The chapter interpretations and repairs documented in the main README remain in effect: queue refinement, '
              'stopping-condition repairs, bounded search, and the printed unweighted consensus update. '
              'No clustering algorithm was changed to improve these results.',
              '- Five sampling seeds describe variability conditional on these same source sequences. '
              'They are not a confidence interval for other patients, studies or laboratory workflows.',
              '- This exposes a limitation of the current implementation under a more difficult input distribution. '
              'It neither reproduces nor disproves the chapter’s published experimental results.', '',
              '## Reproduce', '', 'From the repository root, with the committed aligned data (no network or aligner required):', '',
              '```sh', 'python3 -m viral_pooling.benchmark', 'python3 -m unittest discover -s tests -v', '```', '',
              'For source download and alignment regeneration, see data/public_hcv/README.md. '
              'Machine-readable metrics, confusion matrices and convergence traces are in comparison.json. '
              'Timing is diagnostic only and depends on hardware and interpreter.', '']
    path.write_text('\n'.join(lines))


def compare(synthetic, public, output, seeds=(91, 92, 93, 94, 95), depth=80, sensitivity=True):
    output.mkdir(parents=True, exist_ok=True)
    report = {'protocol': {'seeds': list(seeds), 'depth': depth, 'pool_size_limit': 4,
                           'solver': 'tabu', 'w': 2, 'max_w': 4,
                           'compatibility': 'All pairs assumed compatible for BOTH datasets; no measured public titers or epidemiological exclusions used.',
                           'pooling': 'Independently simulated from source haplotype weights; no added sequence errors.'},
              'datasets': {}}
    for name, folder, is_public in [('synthetic', synthetic, False), ('public_hcv', public, True)]:
        sequences, samples = load_dataset(folder, is_public)
        trials = []
        for seed in seeds:
            trial = run_trial(sequences, samples, seed, depth)
            trials.append(trial)
            print(f'{name} seed={seed} accuracy={trial["accuracy"]:.4f} recall={trial["recall"]:.4f}', flush=True)
        report['datasets'][name] = {'statistics': dataset_statistics(sequences, samples),
                                    'trials': trials, 'aggregate': aggregate(trials)}
        # Checkpoint so completed trials survive an interrupted later experiment.
        (output / 'comparison.json').write_text(json.dumps(report, indent=2) + '\n')
    if sensitivity:
        sequences, samples = load_dataset(public, True)
        report['public_w_sensitivity'] = []
        for w in (1, 2, 3):
            trial = run_trial(sequences, samples, seeds[0], depth, w=w, max_w=w)
            report['public_w_sensitivity'].append(trial)
            print(f'public_hcv fixed W={w} accuracy={trial["accuracy"]:.4f} recall={trial["recall"]:.4f}', flush=True)
    (output / 'comparison.json').write_text(json.dumps(report, indent=2) + '\n')
    write_report(report, output / 'README.md')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--synthetic', type=Path, default=Path('data'))
    parser.add_argument('--public', type=Path, default=Path('data/public_hcv'))
    parser.add_argument('--output', type=Path, default=Path('results/comparison'))
    parser.add_argument('--seeds', nargs='+', type=int, default=[91, 92, 93, 94, 95])
    parser.add_argument('--depth', type=int, default=80)
    parser.add_argument('--skip-sensitivity', action='store_true')
    args = parser.parse_args()
    compare(args.synthetic, args.public, args.output, tuple(args.seeds), args.depth, not args.skip_sensitivity)


if __name__ == '__main__':
    main()
