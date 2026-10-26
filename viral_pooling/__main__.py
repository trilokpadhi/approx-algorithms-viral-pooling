import argparse
from collections import Counter
import json
from pathlib import Path
import random
from .design import compatibility, gpda, validate_design
from .sequence import Pool, read_fasta, write_fasta, deconvolve


def generate_data(folder, seed=37, n=8, length=306):
    """Artificial HCV-HVR1-length amplicon model; no biological provenance."""
    folder.mkdir(parents=True, exist_ok=True)
    rng = random.Random(seed)
    ancestor = ''.join(rng.choices('ACGT', k=length))
    records, manifest = {}, {'seed': seed, 'synthetic': True, 'samples': []}
    def mutate(seq, positions):
        out = list(seq)
        for pos in positions:
            out[pos] = rng.choice([a for a in 'ACGT' if a != out[pos]])
        return ''.join(out)
    for i in range(n):
        founder = mutate(ancestor, rng.sample(range(length), 45))
        variants = [founder] + [mutate(founder, rng.sample(range(length), j)) for j in (2, 4, 6)]
        ids = []
        for j, seq in enumerate(variants):
            key = f's{i}_h{j}'
            ids.append(key)
            records[key] = seq
        manifest['samples'].append({'id': i, 'titer': 10 + (i % 3),
                                    'haplotypes': ids, 'frequencies': [0.65, 0.2, 0.1, 0.05]})
    write_fasta(folder / 'synthetic_haplotypes.fasta', records)
    (folder / 'samples.json').write_text(json.dumps(manifest, indent=2) + '\n')


def demo(data, output, solver='tabu', seed=91, depth=80):
    if depth < 1:
        raise ValueError('depth must be positive')
    manifest = json.loads((data / 'samples.json').read_text())
    records = read_fasta(data / 'synthetic_haplotypes.fasta')
    samples = manifest['samples']
    g = compatibility([s['titer'] for s in samples], 3)
    design, trace = gpda(len(samples), g, 4, solver=solver)
    signatures = validate_design(len(samples), g, design, 4)
    rng = random.Random(seed)
    pools, truth = [], {}
    output.mkdir(parents=True, exist_ok=True)
    for i, members in enumerate(design):
        observations = {}
        for sample in sorted(members):
            s = samples[sample]
            for j, hap in enumerate(rng.choices(s['haplotypes'], weights=s['frequencies'], k=depth)):
                key = f'p{i}_r{len(observations)}'
                observations[key] = records[hap]
                truth[key] = sample
        write_fasta(output / f'pool_{i}.fasta', observations)
        pools.append(Pool(members, observations))
    # Ground-truth labels are withheld from deconvolution.
    recovered, info = deconvolve(pools)
    for sample, reads in recovered.items():
        write_fasta(output / f'recovered_{sample}.fasta', reads)
    assignments = [(key, sample) for sample, reads in recovered.items() for key in reads]
    correct = sum(truth[key] == sample for key, sample in assignments)
    assigned_ids = {key for key, _ in assignments}
    summary = {'synthetic': True, 'seed': seed, 'samples': len(samples), 'pools': len(design),
               'pool_memberships': [sorted(p) for p in design], 'signatures': signatures,
               'pool_to_sample_ratio': len(design) / len(samples),
               'observations': len(truth), 'assigned_observations': len(assigned_ids),
               'assignment_entries': len(assignments), 'correct_assignments': correct,
               'accuracy': correct / len(assignments) if assignments else 0,
               'recall': correct / len(truth) if truth else 0,
               'ambiguous_observations': sum(c > 1 for c in Counter(k for k, _ in assignments).values()),
               'deconvolution': info, 'design_trace': trace,
               'recovered_counts': {str(s): len(reads) for s, reads in recovered.items()}}
    (output / 'result.json').write_text(json.dumps(summary, indent=2) + '\n')
    (output / 'truth.json').write_text(json.dumps(truth, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


def main():
    parser = argparse.ArgumentParser(description='Chapter-based single-process viral pooling')
    sub = parser.add_subparsers(dest='command', required=True)
    gen = sub.add_parser('generate-data')
    gen.add_argument('--output', type=Path, default=Path('data'))
    gen.add_argument('--seed', type=int, default=37)
    run = sub.add_parser('demo')
    run.add_argument('--data', type=Path, default=Path('data'))
    run.add_argument('--output', type=Path, default=Path('results/demo'))
    run.add_argument('--solver', choices=['tabu', 'exact'], default='tabu')
    run.add_argument('--seed', type=int, default=91)
    run.add_argument('--depth', type=int, default=80)
    args = parser.parse_args()
    if args.command == 'generate-data':
        generate_data(args.output, args.seed)
    else:
        demo(args.data, args.output, args.solver, args.seed, args.depth)


if __name__ == '__main__':
    main()
