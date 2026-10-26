"""Aligned FASTA input, chapter kGEM, GI and queue-based deconvolution."""
from dataclasses import dataclass
from collections import Counter
import math

ALPHABET = 'ACTG-'


def read_fasta(path):
    records, name, parts = {}, None, []
    def flush():
        if name is not None:
            seq = ''.join(parts).upper()
            if name in records or not seq or set(seq) - set(ALPHABET):
                raise ValueError('duplicate ID, empty sequence, or unsupported allele')
            records[name] = seq
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                flush()
                name, parts = line[1:].split()[0], []
            elif name is None:
                raise ValueError('sequence before FASTA header')
            else:
                parts.append(line)
    flush()
    if not records or len({len(s) for s in records.values()}) != 1:
        raise ValueError('nonempty, equal-length aligned reads required')
    return records


def write_fasta(path, records):
    with open(path, 'w') as f:
        for name, seq in records.items():
            f.write(f'>{name}\n{seq}\n')


def hamming(a, b):
    if len(a) != len(b):
        raise ValueError('aligned lengths differ')
    return sum(x != y for x, y in zip(a, b))


def kgem(reads, k, epsilon=0.01, tolerance=1e-9, max_outer=100, max_inner=500):
    """Chapter kGeEM, including its discrete consensus allele update.

    Reads are observations (duplicates represent observed counts). E steps are
    computed for unique strings and weighted by multiplicity for mixture fitting.
    Consensus update uses p_i,r without count weights, as printed on p.72.
    """
    if not reads or k < 1 or not 0 < epsilon < 0.2 or tolerance <= 0 or min(max_outer, max_inner) < 1:
        raise ValueError('invalid clustering input')
    if len({len(r) for r in reads}) != 1 or any(not r or set(r) - set(ALPHABET) for r in reads):
        raise ValueError('aligned ACTG- reads required')
    counts = Counter(reads)
    unique = sorted(counts)
    k = min(k, len(unique))
    centers = [min(unique, key=lambda r: (-counts[r], r))]
    while len(centers) < k:
        centers.append(max((r for r in unique if r not in centers),
                           key=lambda r: (min(hamming(r, c) for c in centers), counts[r], r)))
    log_match, log_other = math.log(1 - 4 * epsilon), math.log(epsilon)
    seen, history = set(), []
    converged = False
    inner_converged = True
    def expectation(centroids):
        nonlocal inner_converged
        emissions = [[sum(log_match if a == b else log_other for a, b in zip(r, c))
                      for c in centroids] for r in unique]
        weights = [1 / k] * k
        def posterior():
            result, ll = [], 0.0
            for r, emission in zip(unique, emissions):
                logs = [math.log(max(w, 1e-300)) + h for w, h in zip(weights, emission)]
                peak = max(logs)
                norm = sum(math.exp(x - peak) for x in logs)
                result.append([math.exp(x - peak) / norm for x in logs])
                ll += counts[r] * (peak + math.log(norm))
            return result, ll
        for _ in range(max_inner):
            p, _ = posterior()
            updated = [sum(counts[r] * row[i] for r, row in zip(unique, p)) / len(reads) for i in range(k)]
            change = sum((x - y) ** 2 for x, y in zip(weights, updated))
            weights = updated
            if change < tolerance:
                break
        else:
            inner_converged = False
        p, ll = posterior()
        return p, weights, ll
    for _ in range(max_outer):
        p, weights, ll = expectation(centers)
        history.append(ll)
        updated = [''.join(max(ALPHABET, key=lambda allele: (
            sum(row[i] for r, row in zip(unique, p) if r[pos] == allele),
            -ALPHABET.index(allele))) for pos in range(len(unique[0]))) for i in range(k)]
        if updated == centers:
            converged = True
            break
        state = tuple(centers)
        seen.add(state)
        centers = updated
        if tuple(centers) in seen:
            break
    # Recompute after final centroid change; don't return stale posteriors.
    p, weights, ll = expectation(centers)
    lookup = {r: max(range(k), key=lambda i: row[i]) for r, row in zip(unique, p)}
    return {'labels': [lookup[r] for r in reads], 'centers': centers, 'weights': weights,
            'log_likelihood': ll, 'history': history, 'converged': converged,
            'inner_converged': inner_converged, 'effective_k': k}


@dataclass
class Pool:
    samples: frozenset
    reads: dict  # observation ID -> aligned sequence; identity survives GI/difference


def generalized_intersection(p, q, w=2, diagnostics=None):
    if w < 1:
        raise ValueError('W >= 1 required')
    shared = p.samples & q.samples
    if not shared:
        return Pool(frozenset(), {})
    merged = dict(p.reads)
    for key, seq in q.reads.items():
        if key in merged and merged[key] != seq:
            raise ValueError('observation ID has conflicting sequences')
        merged[key] = seq
    if not merged:
        return Pool(shared, {})
    keys = sorted(merged)
    fit = kgem([merged[key] for key in keys], w * len(p.samples | q.samples))
    if diagnostics is not None:
        diagnostics.append({'w': w, 'requested_k': w * len(p.samples | q.samples),
                            **{key: fit[key] for key in ('effective_k', 'converged', 'inner_converged')}})
    left = {label for key, label in zip(keys, fit['labels']) if key in p.reads}
    right = {label for key, label in zip(keys, fit['labels']) if key in q.reads}
    keep = left & right
    return Pool(shared, {key: merged[key] for key, label in zip(keys, fit['labels']) if label in keep})


def deconvolve(pools, w=2, max_w=4):
    """PD/IS via queue refinement. Replace both operands to prevent duplicates.

    On empty recovered samples repeat from original data with increased W.
    Return all sample IDs, including empty failures, and retry diagnostics.
    """
    if w < 1 or max_w < w:
        raise ValueError('invalid W range')
    expected = set().union(*(p.samples for p in pools)) if pools else set()
    clustering_trace = []
    for current_w in range(w, max_w + 1):
        queue = [Pool(p.samples, dict(p.reads)) for p in pools]
        result, operations = {}, 0
        while queue:
            p = queue.pop(0)
            partner = next((i for i, q in enumerate(queue) if p.samples & q.samples), None)
            if partner is None:
                if len(p.samples) != 1:
                    raise ValueError('pool memberships do not separate samples')
                sample = next(iter(p.samples))
                result.setdefault(sample, {}).update(p.reads)
                continue
            q = queue.pop(partner)
            inter = generalized_intersection(p, q, current_w, clustering_trace)
            queue.append(inter)
            if p.samples - q.samples:
                queue.append(Pool(p.samples - q.samples, {k: v for k, v in p.reads.items() if k not in inter.reads}))
            if q.samples - p.samples:
                queue.append(Pool(q.samples - p.samples, {k: v for k, v in q.reads.items() if k not in inter.reads}))
            operations += 1
        missing = sorted(s for s in expected if not result.get(s))
        if not missing:
            break
    return result, {'w': current_w if expected else w, 'empty_samples': missing if expected else [], 'operations': operations if expected else 0, 'clustering_trace': clustering_trace}
