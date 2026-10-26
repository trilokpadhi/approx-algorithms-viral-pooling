"""Chapter 3 §§3.2.1–3.2.3: compatibility, GPDA, and OCBGRP."""
from itertools import combinations


def compatibility(titers, max_difference, forbidden=()):
    """Edges iff titer difference <= threshold and pair is not forbidden.

    Titer units (raw or log) must be chosen consistently by the caller.
    Forbidden pairs encode externally known epidemiological incompatibility.
    """
    if max_difference < 0 or any(not __import__('math').isfinite(t) for t in titers):
        raise ValueError('finite titers and nonnegative threshold required')
    banned = {frozenset(p) for p in forbidden}
    return {frozenset((i, j)) for i, j in combinations(range(len(titers)), 2)
            if abs(titers[i] - titers[j]) <= max_difference
            and frozenset((i, j)) not in banned}


def cut(a, h):
    return sum(len(e & a) == 1 for e in h)


def is_clique(a, g):
    return all(frozenset(p) in g for p in combinations(a, 2))


def objective(a, g, h):
    """Eq. 3.5 in set form; M=|E(H)|+1."""
    return cut(a, h) - (len(h) + 1) * sum(
        frozenset(p) not in g for p in combinations(a, 2))


def gain(v, a, n, g, h):
    """Eqs. 3.8/3.9, without numerical relaxation."""
    inside = sum(frozenset((v, u)) in h for u in a if u != v)
    outside = sum(frozenset((v, u)) in h for u in range(n) if u not in a and u != v)
    missing = sum(frozenset((v, u)) not in g for u in a if u != v)
    return inside - outside + (len(h) + 1) * missing if v in a else outside - inside - (len(h) + 1) * missing


def exact_clique_cut(n, g, h, limit):
    """Exponential reference oracle; not the chapter's heuristic."""
    best = frozenset()
    for size in range(1, min(limit, n) + 1):
        for vertices in combinations(range(n), size):
            a = frozenset(vertices)
            if is_clique(a, g) and cut(a, h) > cut(best, h):
                best = a
    return best


def maxcut_half(n, h):
    """Greedy conditional-expectation 1/2 approximation, then flip improvement."""
    a, b = set(), set()
    for v in range(n):
        to_a = sum(frozenset((v, u)) in h for u in a)
        to_b = sum(frozenset((v, u)) in h for u in b)
        (a if to_b >= to_a else b).add(v)
    while True:
        moves = [(cut(a ^ {v}, h) - cut(a, h), -v, v) for v in range(n)]
        delta, _, v = max(moves, default=(0, 0, 0))
        if delta <= 0:
            return frozenset(a)
        a ^= {v}


def tabu_clique_cut(n, g, h, limit, tenure=1, max_steps=10000, max_retries=4):
    """OCBGRP: both MaxCut orientations, improving flips, stack undo, tabu.

    Explicit executable interpretations and bounded-search safeguards are in README.
    Returns (clique, diagnostics).
    """
    if limit < 1 or tenure < 1 or max_steps < 1 or max_retries < 1:
        raise ValueError('positive search parameters required')
    if n == 0:
        return frozenset(), {'steps': 0, 'fallback': False, 'bounded': False}
    x = maxcut_half(n, h)
    best = frozenset()
    diagnostics = {'steps': 0, 'fallback': False, 'bounded': False, 'tenure': tenure}
    for kt in range(tenure, tenure + max_retries):
        for initial in (x, frozenset(range(n)) - x):
            a = initial
            tabu, stack, states = [0] * n, [], set()
            record = frozenset()
            for _ in range(max_steps):
                diagnostics['steps'] += 1
                # Blocked gains are zero in the printed algorithm. Clamp negative
                # best gains too, so a strict local maximum enters the undo branch.
                options = [(gain(v, a, n, g, h), -v, v) for v in range(n) if tabu[v] == 0]
                delta, _, v = max(options, default=(0, 0, 0))
                tabu = [max(0, t - 1) for t in tabu]
                if delta > 0:
                    a = a ^ {v}
                    stack.append(v)
                    continue
                # A tabu-restricted stop need not be a clique. Record only feasible
                # states; this corrects the chapter's local-optimum assumption.
                if is_clique(a, g) and cut(a, h) > cut(record, h):
                    record = a
                if stack:
                    v = stack.pop()
                    a = a ^ {v}
                    tabu[v] = kt
                state = (a, tuple(tabu))
                if state in states:
                    break
                states.add(state)
            else:
                diagnostics['bounded'] = True
                if is_clique(a, g) and cut(a, h) > cut(record, h):
                    record = a
            while len(record) > limit:
                # Eq. at line 21: remove vertex with smallest cut loss.
                v = min(record, key=lambda u: (cut(record, h) - cut(record - {u}, h), u))
                record = record - {v}
            if cut(record, h) > cut(best, h):
                best = record
        diagnostics['tenure'] = kt
        if cut(best, h) > 0 or not h:
            return best, diagnostics
    # Bounded retry addition: every endpoint of H has a positive singleton cut.
    best = max((frozenset((v,)) for v in range(n)), key=lambda a: cut(a, h))
    diagnostics['fallback'] = True
    return best, diagnostics


def validate_design(n, g, pools, limit):
    if any(not p or len(p) > limit or not set(p) <= set(range(n)) or not is_clique(p, g) for p in pools):
        raise ValueError('infeasible pool')
    signatures = [tuple(i for i, p in enumerate(pools) if v in p) for v in range(n)]
    if any(not s for s in signatures) or len(set(signatures)) != n:
        raise ValueError('samples must have distinct nonempty signatures')
    return signatures


def gpda(n, g, limit, solver='tabu', **kwargs):
    if n < 0 or limit < 1 or solver not in ('tabu', 'exact'):
        raise ValueError('invalid design parameters')
    h = {frozenset(e) for e in combinations(range(n), 2)}
    pools, trace = [], []
    while h:
        if solver == 'exact':
            a, info = exact_clique_cut(n, g, h, limit), {}
        else:
            a, info = tabu_clique_cut(n, g, h, limit, **kwargs)
        score = cut(a, h)
        if score <= 0:
            raise RuntimeError('clique-cut solver made no progress')
        pools.append(a)
        trace.append(dict(info, separated=score, remaining=len(h) - score))
        h = {e for e in h if len(e & a) != 1}
    covered = set().union(*pools) if pools else set()
    # Chapter p.63: restore formal condition (1) with at most one singleton.
    pools.extend(frozenset((v,)) for v in range(n) if v not in covered)
    validate_design(n, g, pools, limit)
    return pools, trace
