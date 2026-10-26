# Viral sequencing pooling — chapter implementation

A runnable, single-process Python implementation of the computational methods in
Skums et al., **“Pooling Strategy for Massive Viral Sequencing,” Chapter 3** of
*Computational Methods for Next Generation Sequencing Data Analysis* (2016),
printed pages 57–83.

The chapter is the primary algorithm specification. See the source mapping and
explicit reproduction notes below.

## Run

Requires Python 3.10 or later; no runtime packages or network access required.
From this folder:

```sh
python3 -m unittest discover -s tests -v
python3 -m viral_pooling demo
python3 -m viral_pooling demo --solver exact --output results/exact
```

The bundled FASTA is already present. To regenerate it deterministically:

```sh
python3 -m viral_pooling generate-data --output data --seed 37
```

Optional installation: `python3 -m pip install -e .`; then use `viral-pooling demo`.
Run the module commands from the project folder, so the default data path resolves.
`--data`, `--output`, `--seed`, and `--depth` configure the demo. Output folders
should be dedicated to each run; matching output files are overwritten.

`results/demo/` contains the actual saved run: input pool FASTAs, recovered sample
FASTAs, ground truth for evaluation only, and `result.json` with memberships,
nonempty sample signatures, greedy progress, tabu diagnostics, and accuracy/recall.
Recovered observation frequencies can be obtained by counting identical strings in
a sample FASTA; this implementation does not estimate original specimen titers or
correct pooled counts for experimental bias.

## Implemented methods and source mapping

| Component | Implementation | Chapter reference |
|---|---|---|
| Compatibility graph | Titer-difference threshold; optional forbidden epidemiological pairs | §3.2.1 p.62; titer model §3.4.1 p.73 |
| GPDA | H starts complete, clique cuts remove separated pairs; enforce size cap | §3.2.2 p.64 |
| Coverage | Add uncovered singleton after separation | p.63 discussion of condition (1) |
| OCBGRP tabu search | Both MaxCut orientations, penalty objective, flip gains, tabu countdown, undo stack, repeated-state detection, final size trimming, tenure retries | §3.2.3 pp.65–68; Eqs.3.5–3.9 |
| Exact clique-cut oracle | Enumerate all feasible cliques up to T | Added verification tool, not a chapter algorithm |
| GI | k=W times union membership size; keep clusters with observations from both pools | §3.3.1 p.70 |
| PD/IS | Queue refinement using generalized intersection and two directed differences | §3.3.1 p.69; executable clarification below |
| kGEM | Most-frequent/farthest-first seeds, five alleles, likelihood products, nested mixture EM, discrete smoothed consensus update, posterior assignment | §3.3.2 pp.70–72, Eq.3.10 |

The exact solver optimizes **each clique cut**, not the global minimum number of
pools. GPDA remains greedy with either solver. No global approximation ratio or
pool-count optimality is claimed. MaxCut initialization uses a deterministic greedy
conditional-expectation 1/2 approximation followed by improving flips. All ties are
deterministic. Compatibility titers must be supplied in consistent units; whether
to use raw or logarithmic titers is a caller choice, not an inferred clinical rule.

## Explicit departures and executable interpretations

This is a chapter-based reproduction with the following disclosed repairs and
bounds, **not a literal execution of all printed pseudocode or a reproduction of
the published experiments**.

1. **Tabu stopping.** The printed OCBGRP has branches for positive and zero best
   gain, but none for a negative best gain. We treat a nonpositive best gain as
   the stopping/undo branch. Tabu-restricted stops need not be true local optima,
   so only actual cliques may update the record. The initial record is empty.
   Both orientations are searched and their final feasible cuts are compared.
2. **Search limits.** Each orientation has a default 10,000-step cap. Zero-cut
   results retry with tenure 1 through 4, following the chapter's increasing-tenure
   idea. A positive singleton cut is the final guaranteed-progress fallback.
   `design_trace` reports caps and fallbacks. Final truncation removes the vertex
   with smallest cut loss exactly as line 21 prescribes. The demo hits no cap or
   fallback. This Python reference implementation prioritizes readability over
   massive-scale performance; it recomputes gains rather than maintaining arrays.
3. **PD queue ambiguity.** Printed PD removes P′ but retains P′′ while adding their
   intersection and differences, then returns the “last n” singleton entries.
   This can repeatedly reintroduce equivalent memberships. We replace **both**
   operands by the intersection and nonempty differences, and accumulate singleton
   results by sample identity. Repeated memberships therefore merge and overlapping
   memberships split. Original observation IDs survive these operations. Differences
   remove observation IDs selected by GI, not all reads with identical strings.
4. **kGEM convergence typos.** Printed lines 8 and 11 say to stop when frequency
   change is ≥δ and when a genotype changes. We instead stop inner EM when squared
   frequency change is <δ and outer iteration when all centroids are unchanged.
   Both loops are bounded; outer cycles also stop. Final posteriors are recomputed.
   Convergence flags are returned by `kgem`; GI currently uses the bounded result.
5. **kGEM model fidelity.** We use the chapter's discrete 1−4ε/ε consensus update,
   not a replacement soft allele-frequency M-step. Default ε=0.01 and δ=1e−9 are
   implementation settings (the chapter supplies no numerical defaults here).
   The printed consensus update sums posterior p over distinct read strings,
   without observed multiplicity; we follow it, while weighting mixture EM by
   observed counts. This update is not guaranteed to increase count-weighted
   likelihood; `history` exposes values. Duplicate centroids are retained. If W·K
   exceeds distinct observed strings, k is capped; it is impossible to seed more
   distinct observed genotypes. Likelihood calculations use log-sum-exp to avoid
   underflow, preserving Eq.3.10 mathematically.
6. **Input and lab boundaries.** Reads must already be aligned, equal length, and
   use A/C/T/G/`-` (the last represents the chapter's deletion allele d). There is
   no alignment, barcode decoding, error correction, PCR, assembly, or wet-lab
   sequencing. FASTA pooling simulates independently sampled mixtures. API methods
   accept arbitrary aligned reads and explicit memberships; the CLI is a demo,
   not a full real-data pipeline. Known epidemiological exclusions must be given
   externally. GI retries from W=2 to 4 if any sample is empty and reports failures;
   nonempty results can still be wrong and do not trigger a retry. Increasing W
   can also fragment a shared sample: for the four-read toy example below, W=2
   seeds each distinct string separately and GI is empty; W=1 groups the close
   shared variants correctly. The retry strategy is therefore not a recovery guarantee.
7. **Data and experiments.** The bundled data is clearly labeled synthetic, with
   artificial HVR1-length DNA and deliberately separated host populations. It is
   not the chapter's HCV cohort and does not reproduce its empirical claims,
   simulation breadth, real read error model, or reported classification accuracy.
   See `data/README.md` for the exact generator and missing biological effects.

## Actual verification

The saved default run (Python 3, local execution) produced 4 pools for 8 synthetic
samples, pool/sample ratio 0.5, and 1,040 reads. It recovered 1,040/1,040 correctly
with no ambiguous or unassigned reads, W=2, and no empty samples. These numbers
measure this easy seeded synthetic demonstration only.

Actual result: **14 tests passed** (exit code 0). Full output and Python version
are saved in `results/verification.txt`.

The test suite checks:

- Every possible compatibility graph on four vertices, each T=1..4, with both
  exact-cut and tabu GPDA: clique legality, capacity, coverage, unique signatures,
  termination, and strictly positive greedy progress.
- Gain equations against direct objective differences for every subset of every
  four-vertex compatibility graph, plus the penalty/clique equivalence.
- Known exact cut optimum, random tabu cuts against the oracle upper bound,
  MaxCut's half-edge guarantee, and forced bounded-search fallback.
- Known sequence mixtures and expected weights, long-read numerical stability,
  generalized intersection of nonidentical sampled variants, chapter Example 1,
  malformed FASTA, too-large k, empty-sample retries, and nonseparating designs.
- Full deterministic data generation, pooling, deconvolution, and expected perfect
  recovery with both solvers. This is a regression on one synthetic dataset,
  not a general accuracy guarantee.

## Code

- `viral_pooling/design.py`: graph construction, objectives, solvers, validation.
- `viral_pooling/sequence.py`: FASTA, kGEM, GI, PD.
- `viral_pooling/__main__.py`: synthetic generator and runnable experiment.
- `tests/`: standard-library unittest suite.

For custom data, construct `Pool(frozenset(sample_ids), {observation_id: sequence})`
objects and call `deconvolve(pools)`. Observation IDs must be globally unique for
independent reads; retain the same ID only for the same observation propagated
through intermediate pools. Design memberships must cover and separate samples.

## Small worked example

Take three samples S1, S2, S3 with titers 10, 11, 12, threshold 2,
and maximum pool size T=2. All pairs are compatible, so G is a triangle.
H initially contains all three unresolved pairs.

The current deterministic tabu implementation returns:

1. P1={S1,S3}, separating pairs (S1,S2) and (S2,S3). Only (S1,S3) remains in H.
2. P2={S1,S2}, separating the remaining pair. H is now empty.

| Sample | In P1? | In P2? | Signature |
|---|---|---|---|
| S1 | yes | yes | 11 |
| S2 | no | yes | 01 |
| S3 | yes | no | 10 |

Each sample has a different nonzero signature, so the two pools identify all three.
The tabu search's scores are 2, then 1. This example requires no fallback.

For sequence deconvolution, suppose the observed reads are:

- P1: `CCCCCC` (S1) and `AAAAAA` (S3).
- P2: `CCCCCT` (S1) and `GGGGGG` (S2).

The parenthetical sample labels are ground truth for this explanation; the
clustering does not receive them. Literal sequence intersection is empty.
Using W=1 gives k=3. kGEM groups `CCCCCC` and `CCCCCT` together, with the other
two strings in separate clusters. The close cluster contains reads from both
pools, so GI identifies it as S1. The directed differences recover S3 from P1
and S2 from P2. This illustrates the chapter's generalized intersection idea.
W=1 is an explicit toy setting, not the demo default of W=2.

With epsilon=0.01 each centroid has probability 0.96 for its consensus allele
and 0.01 for each of the other four alleles. A six-base read matching a centroid
has emission likelihood 0.96^6; a one-mismatch read has 0.96^5*0.01. Posterior
membership also multiplies these values by the current mixture weights and
normalizes across centroids. kGEM alternates mixture fitting and consensus
updates, then assigns each observation to its highest-posterior cluster.
