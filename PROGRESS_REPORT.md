# Project Progress Report

**Project:** Viral Sample Pool Design and Deconvolution
**Students:** Trilok Padhi (Panther ID: 02716700); Mark Grinshpon (Panther ID: TBD)
**Course:** Approximation Algorithms
**Date:** October 26, 2026

## Project objective and current progress

My goal is to implement the computational methods in *Pooling Strategy for Massive
Viral Sequencing* and evaluate sample recovery using public HCV HVR1 sequences.
The problem is to design overlapping pools that identify individual specimens,
then assign the sequences observed in those pools back to their source specimens.

I have implemented the compatibility graph, Greedy Pool Design Algorithm (GPDA),
tabu-search clique-cut selection, kGEM clustering, and generalized pool
intersections and differences. I have also completed an exploratory evaluation
using eight public specimens. The implementation is a chapter-based reproduction
with documented interpretations of the printed pseudocode, rather than a claim
of exact reproduction of the published experiments.

## Methodology: public HCV dataset

### Data source and specimen selection

I used the HCV HVR1 sequence deposit associated with Janiak et al. (2019), available
on Zenodo under DOI 10.5281/zenodo.3266795 and the CC BY 4.0 license. The source
contains reconstructed, error-corrected haplotypes and estimated frequencies in
FASTA headers. These are not raw sequencing reads.

Before evaluating recovery, I selected the first eight files labeled chronic
genotype 1 in numeric specimen-ID order: **53, 80, 142, 153, 157, 165, 200, and 201**.
The rule is deterministic but does not produce a random or representative cohort.
The filenames use genotype 1; the study describes its genotype-1 cohort as subtype
1b. I did not independently re-genotype the sequences.

| Public specimen | Source haplotype records |
|---|---:|
| 53 | 16 |
| 80 | 70 |
| 142 | 35 |
| 153 | 24 |
| 157 | 20 |
| 165 | 5 |
| 200 | 54 |
| 201 | 114 |
| **Total** | **338** |

### Preparation and pooling protocol

I retained all listed haplotypes without an additional frequency cutoff and
normalized the header frequency estimates within each specimen. I removed the
source gap symbols and aligned all sequences together using MAFFT 7.526 with
`--auto --thread 1`. The alignment contains 180 columns and 333 distinct strings.
Checksums verify the downloaded files, and sequence-content checks confirm that
alignment did not change the ungapped nucleotide strings.

For this evaluation, I assumed all eight specimens were mutually compatible.
Measured titers and epidemiological exclusions were not used. Thus, the experiment
tests deconvolution under a complete compatibility graph, rather than validating
the biological compatibility model.

I used a maximum pool size of four and independently sampled 80 observations per
included specimen per pool using the normalized source frequencies. Each replicate
has 1,040 observations. I ran five replicates using seeds 91 through 95. The initial
clustering multiplier was W=2; deconvolution was repeated from the original pools
with increased W, up to W=4, if any recovered specimen was empty.

**The input haplotypes are real public HCV sequences, but pooling and read sampling
are simulated.** No laboratory pooling, new sequencing, or added sequencing-error
model was used. Ground-truth source labels were withheld from clustering and used
only to evaluate assignments. Independent observations retain distinct identifiers,
even if their nucleotide strings are identical.

## Algorithm

### Pool design with GPDA and tabu search

Let G represent sample compatibility and H represent unresolved sample pairs.
Initially, H contains every pair of specimens. A candidate pool A is valid if it
is a clique in G and its size does not exceed the limit T=4. It separates a pair
when exactly one endpoint belongs to A.

The chapter's clique-cut objective can be written as

$$
f(A)=c_H(A,V\setminus A)-(|E(H)|+1)q_G(A),
$$

where c_H counts unresolved pairs separated by A and q_G counts incompatible pairs
inside A. The penalty discourages invalid pools. The tabu heuristic searches both
orientations of a MaxCut initialization, uses the chapter's single-vertex gain
formulas to make improving moves, and uses undo operations and temporarily forbidden
moves to explore alternatives. It trims an oversized record to satisfy T. The
implementation also uses documented iteration limits and a singleton fallback.

GPDA adds the selected pool and removes the pairs it separates from H. It repeats
until H has no edges, then adds an uncovered singleton if needed. Finally, it checks
that every specimen is covered and has a distinct pool-membership signature.

For the public specimens, the resulting design was:

| Pool | Specimens included |
|---|---|
| P1 | 53, 142, 157, 200 |
| P2 | 53, 80, 157, 165 |
| P3 | 53, 80, 142, 153 |
| P4 | 201 |

This gives four pools for eight specimens. The same design is used across the five
sampling replicates. GPDA and tabu search are heuristics; I do not claim that their
output is globally optimal. The additional exact clique-cut solver is used only
for small-instance verification.

### kGEM clustering

For two pools, I set the requested cluster count to k=W times the number of
specimens in their membership union. If k exceeds the number of distinct observed
strings, the implementation caps it. The initial centroid is the most frequent
read, and subsequent seeds maximize their minimum Hamming distance to the existing
seeds.

For centroid g_i and aligned read r of length L, the emission likelihood is

$$
h_{i,r}=\prod_{m=1}^{L}f_{i,m}(r_m).
$$

The consensus allele at a position has probability 1-4 epsilon and each of the other
four alleles has probability epsilon; the alphabet includes a deletion symbol.
The current implementation uses epsilon=0.01. Given mixture proportions pi_i,
posterior cluster membership is

$$
p_{i,r}=\frac{\pi_i h_{i,r}}{\sum_j\pi_j h_{j,r}}.
$$

The inner iteration estimates mixture proportions using observed read counts. The
outer iteration updates the consensus alleles using the chapter's printed posterior
sum over distinct read strings. This is a discrete consensus update, not a full
soft allele-frequency update. Likelihoods are computed in logarithmic form for
numerical stability. The iterations stop when estimates stabilize or a documented
iteration or cycle limit is reached. Each observation is assigned to its
highest-posterior cluster.

### Generalized intersections and deconvolution

For pools P and Q, I cluster their combined observations and retain clusters that
contain reads from both pools. These clusters form the generalized intersection.
The known pool memberships determine which specimens this intersection represents;
the nucleotide identities determine which reads enter it. Directed differences
remove the selected observation IDs from the corresponding original pool.

The queue procedure repeatedly replaces a pair of overlapping pools by its
intersection and nonempty directed differences until individual specimen sets
remain. Replacing both operands is an explicit interpretation of the chapter's
queue pseudocode. Increasing W is triggered only by an empty recovered specimen;
a nonempty but incorrect result does not trigger another attempt.

## Results on the public HCV dataset

Accuracy is the number of correct assignment entries divided by all assignment
entries. Read recall is the number of correctly recovered observation IDs divided
by all input observation IDs. Every observation was assigned exactly once in these
runs, so accuracy and recall coincide.

| Measure | Public HCV result |
|---|---:|
| Specimens | 8 |
| Source haplotype records | 338 |
| Alignment columns | 180 |
| Pools | 4 |
| Observations per replicate | 1,040 |
| Mean assignment accuracy and read recall | 71.81% |
| Accuracy range across five seeds | 70.48–72.88% |
| Mean per-specimen recall (macro average) | 72.97% |
| Mean recall excluding the singleton specimen | 69.46% |

### Replicate results

| Sampling seed | Correct assignments / observations | Accuracy | Final W |
|---|---:|---:|---:|
| 91 | 746 / 1,040 | 71.73% | 3 |
| 92 | 758 / 1,040 | 72.88% | 3 |
| 93 | 747 / 1,040 | 71.83% | 4 |
| 94 | 750 / 1,040 | 72.12% | 4 |
| 95 | 733 / 1,040 | 70.48% | 3 |

The variation here comes from sampling the same underlying specimens repeatedly.
It does not measure uncertainty across independent cohorts. All 48 kGEM fits in
the main runs, including retries, reported inner and outer convergence. This shows
numerical stabilization, not correctness of the resulting specimen assignments.

### Specimen-level recovery

For seed 91, recovery differed substantially between specimens:

| Public specimen | Correct / true observations | Recall |
|---|---:|---:|
| 53 | 240 / 240 | 100.00% |
| 80 | 6 / 160 | 3.75% |
| 142 | 70 / 160 | 43.75% |
| 153 | 80 / 80 | 100.00% |
| 157 | 160 / 160 | 100.00% |
| 165 | 50 / 80 | 62.50% |
| 200 | 60 / 80 | 75.00% |
| 201 | 80 / 80 | 100.00% |

Specimen 201 occupies a singleton pool, so its complete recovery does not test
separation of mixed populations. Excluding that specimen gives a more informative
additional measure of performance on pooled specimens. Specimen 80 remains poorly
recovered despite a nonempty output, showing why the empty-sample retry rule is
not sufficient to detect failure.

### Sensitivity to W

I also used seed 91 with retries disabled to examine fixed W values:

| Fixed W | Correct / total observations | Accuracy | Empty recovered specimens |
|---|---:|---:|---|
| 1 | 499 / 1,040 | 47.98% | 80, 157 |
| 2 | 710 / 1,040 | 68.27% | 80 |
| 3 | 746 / 1,040 | 71.73% | None |

The results depend on W, but these three runs do not establish that increasing W
always improves recovery. They are a sensitivity analysis on one sampling replicate,
not a separate validation set or a parameter-optimization study.

## Verification and limitations

The saved full test run contains 19 passing tests. Relevant checks include pool
feasibility and separation, clique-cut gain formulas, source-file checksums,
frequency parsing, preservation of ungapped sequence content, and metric accounting
under deliberately incorrect assignments.

The main limitations of the present evaluation are:

- The dataset contains reconstructed haplotypes and frequency estimates, not raw reads or experimentally pooled specimens. The simulation does not include PCR bias, contamination, unequal titers, or sequencing errors.
- Only eight chronic genotype-1 specimens from one study were evaluated. All samples were assumed compatible, so the experiment does not test measured biological compatibility restrictions.
- One aligned sequence occurs in more than one specimen. Identical strings alone cannot uniquely establish the original source of an observation.
- Numerical convergence and nonempty recovered sets do not guarantee correct biological assignments. Poor recovery can remain undetected by the retry condition.
- The implementation includes documented repairs to stopping conditions, bounded search, and queue refinement. It is not a literal execution of every printed pseudocode line or a reproduction of the chapter's published empirical accuracy.
- Four pools instead of eight individually processed specimens demonstrates a reduction in pool count in this setup. It does not establish a 50% monetary saving, because total sequencing depth and laboratory costs were not held fixed.

## Next steps

My next priority is to trace the generalized intersections that misassign
observations from specimens 80 and 142. I plan to inspect the saved confusion
matrices and clustering diagnostics, then vary W, sampling depth, and pool-size
limits separately. I also plan to evaluate additional public specimens using
selection rules set before inspecting their results. Any alternative consensus
update or queue strategy will be identified as a departure from the current
chapter-based implementation. These are planned investigations, not completed
improvements.

## References and supporting results

1. Skums et al. (2016). *Pooling Strategy for Massive Viral Sequencing*. Chapter 3 in *Computational Methods for Next Generation Sequencing Data Analysis*, pp. 57–83.
2. Janiak et al. (2019). *Hepatitis C virus (HCV) genotype 1b displays higher genetic variability of hypervariable region 1 (HVR1) than genotype 3*. [Scientific Reports](https://doi.org/10.1038/s41598-019-49258-y).
3. Caraballo Cortes, K. (2019). Public HCV HVR1 sequence dataset. [Zenodo, DOI 10.5281/zenodo.3266795](https://doi.org/10.5281/zenodo.3266795), CC BY 4.0.

- [Public-data provenance and preparation](data/public_hcv/README.md)
- [Public HCV results and diagnostics](results/public_hcv/report-results.json)
- [Saved full test output](results/verification.txt)
