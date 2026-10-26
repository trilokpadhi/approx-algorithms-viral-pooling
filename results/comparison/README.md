# Public HCV versus synthetic pooling benchmark

Real source sequences; simulated pools and source-labelled read observations.

Source: [Janiak et al. (2019)](https://doi.org/10.1038/s41598-019-49258-y); [Zenodo dataset, DOI 10.5281/zenodo.3266795](https://doi.org/10.5281/zenodo.3266795). Dataset creator: Kamila Caraballo Cortes; CC BY 4.0.

## Protocol

Selected the first eight chronic genotype-1 files in numeric specimen-ID order before measuring accuracy: 53, 80, 142, 153, 157, 165, 200, 201. All source haplotypes were retained, with header frequencies normalized within specimen. No additional frequency cutoff was applied. Source gaps were stripped; all sequences were realigned together with MAFFT 7.526, --auto --thread 1.

Both datasets: eight samples, T=4, tabu GPDA, all pairs assumed compatible, 80 independently drawn observations per included specimen per pool, seeds 91–95, initial W=2, retry through W=4 only on an empty recovered sample. The complete-graph assumption controls the design comparison; measured titers and epidemiological exclusions were not supplied. Both designs have four pools and 1,040 observations per run. Ground-truth source IDs are used only for evaluation.

## Results

| Measure | Synthetic | Public HCV |
|---|---:|---:|
| Haplotype records | 32 | 338 |
| Distinct strings | 32 | 333 |
| Alignment columns | 306 | 180 |
| Mean read-assignment accuracy | 100.00% | 71.81% |
| Mean read recall | 100.00% | 71.81% |
| Mean per-specimen recall (macro) | 100.00% | 72.97% |
| Mean recall excluding singleton specimen | 100.00% | 69.46% |

| Seed | Synthetic correct / total | HCV correct / total | Final HCV W |
|---|---:|---:|---:|
| 91 | 1040/1040 | 746/1040 | 3 |
| 92 | 1040/1040 | 758/1040 | 3 |
| 93 | 1040/1040 | 747/1040 | 4 |
| 94 | 1040/1040 | 750/1040 | 4 |
| 95 | 1040/1040 | 733/1040 | 3 |

Every observation was assigned exactly once in these runs. Accuracy is correct assignment entries / all assignment entries; recall is correctly recovered observation IDs / all input IDs. They coincide here. Repeated seeds are simulation replicates on the same specimens, not independent cohorts.

All 48 public-data clustering fits across the main runs and retries reported outer convergence: True; inner convergence: True. Convergence does not imply biologically correct clustering.

### Per-specimen recall, seed 91

| Public specimen | Correct / true observations | Recall |
|---|---:|---:|
| 53 | 240/240 | 100.00% |
| 80 | 6/160 | 3.75% |
| 142 | 70/160 | 43.75% |
| 153 | 80/80 | 100.00% |
| 157 | 160/160 | 100.00% |
| 165 | 50/80 | 62.50% |
| 200 | 60/80 | 75.00% |
| 201 | 80/80 | 100.00% |

### Fixed-W sensitivity, seed 91

These runs disable retries to isolate W. The main protocol above retains retries.

| Fixed W | Correct / total | Accuracy | Empty recovered specimen IDs |
|---|---:|---:|---|
| 1 | 499/1040 | 47.98% | 80, 157 |
| 2 | 710/1040 | 68.27% | 80 |
| 3 | 746/1040 | 71.73% | none |

## Limitations and interpretation

- These are reconstructed, error-corrected HVR1 haplotypes, not raw FASTQ reads or experimentally pooled specimens. Public specimen membership is a simulation label, not independently validated laboratory truth.
- Frequency estimates come from source FASTA headers. Independent resampling does not reproduce PCR amplification, titer imbalance, contamination, dropout, read errors, or pooling-dependent bias.
- Only eight chronic genotype-1 specimens from one study are included. Numeric-ID selection is deterministic but is not random or representative. No parameters or specimens were selected to maximize the measured score.
- The public alignment has 180 columns versus 306 in the synthetic demo. Diversity, length, haplotype counts and frequency profiles all differ; the accuracy gap cannot be attributed to one factor.
- One distinct aligned string occurs in multiple public specimens. Identical strings alone cannot uniquely establish an originating specimen; the simulated observation labels still retain their original source.
- Specimen 201 occupies a singleton pool and is recovered trivially. The table therefore also reports recall excluding that specimen (synthetic sample 7 is excluded in the matching metric). The 4/8 pool ratio is a count reduction, not a demonstrated 50% monetary saving: total sequencing depth is not fixed.
- Empty-sample retries do not detect nonempty but incorrect assignments. In seed 91, specimen 80 has very poor recall despite a nonempty result. Fixed-W sensitivity also shows parameter dependence.
- The chapter interpretations and repairs documented in the main README remain in effect: queue refinement, stopping-condition repairs, bounded search, and the printed unweighted consensus update. No clustering algorithm was changed to improve these results.
- Five sampling seeds describe variability conditional on these same source sequences. They are not a confidence interval for other patients, studies or laboratory workflows.
- This exposes a limitation of the current implementation under a more difficult input distribution. It neither reproduces nor disproves the chapter’s published experimental results.

## Reproduce

From the repository root, with the committed aligned data (no network or aligner required):

```sh
python3 -m viral_pooling.benchmark
python3 -m unittest discover -s tests -v
```

For source download and alignment regeneration, see data/public_hcv/README.md. Machine-readable metrics, confusion matrices and convergence traces are in comparison.json. Timing is diagnostic only and depends on hardware and interpreter.
