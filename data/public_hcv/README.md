# Public HCV HVR1 source subset

Source dataset: Kamila Caraballo Cortes (2019),
[Zenodo DOI 10.5281/zenodo.3266795](https://doi.org/10.5281/zenodo.3266795).
Related study: Janiak et al., [Scientific Reports 9, 12846 (2019)](https://doi.org/10.1038/s41598-019-49258-y).

The source dataset is licensed [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
The raw files are redistributed with this attribution. They contain public specimen
identifiers and viral sequences; no attempt was made to identify individuals.
The source deposit describes error-corrected, reconstructed intrahost HVR1
populations. Its filenames label genotype 1; the associated paper describes the
study's genotype-1 cohort as subtype 1b. We did not independently re-genotype them.

## Selection and contents

Before evaluating the algorithm, select the first eight files matching
`<numeric ID> chronic genotype 1 HCV infection.fasta`, sorted by numeric ID:
53, 80, 142, 153, 157, 165, 200, 201. This rule avoids choosing specimens based on
measured performance; it is not a random or representative cohort sample.

- `raw/`: eight unmodified source files (338 haplotype records, 5–114 per specimen).
- `provenance.json`: original filenames, download URLs, source MD5 checksums and SHA-256 checksums.
- `unaligned.fasta`: original gap symbols removed, identifiers changed to unique internal IDs.
- `aligned.fasta`: all 338 sequences aligned together using MAFFT 7.526, `--auto --thread 1`;
  180 columns, 333 distinct strings. Ungapped bases are unchanged.
- `samples.json`: normalized source header frequencies, source/internal ID mapping,
  alignment version and checksum. All listed variants are retained, without an
  additional 0.5% cutoff. No measured titers are supplied or invented.

Header `freq=` values are the source's estimated percentages. They are normalized
to sum to one within each specimen for simulation. Original values remain in the
raw files and mapping. Source specimen identity is evaluation truth only; the
clustering receives anonymous observation IDs, aligned sequences, and known
pool-membership sets. Mixtures are independently resampled, not experimental pools.

## Reproduce

With the committed alignment, the comparison requires only Python:

```sh
python3 -m viral_pooling.benchmark
```

To realign the committed source files, install [MAFFT](https://mafft.cbrc.jp/alignment/software/)
(version 7.526 used for the saved results), then run:

```sh
python3 -m viral_pooling.public_data
```

To re-download the same record subset and realign (network and MAFFT required):

```sh
python3 -m viral_pooling.public_data --download
```

Downloads are checked against the record's MD5 values. Alignment preparation checks
local source SHA-256 values and verifies that alignment preserves every ungapped
sequence. Regenerating with another MAFFT version may change alignment and metrics.
Retain the committed provenance for comparison if the source record is revised.

See [the comparison report](../../results/comparison/README.md) for results and limitations.
