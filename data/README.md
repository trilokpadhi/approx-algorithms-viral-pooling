# Synthetic viral amplicon demo

`synthetic_haplotypes.fasta` contains 32 **artificial DNA strings**, representing
8 artificial hosts with 4 variants each. These are NOT patient sequences, public
HCV isolates, a functional viral construct, or sequences validated as biological.
The 306-base length follows the chapter's HVR1 amplicon experiment scale. A seeded
random ancestor is mutated at 45 positions per host; each host has its founder and
variants with 2, 4, and 6 further substitutions. This yields deliberately easy,
well-separated artificial quasispecies, with frequencies 65%, 20%, 10%, and 5%.

`samples.json` defines haplotype membership, frequencies, titers in arbitrary units,
and the generation seed (37). Regenerate with:

    python3 -m viral_pooling generate-data --output data --seed 37

The demo independently samples 80 reads per included host per pool using seed 91.
This allows sampling differences across pools. There are no PCR effects, indels,
sequencing errors, recombination, contamination, or low-titer dropout. Identical
strings from independent observations retain distinct read identifiers. Read IDs
contain only pool and observation indices; the clustering receives no host truth.
