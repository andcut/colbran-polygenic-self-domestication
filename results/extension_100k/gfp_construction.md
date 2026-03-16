# GFP Construction

This repository uses a literature-weighted GFP construction.

Following the phenotypic GFP literature, neuroticism is reversed before combination.
Following Colbran et al.'s downstream polygenic-selection pipeline, only 1240k variants passing p < 1e-8 are retained before PLINK clumping.

- Shared AADR-intersecting Big Five rsIDs: `580,254`
- Harmonized shared Big Five rsIDs: `580,254`
- Alignment failures dropped: `0`
- Stouffer denominator using empirical Big Five GWAS correlation: `1.515412`
- Pre-clump GFP-significant variants: `163`
- Post-clump GFP lead loci: `39`

## Big Five Loadings

- `agreeableness`: `0.57`
- `conscientiousness`: `0.63`
- `extraversion`: `0.57`
- `openness`: `0.42`
- `neuroticism`: `0.62`

## Big Five support among GFP-significant variants

- `4/5` concordant traits: `30`
- `5/5` concordant traits: `133`

This is the main practical defense of the public GFP construction.
Although the standalone Big Five GWAS differ a lot in raw signal density, the final GFP-significant loci are mostly not one-trait artifacts: `133/163` are concordant across all `5/5` Big Five traits, and the remaining `30/163` are concordant across `4/5`.

## Standalone Big Five genome-wide-significant counts on the shared SNP set

- `agreeableness`: `0`
- `conscientiousness`: `0`
- `extraversion`: `7`
- `openness`: `17`
- `neuroticism`: `886`
