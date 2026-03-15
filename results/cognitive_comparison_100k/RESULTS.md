# Cognitive Comparison Results

This run applies the same public Colbran-style workflow used elsewhere in this repository to two cognitive comparison traits:

- `iq`: Savage et al. 2018 intelligence GWAS (full summary statistics)
- `ea4`: Lee et al. 2018 educational attainment GWAS lead-hit supplementary table

## Lead-locus counts

| Trait | Candidate 1240k hits pre-clump | Clumped lead loci |
| --- | ---: | ---: |
| `iq` | `1,401` | `304` |
| `ea4` | `200` | `200` |

## Headline sign-only results

- `iq joint_non_eur = 0.032399676`
- `ea4 joint_non_eur = 0.010149899`
- `iq joint_non_eur_sensitivity = 0.0075799242`
- `ea4 joint_non_eur_sensitivity = 0.088069119`

## Weighted comparison

- `iq weighted joint_non_eur = 0.017859821`
- `ea4 weighted joint_non_eur = 0.012089879`

## Output

- Main sign summary: `polygenic_summary.tsv`
- Main weighted summary: `weighted_summary.tsv`
- Figure 6-style overlay can be generated with `scripts/plot_figure6_comparison.py`.
