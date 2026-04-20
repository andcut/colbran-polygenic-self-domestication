# Cognitive Comparison Results

This run applies the same Colbran-style workflow used elsewhere in this repository to two cognitive comparison traits:

- `iq`: Savage et al. 2018 intelligence GWAS (full summary statistics)
- `ea4`: Lee et al. 2018 educational attainment GWAS lead-hit supplementary table

## Lead-locus counts

| Trait | Candidate 1240k hits pre-clump | Clumped lead loci |
| --- | ---: | ---: |
| `iq` | `1,401` | `295` |
| `ea4` | `200` | `200` |

## Headline sign-only results

- `iq joint_non_eur = 0.025579744`
- `ea4 joint_non_eur = 0.0058699413`
- `iq joint_non_eur_sensitivity = 0.01704983`
- `ea4 joint_non_eur_sensitivity = 0.0069199308`

## Output

- Main summary: `polygenic_summary.tsv`
- Figure 6-style overlay can be generated with `scripts/plot_figure6_comparison.py`.
