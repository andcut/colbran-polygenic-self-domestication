# Extension Results

This run extends Colbran et al. with two traits often placed near the self-domestication literature:

- `gfp`: a literature-weighted General Factor of Personality derived from Big Five GWAS
- `schizophrenia`: European schizophrenia summary statistics processed under the same downstream pipeline

## Lead-locus counts

| Trait | Candidate 1240k hits pre-clump | Clumped lead loci |
| --- | ---: | ---: |
| `gfp` | `163` | `39` |
| `schizophrenia` | `2,553` | `523` |

## Headline sign-only results

- `gfp joint_non_eur = 0.0051399486`
- `schizophrenia joint_non_eur = 9.9999e-06`
- `gfp joint_non_eur_sensitivity = 0.0057399426`
- `schizophrenia joint_non_eur_sensitivity = 9.9999e-06`

## Output

- Main summary: `polygenic_summary.tsv`
- Figure 6-style overlay can be generated with `scripts/plot_figure6_comparison.py`.
