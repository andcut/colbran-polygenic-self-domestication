# Deep Schizophrenia P-Value

The standard result tables use `100,000` sign permutations.
At that depth, schizophrenia reaches the empirical floor in the `joint_non_eur` test.

For the Figure 6-style plot, schizophrenia was rerun with `10,000,000` permutations from the same `counts_joint_non_eur.tsv.gz` and `admixture_joint_non_eur.tsv` files produced by `scripts/run_extension.py`.
That rerun gives:

- `schizophrenia joint_non_eur = 5.9999994e-07`

To regenerate the table:

```bash
python3 scripts/run_deep_schizophrenia.py
```
