# Data Recipe

This repository commits only the small metadata files needed to define the comparison panel:

- [`colbran_eur_polygenic_table3.tsv`](colbran_eur_polygenic_table3.tsv)
- [`sample_info`](sample_info)

The repository also commits the small result artifacts in `results/`, including the exact pre-clump and clumped lead-locus tables for `gfp` and `schizophrenia`.
Those committed result files are meant to make inspection and write-up easier.
They are not substitutes for the raw AADR and schizophrenia inputs below.

The large raw inputs are expected under `data/raw/`.
This note is written as a recipe for a human or agent to reconstruct that tree without needing to infer filenames from the code.

## Expected on-disk layout

```text
data/
  raw/
    aadr/
      10537414
      10537415
      10537126
    full_sumstats/
      PGC3_SCZ_wave3_european.tsv.gz
    cache/
      bigfive_aadr/
```

If you place files elsewhere, the runner can still be used with explicit command-line overrides.
The defaults assume the layout above.

## 1. AADR `1240k` genotype files

Needed files:

- `data/raw/aadr/10537414`
- `data/raw/aadr/10537415`
- `data/raw/aadr/10537126`

These correspond to the `.ind`, `.snp`, and `.geno` EIGENSTRAT files used throughout the project.

Upstream source:

- David Reich Lab, "Allen Ancient DNA Resource (AADR): Downloadable genotypes of present-day and ancient DNA data"
- Search for the public downloadable-genotypes page if the exact URL changes on the Reich Lab site.

What to download:

- the public `1240k` AADR release used by the Colbran-style workflow
- the three EIGENSTRAT files for the same release

Important:

- the runner expects the files exactly as plain files named `10537414`, `10537415`, and `10537126`
- if the upstream download arrives with suffixes or versioned filenames, rename or symlink them into those names

Example shell recipe after manual download:

```bash
mkdir -p data/raw/aadr

mv ~/Downloads/10537414 data/raw/aadr/10537414
mv ~/Downloads/10537415 data/raw/aadr/10537415
mv ~/Downloads/10537126 data/raw/aadr/10537126
```

If your downloads have different names:

```bash
mkdir -p data/raw/aadr

ln -sf /absolute/path/to/AADR_public_1240k.ind  data/raw/aadr/10537414
ln -sf /absolute/path/to/AADR_public_1240k.snp  data/raw/aadr/10537415
ln -sf /absolute/path/to/AADR_public_1240k.geno data/raw/aadr/10537126
```

## 2. Schizophrenia summary statistics

Needed file:

- `data/raw/full_sumstats/PGC3_SCZ_wave3_european.tsv.gz`

Upstream source:

- Figshare dataset: [`scz2022`](https://figshare.com/articles/dataset/scz2022/19426775)
- This is the public Psychiatric Genomics Consortium schizophrenia summary-stat release associated with Trubetskoy et al. 2022 / PGC3

What to download:

- the European schizophrenia summary-stat file from the `scz2022` dataset
- the exact upstream filename may vary by dataset version

What the runner expects:

- a gzipped tab-separated file with columns including `ID`, `CHROM`, `POS`, `A1`, `A2`, `BETA`, and `PVAL`

Recommended procedure:

1. Open the Figshare landing page.
2. Download the European-only summary-statistics file from the dataset.
3. Inspect the header once.
4. Rename or symlink it to `data/raw/full_sumstats/PGC3_SCZ_wave3_european.tsv.gz`.

Example shell recipe:

```bash
mkdir -p data/raw/full_sumstats

mv ~/Downloads/<downloaded_scz_file>.tsv.gz \
  data/raw/full_sumstats/PGC3_SCZ_wave3_european.tsv.gz
```

Or with a symlink:

```bash
mkdir -p data/raw/full_sumstats

ln -sf /absolute/path/to/<downloaded_scz_file>.tsv.gz \
  data/raw/full_sumstats/PGC3_SCZ_wave3_european.tsv.gz
```

Quick sanity check:

```bash
python3 - <<'PY'
import gzip
from pathlib import Path
path = Path("data/raw/full_sumstats/PGC3_SCZ_wave3_european.tsv.gz")
with gzip.open(path, "rt") as handle:
    for line in handle:
        if line.startswith("##") or not line.strip():
            continue
        print(line.rstrip())
        break
PY
```

You should see a header row containing the fields listed above.

## 3. Big Five GWAS

The runner can fetch these automatically on first run, so manual download is optional.

Upstream source page:

- Yale Levey Lab data page: [https://medicine.yale.edu/lab/leveylab/data/](https://medicine.yale.edu/lab/leveylab/data/)

Direct URLs used by the runner:

- Neuroticism: `https://hugesumstats.yale.edu/dl/Neuro_MVP_UKB_sumstat_file`
- Agreeableness: `https://hugesumstats.yale.edu/dl/agree_MVP_GPC1_sumstat_file`
- Conscientiousness: `https://hugesumstats.yale.edu/dl/consc_MVP_GPC1_sumstat_file`
- Openness: `https://hugesumstats.yale.edu/dl/open_MVP_GPC1_sumstat_file`
- Extraversion: `https://hugesumstats.yale.edu/dl/extra_MVP_GPC1_sumstat_file`

What happens on first run:

- each Big Five file is streamed
- only rsIDs present on the AADR SNP panel are kept
- a small gzipped cache is written under `data/raw/cache/bigfive_aadr/`

Expected cache files after the first successful run:

- `data/raw/cache/bigfive_aadr/agreeableness_aadr_intersection.tsv.gz`
- `data/raw/cache/bigfive_aadr/conscientiousness_aadr_intersection.tsv.gz`
- `data/raw/cache/bigfive_aadr/extraversion_aadr_intersection.tsv.gz`
- `data/raw/cache/bigfive_aadr/openness_aadr_intersection.tsv.gz`
- `data/raw/cache/bigfive_aadr/neuroticism_aadr_intersection.tsv.gz`

If you want to prepopulate them manually, you can run the same URLs yourself and then let the script reuse the cached files.
Most users do not need to do that.

## 4. PLINK

The main runner uses PLINK clumping because Colbran et al. explicitly describe PLINK-based clumping in the polygenic methods.

Requirement:

- `plink` must be on `PATH`, or you can pass `--plink-path /absolute/path/to/plink`

Quick check:

```bash
plink --help >/dev/null
```

## 5. Minimal end-to-end recipe

Once the AADR files and schizophrenia file are in place:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

python3 scripts/run_extension.py
python3 scripts/plot_figure6_comparison.py
```

Main outputs:

- `results/extension_100k/polygenic_summary.tsv`
- `results/extension_100k/gfp_construction.md`
- `results/extension_100k/RESULTS.md`
- `results/figure6_comparison_100k/figure6_comparison.png`
