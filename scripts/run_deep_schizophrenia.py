#!/usr/bin/env python3
"""Run the deeper schizophrenia permutation used by the Figure 6-style plot."""

from __future__ import annotations

import argparse
import csv
import gzip
import subprocess
import sys
from pathlib import Path
from typing import TextIO


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT_DIR = ROOT / "results" / "extension_100k"
DEFAULT_OUTPUT_DIR = ROOT / "results" / "deep_schizophrenia_10m"
DEFAULT_POLYGENIC_SCRIPT = ROOT / "polygenic_selection.py"
DEFAULT_COUNTS = DEFAULT_INPUT_DIR / "counts_joint_non_eur.tsv.gz"
DEFAULT_ADMIXTURE = DEFAULT_INPUT_DIR / "admixture_joint_non_eur.tsv"
DEFAULT_OUTPUT = DEFAULT_OUTPUT_DIR / "polygenic_joint_non_eur.tsv"
DEFAULT_P_VALUE_TABLE = DEFAULT_OUTPUT_DIR / "deep_schizophrenia_p_value.tsv"
SCHIZOPHRENIA_ALIASES = ["schizophrenia", "schiz"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Recompute schizophrenia joint_non_eur with more sign permutations. "
            "This leaves the locus set and count/admixture inputs unchanged."
        )
    )
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--counts", type=Path, default=None)
    parser.add_argument("--admixture", type=Path, default=None)
    parser.add_argument("--polygenic-script", type=Path, default=DEFAULT_POLYGENIC_SCRIPT)
    parser.add_argument("--n-permutations", type=int, default=10_000_000)
    parser.add_argument("--block-size", type=int, default=250)
    parser.add_argument("--seed", type=int, default=1)
    return parser.parse_args()


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="")
    return path.open(newline="")


def display_path(path: Path) -> str:
    try:
        return path.resolve().relative_to(ROOT).as_posix()
    except ValueError:
        return path.as_posix()


def select_schizophrenia_trait(counts_path: Path) -> str:
    traits: set[str] = set()
    with open_text(counts_path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "trait" not in (reader.fieldnames or []):
            raise ValueError(f"Counts file must include a trait column: {counts_path}")
        for row in reader:
            trait = row["trait"].strip()
            if trait in SCHIZOPHRENIA_ALIASES:
                return trait
            traits.add(trait)
    raise ValueError(f"No schizophrenia trait rows found in {counts_path}. Saw traits: {sorted(traits)}")


def run_polygenic(args: argparse.Namespace, counts_path: Path, admixture_path: Path, output_path: Path, trait: str) -> None:
    cmd = [
        sys.executable,
        str(args.polygenic_script),
        "--counts",
        str(counts_path),
        "--admixture",
        str(admixture_path),
        "--traits",
        trait,
        "--n-permutations",
        str(args.n_permutations),
        "--block-size",
        str(args.block_size),
        "--seed",
        str(args.seed),
        "--output",
        str(output_path),
    ]
    subprocess.run(cmd, check=True)


def load_schizophrenia_row(path: Path) -> dict[str, str]:
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["trait"] in SCHIZOPHRENIA_ALIASES:
                return row
    raise RuntimeError(f"No schizophrenia row found in {path}")


def write_p_value_table(
    path: Path,
    row: dict[str, str],
    args: argparse.Namespace,
    counts_path: Path,
    admixture_path: Path,
) -> None:
    fieldnames = [
        "analysis",
        "trait",
        "directional_p",
        "n_permutations",
        "seed",
        "block_size",
        "input_counts",
        "input_admixture",
        "source",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(
            {
                "analysis": "joint_non_eur",
                "trait": "schizophrenia",
                "directional_p": row["directional_p"],
                "n_permutations": str(args.n_permutations),
                "seed": str(args.seed),
                "block_size": str(args.block_size),
                "input_counts": display_path(counts_path),
                "input_admixture": display_path(admixture_path),
                "source": display_path(path),
            }
        )


def main() -> int:
    args = parse_args()
    counts_path = args.counts or args.input_dir / DEFAULT_COUNTS.name
    admixture_path = args.admixture or args.input_dir / DEFAULT_ADMIXTURE.name
    output_path = args.output_dir / DEFAULT_OUTPUT.name
    p_value_path = args.output_dir / DEFAULT_P_VALUE_TABLE.name

    for path in [counts_path, admixture_path, args.polygenic_script]:
        if not path.exists():
            raise FileNotFoundError(f"Missing required input: {path}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    trait = select_schizophrenia_trait(counts_path)
    run_polygenic(args, counts_path, admixture_path, output_path, trait)
    row = load_schizophrenia_row(output_path)
    write_p_value_table(p_value_path, row, args, counts_path, admixture_path)
    print(f"Wrote {p_value_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
