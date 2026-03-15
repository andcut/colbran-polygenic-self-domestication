#!/usr/bin/env python3
"""Plot ancient-sample coverage by region and millennium.

This figure is meant to orient readers to the actual temporal lens available
to the Colbran-style polygenic test in this repository. The p-values are driven
by allele frequencies in the ancient-data panels, so it is helpful to show how
many ancient individuals exist in each region across time.

We exclude modern samples (`Date == 0`) here on purpose. The point of the
figure is not total sample size in the sample-info tables; it is the dated
ancient coverage that anchors the ancestry model historically.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SAMPLE_INFO_DIR = ROOT / "data" / "sample_info"
OUTPUT_DIR = ROOT / "results" / "sample_coverage"

REGIONS = [
    ("afr", "Africa", "#D55E00"),
    ("eur", "Europe", "#0072B2"),
    ("eas", "East Asia", "#009E73"),
    ("sas", "South Asia", "#CC79A7"),
    ("sam", "Americas", "#E69F00"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot ancient sample counts by region and millennium.")
    parser.add_argument("--sample-info-dir", type=Path, default=SAMPLE_INFO_DIR)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--bin-years", type=int, default=1000)
    return parser.parse_args()


def load_counts(sample_info_dir: Path, bin_years: int) -> tuple[list[int], dict[str, dict[int, int]], dict[str, int]]:
    counts_by_region: dict[str, dict[int, int]] = {}
    totals: dict[str, int] = {}
    all_bins: set[int] = set()
    for region_key, _, _ in REGIONS:
        path = sample_info_dir / f"{region_key}_sample_info.txt"
        bin_counts: dict[int, int] = defaultdict(int)
        total = 0
        with path.open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                date = float(row["Date"])
                if date <= 0:
                    continue
                bin_start = int(date // bin_years) * bin_years
                bin_counts[bin_start] += 1
                all_bins.add(bin_start)
                total += 1
        counts_by_region[region_key] = dict(bin_counts)
        totals[region_key] = total
    return sorted(all_bins), counts_by_region, totals


def write_table(path: Path, bins: list[int], counts_by_region: dict[str, dict[int, int]]) -> None:
    fieldnames = ["bin_start_bp", "bin_end_bp"] + [region_key for region_key, _, _ in REGIONS]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for bin_start in bins:
            row = {
                "bin_start_bp": bin_start,
                "bin_end_bp": bin_start + 1000,
            }
            for region_key, _, _ in REGIONS:
                row[region_key] = counts_by_region[region_key].get(bin_start, 0)
            writer.writerow(row)


def plot_counts(path: Path, bins: list[int], counts_by_region: dict[str, dict[int, int]], totals: dict[str, int], bin_years: int) -> None:
    fig, ax = plt.subplots(figsize=(13.5, 6.8))
    x = np.arange(len(bins))
    width = 0.16
    offsets = np.linspace(-2, 2, num=len(REGIONS)) * width

    for offset, (region_key, region_label, color) in zip(offsets, REGIONS):
        heights = [counts_by_region[region_key].get(bin_start, 0) for bin_start in bins]
        ax.bar(
            x + offset,
            heights,
            width=width,
            color=color,
            label=f"{region_label} (n={totals[region_key]})",
            edgecolor="white",
            linewidth=0.4,
        )

    ax.set_yscale("log")
    ax.set_ylabel("Ancient individuals per 1,000-year bin", fontsize=11)
    ax.set_xlabel("Years before present (bin start)", fontsize=11)
    ax.set_title("Ancient sample coverage by region across time", fontsize=13)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{bin_start:,}" for bin_start in bins], rotation=45, ha="right")
    ax.grid(axis="y", color="0.9", linewidth=0.8)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, ncol=2, fontsize=9)
    ax.text(
        0.99,
        0.98,
        "Modern samples excluded (Date > 0 only)\nLog y-axis to keep Europe from flattening other regions",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=9,
        color="0.35",
    )
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=220, bbox_inches="tight")
    fig.savefig(path.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    bins, counts_by_region, totals = load_counts(args.sample_info_dir, args.bin_years)
    if not bins:
        raise RuntimeError("No ancient samples were found in the sample-info tables.")
    write_table(args.output_dir / "ancient_sample_coverage.tsv", bins, counts_by_region)
    plot_counts(args.output_dir / "ancient_sample_coverage", bins, counts_by_region, totals, args.bin_years)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
