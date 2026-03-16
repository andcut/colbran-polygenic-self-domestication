#!/usr/bin/env python3
"""Plot a Figure 6-style comparison with the public extension traits.

This plot intentionally mirrors the sign-only non-European comparison in
Colbran et al. Figure 6A. We therefore plot the empirical sign-permutation
`joint_non_eur` p-values for the public extension traits rather than any
weighted or exploratory variants from the parent research repository.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_COLBRAN_TABLE = ROOT / "data" / "colbran_eur_polygenic_table3.tsv"
DEFAULT_OUTPUT_DIR = ROOT / "results" / "figure6_comparison_100k"
DEFAULT_SUMMARIES = [
    ROOT / "results" / "extension_100k" / "polygenic_summary.tsv",
    ROOT / "results" / "cognitive_comparison_100k" / "polygenic_summary.tsv",
]

OUR_LABELS = {
    "gfp": "GFP",
    "schizophrenia": "Schizophrenia",
    "iq": "IQ",
    "ea4": "Educational attainment",
}

# The public extension runs default to 100,000 permutations, matching the
# repository's standard release artifacts. Schizophrenia hits the empirical
# floor at that depth in `joint_non_eur`, so for the headline Figure 6-style
# comparison we use the deeper 10,000,000-permutation estimate for that one
# plotted point. This only changes the displayed p-value, not the trait
# construction, locus set, or test statistic.
P_VALUE_OVERRIDES = {
    ("schizophrenia", "joint_non_eur"): 5.9999994e-07,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create the Figure 6-style panel for this repository.")
    parser.add_argument("--colbran-table", type=Path, default=DEFAULT_COLBRAN_TABLE)
    parser.add_argument("--our-summary", dest="our_summaries", action="append", type=Path)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--analysis", default="joint_non_eur")
    return parser.parse_args()


def figure6_transform(p_value: float) -> float:
    p = min(max(float(p_value), 1e-8), 1.0 - 1e-8)
    if p < 0.5:
        return math.log10(0.5) - math.log10(p)
    if p > 0.5:
        return math.log10(1.0 - p) - math.log10(0.5)
    return 0.0


def load_colbran_rows(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["figure6_nonredundant"].strip().lower() != "yes":
                continue
            p_value = float(row["non_eur_p"])
            rows.append(
                {
                    "label": row["plot_label"],
                    "source": "Colbran",
                    "p_value": p_value,
                    "x": figure6_transform(p_value),
                    "color": "black",
                    "marker": "o",
                    "size": 50,
                }
            )
    return rows


def load_our_rows(paths: list[Path], analysis: str) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for path in paths:
        with path.open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                if row["analysis"] != analysis:
                    continue
                trait = row["trait"]
                if trait not in OUR_LABELS:
                    continue
                p_value = P_VALUE_OVERRIDES.get((trait, analysis), float(row["directional_p"]))
                rows.append(
                    {
                        "label": OUR_LABELS[trait],
                        "source": "Extension",
                        "p_value": p_value,
                        "x": figure6_transform(p_value),
                        "color": "#0d6efd",
                        "marker": "D",
                        "size": 66,
                    }
                )
    return rows


def figure6_ticks() -> tuple[list[float], list[str]]:
    tick_ps = [1 - 1e-8, 1 - 1e-6, 1 - 1e-4, 1 - 1e-2, 0.9, 0.5, 0.1, 1e-2, 1e-4, 1e-6, 1e-8]
    tick_labels = ["1-10^-8", "1-10^-6", "1-10^-4", "1-10^-2", "0.9", "0.5", "0.1", "10^-2", "10^-4", "10^-6", "10^-8"]
    return [figure6_transform(p) for p in tick_ps], tick_labels


def save_plot(rows: list[dict[str, object]], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    rows = sorted(rows, key=lambda row: (row["x"], row["label"]))
    ticks, tick_labels = figure6_ticks()
    bonferroni_count = sum(1 for row in rows if row["source"] == "Colbran")
    directional_threshold = 0.025 / bonferroni_count
    stabilizing_threshold = 1.0 - directional_threshold

    fig_height = max(7.0, 0.35 * len(rows) + 1.8)
    fig, ax = plt.subplots(figsize=(11.5, fig_height))

    for y, row in enumerate(rows):
        ax.scatter(
            row["x"],
            y,
            s=row["size"],
            c=row["color"],
            marker=row["marker"],
            edgecolors="black" if row["source"] == "Extension" else row["color"],
            linewidths=0.6 if row["source"] == "Extension" else 0.0,
            zorder=3,
        )

    ax.axvline(0.0, color="0.75", linewidth=1.5, zorder=1)
    ax.axvline(figure6_transform(stabilizing_threshold), color="#e36b6b", linewidth=1.5, zorder=1)
    ax.axvline(figure6_transform(directional_threshold), color="#e36b6b", linewidth=1.5, zorder=1)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([row["label"] for row in rows], fontsize=10)
    ax.set_xticks(ticks)
    ax.set_xticklabels(tick_labels, fontsize=10)
    ax.set_xlabel("Empirical p-value", fontsize=11)
    ax.set_title("Several self-domestication traits show evidence for polygenic selection", fontsize=13)
    ax.set_xlim(max(ticks) + 0.2, min(ticks) - 0.2)
    ax.set_ylim(-0.8, len(rows) - 0.2)
    ax.grid(axis="x", color="0.92", linewidth=0.8)
    ax.set_axisbelow(True)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    legend_handles = [
        plt.Line2D([0], [0], marker="o", color="black", linestyle="", markersize=7, label="Colbran traits"),
        plt.Line2D([0], [0], marker="D", markerfacecolor="#0d6efd", markeredgecolor="black", color="#0d6efd", linestyle="", markersize=7, label="Self-domestication traits"),
    ]
    ax.legend(handles=legend_handles, loc="upper right", frameon=False)

    fig.tight_layout()
    fig.savefig(output_dir / "figure6_comparison.png", dpi=220, bbox_inches="tight")
    fig.savefig(output_dir / "figure6_comparison.svg", bbox_inches="tight")
    plt.close(fig)


def save_table(rows: list[dict[str, object]], output_dir: Path) -> None:
    path = output_dir / "figure6_comparison.tsv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["label", "source", "p_value", "figure6_x"], delimiter="\t")
        writer.writeheader()
        for row in sorted(rows, key=lambda row: (row["x"], row["label"])):
            writer.writerow(
                {
                    "label": row["label"],
                    "source": row["source"],
                    "p_value": f"{row['p_value']:.8g}",
                    "figure6_x": f"{row['x']:.8f}",
                }
            )


def main() -> int:
    args = parse_args()
    our_summaries = args.our_summaries or DEFAULT_SUMMARIES
    rows = load_colbran_rows(args.colbran_table) + load_our_rows(our_summaries, args.analysis)
    if not rows:
        raise RuntimeError("No rows were loaded for the Figure 6-style plot.")
    save_plot(rows, args.output_dir)
    save_table(rows, args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
