#!/usr/bin/env python3
"""Plot a Figure 6-style comparison with the extension traits.

This plot intentionally mirrors the sign-only non-European comparison in
Colbran et al. Figure 6A. We therefore plot the empirical sign-permutation
`joint_non_eur` p-values for GFP, schizophrenia, IQ, and EA4.
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
DEFAULT_DEEP_SCHIZ_P_VALUES = ROOT / "results" / "deep_schizophrenia_10m" / "deep_schizophrenia_p_value.tsv"

OUR_LABELS = {
    "gfp": "GFP",
    "schizophrenia": "Schizophrenia",
    "iq": "IQ",
    "ea4": "Educational attainment",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create the Figure 6-style panel for this repository.")
    parser.add_argument("--colbran-table", type=Path, default=DEFAULT_COLBRAN_TABLE)
    parser.add_argument("--our-summary", dest="our_summaries", action="append", type=Path)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--analysis", default="joint_non_eur")
    parser.add_argument(
        "--deep-schiz-p-values",
        type=Path,
        default=DEFAULT_DEEP_SCHIZ_P_VALUES,
        help=(
            "Optional p-value table produced by scripts/run_deep_schizophrenia.py. "
            "When present, it supplies the 10,000,000-permutation schizophrenia "
            "joint_non_eur value used in the headline plot."
        ),
    )
    parser.add_argument(
        "--fast",
        action="store_true",
        help="Use the standard 100,000-permutation summaries only.",
    )
    parser.add_argument(
        "--require-deep-schiz",
        action="store_true",
        help="Fail if the deep schizophrenia p-value table is missing.",
    )
    parser.add_argument(
        "--allow-missing-deep-schiz",
        action="store_true",
        help="If the deep schizophrenia p-value table is missing, fall back to the standard 100,000-permutation value.",
    )
    return parser.parse_args()


def figure6_transform(p_value: float) -> float:
    p = min(max(float(p_value), 1e-8), 1.0 - 1e-8)
    if p < 0.5:
        return math.log10(0.5) - math.log10(p)
    if p > 0.5:
        return math.log10(1.0 - p) - math.log10(0.5)
    return 0.0


def display_path(path: Path) -> str:
    try:
        return path.resolve().relative_to(ROOT).as_posix()
    except ValueError:
        return path.as_posix()


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
                    "p_value_source": "Colbran Table 3 transcription",
                    "n_permutations": "NA",
                    "color": "black",
                    "marker": "o",
                    "size": 50,
                }
            )
    return rows


def load_deep_p_values(path: Path | None, require: bool) -> dict[tuple[str, str], dict[str, str]]:
    if path is None or not path.exists():
        if require:
            raise FileNotFoundError(f"Missing deep schizophrenia p-value table: {path}")
        return {}

    values: dict[tuple[str, str], dict[str, str]] = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"analysis", "trait", "directional_p"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing columns in {path}: {sorted(missing)}")
        for row in reader:
            analysis = row["analysis"].strip()
            trait = row["trait"].strip()
            values[(trait, analysis)] = row
    return values


def load_our_rows(paths: list[Path], analysis: str, deep_p_values: dict[tuple[str, str], dict[str, str]]) -> list[dict[str, object]]:
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
                deep_row = deep_p_values.get((trait, analysis))
                if deep_row is None:
                    p_value = float(row["directional_p"])
                    p_value_source = display_path(path)
                    n_permutations = "100000"
                else:
                    p_value = float(deep_row["directional_p"])
                    p_value_source = deep_row.get("source") or deep_row.get("source_path") or display_path(path)
                    n_permutations = deep_row.get("n_permutations", "")
                rows.append(
                    {
                        "label": OUR_LABELS[trait],
                        "source": "Extension",
                        "p_value": p_value,
                        "x": figure6_transform(p_value),
                        "p_value_source": p_value_source,
                        "n_permutations": n_permutations,
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
        writer = csv.DictWriter(
            handle,
            fieldnames=["label", "source", "p_value", "figure6_x", "p_value_source", "n_permutations"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in sorted(rows, key=lambda row: (row["x"], row["label"])):
            writer.writerow(
                {
                    "label": row["label"],
                    "source": row["source"],
                    "p_value": f"{row['p_value']:.8g}",
                    "figure6_x": f"{row['x']:.8f}",
                    "p_value_source": row["p_value_source"],
                    "n_permutations": row["n_permutations"],
                }
            )


def main() -> int:
    args = parse_args()
    our_summaries = args.our_summaries or DEFAULT_SUMMARIES
    deep_path = None if args.fast else args.deep_schiz_p_values
    require_deep = not args.fast and not args.allow_missing_deep_schiz
    deep_p_values = load_deep_p_values(deep_path, require=require_deep or args.require_deep_schiz)
    rows = load_colbran_rows(args.colbran_table) + load_our_rows(our_summaries, args.analysis, deep_p_values)
    if not rows:
        raise RuntimeError("No rows were loaded for the Figure 6-style plot.")
    save_plot(rows, args.output_dir)
    save_table(rows, args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
