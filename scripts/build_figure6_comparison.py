#!/usr/bin/env python3
"""Build the Figure 6-style comparison from the repository inputs."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_EXTENSION_DIR = ROOT / "results" / "extension_100k"
DEFAULT_COGNITIVE_DIR = ROOT / "results" / "cognitive_comparison_100k"
DEFAULT_DEEP_DIR = ROOT / "results" / "deep_schizophrenia_10m"
DEFAULT_FIGURE_DIR = ROOT / "results" / "figure6_comparison_100k"
DEEP_P_VALUE_NAME = "deep_schizophrenia_p_value.tsv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run the standard 100,000-permutation analyses and make the "
            "Figure 6-style comparison plot."
        )
    )
    parser.add_argument("--mode", choices=["exact", "fast"], default="exact")
    parser.add_argument(
        "--fast",
        action="store_true",
        help="Alias for --mode fast.",
    )
    parser.add_argument(
        "--use-existing",
        action="store_true",
        help="Use existing result tables instead of rerunning the source-to-result analyses.",
    )
    parser.add_argument("--extension-dir", type=Path, default=DEFAULT_EXTENSION_DIR)
    parser.add_argument("--cognitive-dir", type=Path, default=DEFAULT_COGNITIVE_DIR)
    parser.add_argument("--deep-schiz-dir", type=Path, default=DEFAULT_DEEP_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    parser.add_argument("--standard-permutations", type=int, default=100_000)
    parser.add_argument("--deep-schiz-permutations", type=int, default=10_000_000)
    parser.add_argument("--block-size", type=int, default=250)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--plink-path", type=Path, default=None)
    parser.add_argument(
        "--force-deep",
        action="store_true",
        help="Rerun the deep schizophrenia step even if its p-value table exists.",
    )
    return parser.parse_args()


def run(cmd: list[str]) -> None:
    print(" ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def append_plink(cmd: list[str], plink_path: Path | None) -> list[str]:
    if plink_path is not None:
        cmd.extend(["--plink-path", str(plink_path)])
    return cmd


def run_standard_analyses(args: argparse.Namespace) -> None:
    extension_cmd = [
        sys.executable,
        str(ROOT / "scripts" / "run_extension.py"),
        "--output-dir",
        str(args.extension_dir),
        "--n-permutations",
        str(args.standard_permutations),
        "--block-size",
        str(args.block_size),
        "--seed",
        str(args.seed),
        "--jobs",
        str(args.jobs),
    ]
    run(append_plink(extension_cmd, args.plink_path))

    cognitive_cmd = [
        sys.executable,
        str(ROOT / "scripts" / "run_cognitive_comparison.py"),
        "--output-dir",
        str(args.cognitive_dir),
        "--n-permutations",
        str(args.standard_permutations),
        "--block-size",
        str(args.block_size),
        "--seed",
        str(args.seed),
        "--jobs",
        str(args.jobs),
    ]
    run(append_plink(cognitive_cmd, args.plink_path))


def ensure_deep_schizophrenia(args: argparse.Namespace) -> Path:
    p_value_path = args.deep_schiz_dir / DEEP_P_VALUE_NAME
    if args.use_existing and p_value_path.exists() and not args.force_deep:
        return p_value_path
    if p_value_path.exists() and not args.force_deep:
        return p_value_path

    cmd = [
        sys.executable,
        str(ROOT / "scripts" / "run_deep_schizophrenia.py"),
        "--input-dir",
        str(args.extension_dir),
        "--output-dir",
        str(args.deep_schiz_dir),
        "--n-permutations",
        str(args.deep_schiz_permutations),
        "--block-size",
        str(args.block_size),
        "--seed",
        str(args.seed),
    ]
    run(cmd)
    return p_value_path


def plot(args: argparse.Namespace, deep_p_value_path: Path | None) -> None:
    cmd = [
        sys.executable,
        str(ROOT / "scripts" / "plot_figure6_comparison.py"),
        "--our-summary",
        str(args.extension_dir / "polygenic_summary.tsv"),
        "--our-summary",
        str(args.cognitive_dir / "polygenic_summary.tsv"),
        "--output-dir",
        str(args.output_dir),
    ]
    if args.mode == "fast":
        cmd.append("--fast")
    else:
        assert deep_p_value_path is not None
        cmd.extend(["--deep-schiz-p-values", str(deep_p_value_path), "--require-deep-schiz"])
    run(cmd)


def main() -> int:
    args = parse_args()
    if args.fast:
        args.mode = "fast"
    if not args.use_existing:
        run_standard_analyses(args)

    deep_p_value_path = None
    if args.mode == "exact":
        deep_p_value_path = ensure_deep_schizophrenia(args)

    plot(args, deep_p_value_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
