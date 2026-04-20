#!/usr/bin/env python3
"""Run the admixture-aware polygenic sign-permutation test from Colbran et al.

This file keeps the sign-only empirical test close to the paper:

- per-population trait-increasing allele counts are compared against an
  admixture-based null
- the test statistic is a likelihood-ratio comparing the admixture null to an
  unconstrained per-population alternative
- empirical p-values come from random sign flips across loci, mirroring the
  paper's permutation logic for polygenic selection

This script is used as-is by the trait runners. The extension work happens in
trait construction and data preparation, not by changing this core sign test.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

import numpy as np
from scipy.optimize import minimize


EPS = 1e-9


def open_text(path: Path, mode: str) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, mode + "t", newline="")
    return path.open(mode, newline="")


def detect_delimiter(path: Path) -> str:
    if path.suffix == ".csv" or path.name.endswith(".csv.gz"):
        return ","
    return "\t"


def clip_prob(values: np.ndarray) -> np.ndarray:
    return np.clip(values, EPS, 1.0 - EPS)


def loglik_binomial(x: np.ndarray, n: np.ndarray, p: np.ndarray) -> float:
    mask = n > 0
    if not np.any(mask):
        return float("nan")
    x_m = x[mask]
    n_m = n[mask]
    p_m = clip_prob(p[mask])
    return float(np.sum(x_m * np.log(p_m) + (n_m - x_m) * np.log1p(-p_m)))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run the paper's admixture-aware polygenic test using trait-oriented "
            "population allele counts."
        )
    )
    parser.add_argument("--counts", type=Path, required=True)
    parser.add_argument("--admixture", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--traits", nargs="*", default=None)
    parser.add_argument("--exclude-populations", nargs="*", default=[])
    parser.add_argument("--counts-delimiter", default=None)
    parser.add_argument("--admixture-delimiter", default=None)
    parser.add_argument("--n-permutations", type=int, default=10000)
    parser.add_argument("--block-size", type=int, default=128)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--population-col", default="population")
    parser.add_argument("--trait-col", default="trait")
    parser.add_argument("--variant-col", default="variant_id")
    parser.add_argument("--effect-count-col", default="effect_allele_count")
    parser.add_argument("--total-count-col", default="non_missing_allele_count")
    parser.add_argument("--component-cols", nargs="*", default=None)
    return parser.parse_args()


@dataclass
class TraitMatrix:
    variants: list[str]
    populations: list[str]
    effect_counts: np.ndarray
    total_counts: np.ndarray


def load_admixture(
    path: Path,
    delimiter: str,
    population_col: str,
    component_cols: list[str] | None,
    excluded_populations: set[str],
) -> tuple[list[str], np.ndarray]:
    rows: list[tuple[str, list[float]]] = []
    with open_text(path, "r") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if not reader.fieldnames or population_col not in reader.fieldnames:
            raise ValueError(f"Admixture file must include '{population_col}'.")
        if component_cols is None:
            component_cols = [c for c in reader.fieldnames if c != population_col]
        missing = set(component_cols).difference(reader.fieldnames)
        if missing:
            raise ValueError(f"Missing admixture columns: {sorted(missing)}")

        for row in reader:
            population = row[population_col].strip()
            if population in excluded_populations:
                continue
            components = [float(row[col]) for col in component_cols]
            rows.append((population, components))

    if not rows:
        raise ValueError("No admixture rows remained after population filtering.")

    populations = [population for population, _ in rows]
    matrix = np.asarray([components for _, components in rows], dtype=float)
    return populations, matrix


def load_counts(
    path: Path,
    delimiter: str,
    trait_col: str,
    variant_col: str,
    population_col: str,
    effect_count_col: str,
    total_count_col: str,
    selected_traits: set[str] | None,
    excluded_populations: set[str],
) -> dict[str, dict[tuple[str, str], tuple[float, float]]]:
    by_trait: dict[str, dict[tuple[str, str], tuple[float, float]]] = defaultdict(dict)
    with open_text(path, "r") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        required = {
            trait_col,
            variant_col,
            population_col,
            effect_count_col,
            total_count_col,
        }
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing count columns: {sorted(missing)}")

        for row in reader:
            trait = row[trait_col].strip()
            if selected_traits and trait not in selected_traits:
                continue

            population = row[population_col].strip()
            if population in excluded_populations:
                continue

            variant_id = row[variant_col].strip()
            effect_count = float(row[effect_count_col])
            total_count = float(row[total_count_col])
            key = (variant_id, population)

            prev = by_trait[trait].get(key)
            if prev is None:
                by_trait[trait][key] = (effect_count, total_count)
            else:
                by_trait[trait][key] = (
                    prev[0] + effect_count,
                    prev[1] + total_count,
                )

    return by_trait


def build_trait_matrix(
    trait_counts: dict[tuple[str, str], tuple[float, float]],
    populations: list[str],
) -> TraitMatrix:
    population_index = {population: idx for idx, population in enumerate(populations)}
    variants = sorted({variant for variant, _ in trait_counts})
    effect_counts = np.zeros((len(variants), len(populations)), dtype=float)
    total_counts = np.zeros((len(variants), len(populations)), dtype=float)

    for variant_idx, variant in enumerate(variants):
        for population in populations:
            value = trait_counts.get((variant, population))
            if value is None:
                continue
            pop_idx = population_index[population]
            effect_counts[variant_idx, pop_idx] = value[0]
            total_counts[variant_idx, pop_idx] = value[1]

    keep = np.any(total_counts > 0, axis=1)
    variants = [variant for variant, flag in zip(variants, keep) if flag]
    effect_counts = effect_counts[keep]
    total_counts = total_counts[keep]

    if not variants:
        raise ValueError("Trait has no non-missing variants after filtering.")

    return TraitMatrix(
        variants=variants,
        populations=populations,
        effect_counts=effect_counts,
        total_counts=total_counts,
    )


def fit_null_q(x: np.ndarray, n: np.ndarray, admixture: np.ndarray) -> tuple[np.ndarray, float]:
    mask = n > 0
    if not np.any(mask):
        raise ValueError("Cannot fit null model with zero total counts in every population.")

    x_m = x[mask]
    n_m = n[mask]
    admixture_m = admixture[mask]
    k = admixture.shape[1]
    pooled = float(np.sum(x_m) / np.sum(n_m))
    q0 = np.full(k, min(max(pooled, EPS), 1.0 - EPS))

    def objective(q: np.ndarray) -> float:
        f = clip_prob(admixture_m @ q)
        return -loglik_binomial(x_m, n_m, f)

    def gradient(q: np.ndarray) -> np.ndarray:
        f = clip_prob(admixture_m @ q)
        grad = admixture_m.T @ (x_m / f - (n_m - x_m) / (1.0 - f))
        return -grad

    bounds = [(EPS, 1.0 - EPS)] * k
    result = minimize(objective, q0, jac=gradient, bounds=bounds, method="L-BFGS-B")
    if not result.success:
        fallback = minimize(objective, q0, bounds=bounds, method="SLSQP")
        if not fallback.success:
            raise RuntimeError(
                "Null optimization failed: "
                f"L-BFGS-B={result.message!r}, SLSQP={fallback.message!r}"
            )
        result = fallback

    q_hat = clip_prob(result.x)
    expected = clip_prob(admixture_m @ q_hat)
    loglik = loglik_binomial(x_m, n_m, expected)
    return q_hat, loglik


def likelihood_ratio_test(
    x: np.ndarray,
    n: np.ndarray,
    admixture: np.ndarray,
) -> tuple[float, float, float, np.ndarray, np.ndarray, np.ndarray]:
    observed_freq = np.divide(
        x,
        n,
        out=np.full_like(x, 0.5, dtype=float),
        where=n > 0,
    )
    alt_loglik = loglik_binomial(x, n, observed_freq)
    q_hat, null_loglik = fit_null_q(x, n, admixture)
    expected_freq = clip_prob(admixture @ q_hat)
    statistic = max(0.0, 2.0 * (alt_loglik - null_loglik))
    return statistic, alt_loglik, null_loglik, q_hat, expected_freq, observed_freq


def run_permutations(
    effect_counts: np.ndarray,
    total_counts: np.ndarray,
    admixture: np.ndarray,
    observed_statistic: float,
    n_permutations: int,
    block_size: int,
    seed: int,
) -> tuple[float, float]:
    if n_permutations <= 0:
        return float("nan"), float("nan")

    rng = np.random.default_rng(seed)
    baseline = np.sum(total_counts - effect_counts, axis=0)
    delta = 2.0 * effect_counts - total_counts
    total_by_population = np.sum(total_counts, axis=0)
    ge_count = 0
    le_count = 0

    for start in range(0, n_permutations, block_size):
        block = min(block_size, n_permutations - start)
        signs = rng.integers(0, 2, size=(block, effect_counts.shape[0]), dtype=np.int8)
        aggregated = baseline + signs.astype(np.float64) @ delta
        for permuted_x in aggregated:
            statistic = likelihood_ratio_test(permuted_x, total_by_population, admixture)[0]
            if statistic >= observed_statistic - 1e-12:
                ge_count += 1
            if statistic <= observed_statistic + 1e-12:
                le_count += 1

    upper_tail_p = (ge_count + 1.0) / (n_permutations + 1.0)
    lower_tail_p = (le_count + 1.0) / (n_permutations + 1.0)
    return upper_tail_p, lower_tail_p


def json_vector(values: np.ndarray) -> str:
    return json.dumps([round(float(value), 8) for value in values.tolist()])


def main() -> int:
    args = parse_args()
    selected_traits = set(args.traits) if args.traits else None
    excluded_populations = set(args.exclude_populations)
    counts_delimiter = args.counts_delimiter or detect_delimiter(args.counts)
    admixture_delimiter = args.admixture_delimiter or detect_delimiter(args.admixture)

    populations, admixture = load_admixture(
        path=args.admixture,
        delimiter=admixture_delimiter,
        population_col=args.population_col,
        component_cols=args.component_cols,
        excluded_populations=excluded_populations,
    )
    by_trait = load_counts(
        path=args.counts,
        delimiter=counts_delimiter,
        trait_col=args.trait_col,
        variant_col=args.variant_col,
        population_col=args.population_col,
        effect_count_col=args.effect_count_col,
        total_count_col=args.total_count_col,
        selected_traits=selected_traits,
        excluded_populations=excluded_populations,
    )
    if not by_trait:
        raise ValueError("No trait rows remained after filtering.")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open_text(args.output, "w") as handle:
        fieldnames = [
            "trait",
            "n_variants",
            "n_populations",
            "observed_statistic",
            "directional_p",
            "lower_tail_p",
            "alt_loglik",
            "null_loglik",
            "populations_json",
            "source_frequency_json",
            "observed_frequency_json",
            "expected_frequency_json",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for trait in sorted(by_trait):
            matrix = build_trait_matrix(by_trait[trait], populations)
            total_effect = np.sum(matrix.effect_counts, axis=0)
            total_counts = np.sum(matrix.total_counts, axis=0)

            statistic, alt_loglik, null_loglik, q_hat, expected_freq, observed_freq = (
                likelihood_ratio_test(total_effect, total_counts, admixture)
            )
            directional_p, lower_tail_p = run_permutations(
                effect_counts=matrix.effect_counts,
                total_counts=matrix.total_counts,
                admixture=admixture,
                observed_statistic=statistic,
                n_permutations=args.n_permutations,
                block_size=args.block_size,
                seed=args.seed,
            )

            writer.writerow(
                {
                    "trait": trait,
                    "n_variants": len(matrix.variants),
                    "n_populations": len(matrix.populations),
                    "observed_statistic": f"{statistic:.8f}",
                    "directional_p": (
                        "" if math.isnan(directional_p) else f"{directional_p:.8g}"
                    ),
                    "lower_tail_p": (
                        "" if math.isnan(lower_tail_p) else f"{lower_tail_p:.8g}"
                    ),
                    "alt_loglik": f"{alt_loglik:.8f}",
                    "null_loglik": f"{null_loglik:.8f}",
                    "populations_json": json.dumps(matrix.populations),
                    "source_frequency_json": json_vector(q_hat),
                    "observed_frequency_json": json_vector(observed_freq),
                    "expected_frequency_json": json_vector(expected_freq),
                }
            )

            print(
                f"{trait}: n_variants={len(matrix.variants)} "
                f"stat={statistic:.4f} directional_p={directional_p:.6g} "
                f"lower_tail_p={lower_tail_p:.6g}",
                file=sys.stderr,
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
