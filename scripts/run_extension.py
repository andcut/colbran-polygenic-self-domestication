#!/usr/bin/env python3
"""Build the GFP extension and compare it to schizophrenia.

Method choices in this script follow Colbran, Terhorst, and Mathieson,
"Global patterns of natural selection inferred using ancient DNA",
especially the paper's polygenic-selection methods:

- filter to variants present on the 1240k panel
- use PLINK clumping with a 200 kb window and r^2 > 0.4
- orient to the trait-increasing allele
- run the sign-permutation test on non-European target panels

The one substantive extension is trait construction for GFP. We build a
GWAS-like GFP summary statistic from the Big Five using published GFP-style
loadings, then feed that synthetic trait through the same downstream pipeline
as schizophrenia.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import math
import subprocess
import sys
import urllib.request
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
RESULTS_DIR = ROOT / "results" / "extension_100k"
POLYGENIC_SCRIPT = ROOT / "polygenic_selection.py"

# These small sample-info tables are committed because they are lightweight,
# stable, and central to mirroring Colbran's regional analyses. The large
# genotype and GWAS inputs are not committed.
SAMPLE_INFO_DIR = DATA_DIR / "sample_info"

AADR_IND = RAW_DIR / "aadr" / "10537414"
AADR_SNP = RAW_DIR / "aadr" / "10537415"
AADR_GENO = RAW_DIR / "aadr" / "10537126"
SCZ_SUMSTATS = RAW_DIR / "full_sumstats" / "PGC3_SCZ_wave3_european.tsv.gz"
BIG5_CACHE_DIR = RAW_DIR / "cache" / "bigfive_aadr"

BIG5_URLS = {
    "agreeableness": "https://hugesumstats.yale.edu/dl/agree_MVP_GPC1_sumstat_file",
    "conscientiousness": "https://hugesumstats.yale.edu/dl/consc_MVP_GPC1_sumstat_file",
    "extraversion": "https://hugesumstats.yale.edu/dl/extra_MVP_GPC1_sumstat_file",
    "openness": "https://hugesumstats.yale.edu/dl/open_MVP_GPC1_sumstat_file",
    "neuroticism": "https://hugesumstats.yale.edu/dl/Neuro_MVP_UKB_sumstat_file",
}
TRAIT_ORDER = [
    "agreeableness",
    "conscientiousness",
    "extraversion",
    "openness",
    "neuroticism",
]

# The GFP construction follows the standard phenotypic interpretation: combine
# the Big Five using literature-informed loadings, with neuroticism reversed so
# higher scores point in the same broad direction as the other components.
LITERATURE_GFP_LOADINGS = {
    "agreeableness": 0.57,
    "conscientiousness": 0.63,
    "extraversion": 0.57,
    "openness": 0.42,
    "neuroticism": 0.62,
}

REGION_SPECS = {
    "afr": {"label": "AFR", "sample_info": "afr_sample_info.txt", "component_cols": ["k1", "k2", "k3"]},
    "eas": {"label": "EAS", "sample_info": "eas_sample_info.txt", "component_cols": ["k1", "k2", "k3"]},
    "sas": {"label": "SAS", "sample_info": "sas_sample_info.txt", "component_cols": ["k1", "k2"]},
    "sam": {"label": "SAM", "sample_info": "sam_sample_info.txt", "component_cols": ["k1", "k2"]},
}
ALL_ANALYSES = ["afr", "eas", "sas", "sam", "joint_non_eur"]
SAMPLE_ID_SUFFIXES = [
    ".AG",
    ".SG",
    ".DG",
    ".HO",
    "_genotyping_noUDG",
    "_genotyping",
    "_in.preparation",
    "_noUDG",
    "_d",
]
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


@dataclass(frozen=True)
class Hit:
    trait: str
    variant_id: str
    chrom: str
    position: int
    effect_allele: str
    other_allele: str
    effect_value: float
    p_value: float
    source_gwas: str


@dataclass(frozen=True)
class VariantRecord:
    variant_id: str
    chrom: str
    position: int
    ref_allele: str
    alt_allele: str
    snp_index: int


@dataclass(frozen=True)
class MatchedHit:
    hit: Hit
    variant: VariantRecord
    effect_allele_aadr: str
    harmonization: str


@dataclass
class GroupData:
    region: str
    population: str
    component_cols: list[str]
    admixture: np.ndarray
    union_positions: np.ndarray
    modern_mask: np.ndarray
    sample_count: int


@dataclass(frozen=True)
class VariantRow:
    rsid: str
    chrom: str
    position: int
    a1: str
    a2: str
    beta: float
    se: float
    p_value: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Reproduce the Colbran-style extension using a "
            "literature-weighted GFP and European schizophrenia GWAS."
        )
    )
    parser.add_argument("--output-dir", type=Path, default=RESULTS_DIR)
    parser.add_argument("--sample-info-dir", type=Path, default=SAMPLE_INFO_DIR)
    parser.add_argument("--aadr-ind", type=Path, default=AADR_IND)
    parser.add_argument("--aadr-snp", type=Path, default=AADR_SNP)
    parser.add_argument("--aadr-geno", type=Path, default=AADR_GENO)
    parser.add_argument("--schizophrenia-sumstats", type=Path, default=SCZ_SUMSTATS)
    parser.add_argument("--big5-cache-dir", type=Path, default=BIG5_CACHE_DIR)
    parser.add_argument("--gwas-p-threshold", type=float, default=1e-8)
    parser.add_argument("--clump-kb", type=int, default=200)
    parser.add_argument("--clump-r2", type=float, default=0.4)
    parser.add_argument("--plink-path", type=Path, default=Path("plink"))
    parser.add_argument("--polygenic-script", type=Path, default=POLYGENIC_SCRIPT)
    parser.add_argument("--n-permutations", type=int, default=100000)
    parser.add_argument("--block-size", type=int, default=250)
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--seed", type=int, default=1)
    return parser.parse_args()


def clean_value(value: object) -> str:
    return "" if value is None else str(value).strip()


def normalize_allele(value: object) -> str:
    return clean_value(value).upper()


def float_value(value: object) -> float:
    return float(value)


def canonical_sample_id(value: str) -> str:
    sample_id = clean_value(value)
    changed = True
    while changed:
        changed = False
        for suffix in SAMPLE_ID_SUFFIXES:
            if sample_id.endswith(suffix):
                sample_id = sample_id[: -len(suffix)]
                changed = True
                break
    return sample_id.casefold()


def sample_candidate_score(original: str, candidate: str) -> int:
    score = 0
    if original.endswith(".SG") and candidate.endswith(".SG"):
        score += 100
    if original.endswith(".DG") and candidate.endswith(".DG"):
        score += 50
    if original.endswith(".AG") and candidate.endswith(".AG"):
        score += 50
    if "_d" in original and ("_d." in candidate or candidate.endswith(".DG")):
        score += 25
    if "_d" not in original and "_d." not in candidate:
        score += 5
    return score


def load_aadr_variant_ids(path: Path) -> set[str]:
    variant_ids: set[str] = set()
    with path.open() as handle:
        for line in handle:
            fields = line.strip().split()
            if fields:
                variant_ids.add(fields[0])
    return variant_ids


def load_aadr_ind(path: Path) -> tuple[list[str], dict[str, int], dict[str, list[str]]]:
    sample_ids: list[str] = []
    iid_to_index: dict[str, int] = {}
    canonical_to_ids: dict[str, list[str]] = defaultdict(list)
    with path.open() as handle:
        for idx, line in enumerate(handle):
            fields = line.strip().split()
            if not fields:
                continue
            iid = fields[0]
            sample_ids.append(iid)
            iid_to_index[iid] = idx
            canonical_to_ids[canonical_sample_id(iid)].append(iid)
    return sample_ids, iid_to_index, canonical_to_ids


def resolve_aadr_iid(
    sample_iid: str,
    iid_to_index: dict[str, int],
    canonical_to_ids: dict[str, list[str]],
) -> tuple[int | None, str | None]:
    if sample_iid in iid_to_index:
        return iid_to_index[sample_iid], sample_iid
    candidates = canonical_to_ids.get(canonical_sample_id(sample_iid), [])
    if not candidates:
        return None, None
    if len(candidates) == 1:
        candidate = candidates[0]
        return iid_to_index[candidate], candidate
    scored = sorted(
        ((sample_candidate_score(sample_iid, candidate), candidate) for candidate in candidates),
        reverse=True,
    )
    if len(scored) >= 2 and scored[0][0] == scored[1][0]:
        return None, None
    candidate = scored[0][1]
    return iid_to_index[candidate], candidate


def parse_region_groups(
    sample_info_dir: Path,
    iid_to_index: dict[str, int],
    canonical_to_ids: dict[str, list[str]],
) -> tuple[dict[str, list[GroupData]], list[int]]:
    union_entries: dict[int, dict[str, object]] = {}
    ordered_union_indices: list[int] = []
    group_specs: dict[str, list[dict[str, object]]] = {}

    for region, spec in REGION_SPECS.items():
        path = sample_info_dir / spec["sample_info"]
        groups_in_region: dict[str, list[dict[str, object]]] = defaultdict(list)
        group_order: list[str] = []
        with path.open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                group_name = clean_value(row["Group"])
                iid = clean_value(row["IID"])
                ind_index, matched_iid = resolve_aadr_iid(iid, iid_to_index, canonical_to_ids)
                if ind_index is None or matched_iid is None:
                    continue
                if group_name not in groups_in_region:
                    group_order.append(group_name)
                row = dict(row)
                row["_matched_aadr_iid"] = matched_iid
                groups_in_region[group_name].append(row)
                if ind_index not in union_entries:
                    union_entries[ind_index] = {
                        "iid": matched_iid,
                        "is_modern": float_value(row["Date"]) == 0.0,
                    }
                    ordered_union_indices.append(ind_index)

        region_group_specs: list[dict[str, object]] = []
        for group_name in group_order:
            rows = groups_in_region[group_name]
            weights = np.array([float_value(row["NumSNPs"]) for row in rows], dtype=float)
            admixture = np.array(
                [
                    np.average(
                        np.array([float_value(row[col]) for row in rows], dtype=float),
                        weights=weights,
                    )
                    for col in spec["component_cols"]
                ],
                dtype=float,
            )
            if admixture.sum() > 0:
                admixture = admixture / admixture.sum()
            region_group_specs.append(
                {
                    "population": group_name,
                    "sample_iids": [clean_value(row["_matched_aadr_iid"]) for row in rows],
                    "component_cols": list(spec["component_cols"]),
                    "admixture": admixture,
                }
            )
        group_specs[region] = region_group_specs

    ordered_union_indices.sort()
    union_index_to_pos = {ind_index: pos for pos, ind_index in enumerate(ordered_union_indices)}
    modern_mask = np.array(
        [bool(union_entries[ind_index]["is_modern"]) for ind_index in ordered_union_indices],
        dtype=bool,
    )

    groups_by_region: dict[str, list[GroupData]] = {}
    for region in REGION_SPECS:
        groups: list[GroupData] = []
        for spec in group_specs[region]:
            sample_positions = []
            group_modern = []
            for iid in spec["sample_iids"]:
                pos = union_index_to_pos[iid_to_index[iid]]
                sample_positions.append(pos)
                group_modern.append(modern_mask[pos])
            groups.append(
                GroupData(
                    region=region,
                    population=str(spec["population"]),
                    component_cols=list(spec["component_cols"]),
                    admixture=np.array(spec["admixture"], dtype=float),
                    union_positions=np.array(sample_positions, dtype=int),
                    modern_mask=np.array(group_modern, dtype=bool),
                    sample_count=len(sample_positions),
                )
            )
        groups_by_region[region] = groups
    return groups_by_region, ordered_union_indices


def open_text(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode + "t", newline="")
    return path.open(mode, newline="")


def load_scz_full_hits(path: Path, aadr_variant_ids: set[str], p_threshold: float) -> list[Hit]:
    hits: list[Hit] = []
    # This loader is intentionally narrow: it supports the PGC3 European
    # schizophrenia summary-stat schema used here.
    with open_text(path, "r") as handle:
        header: list[str] | None = None
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line or line.startswith("##"):
                continue
            if header is None:
                header = line.lstrip("#").split("\t")
                continue
            assert header is not None
            values = line.split("\t")
            row = {header[idx]: values[idx] if idx < len(values) else "" for idx in range(len(header))}
            variant_id = clean_value(row.get("ID"))
            if not variant_id or variant_id not in aadr_variant_ids:
                continue
            try:
                p_value = float_value(row["PVAL"])
                beta = float_value(row["BETA"])
            except (KeyError, TypeError, ValueError):
                continue
            if not math.isfinite(p_value) or p_value > p_threshold:
                continue
            a1 = normalize_allele(row["A1"])
            a2 = normalize_allele(row["A2"])
            if not a1 or not a2:
                continue
            hits.append(
                Hit(
                    trait="schizophrenia",
                    variant_id=variant_id,
                    chrom=clean_value(row["CHROM"]),
                    position=int(round(float_value(row["POS"]))),
                    effect_allele=a1 if beta >= 0 else a2,
                    other_allele=a2 if beta >= 0 else a1,
                    effect_value=beta,
                    p_value=p_value,
                    source_gwas="Trubetskoy et al. schizophrenia GWAS (European summary statistics)",
                )
            )
    return hits


def load_aadr_rsids(path: Path) -> set[str]:
    rsids: set[str] = set()
    with path.open() as handle:
        for line in handle:
            fields = line.split()
            if fields:
                rsids.add(fields[0])
    return rsids


def stream_big5_rows(url: str, target_rsids: set[str]) -> dict[str, VariantRow]:
    rows: dict[str, VariantRow] = {}
    with urllib.request.urlopen(url, timeout=300) as response:
        text_handle = io.TextIOWrapper(response, encoding="utf-8")
        reader = csv.DictReader(text_handle, delimiter="\t")
        required = {"SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"Unexpected Big Five header from {url}: {reader.fieldnames}")
        for row in reader:
            rsid = row["SNP"]
            if rsid not in target_rsids:
                continue
            rows[rsid] = VariantRow(
                rsid=rsid,
                chrom=row["CHR"],
                position=int(row["BP"]),
                a1=row["A1"].upper(),
                a2=row["A2"].upper(),
                beta=float(row["BETA"]),
                se=float(row["SE"]),
                p_value=float(row["P"]),
            )
    return rows


def write_big5_cache(path: Path, rows: dict[str, VariantRow]) -> None:
    fieldnames = ["SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P"]
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows.values():
            writer.writerow(
                {
                    "SNP": row.rsid,
                    "CHR": row.chrom,
                    "BP": row.position,
                    "A1": row.a1,
                    "A2": row.a2,
                    "BETA": row.beta,
                    "SE": row.se,
                    "P": row.p_value,
                }
            )


def read_big5_cache(path: Path) -> dict[str, VariantRow]:
    rows: dict[str, VariantRow] = {}
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows[row["SNP"]] = VariantRow(
                rsid=row["SNP"],
                chrom=row["CHR"],
                position=int(row["BP"]),
                a1=row["A1"].upper(),
                a2=row["A2"].upper(),
                beta=float(row["BETA"]),
                se=float(row["SE"]),
                p_value=float(row["P"]),
            )
    return rows


def complement(base: str) -> str:
    return COMPLEMENT[base]


def align_to_reference(reference: VariantRow, row: VariantRow) -> tuple[float, bool]:
    if row.se == 0 or not np.isfinite(row.se):
        return 0.0, False
    z_score = row.beta / row.se
    if (row.a1, row.a2) == (reference.a1, reference.a2):
        return z_score, True
    if (row.a1, row.a2) == (reference.a2, reference.a1):
        return -z_score, True
    if len(row.a1) != 1 or len(row.a2) != 1:
        return 0.0, False
    if row.a1 not in COMPLEMENT or row.a2 not in COMPLEMENT:
        return 0.0, False
    comp = (complement(row.a1), complement(row.a2))
    if comp == (reference.a1, reference.a2):
        return z_score, True
    if comp == (reference.a2, reference.a1):
        return -z_score, True
    return 0.0, False


def two_sided_p_from_z(z_score: float) -> float:
    return math.erfc(abs(z_score) / math.sqrt(2.0))


def ensure_big5_cache(cache_dir: Path, aadr_snp: Path) -> dict[str, dict[str, VariantRow]]:
    cache_dir.mkdir(parents=True, exist_ok=True)
    aadr_rsids = load_aadr_rsids(aadr_snp)
    out: dict[str, dict[str, VariantRow]] = {}
    for trait in TRAIT_ORDER:
        cache_path = cache_dir / f"{trait}_aadr_intersection.tsv.gz"
        if not cache_path.exists():
            rows = stream_big5_rows(BIG5_URLS[trait], aadr_rsids)
            write_big5_cache(cache_path, rows)
        out[trait] = read_big5_cache(cache_path)
    return out


def construct_gfp_hits(
    trait_rows: dict[str, dict[str, VariantRow]],
    p_threshold: float,
) -> tuple[list[Hit], dict[str, object]]:
    shared_rsids = set.intersection(*(set(trait_rows[trait]) for trait in TRAIT_ORDER))
    if not shared_rsids:
        raise RuntimeError("No shared AADR rsIDs were found across the Big Five caches.")

    matrix_rows: list[list[float]] = []
    variant_rows: list[VariantRow] = []
    skipped_alignment = 0
    for rsid in shared_rsids:
        reference = trait_rows[TRAIT_ORDER[0]][rsid]
        z_values: list[float] = []
        ok = True
        for trait in TRAIT_ORDER:
            aligned_z, matched = align_to_reference(reference, trait_rows[trait][rsid])
            if not matched or not np.isfinite(aligned_z):
                ok = False
                break
            # Following the standard phenotypic GFP interpretation, we reverse
            # neuroticism so that higher values move in the "socially stable"
            # direction with the other Big Five components.
            if trait == "neuroticism":
                aligned_z *= -1.0
            z_values.append(aligned_z)
        if not ok:
            skipped_alignment += 1
            continue
        matrix_rows.append(z_values)
        variant_rows.append(reference)
    if not matrix_rows:
        raise RuntimeError("No Big Five variants survived allele harmonization.")

    z_matrix = np.asarray(matrix_rows, dtype=float)
    loadings = np.asarray([LITERATURE_GFP_LOADINGS[trait] for trait in TRAIT_ORDER], dtype=float)
    corr = np.corrcoef(z_matrix, rowvar=False)
    denom = math.sqrt(float(loadings @ corr @ loadings))
    score_z = (z_matrix @ loadings) / denom
    score_direction = np.where(score_z >= 0, 1.0, -1.0)
    support_counts = np.sum(np.sign(z_matrix) == score_direction[:, None], axis=1)

    hits: list[Hit] = []
    support_hist: Counter[int] = Counter()
    for row, z_score, support in zip(variant_rows, score_z, support_counts, strict=True):
        p_value = two_sided_p_from_z(float(z_score))
        if p_value > p_threshold:
            continue
        support_hist[int(support)] += 1
        effect_allele, other_allele = row.a1, row.a2
        if z_score < 0:
            effect_allele, other_allele = other_allele, effect_allele
        hits.append(
            Hit(
                trait="gfp",
                variant_id=row.rsid,
                chrom=row.chrom,
                position=row.position,
                effect_allele=effect_allele,
                other_allele=other_allele,
                effect_value=abs(float(z_score)),
                p_value=p_value,
                source_gwas="GFP built from Yale/Levey Big Five GWAS with literature-weighted Stouffer combination",
            )
        )

    diagnostics = {
        "shared_rsids": len(shared_rsids),
        "harmonized_rows": len(matrix_rows),
        "skipped_alignment": skipped_alignment,
        "stouffer_denom": denom,
        "support_hist": dict(sorted(support_hist.items())),
        "loadings": LITERATURE_GFP_LOADINGS,
        "per_trait_gws": {
            trait: int(np.sum(np.array([two_sided_p_from_z(float(z)) < p_threshold for z in z_matrix[:, idx]])))
            for idx, trait in enumerate(TRAIT_ORDER)
        },
    }
    return hits, diagnostics


def deduplicate_hits(hits: list[Hit]) -> list[Hit]:
    kept: dict[tuple[str, str], Hit] = {}
    for hit in hits:
        key = (hit.trait, hit.variant_id)
        prev = kept.get(key)
        if prev is None or hit.p_value < prev.p_value or (hit.p_value == prev.p_value and abs(hit.effect_value) > abs(prev.effect_value)):
            kept[key] = hit
    return sorted(kept.values(), key=lambda hit: (hit.trait, hit.chrom, hit.position, hit.variant_id))


def load_target_variants(path: Path, target_ids: set[str]) -> dict[str, VariantRecord]:
    matches: dict[str, VariantRecord] = {}
    with path.open() as handle:
        for snp_index, line in enumerate(handle):
            fields = line.strip().split()
            if len(fields) < 6:
                continue
            variant_id = fields[0]
            if variant_id not in target_ids or variant_id in matches:
                continue
            matches[variant_id] = VariantRecord(
                variant_id=variant_id,
                chrom=fields[1],
                position=int(fields[3]),
                ref_allele=normalize_allele(fields[4]),
                alt_allele=normalize_allele(fields[5]),
                snp_index=snp_index,
            )
    return matches


def complement_allele(value: str) -> str | None:
    if not value or any(base not in COMPLEMENT for base in value):
        return None
    return "".join(COMPLEMENT[base] for base in value)


def harmonize_hit(hit: Hit, variant: VariantRecord) -> tuple[str, str] | None:
    direct = {hit.effect_allele, hit.other_allele}
    target = {variant.ref_allele, variant.alt_allele}
    if direct == target and hit.effect_allele in target:
        return hit.effect_allele, "direct"
    comp_effect = complement_allele(hit.effect_allele)
    comp_other = complement_allele(hit.other_allele)
    if comp_effect and comp_other and {comp_effect, comp_other} == target:
        return comp_effect, "complement"
    if hit.effect_allele in target:
        return hit.effect_allele, "effect_only"
    if comp_effect and comp_effect in target:
        return comp_effect, "effect_only_complement"
    return None


def build_matched_hits(hits: list[Hit], aadr_matches: dict[str, VariantRecord]) -> list[MatchedHit]:
    matched: list[MatchedHit] = []
    for hit in hits:
        variant = aadr_matches.get(hit.variant_id)
        if variant is None:
            continue
        harmonized = harmonize_hit(hit, variant)
        if harmonized is None:
            continue
        effect_allele_aadr, mode = harmonized
        matched.append(
            MatchedHit(
                hit=hit,
                variant=variant,
                effect_allele_aadr=effect_allele_aadr,
                harmonization=mode,
            )
        )
    kept: dict[tuple[str, str], MatchedHit] = {}
    for match in matched:
        key = (match.hit.trait, match.variant.variant_id)
        prev = kept.get(key)
        if prev is None or match.hit.p_value < prev.hit.p_value:
            kept[key] = match
    return sorted(
        kept.values(),
        key=lambda match: (match.hit.trait, match.variant.chrom, match.variant.position, match.variant.variant_id),
    )


def write_hits(path: Path, hits: list[Hit]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "trait",
                "variant_id",
                "chrom",
                "position",
                "effect_allele",
                "other_allele",
                "effect_value",
                "p_value",
                "source_gwas",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for hit in hits:
            writer.writerow(
                {
                    "trait": hit.trait,
                    "variant_id": hit.variant_id,
                    "chrom": hit.chrom,
                    "position": hit.position,
                    "effect_allele": hit.effect_allele,
                    "other_allele": hit.other_allele,
                    "effect_value": f"{hit.effect_value:.8g}",
                    "p_value": f"{hit.p_value:.8g}",
                    "source_gwas": hit.source_gwas,
                }
            )


def write_match_summary(path: Path, hits: list[Hit], matched_hits: list[MatchedHit], aadr_matches: dict[str, VariantRecord]) -> None:
    by_trait = Counter(hit.trait for hit in hits)
    by_trait_rsid = Counter(hit.trait for hit in hits if hit.variant_id in aadr_matches)
    by_trait_harmonized = Counter(match.hit.trait for match in matched_hits)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["trait", "gwas_hits", "matched_rsid", "allele_harmonized"], delimiter="\t")
        writer.writeheader()
        for trait in sorted(by_trait):
            writer.writerow(
                {
                    "trait": trait,
                    "gwas_hits": by_trait[trait],
                    "matched_rsid": by_trait_rsid[trait],
                    "allele_harmonized": by_trait_harmonized[trait],
                }
            )


def decode_packed_genotypes(raw_bytes: bytes, n_ind: int, sample_indices: list[int]) -> np.ndarray:
    genos = np.zeros(len(sample_indices), dtype=np.int8)
    sample_index_set = {sample_idx: out_idx for out_idx, sample_idx in enumerate(sample_indices)}
    for byte_idx, byte in enumerate(raw_bytes):
        for bit_pos in range(4):
            ind_idx = byte_idx * 4 + bit_pos
            if ind_idx >= n_ind:
                break
            out_idx = sample_index_set.get(ind_idx)
            if out_idx is None:
                continue
            value = (byte >> (6 - 2 * bit_pos)) & 0x03
            genos[out_idx] = 9 if value == 3 else 2 - value
    return genos


def _genotype_pair(genotype: int, ref_allele: str, alt_allele: str) -> tuple[str, str]:
    if genotype == 9:
        return ("0", "0")
    if genotype == 2:
        return (ref_allele, ref_allele)
    if genotype == 1:
        return (ref_allele, alt_allele)
    if genotype == 0:
        return (alt_allele, alt_allele)
    return ("0", "0")


def alt_total_from_genotypes(genotypes: np.ndarray, modern_mask: np.ndarray) -> tuple[float, float]:
    called = genotypes != 9
    if not np.any(called):
        return 0.0, 0.0
    dosage = 2.0 - genotypes.astype(float)
    dosage[~called] = 0.0
    alt = np.where(modern_mask, dosage, dosage / 2.0)
    total = np.where(modern_mask, 2.0, 1.0)
    total[~called] = 0.0
    return float(np.sum(alt)), float(np.sum(total))


def aggregate_variant_counts(
    geno_path: Path,
    n_ind: int,
    union_indices: list[int],
    matched_hits: list[MatchedHit],
    groups_by_region: dict[str, list[GroupData]],
) -> dict[tuple[str, str, str], tuple[float, float]]:
    bytes_per_snp = (n_ind + 3) // 4
    key_order = sorted(
        {(match.variant.variant_id, match.variant.snp_index) for match in matched_hits},
        key=lambda item: item[1],
    )
    counts: dict[tuple[str, str, str], tuple[float, float]] = {}
    with geno_path.open("rb") as handle:
        header_size = len(handle.readline())
        for count, (variant_id, snp_index) in enumerate(key_order, start=1):
            offset = header_size + snp_index * bytes_per_snp
            handle.seek(offset)
            raw = handle.read(bytes_per_snp)
            genotypes = decode_packed_genotypes(raw, n_ind=n_ind, sample_indices=union_indices)
            for region, groups in groups_by_region.items():
                for group in groups:
                    subset = genotypes[group.union_positions]
                    alt_count, total_count = alt_total_from_genotypes(subset, group.modern_mask)
                    counts[(region, group.population, variant_id)] = (alt_count, total_count)
            if count % 250 == 0:
                print(f"Read {count:,}/{len(key_order):,} matched AADR variants", file=sys.stderr)
    return counts


def write_admixture_files(output_dir: Path, groups_by_region: dict[str, list[GroupData]]) -> dict[str, Path]:
    paths: dict[str, Path] = {}
    for region, groups in groups_by_region.items():
        region_label = REGION_SPECS[region]["label"].lower()
        path = output_dir / f"admixture_{region_label}.tsv"
        with path.open("w", newline="") as handle:
            fieldnames = ["population"] + groups[0].component_cols
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for group in groups:
                row = {"population": group.population}
                for col, value in zip(group.component_cols, group.admixture):
                    row[col] = f"{value:.8f}"
                writer.writerow(row)
        paths[region] = path

    joint_path = output_dir / "admixture_joint_non_eur.tsv"
    with joint_path.open("w", newline="") as handle:
        fieldnames = ["population", "k1", "k2", "k3"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for region in ["afr", "eas"]:
            for group in groups_by_region[region]:
                writer.writerow(
                    {
                        "population": f"{REGION_SPECS[region]['label']}__{group.population}",
                        "k1": f"{group.admixture[0]:.8f}",
                        "k2": f"{group.admixture[1]:.8f}",
                        "k3": f"{group.admixture[2]:.8f}",
                    }
                )
        for region in ["sas", "sam"]:
            for group in groups_by_region[region]:
                row = {"population": f"{REGION_SPECS[region]['label']}__{group.population}", "k1": f"{group.admixture[0]:.8f}", "k2": f"{group.admixture[1]:.8f}", "k3": "0.00000000"}
                writer.writerow(row)
    paths["joint_non_eur"] = joint_path
    return paths


def write_counts_files(
    output_dir: Path,
    matched_hits: list[MatchedHit],
    counts: dict[tuple[str, str, str], tuple[float, float]],
    groups_by_region: dict[str, list[GroupData]],
) -> dict[str, Path]:
    counts_by_region: dict[str, list[dict[str, str]]] = defaultdict(list)
    joint_rows: list[dict[str, str]] = []
    for match in matched_hits:
        variant = match.variant
        for region, groups in groups_by_region.items():
            label = REGION_SPECS[region]["label"]
            for group in groups:
                alt_count, total_count = counts[(region, group.population, variant.variant_id)]
                effect_count = alt_count if match.effect_allele_aadr == variant.alt_allele else total_count - alt_count
                row = {
                    "trait": match.hit.trait,
                    "variant_id": match.hit.variant_id,
                    "population": group.population,
                    "effect_allele": match.effect_allele_aadr,
                    "ref_allele": variant.ref_allele,
                    "alt_allele": variant.alt_allele,
                    "effect_allele_count": f"{effect_count:.6f}",
                    "non_missing_allele_count": f"{total_count:.6f}",
                }
                counts_by_region[region].append(row)
                joint_rows.append({**row, "population": f"{label}__{group.population}"})

    paths: dict[str, Path] = {}
    fieldnames = [
        "trait",
        "variant_id",
        "population",
        "effect_allele",
        "ref_allele",
        "alt_allele",
        "effect_allele_count",
        "non_missing_allele_count",
    ]
    for region, rows in counts_by_region.items():
        path = output_dir / f"counts_{REGION_SPECS[region]['label'].lower()}.tsv.gz"
        with gzip.open(path, "wt", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        paths[region] = path
    joint_path = output_dir / "counts_joint_non_eur.tsv.gz"
    with gzip.open(joint_path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(joint_rows)
    paths["joint_non_eur"] = joint_path
    return paths


def write_tfam(path: Path, sample_ids: list[str]) -> None:
    with path.open("w") as handle:
        for iid in sample_ids:
            handle.write(f"{iid}\t{iid}\t0\t0\t0\t-9\n")


def write_tped(
    path: Path,
    geno_path: Path,
    n_ind: int,
    sample_indices: list[int],
    sample_ids: list[str],
    variants: list[VariantRecord],
    plink_ids: dict[int, str],
) -> None:
    del sample_ids
    bytes_per_snp = (n_ind + 3) // 4
    variants = sorted(variants, key=lambda item: (item.chrom, item.position, item.variant_id))
    with geno_path.open("rb") as geno_handle, path.open("w") as handle:
        header_size = len(geno_handle.readline())
        for variant in variants:
            offset = header_size + variant.snp_index * bytes_per_snp
            geno_handle.seek(offset)
            raw = geno_handle.read(bytes_per_snp)
            genotypes = decode_packed_genotypes(raw, n_ind=n_ind, sample_indices=sample_indices)
            row = [variant.chrom, plink_ids[variant.snp_index], "0", str(variant.position)]
            for genotype in genotypes:
                a1, a2 = _genotype_pair(int(genotype), variant.ref_allele, variant.alt_allele)
                row.extend([a1, a2])
            handle.write("\t".join(row) + "\n")


def build_plink_reference(
    output_prefix: Path,
    geno_path: Path,
    n_ind: int,
    sample_indices: list[int],
    sample_ids: list[str],
    variants: list[VariantRecord],
    plink_path: Path,
) -> tuple[Path, dict[int, str]]:
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    plink_ids = {
        variant.snp_index: f"{variant.variant_id}__{variant.chrom}:{variant.position}__{variant.snp_index}"
        for variant in variants
    }
    write_tped(output_prefix.with_suffix(".tped"), geno_path, n_ind, sample_indices, sample_ids, variants, plink_ids)
    write_tfam(output_prefix.with_suffix(".tfam"), sample_ids)
    subprocess.run(
        [
            str(plink_path),
            "--tfile",
            str(output_prefix),
            "--make-bed",
            "--out",
            str(output_prefix),
            "--allow-extra-chr",
            "--silent",
        ],
        check=True,
    )
    return output_prefix, plink_ids


def write_assoc_for_clump(path: Path, matches: list[MatchedHit], plink_ids: dict[int, str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["SNP", "P"], delimiter="\t")
        writer.writeheader()
        for match in sorted(matches, key=lambda item: (item.hit.p_value, item.variant.chrom, item.variant.position)):
            writer.writerow({"SNP": plink_ids[match.variant.snp_index], "P": f"{match.hit.p_value:.8g}"})


def run_plink_clump(
    matches: list[MatchedHit],
    reference_prefix: Path,
    plink_ids: dict[int, str],
    output_prefix: Path,
    clump_kb: int,
    clump_r2: float,
    plink_path: Path,
) -> tuple[set[int], Path]:
    assoc_path = output_prefix.with_suffix(".assoc.tsv")
    write_assoc_for_clump(assoc_path, matches, plink_ids)
    id_to_snp_index = {plink_id: snp_index for snp_index, plink_id in plink_ids.items()}
    subprocess.run(
        [
            str(plink_path),
            "--bfile",
            str(reference_prefix),
            "--clump",
            str(assoc_path),
            "--clump-snp-field",
            "SNP",
            "--clump-field",
            "P",
            "--clump-p1",
            "1",
            "--clump-r2",
            str(clump_r2),
            "--clump-kb",
            str(clump_kb),
            "--allow-extra-chr",
            "--out",
            str(output_prefix),
            "--silent",
        ],
        check=True,
    )
    clumped_path = output_prefix.with_suffix(".clumped")
    selected: set[int] = set()
    with clumped_path.open() as handle:
        header = handle.readline().strip().split()
        index_idx = header.index("SNP")
        for line in handle:
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) <= index_idx:
                continue
            plink_id = fields[index_idx]
            if plink_id in id_to_snp_index:
                selected.add(id_to_snp_index[plink_id])
    return selected, clumped_path


def load_eur_reference_indices(sample_info_path: Path, iid_to_index: dict[str, int], canonical_to_ids: dict[str, list[str]]) -> list[int]:
    indices: list[int] = []
    with sample_info_path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            try:
                if float_value(row["Date"]) != 0.0:
                    continue
            except (KeyError, TypeError, ValueError):
                continue
            ind_index, _ = resolve_aadr_iid(clean_value(row["IID"]), iid_to_index, canonical_to_ids)
            if ind_index is not None:
                indices.append(ind_index)
    unique = sorted(set(indices))
    if not unique:
        raise RuntimeError("No modern European LD-reference samples were resolved.")
    return unique


def run_one_polygenic_test(
    script_path: Path,
    counts_path: Path,
    admixture_path: Path,
    output_path: Path,
    traits: list[str],
    n_permutations: int,
    block_size: int,
    seed: int,
    exclude_populations: list[str] | None = None,
) -> None:
    cmd = [
        sys.executable,
        str(script_path),
        "--counts",
        str(counts_path),
        "--admixture",
        str(admixture_path),
        "--traits",
        *traits,
        "--n-permutations",
        str(n_permutations),
        "--block-size",
        str(block_size),
        "--seed",
        str(seed),
        "--output",
        str(output_path),
    ]
    if exclude_populations:
        cmd.extend(["--exclude-populations", *exclude_populations])
    subprocess.run(cmd, check=True)


def run_and_collect_rows(
    analysis: str,
    output_dir: Path,
    counts_paths: dict[str, Path],
    admixture_paths: dict[str, Path],
    traits: list[str],
    script_path: Path,
    n_permutations: int,
    block_size: int,
    seed: int,
) -> tuple[str, list[dict[str, str]]]:
    out_path = output_dir / f"polygenic_{analysis}.tsv"
    run_one_polygenic_test(
        script_path,
        counts_paths[analysis],
        admixture_paths[analysis],
        out_path,
        traits,
        n_permutations,
        block_size,
        seed,
    )
    return analysis, load_rows(out_path, analysis)


def load_rows(path: Path, analysis: str) -> list[dict[str, str]]:
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = []
        for row in reader:
            row["analysis"] = analysis
            rows.append(row)
        return rows


def write_summary(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as handle:
        fieldnames = ["analysis"] + [name for name in rows[0] if name != "analysis"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def run_sign_suite(
    output_dir: Path,
    counts_paths: dict[str, Path],
    admixture_paths: dict[str, Path],
    traits: list[str],
    script_path: Path,
    n_permutations: int,
    block_size: int,
    seed: int,
    jobs: int = 1,
) -> Path:
    rows_by_analysis: dict[str, list[dict[str, str]]] = {}
    if jobs <= 1:
        for analysis in ALL_ANALYSES:
            _, analysis_rows = run_and_collect_rows(
                analysis,
                output_dir,
                counts_paths,
                admixture_paths,
                traits,
                script_path,
                n_permutations,
                block_size,
                seed,
            )
            rows_by_analysis[analysis] = analysis_rows
    else:
        with ThreadPoolExecutor(max_workers=min(jobs, len(ALL_ANALYSES))) as executor:
            futures = [
                executor.submit(
                    run_and_collect_rows,
                    analysis,
                    output_dir,
                    counts_paths,
                    admixture_paths,
                    traits,
                    script_path,
                    n_permutations,
                    block_size,
                    seed,
                )
                for analysis in ALL_ANALYSES
            ]
            for future in futures:
                analysis, analysis_rows = future.result()
                rows_by_analysis[analysis] = analysis_rows
    rows: list[dict[str, str]] = []
    for analysis in ALL_ANALYSES:
        rows.extend(rows_by_analysis[analysis])
    summary_path = output_dir / "polygenic_summary.tsv"
    write_summary(summary_path, rows)
    return summary_path


def run_joint_sensitivity(
    output_dir: Path,
    counts_paths: dict[str, Path],
    admixture_paths: dict[str, Path],
    traits: list[str],
    script_path: Path,
    n_permutations: int,
    block_size: int,
    seed: int,
    jobs: int = 1,
) -> None:
    del jobs
    out_path = output_dir / "polygenic_joint_non_eur_sensitivity.tsv"
    run_one_polygenic_test(
        script_path,
        counts_paths["joint_non_eur"],
        admixture_paths["joint_non_eur"],
        out_path,
        traits,
        n_permutations,
        block_size,
        seed,
        exclude_populations=["SAM__IBS.SG", "SAM__BA-S"],
    )


def write_gfp_notes(path: Path, diagnostics: dict[str, object], clumped_hits: list[MatchedHit]) -> None:
    lines = [
        "# GFP Construction",
        "",
        "This repository uses a literature-weighted GFP construction.",
        "",
        "Following the phenotypic GFP literature, neuroticism is reversed before combination.",
        "Following Colbran et al.'s downstream polygenic-selection pipeline, only 1240k variants passing p < 1e-8 are retained before PLINK clumping.",
        "",
        f"- Shared AADR-intersecting Big Five rsIDs: `{diagnostics['shared_rsids']:,}`",
        f"- Harmonized shared Big Five rsIDs: `{diagnostics['harmonized_rows']:,}`",
        f"- Alignment failures dropped: `{diagnostics['skipped_alignment']:,}`",
        f"- Stouffer denominator using empirical Big Five GWAS correlation: `{float(diagnostics['stouffer_denom']):.6f}`",
        f"- Pre-clump GFP-significant variants: `{sum(int(v) for v in diagnostics['support_hist'].values()):,}`",
        f"- Post-clump GFP lead loci: `{sum(1 for hit in clumped_hits if hit.hit.trait == 'gfp'):,}`",
        "",
        "## Big Five Loadings",
        "",
    ]
    for trait in TRAIT_ORDER:
        lines.append(f"- `{trait}`: `{float(diagnostics['loadings'][trait]):.2f}`")
    lines.extend(["", "## Big Five support among GFP-significant variants", ""])
    for support, count in diagnostics["support_hist"].items():
        lines.append(f"- `{support}/5` concordant traits: `{count:,}`")
    lines.extend(["", "## Standalone Big Five genome-wide-significant counts on the shared SNP set", ""])
    for trait in TRAIT_ORDER:
        lines.append(f"- `{trait}`: `{int(diagnostics['per_trait_gws'][trait]):,}`")
    path.write_text("\n".join(lines) + "\n")


def load_summary_lookup(path: Path) -> dict[tuple[str, str], dict[str, str]]:
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {(row["analysis"], row["trait"]): row for row in reader}


def load_trait_lookup(path: Path) -> dict[str, dict[str, str]]:
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {row["trait"]: row for row in reader}


def write_results_md(
    path: Path,
    candidate_hits: list[Hit],
    clumped_hits: list[MatchedHit],
    summary_path: Path,
    sensitivity_path: Path,
) -> None:
    candidate_counts = Counter(hit.trait for hit in candidate_hits)
    clumped_counts = Counter(match.hit.trait for match in clumped_hits)
    summary = load_summary_lookup(summary_path)
    sensitivity = load_trait_lookup(sensitivity_path)
    lines = [
        "# Extension Results",
        "",
        "This run extends Colbran et al. with two traits often placed near the self-domestication literature:",
        "",
        "- `gfp`: a literature-weighted General Factor of Personality derived from Big Five GWAS",
        "- `schizophrenia`: European schizophrenia summary statistics processed under the same downstream pipeline",
        "",
        "## Lead-locus counts",
        "",
        "| Trait | Candidate 1240k hits pre-clump | Clumped lead loci |",
        "| --- | ---: | ---: |",
    ]
    for trait in ["gfp", "schizophrenia"]:
        lines.append(f"| `{trait}` | `{candidate_counts[trait]:,}` | `{clumped_counts[trait]:,}` |")
    lines.extend(
        [
            "",
            "## Headline sign-only results",
            "",
            f"- `gfp joint_non_eur = {summary[('joint_non_eur', 'gfp')]['directional_p']}`",
            f"- `schizophrenia joint_non_eur = {summary[('joint_non_eur', 'schizophrenia')]['directional_p']}`",
            f"- `gfp joint_non_eur_sensitivity = {sensitivity['gfp']['directional_p']}`",
            f"- `schizophrenia joint_non_eur_sensitivity = {sensitivity['schizophrenia']['directional_p']}`",
            "",
            "## Output",
            "",
            f"- Main summary: `{summary_path.name}`",
            "- Figure 6-style overlay can be generated with `scripts/plot_figure6_comparison.py`.",
        ]
    )
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.big5_cache_dir.mkdir(parents=True, exist_ok=True)

    aadr_variant_ids = load_aadr_variant_ids(args.aadr_snp)
    big5_rows = ensure_big5_cache(args.big5_cache_dir, args.aadr_snp)
    gfp_hits, gfp_diagnostics = construct_gfp_hits(big5_rows, args.gwas_p_threshold)
    scz_hits = load_scz_full_hits(args.schizophrenia_sumstats, aadr_variant_ids, args.gwas_p_threshold)
    candidate_hits = deduplicate_hits(gfp_hits + scz_hits)
    if not candidate_hits:
        raise RuntimeError("No candidate GFP or schizophrenia hits survived the genome-wide threshold.")
    write_hits(args.output_dir / "candidate_lead_variants.tsv", candidate_hits)

    sample_ids, iid_to_index, canonical_to_ids = load_aadr_ind(args.aadr_ind)
    groups_by_region, union_indices = parse_region_groups(args.sample_info_dir, iid_to_index, canonical_to_ids)
    aadr_matches = load_target_variants(args.aadr_snp, {hit.variant_id for hit in candidate_hits})
    matched_preclump = build_matched_hits(candidate_hits, aadr_matches)
    if not matched_preclump:
        raise RuntimeError("No candidate hits matched to AADR after allele harmonization.")
    write_match_summary(args.output_dir / "candidate_match_summary.tsv", candidate_hits, matched_preclump, aadr_matches)

    eur_indices = load_eur_reference_indices(args.sample_info_dir / "eur_sample_info.txt", iid_to_index, canonical_to_ids)
    eur_variants = sorted({match.variant for match in matched_preclump}, key=lambda variant: variant.snp_index)
    ref_prefix, plink_ids = build_plink_reference(
        args.output_dir / "eur_reference",
        args.aadr_geno,
        len(sample_ids),
        eur_indices,
        [sample_ids[idx] for idx in eur_indices],
        eur_variants,
        args.plink_path,
    )
    selected_snp_indices, clumped_path = run_plink_clump(
        matched_preclump,
        ref_prefix,
        plink_ids,
        args.output_dir / "plink_clump",
        args.clump_kb,
        args.clump_r2,
        args.plink_path,
    )
    clumped_matches = [match for match in matched_preclump if match.variant.snp_index in selected_snp_indices]
    if not clumped_matches:
        raise RuntimeError("PLINK clumping removed every matched hit.")

    write_hits(args.output_dir / "lead_variants.tsv", [match.hit for match in clumped_matches])
    write_match_summary(args.output_dir / "lead_variant_match_summary.tsv", [match.hit for match in clumped_matches], clumped_matches, aadr_matches)
    (args.output_dir / "plink_clump.note.txt").write_text(f"PLINK clump output: {clumped_path}\n")

    counts = aggregate_variant_counts(args.aadr_geno, len(sample_ids), union_indices, clumped_matches, groups_by_region)
    admixture_paths = write_admixture_files(args.output_dir, groups_by_region)
    counts_paths = write_counts_files(args.output_dir, clumped_matches, counts, groups_by_region)
    summary_path = run_sign_suite(
        args.output_dir,
        counts_paths,
        admixture_paths,
        ["gfp", "schizophrenia"],
        args.polygenic_script,
        args.n_permutations,
        args.block_size,
        args.seed,
        args.jobs,
    )
    run_joint_sensitivity(
        args.output_dir,
        counts_paths,
        admixture_paths,
        ["gfp", "schizophrenia"],
        args.polygenic_script,
        args.n_permutations,
        args.block_size,
        args.seed,
        args.jobs,
    )
    write_gfp_notes(args.output_dir / "gfp_construction.md", gfp_diagnostics, clumped_matches)
    write_results_md(
        args.output_dir / "RESULTS.md",
        candidate_hits,
        clumped_matches,
        summary_path,
        args.output_dir / "polygenic_joint_non_eur_sensitivity.tsv",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
