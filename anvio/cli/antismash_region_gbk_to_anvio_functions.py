#!/usr/bin/env python3
"""
Convert antiSMASH per-region GenBank files into an anvi'o functions-txt table.
"""

from __future__ import annotations

import os
import csv
import re
import sys
if __package__ is None or __package__ == "":
    script_dir = os.path.dirname(__file__)
    parent_dir = os.path.dirname(os.path.dirname(script_dir))
    if sys.path and sys.path[0] == script_dir:
        sys.path.pop(0)
    if parent_dir not in sys.path:
        sys.path.insert(0, parent_dir)

from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import anvio
import anvio.terminal as terminal
from anvio.argparse import ArgumentParser
from anvio.errors import ConfigError, FilesNPathsError


__authors__ = ["kpf734"]
__credits__ = []
__license__ = "GPL 3.0"
__description__ = "Convert antiSMASH per-region GenBank files into an anvi'o functions-txt table"
__requires__ = ["antismash-region-gbk"]
__provides__ = ["functions-txt"]
__anvio_workflows__ = ["contigs"]


run = terminal.Run()
progress = terminal.Progress()


FEATURE_RE = re.compile(r"^     (\S+)\s+(.+)")
QUALIFIER_RE = re.compile(r"^                     /([^=\s]+)(?:=(.*))?$")
CONTINUATION_RE = re.compile(r"^                     (.*)$")
ANVIO_REGION_GBK_RE = re.compile(r"contig\d+\.region\d+\.gbk$", re.IGNORECASE)
ANVIO_LOCUS_TAG_RE = re.compile(r"^anvio_gene_(?P<gene_callers_id>\d+)$")


@dataclass
class Feature:
    kind: str
    qualifiers: dict[str, list[str]] = field(default_factory=lambda: defaultdict(list))


@dataclass(frozen=True)
class RegionContext:
    region_id: str
    region_accession: str
    region_type: str


@dataclass(frozen=True)
class AntiSMASHGene:
    gene_callers_id: str
    gene_kind: str
    region: RegionContext


def parse_args() -> object:
    parser = ArgumentParser(description=__description__)
    parser.add_argument(
        "-i",
        "--antismash-dir",
        required=True,
        type=Path,
        help="antiSMASH output directory containing contig<number>.region<number>.gbk files.",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        type=Path,
        help="Output functions-txt TSV path for anvi-import-functions.",
    )
    return parser.get_args(parser)


def clean_qualifier(parts: Iterable[str]) -> str:
    value = " ".join(part.strip() for part in parts if part.strip()).strip()
    if value.startswith('"'):
        value = value[1:]
    if value.endswith('"'):
        value = value[:-1]
    return value.replace('""', '"')


def parse_genbank(path: Path) -> list[Feature]:
    features: list[Feature] = []
    current_feature: Feature | None = None
    current_qualifier: str | None = None
    current_value_parts: list[str] = []
    in_features = False

    def finish_qualifier() -> None:
        nonlocal current_qualifier, current_value_parts
        if current_feature is not None and current_qualifier is not None:
            current_feature.qualifiers[current_qualifier].append(clean_qualifier(current_value_parts))
        current_qualifier = None
        current_value_parts = []

    def finish_feature() -> None:
        nonlocal current_feature
        finish_qualifier()
        if current_feature is not None:
            features.append(current_feature)
        current_feature = None

    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")

            if line.startswith("FEATURES"):
                in_features = True
                continue

            if not in_features:
                continue

            if line.startswith("ORIGIN") or line.startswith("//"):
                finish_feature()
                break

            feature_match = FEATURE_RE.match(line)
            if feature_match:
                finish_feature()
                current_feature = Feature(feature_match.group(1))
                continue

            if current_feature is None:
                continue

            qualifier_match = QUALIFIER_RE.match(line)
            if qualifier_match:
                finish_qualifier()
                current_qualifier = qualifier_match.group(1)
                current_value_parts = [qualifier_match.group(2) or "true"]
                continue

            continuation_match = CONTINUATION_RE.match(line)
            if continuation_match:
                continuation = continuation_match.group(1).strip()
                if current_qualifier is not None:
                    current_value_parts.append(continuation)

    return features


def one(qualifiers: dict[str, list[str]], key: str) -> str:
    values = qualifiers.get(key) or []
    return values[0] if values else ""


def joined(qualifiers: dict[str, list[str]], key: str) -> str:
    return "; ".join(qualifiers.get(key) or [])


def region_accession(region_id: str, contig_edge: str) -> str:
    edge = contig_edge.strip().upper()
    if edge not in {"TRUE", "FALSE"}:
        edge = "UNKNOWN"
    return f"{region_id}_is_edge_{edge}"


def gene_callers_id_from_locus_tag(locus_tag: str) -> str | None:
    match = ANVIO_LOCUS_TAG_RE.match(locus_tag.strip())
    if match is None:
        return None
    return match.group("gene_callers_id")


def discover_region_gbk_files(antismash_dir: Path) -> list[Path]:
    if not antismash_dir.is_dir():
        raise FilesNPathsError(f"Not a directory: {antismash_dir}")

    files = sorted(
        path
        for path in antismash_dir.glob("*.gbk")
        if path.is_file()
        and not path.name.startswith(".")
        and ANVIO_REGION_GBK_RE.search(path.name)
    )
    if not files:
        raise FilesNPathsError(
            "No antiSMASH per-region GenBank files found in "
            f"{antismash_dir}. Expected names ending in contig<number>.region<number>.gbk"
        )
    return files


def extract_region_context(path: Path, features: list[Feature]) -> RegionContext:
    region_features = [feature for feature in features if feature.kind == "region"]
    if not region_features:
        raise ConfigError(f"No antiSMASH region feature found in: {path}")
    if len(region_features) > 1:
        raise ConfigError(f"Expected one region feature in {path}, found {len(region_features)}")

    region = region_features[0]
    region_id = path.stem
    return RegionContext(
        region_id=region_id,
        region_accession=region_accession(region_id, one(region.qualifiers, "contig_edge")),
        region_type=joined(region.qualifiers, "product"),
    )


def extract_antismash_genes(path: Path) -> tuple[list[AntiSMASHGene], int, int]:
    features = parse_genbank(path)
    region = extract_region_context(path, features)
    genes: list[AntiSMASHGene] = []
    cds_with_gene_kind = 0
    skipped_locus_tags = 0

    for feature in features:
        if feature.kind != "CDS":
            continue

        gene_kind = joined(feature.qualifiers, "gene_kind")
        if not gene_kind:
            continue

        cds_with_gene_kind += 1
        gene_callers_id = gene_callers_id_from_locus_tag(one(feature.qualifiers, "locus_tag"))
        if gene_callers_id is None:
            skipped_locus_tags += 1
            continue

        genes.append(
            AntiSMASHGene(
                gene_callers_id=gene_callers_id,
                gene_kind=gene_kind,
                region=region,
            )
        )

    return genes, cds_with_gene_kind, skipped_locus_tags


def function_rows(genes: list[AntiSMASHGene]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for gene in genes:
        rows.append(
            {
                "gene_callers_id": gene.gene_callers_id,
                "source": "antiSMASH_type",
                "accession": gene.region.region_accession,
                "function": gene.region.region_type,
                "e_value": "",
            }
        )
        rows.append(
            {
                "gene_callers_id": gene.gene_callers_id,
                "source": "antiSMASH_function",
                "accession": gene.region.region_accession,
                "function": gene.gene_kind,
                "e_value": "",
            }
        )
    return rows


def write_functions(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["gene_callers_id", "source", "accession", "function", "e_value"]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


@terminal.time_program
def main() -> None:
    args = parse_args()

    try:
        gbk_files = discover_region_gbk_files(args.antismash_dir)
        genes: list[AntiSMASHGene] = []
        cds_with_gene_kind = 0
        skipped_locus_tags = 0

        progress.new("Reading antiSMASH region GenBank files", progress_total_items=len(gbk_files))
        for gbk_file in gbk_files:
            file_genes, file_cds_with_gene_kind, file_skipped_locus_tags = extract_antismash_genes(gbk_file)
            genes.extend(file_genes)
            cds_with_gene_kind += file_cds_with_gene_kind
            skipped_locus_tags += file_skipped_locus_tags
            progress.increment()
        progress.end()

        rows = function_rows(genes)
        if not rows:
            raise ConfigError(
                "No importable antiSMASH annotations found. The parser only uses CDS features "
                "with /gene_kind and /locus_tag values exactly like anvio_gene_<integer>."
            )

        write_functions(args.output, rows)

        run.info_single(f"Read {len(gbk_files)} antiSMASH region GenBank files")
        run.info_single(f"Found {cds_with_gene_kind} CDS features with gene_kind")
        run.info_single(f"Used {len(genes)} CDS features with anvio_gene_<integer> locus_tag")
        if skipped_locus_tags:
            run.warning(
                f"Skipped {skipped_locus_tags} CDS features whose locus_tag was not "
                "anvio_gene_<integer>"
            )
        run.info_single(f"Wrote {len(rows)} anvi'o function rows: {args.output}")
    except (ConfigError, FilesNPathsError) as error:
        print(error, file=sys.stderr)
        sys.exit(-1 if isinstance(error, ConfigError) else -2)


if __name__ == "__main__":
    main()
