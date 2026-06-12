#!/usr/bin/env python3
"""
Extract antiSMASH regions and gene annotations for anvi'o contigs-db workflows.

Primary outputs:
  PREFIX.regions.tsv   antiSMASH region coordinates and products/types
  PREFIX.genes.tsv     CDS features in antiSMASH regions, with gene_kind when set

Additional output:
  PREFIX.functions.tsv five-column anvi-import-functions-style annotations

The parser supports antiSMASH JSON and anvi'o-derived per-region GenBank files
named like Day17a_QCcontig235.region001.gbk. Coordinates in the TSV outputs are
0-based, end-exclusive unless explicitly named start_1based/end_1based.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable


LOCATION_RE = re.compile(r"\[<?(\d+):>?(\d+)\]\(([+\-?])\)")
COORD_RE = re.compile(r"\[<?(\d+):>?(\d+)\]")
GBK_FEATURE_RE = re.compile(r"^     (\S+)\s+(.+)")
GBK_QUALIFIER_RE = re.compile(r"^                     /([^=\s]+)(?:=(.*))?$")
GBK_CONTINUATION_RE = re.compile(r"^                     (.*)$")
GBK_COORD_RE = re.compile(r"<?(\d+)(?:\.\.>?(\d+))?")
ORIG_START_RE = re.compile(r"Orig\. start\s*::\s*(\d+)")
ORIG_END_RE = re.compile(r"Orig\. end\s*::\s*(\d+)")
ANVIO_REGION_GBK_RE = re.compile(r"contig\d+\.region\d+\.gbk$", re.IGNORECASE)
ANVIO_LOCUS_TAG_RE = re.compile(r"^anvio_gene_(?P<gene_callers_id>\d+)$")


@dataclass
class Region:
    source_json: str
    record_id: str
    contig: str
    region_number: str
    start: int
    stop: int
    strand: str
    region_type: str
    category: str = ""
    candidate_cluster_numbers: str = ""
    protocluster_numbers: str = ""
    core_locations: str = ""
    contig_edge: str = ""
    rules: str = ""


@dataclass
class Gene:
    source_json: str
    record_id: str
    contig: str
    start: int
    stop: int
    strand: str
    locus_tag: str = ""
    protein_id: str = ""
    product: str = ""
    gene_kind: str = ""
    gene_functions: str = ""
    sec_met_domain: str = ""
    region_numbers: str = ""
    region_types: str = ""
    region_categories: str = ""
    region_locations: str = ""
    region_contig_edges: str = ""
    gene_callers_id: str = ""
    gene_call_match: str = ""


@dataclass
class GeneCall:
    gene_callers_id: str
    contig: str
    start: int
    stop: int
    direction: str


@dataclass
class GeneCallIndex:
    exact: dict[tuple[str, int, int, str], list[GeneCall]] = field(default_factory=lambda: defaultdict(list))
    by_contig: dict[str, list[GeneCall]] = field(default_factory=lambda: defaultdict(list))


@dataclass
class GBKFeature:
    type: str
    location: str
    qualifiers: dict[str, list[str]] = field(default_factory=lambda: defaultdict(list))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Extract gene locations, antiSMASH gene_kind, region locations, and "
            "region product/type annotations from antiSMASH GenBank or JSON output."
        )
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        type=Path,
        help=(
            "antiSMASH output directory, antiSMASH JSON file, or per-region GenBank file. "
            "For directories, files ending in contig<number>.region<number>.gbk are preferred."
        ),
    )
    parser.add_argument(
        "-o",
        "--output-prefix",
        required=True,
        type=Path,
        help="Output prefix. Files will be written as PREFIX.regions.tsv, PREFIX.genes.tsv, etc.",
    )
    parser.add_argument(
        "--input-format",
        choices=["auto", "gbk", "json"],
        default="auto",
        help="Input format. Default: auto.",
    )
    parser.add_argument(
        "--gene-calls",
        type=Path,
        help=(
            "Optional anvi'o gene calls file from anvi-export-gene-calls. "
            "When provided, PREFIX.functions.tsv uses mapped anvi'o gene_callers_id values."
        ),
    )
    parser.add_argument(
        "--id-source",
        choices=["auto", "gene_callers_id", "locus_tag", "protein_id"],
        default="auto",
        help=(
            "Identifier to place in the gene_callers_id column of PREFIX.functions.tsv. "
            "auto uses gene_callers_id when --gene-calls is supplied, otherwise the "
            "integer parsed from locus_tag values like anvio_gene_707. "
            "Only gene_callers_id is directly importable by anvi'o."
        ),
    )
    parser.add_argument(
        "--genes",
        choices=["all-region-cds", "kinded"],
        default="all-region-cds",
        help=(
            "Which CDS features to emit. all-region-cds includes every CDS overlapping "
            "an antiSMASH region; kinded includes only CDS features with gene_kind. "
            "Default: all-region-cds."
        ),
    )
    parser.add_argument(
        "--min-reciprocal-overlap",
        type=float,
        default=0.95,
        help=(
            "Fallback coordinate-match threshold for mapping antiSMASH CDS to anvi'o "
            "gene_callers_id when exact coordinates do not match. Both genes must cover "
            "at least this fraction of each other. Default: 0.95."
        ),
    )
    parser.add_argument(
        "--write-nucleotide-misc",
        action="store_true",
        help=(
            "Also write PREFIX.nucleotides.tsv for anvi-import-misc-data -t nucleotides. "
            "This emits one row per nucleotide in each antiSMASH region and can be large."
        ),
    )
    parser.add_argument(
        "--skip-region-gene-tables",
        action="store_true",
        help="Only write PREFIX.functions.tsv, plus nucleotide misc-data if requested.",
    )
    return parser.parse_args()


def discover_input_files(path: Path, input_format: str) -> tuple[str, list[Path]]:
    if path.is_file():
        if input_format == "auto":
            suffix = path.suffix.lower()
            if suffix in {".gbk", ".gbff", ".gb"}:
                if not ANVIO_REGION_GBK_RE.search(path.name):
                    raise ValueError(
                        "GenBank input must be an anvi'o antiSMASH per-region file "
                        "with a name ending in contig<number>.region<number>.gbk"
                    )
                return "gbk", [path]
            return "json", [path]
        if input_format == "gbk" and not ANVIO_REGION_GBK_RE.search(path.name):
            raise ValueError(
                "GenBank input must be an anvi'o antiSMASH per-region file "
                "with a name ending in contig<number>.region<number>.gbk"
            )
        return input_format, [path]

    if not path.is_dir():
        raise FileNotFoundError(f"Input path does not exist: {path}")

    if input_format in {"auto", "gbk"}:
        gbk_files = sorted(
            p
            for p in path.glob("*.region*.gbk")
            if p.is_file()
            and not p.name.startswith(".")
            and ANVIO_REGION_GBK_RE.search(p.name)
        )
        if gbk_files:
            return "gbk", gbk_files
        if input_format == "gbk":
            raise FileNotFoundError(
                "No anvi'o antiSMASH per-region GenBank files found in "
                f"{path}. Expected names ending in contig<number>.region<number>.gbk"
            )

    json_files = sorted(
        p
        for p in path.glob("*.json")
        if p.is_file() and not p.name.startswith(".")
    )
    if not json_files:
        raise FileNotFoundError(f"No antiSMASH JSON or per-region GenBank files found in: {path}")
    return "json", json_files


def as_list(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, list):
        return [str(v) for v in value if v is not None]
    return [str(value)]


def q_join(qualifiers: dict[str, Any], key: str, sep: str = "; ") -> str:
    return sep.join(as_list(qualifiers.get(key)))


def q_one(qualifiers: dict[str, Any], key: str, default: str = "") -> str:
    values = as_list(qualifiers.get(key))
    return values[0] if values else default


def uniq_join(values: Iterable[str], sep: str = "; ") -> str:
    seen: set[str] = set()
    ordered: list[str] = []
    for value in values:
        if value and value not in seen:
            ordered.append(value)
            seen.add(value)
    return sep.join(ordered)


def parse_location(location: str) -> tuple[int, int, str]:
    """Return 0-based start, end-exclusive stop, and +|-|? strand."""
    parts = [
        (int(match.group(1)), int(match.group(2)), match.group(3))
        for match in LOCATION_RE.finditer(location)
    ]

    if not parts:
        parts = [
            (int(match.group(1)), int(match.group(2)), "?")
            for match in COORD_RE.finditer(location)
        ]

    if not parts:
        raise ValueError(f"Could not parse antiSMASH location string: {location!r}")

    start = min(part[0] for part in parts)
    stop = max(part[1] for part in parts)
    strands = {part[2] for part in parts}
    strand = strands.pop() if len(strands) == 1 else "?"
    return start, stop, strand


def parse_gbk_location(location: str) -> tuple[int, int, str]:
    """Return 0-based start, end-exclusive stop, and +|-|? strand from GenBank."""
    coords: list[tuple[int, int]] = []
    for match in GBK_COORD_RE.finditer(location):
        start_1based = int(match.group(1))
        end_1based = int(match.group(2) or match.group(1))
        coords.append((start_1based - 1, end_1based))

    if not coords:
        raise ValueError(f"Could not parse GenBank location string: {location!r}")

    strand = "-" if "complement(" in location else "+"
    return min(start for start, _ in coords), max(stop for _, stop in coords), strand


def offset_core_location(core_location: str, offset: int) -> str:
    """Convert an antiSMASH-style local core location to original coordinates."""
    def replace(match: re.Match[str]) -> str:
        start = int(match.group(1)) + offset
        stop = int(match.group(2)) + offset
        return f"[{start}:{stop}]"

    return re.sub(r"\[(\d+):(\d+)\]", replace, core_location)


def clean_gbk_qualifier_value(value_parts: list[str]) -> str:
    value = " ".join(part.strip() for part in value_parts if part.strip()).strip()
    if value.startswith('"'):
        value = value[1:]
    if value.endswith('"'):
        value = value[:-1]
    return value.replace('""', '"')


def parse_gbk_file(path: Path) -> tuple[str, int | None, int | None, list[GBKFeature]]:
    record_id = ""
    orig_start: int | None = None
    orig_end: int | None = None
    features: list[GBKFeature] = []
    current_feature: GBKFeature | None = None
    current_key: str | None = None
    current_value_parts: list[str] = []
    in_features = False

    def finish_qualifier() -> None:
        nonlocal current_key, current_value_parts
        if current_feature is None or current_key is None:
            return
        current_feature.qualifiers[current_key].append(clean_gbk_qualifier_value(current_value_parts))
        current_key = None
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

            if line.startswith("ACCESSION"):
                parts = line.split()
                if len(parts) > 1:
                    record_id = parts[1]

            start_match = ORIG_START_RE.search(line)
            if start_match:
                orig_start = int(start_match.group(1))

            end_match = ORIG_END_RE.search(line)
            if end_match:
                orig_end = int(end_match.group(1))

            if line.startswith("FEATURES"):
                in_features = True
                continue

            if not in_features:
                continue

            if line.startswith("ORIGIN") or line.startswith("//"):
                finish_feature()
                break

            feature_match = GBK_FEATURE_RE.match(line)
            if feature_match:
                finish_feature()
                current_feature = GBKFeature(feature_match.group(1), feature_match.group(2).strip())
                continue

            if current_feature is None:
                continue

            qualifier_match = GBK_QUALIFIER_RE.match(line)
            if qualifier_match:
                finish_qualifier()
                current_key = qualifier_match.group(1)
                current_value_parts = [qualifier_match.group(2) or "true"]
                continue

            continuation_match = GBK_CONTINUATION_RE.match(line)
            if continuation_match:
                continuation = continuation_match.group(1).strip()
                if current_key is not None:
                    current_value_parts.append(continuation)
                elif continuation:
                    current_feature.location += continuation

    if not record_id:
        record_id = path.stem.split(".region", 1)[0]

    return record_id, orig_start, orig_end, features


def direction_from_strand(strand: str) -> str:
    if strand == "+":
        return "f"
    if strand == "-":
        return "r"
    return "?"


def format_interval(contig: str, start: int, stop: int, strand: str) -> str:
    return f"{contig}:{start}-{stop}({strand})"


def overlaps(a_start: int, a_stop: int, b_start: int, b_stop: int) -> bool:
    return max(a_start, b_start) < min(a_stop, b_stop)


def contained_or_overlapping_regions(gene: Gene, regions: list[Region]) -> list[Region]:
    return [
        region
        for region in regions
        if gene.contig == region.contig and overlaps(gene.start, gene.stop, region.start, region.stop)
    ]


def sanitize_accession(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.:-]+", "_", value.strip())
    cleaned = cleaned.strip("_")
    return cleaned[:200] if cleaned else "unknown"


def gene_callers_id_from_locus_tag(locus_tag: str) -> str | None:
    match = ANVIO_LOCUS_TAG_RE.match(locus_tag.strip())
    if not match:
        return None
    return match.group("gene_callers_id")


def apply_function_id_source(genes: list[Gene], id_source: str) -> None:
    if id_source == "gene_callers_id":
        return

    for gene in genes:
        if id_source == "locus_tag":
            replacement = gene_callers_id_from_locus_tag(gene.locus_tag)
        else:
            replacement = getattr(gene, id_source)
        if replacement:
            gene.gene_callers_id = replacement
            gene.gene_call_match = id_source


def load_antismash_json(path: Path) -> tuple[list[Region], list[Gene]]:
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)

    regions: list[Region] = []
    genes: list[Gene] = []

    for record in payload.get("records", []):
        record_id = str(record.get("id") or record.get("name") or "")
        contig = record_id
        features = record.get("features", [])

        protoclusters: list[dict[str, Any]] = []
        for feature in features:
            if feature.get("type") != "protocluster":
                continue
            start, stop, strand = parse_location(feature.get("location", ""))
            qualifiers = feature.get("qualifiers", {})
            protoclusters.append(
                {
                    "start": start,
                    "stop": stop,
                    "strand": strand,
                    "number": q_one(qualifiers, "protocluster_number"),
                    "category": q_one(qualifiers, "category"),
                    "core_location": q_one(qualifiers, "core_location"),
                }
            )

        for feature in features:
            if feature.get("type") != "region":
                continue
            start, stop, strand = parse_location(feature.get("location", ""))
            qualifiers = feature.get("qualifiers", {})
            overlapping_protoclusters = [
                pc for pc in protoclusters if overlaps(start, stop, pc["start"], pc["stop"])
            ]
            region = Region(
                source_json=str(path),
                record_id=record_id,
                contig=contig,
                region_number=q_one(qualifiers, "region_number", str(len(regions) + 1)),
                start=start,
                stop=stop,
                strand=strand,
                region_type=q_join(qualifiers, "product"),
                category=uniq_join(pc["category"] for pc in overlapping_protoclusters),
                candidate_cluster_numbers=q_join(qualifiers, "candidate_cluster_numbers"),
                protocluster_numbers=uniq_join(pc["number"] for pc in overlapping_protoclusters),
                core_locations=uniq_join(pc["core_location"] for pc in overlapping_protoclusters),
                contig_edge=q_one(qualifiers, "contig_edge"),
                rules=q_join(qualifiers, "rules"),
            )
            regions.append(region)

        record_regions = [region for region in regions if region.record_id == record_id]
        for feature in features:
            if feature.get("type") != "CDS":
                continue

            start, stop, strand = parse_location(feature.get("location", ""))
            qualifiers = feature.get("qualifiers", {})
            gene = Gene(
                source_json=str(path),
                record_id=record_id,
                contig=contig,
                start=start,
                stop=stop,
                strand=strand,
                locus_tag=q_one(qualifiers, "locus_tag") or q_one(qualifiers, "gene"),
                protein_id=q_one(qualifiers, "protein_id"),
                product=q_join(qualifiers, "product"),
                gene_kind=q_join(qualifiers, "gene_kind"),
                gene_functions=q_join(qualifiers, "gene_functions"),
                sec_met_domain=q_join(qualifiers, "sec_met_domain"),
            )

            gene_regions = contained_or_overlapping_regions(gene, record_regions)
            if not gene_regions:
                continue

            gene.region_numbers = uniq_join(region.region_number for region in gene_regions)
            gene.region_types = uniq_join(region.region_type for region in gene_regions)
            gene.region_categories = uniq_join(region.category for region in gene_regions)
            gene.region_locations = uniq_join(
                format_interval(region.contig, region.start, region.stop, region.strand)
                for region in gene_regions
            )
            gene.region_contig_edges = uniq_join(region.contig_edge for region in gene_regions)
            genes.append(gene)

    return regions, genes


def load_antismash_gbk(path: Path) -> tuple[list[Region], list[Gene]]:
    record_id, orig_start, _orig_end, features = parse_gbk_file(path)
    if orig_start is None:
        raise ValueError(
            f"{path} does not look like a per-region antiSMASH GenBank file: missing 'Orig. start'"
        )

    contig = record_id
    regions: list[Region] = []
    genes: list[Gene] = []

    protoclusters: list[dict[str, Any]] = []
    for feature in features:
        if feature.type != "protocluster":
            continue
        local_start, local_stop, strand = parse_gbk_location(feature.location)
        qualifiers = feature.qualifiers
        protoclusters.append(
            {
                "start": orig_start + local_start,
                "stop": orig_start + local_stop,
                "strand": strand,
                "number": q_one(qualifiers, "protocluster_number"),
                "category": q_one(qualifiers, "category"),
                "core_location": offset_core_location(q_one(qualifiers, "core_location"), orig_start)
                if q_one(qualifiers, "core_location")
                else "",
            }
        )

    for feature in features:
        if feature.type != "region":
            continue

        local_start, local_stop, strand = parse_gbk_location(feature.location)
        start = orig_start + local_start
        stop = orig_start + local_stop
        qualifiers = feature.qualifiers
        overlapping_protoclusters = [
            pc for pc in protoclusters if overlaps(start, stop, pc["start"], pc["stop"])
        ]
        regions.append(
            Region(
                source_json=str(path),
                record_id=record_id,
                contig=contig,
                region_number=q_one(qualifiers, "region_number", path.stem.rsplit("region", 1)[-1]),
                start=start,
                stop=stop,
                strand=strand,
                region_type=q_join(qualifiers, "product"),
                category=uniq_join(pc["category"] for pc in overlapping_protoclusters),
                candidate_cluster_numbers=q_join(qualifiers, "candidate_cluster_numbers"),
                protocluster_numbers=uniq_join(pc["number"] for pc in overlapping_protoclusters),
                core_locations=uniq_join(pc["core_location"] for pc in overlapping_protoclusters),
                contig_edge=q_one(qualifiers, "contig_edge"),
                rules=q_join(qualifiers, "rules"),
            )
        )

    for feature in features:
        if feature.type != "CDS":
            continue

        local_start, local_stop, strand = parse_gbk_location(feature.location)
        gene = Gene(
            source_json=str(path),
            record_id=record_id,
            contig=contig,
            start=orig_start + local_start,
            stop=orig_start + local_stop,
            strand=strand,
            locus_tag=q_one(feature.qualifiers, "locus_tag") or q_one(feature.qualifiers, "gene"),
            protein_id=q_one(feature.qualifiers, "protein_id"),
            product=q_join(feature.qualifiers, "product"),
            gene_kind=q_join(feature.qualifiers, "gene_kind"),
            gene_functions=q_join(feature.qualifiers, "gene_functions"),
            sec_met_domain=q_join(feature.qualifiers, "sec_met_domain"),
        )

        gene_regions = contained_or_overlapping_regions(gene, regions)
        if not gene_regions:
            continue

        gene.region_numbers = uniq_join(region.region_number for region in gene_regions)
        gene.region_types = uniq_join(region.region_type for region in gene_regions)
        gene.region_categories = uniq_join(region.category for region in gene_regions)
        gene.region_locations = uniq_join(
            format_interval(region.contig, region.start, region.stop, region.strand)
            for region in gene_regions
        )
        gene.region_contig_edges = uniq_join(region.contig_edge for region in gene_regions)
        genes.append(gene)

    return regions, genes


def read_gene_calls(path: Path) -> GeneCallIndex:
    with path.open("r", encoding="utf-8", newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters="\t,")
        except csv.Error:
            dialect = csv.excel_tab
        reader = csv.DictReader(handle, dialect=dialect)
        required = {"gene_callers_id", "contig", "start", "stop", "direction"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"Gene calls file is missing required column(s): {', '.join(sorted(missing))}"
            )

        index = GeneCallIndex()
        for row in reader:
            call = GeneCall(
                gene_callers_id=str(row["gene_callers_id"]),
                contig=str(row["contig"]),
                start=int(row["start"]),
                stop=int(row["stop"]),
                direction=str(row["direction"]).lower(),
            )
            index.exact[(call.contig, call.start, call.stop, call.direction)].append(call)
            index.by_contig[call.contig].append(call)
    return index


def reciprocal_overlap_fraction(gene: Gene, call: GeneCall) -> float:
    overlap = max(0, min(gene.stop, call.stop) - max(gene.start, call.start))
    gene_len = max(1, gene.stop - gene.start)
    call_len = max(1, call.stop - call.start)
    return min(overlap / gene_len, overlap / call_len)


def map_gene_to_call(gene: Gene, index: GeneCallIndex, min_reciprocal_overlap: float) -> None:
    direction = direction_from_strand(gene.strand)
    exact_matches = index.exact.get((gene.contig, gene.start, gene.stop, direction), [])
    if len(exact_matches) == 1:
        gene.gene_callers_id = exact_matches[0].gene_callers_id
        gene.gene_call_match = "exact"
        return
    if len(exact_matches) > 1:
        gene.gene_call_match = "ambiguous_exact"
        return

    candidates: list[tuple[float, GeneCall]] = []
    for call in index.by_contig.get(gene.contig, []):
        if direction != "?" and call.direction and call.direction != direction:
            continue
        score = reciprocal_overlap_fraction(gene, call)
        if score >= min_reciprocal_overlap:
            candidates.append((score, call))

    if not candidates:
        gene.gene_call_match = "unmapped"
        return

    candidates.sort(key=lambda item: (item[0], -abs(item[1].start - gene.start), -abs(item[1].stop - gene.stop)), reverse=True)
    best_score = candidates[0][0]
    best = [call for score, call in candidates if score == best_score]
    if len(best) == 1:
        gene.gene_callers_id = best[0].gene_callers_id
        gene.gene_call_match = f"overlap:{best_score:.3f}"
    else:
        gene.gene_call_match = "ambiguous_overlap"


def write_regions(path: Path, regions: list[Region]) -> None:
    fields = [
        "source_file",
        "record_id",
        "contig",
        "region_number",
        "start",
        "stop",
        "start_1based",
        "end_1based",
        "strand",
        "region_location",
        "region_type",
        "category",
        "candidate_cluster_numbers",
        "protocluster_numbers",
        "core_locations",
        "contig_edge",
        "rules",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for region in regions:
            writer.writerow(
                {
                    "source_file": region.source_json,
                    "record_id": region.record_id,
                    "contig": region.contig,
                    "region_number": region.region_number,
                    "start": region.start,
                    "stop": region.stop,
                    "start_1based": region.start + 1,
                    "end_1based": region.stop,
                    "strand": region.strand,
                    "region_location": format_interval(region.contig, region.start, region.stop, region.strand),
                    "region_type": region.region_type,
                    "category": region.category,
                    "candidate_cluster_numbers": region.candidate_cluster_numbers,
                    "protocluster_numbers": region.protocluster_numbers,
                    "core_locations": region.core_locations,
                    "contig_edge": region.contig_edge,
                    "rules": region.rules,
                }
            )


def write_genes(path: Path, genes: list[Gene]) -> None:
    fields = [
        "source_file",
        "record_id",
        "contig",
        "locus_tag",
        "protein_id",
        "start",
        "stop",
        "start_1based",
        "end_1based",
        "strand",
        "direction",
        "gene_location",
        "gene_kind",
        "gene_functions",
        "product",
        "sec_met_domain",
        "region_numbers",
        "region_types",
        "region_categories",
        "region_locations",
        "region_contig_edges",
        "gene_callers_id",
        "gene_call_match",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for gene in genes:
            writer.writerow(
                {
                    "source_file": gene.source_json,
                    "record_id": gene.record_id,
                    "contig": gene.contig,
                    "locus_tag": gene.locus_tag,
                    "protein_id": gene.protein_id,
                    "start": gene.start,
                    "stop": gene.stop,
                    "start_1based": gene.start + 1,
                    "end_1based": gene.stop,
                    "strand": gene.strand,
                    "direction": direction_from_strand(gene.strand),
                    "gene_location": format_interval(gene.contig, gene.start, gene.stop, gene.strand),
                    "gene_kind": gene.gene_kind,
                    "gene_functions": gene.gene_functions,
                    "product": gene.product,
                    "sec_met_domain": gene.sec_met_domain,
                    "region_numbers": gene.region_numbers,
                    "region_types": gene.region_types,
                    "region_categories": gene.region_categories,
                    "region_locations": gene.region_locations,
                    "region_contig_edges": gene.region_contig_edges,
                    "gene_callers_id": gene.gene_callers_id,
                    "gene_call_match": gene.gene_call_match,
                }
            )


def function_rows_for_gene(gene: Gene) -> list[dict[str, str]]:
    if not gene.gene_callers_id or not gene.gene_kind:
        return []

    gene_location = format_interval(gene.contig, gene.start, gene.stop, gene.strand)
    gene_kind = gene.gene_kind
    function_parts = [
        f"gene_kind: {gene_kind}",
        f"gene_location: {gene_location}",
        f"region_location: {gene.region_locations}",
        f"region_type: {gene.region_types}",
        f"contig_edge: {gene.region_contig_edges}",
    ]
    if gene.gene_functions:
        function_parts.append(f"gene_functions: {gene.gene_functions}")

    return [
        {
            "gene_callers_id": gene.gene_callers_id,
            "source": "antiSMASH",
            "accession": sanitize_accession(gene_kind),
            "function": "; ".join(function_parts),
            "e_value": "0",
        }
    ]


def write_functions(path: Path, genes: list[Gene]) -> int:
    fields = ["gene_callers_id", "source", "accession", "function", "e_value"]
    seen: set[tuple[str, str, str, str]] = set()
    count = 0
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for gene in genes:
            for row in function_rows_for_gene(gene):
                key = (row["gene_callers_id"], row["source"], row["accession"], row["function"])
                if key in seen:
                    continue
                writer.writerow(row)
                seen.add(key)
                count += 1
    return count


def write_nucleotide_misc(path: Path, regions: list[Region]) -> int:
    fields = [
        "item_name",
        "antismash_region_number",
        "antismash_region_type",
        "antismash_region_category",
        "data_group",
    ]
    count = 0
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for region in regions:
            for pos in range(region.start, region.stop):
                writer.writerow(
                    {
                        "item_name": f"{region.contig}:{pos}",
                        "antismash_region_number": region.region_number,
                        "antismash_region_type": region.region_type,
                        "antismash_region_category": region.category,
                        "data_group": "antiSMASH_regions",
                    }
                )
                count += 1
    return count


def main() -> int:
    args = parse_args()

    try:
        input_format, input_files = discover_input_files(args.input, args.input_format)
        all_regions: list[Region] = []
        all_genes: list[Gene] = []
        for input_file in input_files:
            if input_format == "gbk":
                regions, genes = load_antismash_gbk(input_file)
            else:
                regions, genes = load_antismash_json(input_file)
            all_regions.extend(regions)
            all_genes.extend(genes)

        if args.genes == "kinded":
            all_genes = [gene for gene in all_genes if gene.gene_kind]

        if args.gene_calls:
            index = read_gene_calls(args.gene_calls)
            for gene in all_genes:
                map_gene_to_call(gene, index, args.min_reciprocal_overlap)

        id_source = args.id_source
        if id_source == "auto":
            id_source = "gene_callers_id" if args.gene_calls else "locus_tag"
        apply_function_id_source(all_genes, id_source)

        args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
        regions_path = args.output_prefix.with_suffix(".regions.tsv")
        genes_path = args.output_prefix.with_suffix(".genes.tsv")
        functions_path = args.output_prefix.with_suffix(".functions.tsv")

        if not args.skip_region_gene_tables:
            write_regions(regions_path, all_regions)
            write_genes(genes_path, all_genes)

            print(f"Parsed {len(input_files)} antiSMASH {input_format.upper()} file(s)", file=sys.stderr)
            print(f"Wrote {len(all_regions)} antiSMASH regions: {regions_path}", file=sys.stderr)
            print(f"Wrote {len(all_genes)} CDS rows: {genes_path}", file=sys.stderr)

        function_count = write_functions(functions_path, all_genes)
        mapped = sum(1 for gene in all_genes if gene.gene_callers_id)
        print(
            f"Wrote {function_count} function rows for {mapped} genes using {id_source}: {functions_path}",
            file=sys.stderr,
        )
        if id_source != "gene_callers_id":
            print(
                "Used numeric gene_callers_id values parsed from anvio_gene_* locus_tag qualifiers.",
                file=sys.stderr,
            )

        if args.write_nucleotide_misc:
            nucleotide_path = args.output_prefix.with_suffix(".nucleotides.tsv")
            nt_count = write_nucleotide_misc(nucleotide_path, all_regions)
            print(f"Wrote {nt_count} nucleotide misc-data rows: {nucleotide_path}", file=sys.stderr)

    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
