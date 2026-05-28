#!/usr/bin/env python3
"""Filter OrthoFinder orthogroups by phenotype-specific presence/absence."""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path


ABSENT_VALUES = {"", "na", "nan", "none", "-"}
GROUP_ALIASES = {
    "calcifying": "calcifying",
    "calc": "calcifying",
    "yes": "calcifying",
    "true": "calcifying",
    "1": "calcifying",
    "non_calcifying": "non_calcifying",
    "non-calcifying": "non_calcifying",
    "noncalcifying": "non_calcifying",
    "noncalc": "non_calcifying",
    "no": "non_calcifying",
    "false": "non_calcifying",
    "0": "non_calcifying",
}


@dataclass(frozen=True)
class OrthogroupData:
    rows: list[dict[str, str]]
    species_columns: list[str]


def normalize_group(value: str) -> str:
    key = value.strip().lower().replace(" ", "_")
    if key not in GROUP_ALIASES:
        raise ValueError(
            f"Unknown species group '{value}'. Use calcifying or non_calcifying."
        )
    return GROUP_ALIASES[key]


def is_present(value: str | None) -> bool:
    if value is None:
        return False
    return value.strip().lower() not in ABSENT_VALUES


def split_members(value: str | None) -> list[str]:
    if not is_present(value):
        return []
    return [part.strip() for part in str(value).split(",") if part.strip()]


def read_orthogroups(path: Path) -> OrthogroupData:
    if not path.exists():
        raise ValueError(f"Orthogroups table does not exist: {path}")
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Orthogroups table is empty: {path}")
        fieldnames = [name.strip() for name in reader.fieldnames]
        if "Orthogroup" not in fieldnames:
            raise ValueError("Orthogroups table must contain an Orthogroup column.")
        rows = []
        for row in reader:
            rows.append({(key or "").strip(): value for key, value in row.items()})
    species_columns = [name for name in fieldnames if name != "Orthogroup"]
    if not species_columns:
        raise ValueError("Orthogroups table must contain at least one species column.")
    return OrthogroupData(rows=rows, species_columns=species_columns)


def read_species_metadata(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        raise ValueError(f"Species metadata does not exist: {path}")
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Species metadata is empty: {path}")
        required = {"species_id", "group"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(
                f"Species metadata is missing column(s): {', '.join(sorted(missing))}"
            )
        metadata = {}
        for row in reader:
            species_id = row.get("species_id", "").strip()
            if not species_id:
                continue
            group = normalize_group(row.get("group", ""))
            label = row.get("label", species_id).strip() or species_id
            metadata[species_id] = {
                "species_id": species_id,
                "label": label,
                "group": group,
            }
    if not metadata:
        raise ValueError(f"No species rows found in metadata: {path}")
    return metadata


def validate_inputs(
    og_data: OrthogroupData, metadata: dict[str, dict[str, str]]
) -> list[str]:
    warnings = []
    og_species = set(og_data.species_columns)
    metadata_species = set(metadata)
    missing_metadata = sorted(og_species - metadata_species)
    extra_metadata = sorted(metadata_species - og_species)
    if missing_metadata:
        raise ValueError(
            "Species missing from metadata: " + ", ".join(missing_metadata)
        )
    if extra_metadata:
        warnings.append(
            "Metadata species not found in Orthogroups.tsv: " + ", ".join(extra_metadata)
        )
    group_counts = {
        "calcifying": sum(
            1
            for species_id in og_data.species_columns
            if metadata[species_id]["group"] == "calcifying"
        ),
        "non_calcifying": sum(
            1
            for species_id in og_data.species_columns
            if metadata[species_id]["group"] == "non_calcifying"
        ),
    }
    if group_counts["calcifying"] == 0 or group_counts["non_calcifying"] == 0:
        raise ValueError("Metadata must include both calcifying and non_calcifying species.")
    return warnings


def summarize(
    orthogroups_path: Path,
    metadata_path: Path,
) -> dict[str, object]:
    og_data = read_orthogroups(orthogroups_path)
    metadata = read_species_metadata(metadata_path)
    warnings = validate_inputs(og_data, metadata)
    calcifying = [
        species_id
        for species_id in og_data.species_columns
        if metadata[species_id]["group"] == "calcifying"
    ]
    non_calcifying = [
        species_id
        for species_id in og_data.species_columns
        if metadata[species_id]["group"] == "non_calcifying"
    ]
    return {
        "orthogroup_count": len(og_data.rows),
        "species_count": len(og_data.species_columns),
        "calcifying_count": len(calcifying),
        "non_calcifying_count": len(non_calcifying),
        "calcifying_species": calcifying,
        "non_calcifying_species": non_calcifying,
        "warnings": warnings,
    }


def filter_orthogroups(
    orthogroups_path: Path,
    metadata_path: Path,
    direction: str,
    min_target_present: int,
    max_background_present: int,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    if direction not in {"calcifying", "non_calcifying"}:
        raise ValueError("direction must be calcifying or non_calcifying")
    og_data = read_orthogroups(orthogroups_path)
    metadata = read_species_metadata(metadata_path)
    warnings = validate_inputs(og_data, metadata)
    background_group = (
        "non_calcifying" if direction == "calcifying" else "calcifying"
    )
    target_species = [
        species_id
        for species_id in og_data.species_columns
        if metadata[species_id]["group"] == direction
    ]
    background_species = [
        species_id
        for species_id in og_data.species_columns
        if metadata[species_id]["group"] == background_group
    ]
    if min_target_present < 0 or max_background_present < 0:
        raise ValueError("Thresholds must be zero or greater.")
    if min_target_present > len(target_species):
        raise ValueError("Target threshold exceeds number of target species.")
    if max_background_present > len(background_species):
        raise ValueError("Background threshold exceeds number of background species.")

    results = []
    tier_counts: dict[str, int] = {}
    for row in og_data.rows:
        og = row["Orthogroup"]
        present_target = [sid for sid in target_species if is_present(row.get(sid))]
        present_background = [
            sid for sid in background_species if is_present(row.get(sid))
        ]
        target_count = len(present_target)
        background_count = len(present_background)
        if target_count < min_target_present:
            continue
        if background_count > max_background_present:
            continue
        target_fraction = target_count / len(target_species)
        background_fraction = background_count / len(background_species)
        member_count = sum(len(split_members(row.get(sid))) for sid in og_data.species_columns)
        tier = f"{target_count}/{len(target_species)}"
        tier_counts[tier] = tier_counts.get(tier, 0) + 1
        results.append(
            {
                "Orthogroup": og,
                "target_group": direction,
                "target_present": target_count,
                "target_total": len(target_species),
                "target_fraction": round(target_fraction, 4),
                "background_group": background_group,
                "background_present": background_count,
                "background_total": len(background_species),
                "background_fraction": round(background_fraction, 4),
                "specificity": round(target_fraction - background_fraction, 4),
                "tier": tier,
                "member_count": member_count,
                "target_species_present": ";".join(present_target),
                "background_species_present": ";".join(present_background),
            }
        )

    results.sort(
        key=lambda item: (
            -float(item["specificity"]),
            -int(item["target_present"]),
            int(item["background_present"]),
            str(item["Orthogroup"]),
        )
    )
    summary = {
        "orthogroup_count": len(og_data.rows),
        "species_count": len(og_data.species_columns),
        "target_group": direction,
        "target_total": len(target_species),
        "background_group": background_group,
        "background_total": len(background_species),
        "candidate_count": len(results),
        "tier_counts": tier_counts,
        "target_species": target_species,
        "background_species": background_species,
        "calcifying_species": [
            species_id
            for species_id in og_data.species_columns
            if metadata[species_id]["group"] == "calcifying"
        ],
        "non_calcifying_species": [
            species_id
            for species_id in og_data.species_columns
            if metadata[species_id]["group"] == "non_calcifying"
        ],
        "warnings": warnings,
    }
    return results, summary


def write_results(
    results: list[dict[str, object]],
    summary: dict[str, object],
    output_dir: Path,
) -> dict[str, str]:
    output_dir.mkdir(parents=True, exist_ok=True)
    candidate_path = output_dir / "candidate_orthogroups.tsv"
    summary_path = output_dir / "candidate_orthogroups_summary.tsv"

    fieldnames = [
        "Orthogroup",
        "target_group",
        "target_present",
        "target_total",
        "target_fraction",
        "background_group",
        "background_present",
        "background_total",
        "background_fraction",
        "specificity",
        "tier",
        "member_count",
        "target_species_present",
        "background_species_present",
    ]
    with candidate_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)

    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        for key, value in summary.items():
            writer.writerow([key, value])

    return {
        "candidate_table": str(candidate_path),
        "summary_table": str(summary_path),
    }


def write_run_report(
    results: list[dict[str, object]],
    summary: dict[str, object],
    output_dir: Path,
    orthogroups_path: Path,
    metadata_path: Path,
    direction: str,
    min_target_present: int,
    max_background_present: int,
) -> dict[str, str]:
    paths = write_results(results, summary, output_dir)
    settings_path = output_dir / "settings.json"
    report_path = output_dir / "run_report.md"

    settings = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "orthogroups_path": str(orthogroups_path),
        "metadata_path": str(metadata_path),
        "direction": direction,
        "min_target_present": min_target_present,
        "max_background_present": max_background_present,
        "summary": summary,
    }
    with settings_path.open("w", encoding="utf-8") as handle:
        json.dump(settings, handle, indent=2)

    tier_lines = []
    for tier, count in sorted(
        summary.get("tier_counts", {}).items(),
        key=lambda item: tuple(int(part) for part in str(item[0]).split("/", 1)),
        reverse=True,
    ):
        tier_lines.append(f"- {tier}: {count}")
    if not tier_lines:
        tier_lines.append("- No candidate tiers found")

    warning_lines = [
        f"- {warning}" for warning in summary.get("warnings", [])
    ] or ["- None"]
    top_rows = results[:20]
    top_table = [
        "| Orthogroup | Target | Background | Specificity | Tier | Members |",
        "| --- | ---: | ---: | ---: | --- | ---: |",
    ]
    for row in top_rows:
        top_table.append(
            "| {og} | {target}/{target_total} | {background}/{background_total} | "
            "{specificity} | {tier} | {members} |".format(
                og=row["Orthogroup"],
                target=row["target_present"],
                target_total=row["target_total"],
                background=row["background_present"],
                background_total=row["background_total"],
                specificity=row["specificity"],
                tier=row["tier"],
                members=row["member_count"],
            )
        )
    if not top_rows:
        top_table.append("| No candidates |  |  |  |  |  |")

    target_species = "\n".join(
        f"- {species_id}" for species_id in summary.get("target_species", [])
    ) or "- None"
    background_species = "\n".join(
        f"- {species_id}" for species_id in summary.get("background_species", [])
    ) or "- None"

    report = f"""# Orthogroup Presence/Absence Run Report

Generated: {settings["created_at"]}

## Inputs

- Orthogroups table: `{orthogroups_path}`
- Species metadata: `{metadata_path}`

## Filter

- Direction: `{direction}`
- Target group: `{summary["target_group"]}`
- Background group: `{summary["background_group"]}`
- Minimum target presence: `{min_target_present}/{summary["target_total"]}`
- Maximum background presence: `{max_background_present}/{summary["background_total"]}`

## Results

- Orthogroups evaluated: {summary["orthogroup_count"]}
- Species evaluated: {summary["species_count"]}
- Candidate orthogroups: {summary["candidate_count"]}

## Tier Breakdown

{chr(10).join(tier_lines)}

## Target Species

{target_species}

## Background Species

{background_species}

## Warnings

{chr(10).join(warning_lines)}

## Top Candidate Orthogroups

{chr(10).join(top_table)}

## Output Files

- Candidate table: `{paths["candidate_table"]}`
- Summary table: `{paths["summary_table"]}`
- Settings JSON: `{settings_path}`
"""
    report_path.write_text(report, encoding="utf-8")

    return {
        **paths,
        "settings_json": str(settings_path),
        "run_report": str(report_path),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Filter OrthoFinder Orthogroups.tsv by phenotype presence/absence."
    )
    parser.add_argument("--orthogroups", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument(
        "--direction",
        choices=["calcifying", "non_calcifying"],
        default="calcifying",
    )
    parser.add_argument("--min-target-present", required=True, type=int)
    parser.add_argument("--max-background-present", required=True, type=int)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--report", action="store_true")
    args = parser.parse_args()

    results, summary = filter_orthogroups(
        orthogroups_path=args.orthogroups,
        metadata_path=args.metadata,
        direction=args.direction,
        min_target_present=args.min_target_present,
        max_background_present=args.max_background_present,
    )
    print(f"Candidate orthogroups: {summary['candidate_count']}")
    if args.output_dir:
        if args.report:
            paths = write_run_report(
                results,
                summary,
                args.output_dir,
                args.orthogroups,
                args.metadata,
                args.direction,
                args.min_target_present,
                args.max_background_present,
            )
        else:
            paths = write_results(results, summary, args.output_dir)
        for label, path in paths.items():
            print(f"{label}: {path}")


if __name__ == "__main__":
    main()
