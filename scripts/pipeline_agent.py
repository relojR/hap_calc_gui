#!/usr/bin/env python3
"""Small deterministic copilot for the orthogroup pipeline GUI."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any

from scripts.orthogroup_filter import (
    filter_orthogroups,
    is_present,
    read_orthogroups,
    read_species_metadata,
    split_members,
    summarize,
    validate_inputs,
)


DEFAULT_COMPARISONS = [(7, 0), (6, 0), (4, 2)]


def _path(value: Any) -> Path:
    return Path(str(value or "")).expanduser().resolve()


def _int(value: Any, default: int = 0) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def normalize_context(context: dict[str, Any]) -> dict[str, Any]:
    return {
        "orthogroups_path": _path(context.get("orthogroups_path")),
        "metadata_path": _path(context.get("metadata_path")),
        "output_dir": _path(context.get("output_dir", "results/presence_absence")),
        "direction": str(context.get("direction", "calcifying") or "calcifying"),
        "min_target_present": _int(context.get("min_target_present"), 4),
        "max_background_present": _int(context.get("max_background_present"), 2),
    }


def format_count(value: Any) -> str:
    try:
        return f"{int(value):,}"
    except (TypeError, ValueError):
        return str(value)


def summarize_dataset(context: dict[str, Any]) -> dict[str, Any]:
    ctx = normalize_context(context)
    summary = summarize(ctx["orthogroups_path"], ctx["metadata_path"])
    text = (
        f"Dataset loaded: {format_count(summary['orthogroup_count'])} orthogroups "
        f"across {summary['species_count']} species. "
        f"{summary['calcifying_count']} are calcifying and "
        f"{summary['non_calcifying_count']} are non-calcifying."
    )
    return {"text": text, "summary": summary}


def current_filter(context: dict[str, Any]) -> dict[str, Any]:
    ctx = normalize_context(context)
    rows, summary = filter_orthogroups(
        orthogroups_path=ctx["orthogroups_path"],
        metadata_path=ctx["metadata_path"],
        direction=ctx["direction"],
        min_target_present=ctx["min_target_present"],
        max_background_present=ctx["max_background_present"],
    )
    text = (
        f"The current filter found {format_count(summary['candidate_count'])} "
        f"{summary['target_group']}-specific candidate orthogroups. "
        f"Threshold: at least {ctx['min_target_present']}/{summary['target_total']} "
        f"target species and at most {ctx['max_background_present']}/"
        f"{summary['background_total']} background species."
    )
    return {"text": text, "summary": summary, "rows": rows[:10]}


def compare_thresholds(context: dict[str, Any]) -> dict[str, Any]:
    ctx = normalize_context(context)
    comparisons = list(DEFAULT_COMPARISONS)
    current = (ctx["min_target_present"], ctx["max_background_present"])
    if current not in comparisons:
        comparisons.append(current)

    lines = ["Threshold comparison:"]
    table = []
    for min_target, max_background in comparisons:
        _rows, summary = filter_orthogroups(
            orthogroups_path=ctx["orthogroups_path"],
            metadata_path=ctx["metadata_path"],
            direction=ctx["direction"],
            min_target_present=min_target,
            max_background_present=max_background,
        )
        label = (
            f">={min_target}/{summary['target_total']} {summary['target_group']}, "
            f"<={max_background}/{summary['background_total']} {summary['background_group']}"
        )
        count = int(summary["candidate_count"])
        lines.append(f"- {label}: {format_count(count)} candidates")
        table.append(
            {
                "min_target_present": min_target,
                "max_background_present": max_background,
                "candidate_count": count,
                "target_total": summary["target_total"],
                "background_total": summary["background_total"],
            }
        )
    return {"text": "\n".join(lines), "comparisons": table}


def inspect_orthogroup(context: dict[str, Any], orthogroup_id: str) -> dict[str, Any]:
    ctx = normalize_context(context)
    og_data = read_orthogroups(ctx["orthogroups_path"])
    metadata = read_species_metadata(ctx["metadata_path"])
    validate_inputs(og_data, metadata)
    background_group = (
        "non_calcifying" if ctx["direction"] == "calcifying" else "calcifying"
    )
    target_species = [
        species_id
        for species_id in og_data.species_columns
        if metadata[species_id]["group"] == ctx["direction"]
    ]
    background_species = [
        species_id
        for species_id in og_data.species_columns
        if metadata[species_id]["group"] == background_group
    ]

    row = next(
        (
            candidate
            for candidate in og_data.rows
            if str(candidate.get("Orthogroup", "")).lower() == orthogroup_id.lower()
        ),
        None,
    )
    if not row:
        return {"text": f"I could not find {orthogroup_id} in the orthogroups table."}

    target_present = [sid for sid in target_species if is_present(row.get(sid))]
    background_present = [sid for sid in background_species if is_present(row.get(sid))]
    member_count = sum(len(split_members(row.get(sid))) for sid in og_data.species_columns)
    passes = (
        len(target_present) >= ctx["min_target_present"]
        and len(background_present) <= ctx["max_background_present"]
    )
    text = (
        f"{row['Orthogroup']} has {len(target_present)}/{len(target_species)} "
        f"{ctx['direction']} species present and {len(background_present)}/"
        f"{len(background_species)} {background_group} species present. "
        f"It has {member_count} total member proteins. "
        f"With the current thresholds, it {'passes' if passes else 'does not pass'}."
    )
    if target_present:
        text += "\nTarget species present: " + "; ".join(target_present[:12])
        if len(target_present) > 12:
            text += f"; ... +{len(target_present) - 12} more"
    if background_present:
        text += "\nBackground species present: " + "; ".join(background_present[:12])
        if len(background_present) > 12:
            text += f"; ... +{len(background_present) - 12} more"
    return {
        "text": text,
        "orthogroup": row["Orthogroup"],
        "passes": passes,
        "target_present": len(target_present),
        "target_total": len(target_species),
        "background_present": len(background_present),
        "background_total": len(background_species),
        "member_count": member_count,
    }


def check_metadata(context: dict[str, Any]) -> dict[str, Any]:
    ctx = normalize_context(context)
    og_data = read_orthogroups(ctx["orthogroups_path"])
    metadata = read_species_metadata(ctx["metadata_path"])
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
    text = (
        "Metadata check passed. "
        f"{len(calcifying)} calcifying and {len(non_calcifying)} non-calcifying "
        "species match the Orthogroups.tsv columns."
    )
    if warnings:
        text += "\nWarnings:\n- " + "\n- ".join(warnings)
    return {"text": text, "warnings": warnings}


def methods_summary(context: dict[str, Any]) -> dict[str, Any]:
    ctx = normalize_context(context)
    rows, summary = filter_orthogroups(
        orthogroups_path=ctx["orthogroups_path"],
        metadata_path=ctx["metadata_path"],
        direction=ctx["direction"],
        min_target_present=ctx["min_target_present"],
        max_background_present=ctx["max_background_present"],
    )
    text = (
        "Methods summary:\n"
        f"Orthogroups were loaded from {ctx['orthogroups_path'].name}. Species were "
        f"classified using {ctx['metadata_path'].name}. Orthogroups were treated as "
        "present in a species when the corresponding OrthoFinder cell contained at "
        "least one protein identifier. Candidate orthogroups were selected by "
        f"requiring presence in at least {ctx['min_target_present']} of "
        f"{summary['target_total']} {summary['target_group']} species and presence in "
        f"no more than {ctx['max_background_present']} of {summary['background_total']} "
        f"{summary['background_group']} species. This produced "
        f"{format_count(len(rows))} candidate orthogroups."
    )
    return {"text": text, "summary": summary}


def help_text() -> str:
    return (
        "I can help with the current orthogroup analysis. Try:\n"
        "- summarize dataset\n"
        "- compare thresholds\n"
        "- check metadata\n"
        "- explain OG0000049\n"
        "- write a methods summary"
    )


def chat(message: str, context: dict[str, Any]) -> dict[str, Any]:
    lower = message.strip().lower()
    try:
        if not lower:
            return {"reply": help_text(), "tool": "help"}
        og_match = re.search(r"\bOG\d+\b", message, flags=re.IGNORECASE)
        if og_match:
            result = inspect_orthogroup(context, og_match.group(0).upper())
            return {"reply": result["text"], "tool": "inspect_orthogroup", "data": result}
        if "compare" in lower or "threshold" in lower and "?" in lower:
            result = compare_thresholds(context)
            return {"reply": result["text"], "tool": "compare_thresholds", "data": result}
        if "metadata" in lower or "check" in lower or "validate" in lower:
            result = check_metadata(context)
            return {"reply": result["text"], "tool": "check_metadata", "data": result}
        if "method" in lower or "report" in lower or "summary paragraph" in lower:
            result = methods_summary(context)
            return {"reply": result["text"], "tool": "methods_summary", "data": result}
        if "filter" in lower or "candidate" in lower or "current" in lower:
            result = current_filter(context)
            return {"reply": result["text"], "tool": "current_filter", "data": result}
        if "summar" in lower or "dataset" in lower or "species" in lower:
            result = summarize_dataset(context)
            return {"reply": result["text"], "tool": "summarize_dataset", "data": result}
        return {
            "reply": (
                "I can help with this, but I only have a few pipeline tools wired up so far.\n\n"
                + help_text()
            ),
            "tool": "help",
        }
    except Exception as exc:
        return {
            "reply": f"I hit a pipeline error while handling that request: {exc}",
            "tool": "error",
            "error": str(exc),
        }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--message", required=True)
    parser.add_argument("--context-json", required=True)
    args = parser.parse_args()
    print(json.dumps(chat(args.message, json.loads(args.context_json)), indent=2))


if __name__ == "__main__":
    main()
