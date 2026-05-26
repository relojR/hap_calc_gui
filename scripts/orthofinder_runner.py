#!/usr/bin/env python3
"""Validate FASTA inputs and run OrthoFinder."""

from __future__ import annotations

import argparse
import csv
import shlex
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path


FASTA_EXTENSIONS = {".fa", ".faa", ".fas", ".fasta", ".fna"}
NUCLEOTIDE_EXTENSIONS = {".fna"}
SEARCH_PROGRAMS = {"auto", "diamond", "blast", "mmseqs"}


@dataclass(frozen=True)
class FastaSummary:
    species_id: str
    path: Path
    sequence_count: int
    residue_count: int
    warnings: tuple[str, ...]


def find_fasta_files(input_dir: Path) -> list[Path]:
    if not input_dir.exists():
        raise ValueError(f"Input directory does not exist: {input_dir}")
    if not input_dir.is_dir():
        raise ValueError(f"Input path is not a directory: {input_dir}")
    return sorted(
        path
        for path in input_dir.iterdir()
        if path.is_file() and path.suffix.lower() in FASTA_EXTENSIONS
    )


def summarize_fasta(path: Path) -> FastaSummary:
    sequence_count = 0
    residue_count = 0
    seen_sequence = False
    warnings: list[str] = []
    sequence_ids: set[str] = set()
    duplicate_ids: set[str] = set()

    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                seen_sequence = True
                sequence_count += 1
                sequence_id = line[1:].split(maxsplit=1)[0]
                if not sequence_id:
                    warnings.append(f"line {line_number}: blank sequence id")
                elif sequence_id in sequence_ids:
                    duplicate_ids.add(sequence_id)
                else:
                    sequence_ids.add(sequence_id)
                continue
            if not seen_sequence:
                warnings.append(f"line {line_number}: sequence found before first header")
            residue_count += len(line.replace(" ", ""))

    if sequence_count == 0:
        warnings.append("no FASTA records found")
    if " " in path.name:
        warnings.append("filename contains spaces; use underscores for species names")
    if path.suffix.lower() in NUCLEOTIDE_EXTENSIONS:
        warnings.append("file looks nucleotide-based; OrthoFinder usually expects protein FASTA")
    if duplicate_ids:
        warnings.append(f"{len(duplicate_ids)} duplicate sequence id(s)")

    return FastaSummary(
        species_id=path.stem,
        path=path,
        sequence_count=sequence_count,
        residue_count=residue_count,
        warnings=tuple(warnings),
    )


def validate_fasta_dir(input_dir: Path) -> list[FastaSummary]:
    files = find_fasta_files(input_dir)
    if not files:
        raise ValueError(
            f"No FASTA files found in {input_dir}. Expected one species per .faa/.fa/.fasta file."
        )
    summaries = [summarize_fasta(path) for path in files]
    species_ids = [summary.species_id for summary in summaries]
    duplicates = sorted({sid for sid in species_ids if species_ids.count(sid) > 1})
    if duplicates:
        raise ValueError(f"Duplicate species ids from file names: {', '.join(duplicates)}")
    return summaries


def write_manifest(summaries: list[FastaSummary], output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest = output_dir / "orthofinder_input_manifest.tsv"
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "species_id",
                "fasta_file",
                "sequence_count",
                "residue_count",
                "warnings",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for summary in summaries:
            writer.writerow(
                {
                    "species_id": summary.species_id,
                    "fasta_file": str(summary.path),
                    "sequence_count": summary.sequence_count,
                    "residue_count": summary.residue_count,
                    "warnings": "; ".join(summary.warnings),
                }
            )
    return manifest


def build_orthofinder_command(
    input_dir: Path,
    output_dir: Path,
    orthofinder_bin: str = "orthofinder",
    threads: int = 8,
    analysis_threads: int | None = None,
    search_program: str = "auto",
    orthogroups_only: bool = True,
) -> list[str]:
    if threads < 1:
        raise ValueError("threads must be at least 1")
    if analysis_threads is None:
        analysis_threads = threads
    if analysis_threads < 1:
        raise ValueError("analysis_threads must be at least 1")
    if search_program not in SEARCH_PROGRAMS:
        raise ValueError(
            f"search_program must be one of: {', '.join(sorted(SEARCH_PROGRAMS))}"
        )

    command = [
        orthofinder_bin,
        "-f",
        str(input_dir),
        "-t",
        str(threads),
        "-a",
        str(analysis_threads),
        "-o",
        str(output_dir),
    ]
    if search_program != "auto":
        command.extend(["-S", search_program])
    if orthogroups_only:
        command.append("-og")
    return command


def command_display(command: list[str]) -> str:
    return " ".join(shlex.quote(part) for part in command)


def check_orthofinder_available(orthofinder_bin: str) -> str | None:
    path = Path(orthofinder_bin)
    if path.exists():
        return str(path)
    return shutil.which(orthofinder_bin)


def run_orthofinder(command: list[str]) -> int:
    completed = subprocess.run(command, check=False)
    return completed.returncode


def main() -> None:
    parser = argparse.ArgumentParser(description="Run OrthoFinder on a folder of FASTA files.")
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--orthofinder-bin", default="orthofinder")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--analysis-threads", type=int)
    parser.add_argument(
        "--search-program",
        default="auto",
        choices=sorted(SEARCH_PROGRAMS),
    )
    parser.add_argument(
        "--full-run",
        action="store_true",
        help="Run the full OrthoFinder workflow instead of stopping after orthogroups.",
    )
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()
    summaries = validate_fasta_dir(input_dir)
    manifest = write_manifest(summaries, output_dir)
    command = build_orthofinder_command(
        input_dir=input_dir,
        output_dir=output_dir,
        orthofinder_bin=args.orthofinder_bin,
        threads=args.threads,
        analysis_threads=args.analysis_threads,
        search_program=args.search_program,
        orthogroups_only=not args.full_run,
    )

    print(f"Validated {len(summaries)} FASTA file(s).")
    print(f"Wrote input manifest: {manifest}")
    print(command_display(command))

    if args.dry_run:
        return
    if not check_orthofinder_available(args.orthofinder_bin):
        raise SystemExit(
            "OrthoFinder was not found. Install it or pass --orthofinder-bin with the executable path."
        )
    raise SystemExit(run_orthofinder(command))


if __name__ == "__main__":
    main()
