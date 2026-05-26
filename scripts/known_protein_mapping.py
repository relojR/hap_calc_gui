#!/usr/bin/env python3
"""Map known calcification query proteins to OrthoFinder orthogroups."""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import zipfile
from dataclasses import dataclass
from pathlib import Path


FASTA_EXTENSIONS = {".fa", ".faa", ".fas", ".fasta", ".fna"}
DIAMOND_FIELDS = ["qseqid", "sseqid", "evalue", "bitscore", "pident", "length", "qlen", "slen"]


@dataclass(frozen=True)
class FastaRecord:
    record_id: str
    description: str
    sequence: str


def is_fasta(path: Path) -> bool:
    return path.suffix.lower() in FASTA_EXTENSIONS


def iter_fasta(path: Path):
    current_header = ""
    current_seq: list[str] = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header:
                    record_id = current_header[1:].split(maxsplit=1)[0]
                    yield FastaRecord(record_id, current_header[1:], "".join(current_seq))
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)
        if current_header:
            record_id = current_header[1:].split(maxsplit=1)[0]
            yield FastaRecord(record_id, current_header[1:], "".join(current_seq))


def count_fasta_records(path: Path) -> int:
    return sum(1 for _record in iter_fasta(path))


def list_fasta_files(path: Path, recursive: bool = False) -> list[Path]:
    if not path.exists():
        raise ValueError(f"Path does not exist: {path}")
    if path.is_file():
        if not is_fasta(path):
            raise ValueError(f"Not a FASTA file: {path}")
        return [path]
    globber = path.rglob if recursive else path.glob
    return sorted(file for file in globber("*") if file.is_file() and is_fasta(file))


def extract_query_fastas(query_source: Path, output_dir: Path) -> list[Path]:
    query_dir = output_dir / "query_fastas"
    query_dir.mkdir(parents=True, exist_ok=True)
    extracted: list[Path] = []

    if query_source.is_file() and query_source.suffix.lower() == ".zip":
        with zipfile.ZipFile(query_source) as archive:
            for member in archive.namelist():
                member_path = Path(member)
                if member.endswith("/") or member_path.suffix.lower() not in FASTA_EXTENSIONS:
                    continue
                target = query_dir / member_path.name
                with archive.open(member) as source, target.open("wb") as dest:
                    shutil.copyfileobj(source, dest)
                extracted.append(target)
        return sorted(extracted)

    if query_source.is_file():
        target = query_dir / query_source.name
        shutil.copy2(query_source, target)
        return [target]

    if query_source.is_dir():
        for fasta in list_fasta_files(query_source):
            target = query_dir / fasta.name
            shutil.copy2(fasta, target)
            extracted.append(target)
        return sorted(extracted)

    raise ValueError(f"Query source does not exist: {query_source}")


def read_orthogroup_lookup(orthogroups_tsv: Path) -> tuple[dict[str, str], int, list[str]]:
    if not orthogroups_tsv.exists():
        raise ValueError(f"Orthogroups table does not exist: {orthogroups_tsv}")
    protein_to_og: dict[str, str] = {}
    with orthogroups_tsv.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "Orthogroup" not in reader.fieldnames:
            raise ValueError("Orthogroups table must contain an Orthogroup column.")
        species_columns = [col.strip() for col in reader.fieldnames if col != "Orthogroup"]
        orthogroup_count = 0
        for row in reader:
            orthogroup_count += 1
            og = row.get("Orthogroup", "").strip()
            if not og:
                continue
            for species in species_columns:
                value = row.get(species)
                if not value:
                    continue
                for protein in str(value).split(","):
                    protein_id = protein.strip().split(maxsplit=1)[0]
                    if protein_id:
                        protein_to_og[protein_id] = og
    return protein_to_og, orthogroup_count, species_columns


def summarize_query_source(query_source: Path) -> dict[str, object]:
    files = []
    if query_source.is_file() and query_source.suffix.lower() == ".zip":
        with zipfile.ZipFile(query_source) as archive:
            for member in archive.namelist():
                member_path = Path(member)
                if member.endswith("/") or member_path.suffix.lower() not in FASTA_EXTENSIONS:
                    continue
                lines = archive.read(member).decode("utf-8", errors="replace").splitlines()
                count = sum(1 for line in lines if line.startswith(">"))
                files.append({"name": member_path.name, "path": member, "query_count": count})
    else:
        for fasta in list_fasta_files(query_source):
            files.append({"name": fasta.name, "path": str(fasta), "query_count": count_fasta_records(fasta)})
    return {
        "query_file_count": len(files),
        "query_count": sum(int(item["query_count"]) for item in files),
        "query_files": files,
    }


def summarize_inputs(
    orthogroups_tsv: Path,
    species_fasta_dir: Path,
    query_source: Path,
    diamond_bin: str,
    recursive_species: bool = False,
) -> dict[str, object]:
    protein_to_og, orthogroup_count, species_columns = read_orthogroup_lookup(orthogroups_tsv)
    query_summary = summarize_query_source(query_source)
    target_files = []
    target_sequence_count = 0
    target_error = ""
    try:
        target_paths = list_fasta_files(species_fasta_dir, recursive=recursive_species)
        for fasta in target_paths:
            count = count_fasta_records(fasta)
            target_sequence_count += count
            target_files.append({"name": fasta.name, "path": str(fasta), "sequence_count": count})
    except Exception as exc:
        target_error = str(exc)

    diamond_path = shutil.which(diamond_bin) or (str(Path(diamond_bin)) if Path(diamond_bin).exists() else "")
    return {
        "orthogroup_count": orthogroup_count,
        "orthogroup_species_count": len(species_columns),
        "og_assigned_protein_count": len(protein_to_og),
        "query_file_count": query_summary["query_file_count"],
        "query_count": query_summary["query_count"],
        "query_files": query_summary["query_files"],
        "species_fasta_count": len(target_files),
        "species_sequence_count": target_sequence_count,
        "species_fastas": target_files[:100],
        "species_error": target_error,
        "diamond_found": bool(diamond_path),
        "diamond_path": diamond_path,
    }


def write_og_tagged_fasta(
    protein_to_og: dict[str, str],
    species_fasta_dir: Path,
    output_fasta: Path,
    recursive_species: bool = False,
) -> dict[str, int]:
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    fasta_files = list_fasta_files(species_fasta_dir, recursive=recursive_species)
    if not fasta_files:
        raise ValueError(f"No species FASTA files found in {species_fasta_dir}")

    total_sequences = 0
    written = 0
    seen_assigned: set[str] = set()
    with output_fasta.open("w", encoding="utf-8") as out:
        for fasta in fasta_files:
            for record in iter_fasta(fasta):
                total_sequences += 1
                og = protein_to_og.get(record.record_id)
                if not og:
                    continue
                seen_assigned.add(record.record_id)
                out.write(f">{og}|||{record.record_id}\n")
                sequence = record.sequence
                for index in range(0, len(sequence), 60):
                    out.write(sequence[index : index + 60] + "\n")
                written += 1
    return {
        "species_fasta_count": len(fasta_files),
        "species_sequence_count": total_sequences,
        "og_tagged_sequence_count": written,
        "og_assigned_missing_from_fastas": len(protein_to_og) - len(seen_assigned),
    }


def run_checked(command: list[str]) -> None:
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    if completed.returncode != 0:
        raise RuntimeError(
            "Command failed: "
            + " ".join(command)
            + "\nSTDOUT:\n"
            + completed.stdout
            + "\nSTDERR:\n"
            + completed.stderr
        )


def build_diamond_database(diamond_bin: str, combined_fasta: Path, db_path: Path, threads: int) -> None:
    run_checked(
        [
            diamond_bin,
            "makedb",
            "--in",
            str(combined_fasta),
            "--db",
            str(db_path),
            "--threads",
            str(threads),
            "--quiet",
        ]
    )


def run_diamond_searches(
    diamond_bin: str,
    query_fastas: list[Path],
    db_path: Path,
    diamond_results_dir: Path,
    threads: int,
    evalue: str,
    max_target_seqs: int,
    sensitive: bool,
) -> dict[str, Path]:
    diamond_results_dir.mkdir(parents=True, exist_ok=True)
    outputs = {}
    for query_fasta in query_fastas:
        category = query_fasta.stem[:2] if query_fasta.stem[:2].isdigit() else query_fasta.stem
        out_tsv = diamond_results_dir / f"cat{category}_diamond.tsv"
        command = [
            diamond_bin,
            "blastp",
            "--query",
            str(query_fasta),
            "--db",
            str(db_path),
            "--out",
            str(out_tsv),
            "--outfmt",
            "6",
            *DIAMOND_FIELDS,
            "--evalue",
            str(evalue),
            "--threads",
            str(threads),
            "--max-target-seqs",
            str(max_target_seqs),
            "--quiet",
        ]
        if sensitive:
            command.insert(command.index("--threads"), "--sensitive")
        run_checked(command)
        outputs[category] = out_tsv
    return outputs


def best_hits_by_query(diamond_tsv: Path) -> dict[str, dict[str, str]]:
    best: dict[str, dict[str, str]] = {}
    if not diamond_tsv.exists() or diamond_tsv.stat().st_size == 0:
        return best
    with diamond_tsv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, fieldnames=DIAMOND_FIELDS, delimiter="\t")
        for row in reader:
            qseqid = row["qseqid"]
            bitscore = float(row["bitscore"])
            if qseqid not in best or bitscore > float(best[qseqid]["bitscore"]):
                best[qseqid] = row
    return best


def postprocess_diamond_outputs(
    query_fastas: list[Path],
    diamond_outputs: dict[str, Path],
    output_dir: Path,
) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_rows = []
    aggregate_rows = []
    for query_fasta in query_fastas:
        category = query_fasta.stem[:2] if query_fasta.stem[:2].isdigit() else query_fasta.stem
        query_ids = [record.record_id for record in iter_fasta(query_fasta)]
        best_hits = best_hits_by_query(diamond_outputs.get(category, Path()))
        output_rows = []
        mapped = 0
        unique_ogs = set()
        for query_id in query_ids:
            hit = best_hits.get(query_id)
            if hit:
                orthogroup = hit["sseqid"].split("|||", 1)[0] if "|||" in hit["sseqid"] else ""
                qcoverage = round(float(hit["length"]) / float(hit["qlen"]) * 100, 1)
                mapped += 1
                if orthogroup:
                    unique_ogs.add(orthogroup)
                row = {
                    "protein_id": query_id,
                    "orthogroup": orthogroup,
                    "subject_id": hit["sseqid"].split("|||", 1)[1] if "|||" in hit["sseqid"] else hit["sseqid"],
                    "evalue": hit["evalue"],
                    "bitscore": hit["bitscore"],
                    "pident": hit["pident"],
                    "qcoverage": qcoverage,
                }
            else:
                row = {
                    "protein_id": query_id,
                    "orthogroup": "",
                    "subject_id": "",
                    "evalue": "",
                    "bitscore": "",
                    "pident": "",
                    "qcoverage": "",
                }
            output_rows.append(row)
            aggregate_rows.append({"category": category, **row})

        out_tsv = output_dir / f"{query_fasta.stem}_orthogroup_mapping.tsv"
        out_csv = output_dir / f"{query_fasta.stem}_orthogroup_mapping.csv"
        fieldnames = ["protein_id", "orthogroup", "subject_id", "evalue", "bitscore", "pident", "qcoverage"]
        with out_tsv.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(output_rows)
        with out_csv.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(output_rows)

        summary_rows.append(
            {
                "category": category,
                "query_file": query_fasta.name,
                "query_count": len(query_ids),
                "mapped": mapped,
                "unmapped": len(query_ids) - mapped,
                "unique_orthogroups": len(unique_ogs),
                "mapping_tsv": str(out_tsv),
            }
        )

    summary_tsv = output_dir / "known_protein_mapping_summary.tsv"
    with summary_tsv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["category", "query_file", "query_count", "mapped", "unmapped", "unique_orthogroups", "mapping_tsv"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(summary_rows)

    aggregate_tsv = output_dir / "known_protein_mapping_all.tsv"
    with aggregate_tsv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["category", "protein_id", "orthogroup", "subject_id", "evalue", "bitscore", "pident", "qcoverage"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(aggregate_rows)

    return {
        "summary_rows": summary_rows,
        "summary_tsv": str(summary_tsv),
        "aggregate_tsv": str(aggregate_tsv),
    }


def run_mapping(
    orthogroups_tsv: Path,
    species_fasta_dir: Path,
    query_source: Path,
    output_dir: Path,
    diamond_bin: str = "diamond",
    threads: int = 8,
    evalue: str = "1e-5",
    max_target_seqs: int = 5,
    sensitive: bool = True,
    recursive_species: bool = False,
) -> dict[str, object]:
    diamond_path = shutil.which(diamond_bin) or (str(Path(diamond_bin)) if Path(diamond_bin).exists() else "")
    if not diamond_path:
        raise ValueError("DIAMOND was not found. Install DIAMOND or provide the executable path.")

    output_dir.mkdir(parents=True, exist_ok=True)
    work_dir = output_dir / "work"
    db_dir = work_dir / "diamond_db"
    diamond_results = work_dir / "diamond_results"
    mapping_dir = output_dir / "orthogroup_mapping2"
    for path in [work_dir, db_dir, diamond_results, mapping_dir]:
        path.mkdir(parents=True, exist_ok=True)

    print("Parsing Orthogroups.tsv...")
    protein_to_og, orthogroup_count, species_columns = read_orthogroup_lookup(orthogroups_tsv)
    print(f"Proteins with OG assignments: {len(protein_to_og):,}")
    print(f"Unique species columns: {len(species_columns):,}")
    print(f"Orthogroups: {orthogroup_count:,}")

    print("Extracting query FASTAs...")
    query_fastas = extract_query_fastas(query_source, work_dir)
    if not query_fastas:
        raise ValueError("No query FASTA files found.")
    print(f"Query FASTA files: {len(query_fastas)}")

    combined_fasta = db_dir / "all_orthogroups_combined.fasta"
    print("Writing OG-tagged combined FASTA...")
    combined_stats = write_og_tagged_fasta(
        protein_to_og,
        species_fasta_dir,
        combined_fasta,
        recursive_species=recursive_species,
    )
    print(f"OG-tagged sequences written: {combined_stats['og_tagged_sequence_count']:,}")
    if int(combined_stats["og_tagged_sequence_count"]) == 0:
        raise ValueError(
            "No species FASTA records matched the Orthogroups.tsv protein IDs. "
            "Check that the species FASTA folder belongs to the same OrthoFinder run."
        )

    db_path = db_dir / "all_orthogroups"
    print("Building DIAMOND database...")
    build_diamond_database(diamond_path, combined_fasta, db_path, threads)

    print("Running DIAMOND BLASTP searches...")
    diamond_outputs = run_diamond_searches(
        diamond_path,
        query_fastas,
        db_path,
        diamond_results,
        threads=threads,
        evalue=evalue,
        max_target_seqs=max_target_seqs,
        sensitive=sensitive,
    )

    print("Post-processing DIAMOND outputs...")
    result = postprocess_diamond_outputs(query_fastas, diamond_outputs, mapping_dir)
    print(f"Summary table: {result['summary_tsv']}")
    print(f"Aggregate table: {result['aggregate_tsv']}")
    return {
        "combined_stats": combined_stats,
        **result,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Map known proteins to OrthoFinder orthogroups.")
    parser.add_argument("--orthogroups", required=True, type=Path)
    parser.add_argument("--species-fasta-dir", required=True, type=Path)
    parser.add_argument("--queries", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--diamond-bin", default="diamond")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--evalue", default="1e-5")
    parser.add_argument("--max-target-seqs", type=int, default=5)
    parser.add_argument("--fast", action="store_true", help="Do not pass DIAMOND --sensitive.")
    parser.add_argument("--recursive-species", action="store_true")
    parser.add_argument("--summary-only", action="store_true")
    args = parser.parse_args()

    if args.summary_only:
        summary = summarize_inputs(
            args.orthogroups,
            args.species_fasta_dir,
            args.queries,
            args.diamond_bin,
            recursive_species=args.recursive_species,
        )
        for key, value in summary.items():
            if key not in {"query_files", "species_fastas"}:
                print(f"{key}: {value}")
        return

    run_mapping(
        orthogroups_tsv=args.orthogroups,
        species_fasta_dir=args.species_fasta_dir,
        query_source=args.queries,
        output_dir=args.output_dir,
        diamond_bin=args.diamond_bin,
        threads=args.threads,
        evalue=args.evalue,
        max_target_seqs=args.max_target_seqs,
        sensitive=not args.fast,
        recursive_species=args.recursive_species,
    )


if __name__ == "__main__":
    main()
