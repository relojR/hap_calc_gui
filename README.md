# Haptophyte Orthogroup Pipeline GUI

This project contains a local GUI for the first two pipeline stages:

- running OrthoFinder on one protein FASTA per species
- filtering an existing `Orthogroups.tsv` by calcifying/non-calcifying presence and absence

The presence/absence screen is the main working view right now. It does not include ranking steps.

## Start the GUI

From this folder:

```powershell
.\run_gui.ps1
```

Then open:

```text
http://127.0.0.1:8765
```

The default page is the `Orthogroups.tsv` presence/absence filter. The OrthoFinder launcher is available from the `OrthoFinder` tab.

The `Known Proteins` tab maps known calcification proteins to orthogroups using the V2 DIAMOND-style workflow from the reference notebook. Plotting sections are not included.

## Presence/Absence Filtering

Default demo inputs:

```text
examples/orthogroups/Orthogroups.tsv
examples/species_metadata.tsv
```

The species metadata table requires:

```text
species_id    label    group
```

`species_id` must exactly match the species column names in `Orthogroups.tsv`. Use `calcifying` or `non_calcifying` in the `group` column.

The GUI controls:

- `Direction`: calcifying-specific or non-calcifying-specific
- `Target Min Present`: minimum number of target-group species with the orthogroup
- `Background Max Present`: maximum number of background-group species with the orthogroup

For example, with 7 calcifying and 20 non-calcifying species:

```text
Direction = Calcifying-specific
Target Min Present = 6
Background Max Present = 0
```

This returns orthogroups present in at least 6 of 7 calcifying species and absent from all 20 non-calcifying species.

Exports are written to:

```text
results/presence_absence/candidate_orthogroups.tsv
results/presence_absence/candidate_orthogroups_summary.tsv
```

## Command-Line Filtering

```powershell
python scripts/orthogroup_filter.py `
  --orthogroups examples/orthogroups/Orthogroups.tsv `
  --metadata examples/species_metadata.tsv `
  --direction calcifying `
  --min-target-present 6 `
  --max-background-present 0 `
  --output-dir results/presence_absence
```

## Known Calcification Protein Mapping

The mapping workflow follows these non-plotting steps:

1. Parse `Orthogroups.tsv` into a protein-to-orthogroup lookup.
2. Read the species proteome FASTA files from the same OrthoFinder run.
3. Write an OG-tagged FASTA containing only proteins assigned to orthogroups.
4. Build a DIAMOND database from that OG-tagged FASTA.
5. Search known calcification protein FASTAs against the database.
6. Save per-category mapping tables and a combined summary.

Current default query source:

```text
C:\Users\j2rel\Downloads\known_calc_proteins.zip
```

That zip contains 6 query FASTA files with 2,070 known calcification proteins total.

This step also requires the 27 species proteome FASTA files used in the OrthoFinder run. `Orthogroups.tsv` alone is not enough because it contains protein IDs, not protein sequences.

Expected outputs:

```text
results/known_protein_mapping/orthogroup_mapping2/
  *_orthogroup_mapping.tsv
  *_orthogroup_mapping.csv
  known_protein_mapping_summary.tsv
  known_protein_mapping_all.tsv
```

Command-line version:

```powershell
python scripts/known_protein_mapping.py `
  --orthogroups "C:\Users\j2rel\Downloads\Orthogroups_diamond (1).tsv" `
  --species-fasta-dir path/to/species_fastas `
  --queries "C:\Users\j2rel\Downloads\known_calc_proteins.zip" `
  --output-dir results/known_protein_mapping `
  --diamond-bin diamond `
  --threads 8
```

## OrthoFinder Input Format

Place one protein FASTA file per species in a single folder:

```text
my_fastas/
  Emiliania_huxleyi.faa
  Gephyrocapsa_oceanica.faa
  Phaeocystis_globosa.faa
```

OrthoFinder normally expects protein FASTA files. The GUI accepts `.faa`, `.fa`, `.fas`, `.fasta`, and `.fna`, but it warns when files look nucleotide-based.

The species name used by OrthoFinder comes from the file stem, so prefer concise filenames with underscores and no spaces.

## OrthoFinder Workflow

The OrthoFinder tab defaults to:

```text
examples/proteomes
results/orthofinder
```

1. Set the species FASTA folder.
2. Set the output folder.
3. Set the OrthoFinder executable, usually `orthofinder`.
4. Choose threads and search program.
5. Leave `Stop after orthogroups` checked for this first pipeline stage.
6. Click `Validate`.
7. Click `Build Command` to preview the exact command.
8. Click `Run` when OrthoFinder is installed and available.

Before each command preview or run, the app writes:

```text
results/orthofinder/orthofinder_input_manifest.tsv
```

## Command-Line OrthoFinder Runner

```powershell
python scripts/orthofinder_runner.py `
  --input-dir examples/proteomes `
  --output-dir results/orthofinder `
  --threads 8 `
  --analysis-threads 8 `
  --dry-run
```

Remove `--dry-run` to launch OrthoFinder. Add `--full-run` only if you want OrthoFinder to continue beyond orthogroup inference.
