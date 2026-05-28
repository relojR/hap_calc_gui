#!/usr/bin/env python3
"""Small local web GUI for launching OrthoFinder."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import threading
import time
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlparse

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.orthofinder_runner import (  # noqa: E402
    build_orthofinder_command,
    check_orthofinder_available,
    command_display,
    validate_fasta_dir,
    write_manifest,
)
from scripts.orthogroup_filter import (  # noqa: E402
    filter_orthogroups,
    summarize,
    write_run_report,
    write_results,
)
from scripts.known_protein_mapping import summarize_inputs as summarize_mapping_inputs  # noqa: E402
from scripts.pipeline_agent import chat as agent_chat  # noqa: E402


DEFAULT_INPUT = PROJECT_ROOT / "examples" / "proteomes"
DEFAULT_OUTPUT = PROJECT_ROOT / "results" / "orthofinder"
USER_OG_TSV = Path(r"C:\Users\j2rel\Downloads\Orthogroups_diamond (1).tsv")
DEFAULT_OG_TSV = (
    USER_OG_TSV
    if USER_OG_TSV.exists()
    else PROJECT_ROOT / "examples" / "orthogroups" / "Orthogroups.tsv"
)
DEFAULT_METADATA = PROJECT_ROOT / "config" / "orthogroups_diamond_species_metadata.tsv"
DEFAULT_ANALYSIS_OUTPUT = PROJECT_ROOT / "results" / "presence_absence"
USER_QUERY_SOURCE = Path(r"C:\Users\j2rel\Downloads\known_calc_proteins.zip")
DEFAULT_QUERY_SOURCE = (
    USER_QUERY_SOURCE
    if USER_QUERY_SOURCE.exists()
    else PROJECT_ROOT / "examples" / "known_calc_proteins"
)
DEFAULT_SPECIES_FASTA_DIR = PROJECT_ROOT / "examples" / "proteomes"
DEFAULT_MAPPING_OUTPUT = PROJECT_ROOT / "results" / "known_protein_mapping"


class JobState:
    def __init__(self) -> None:
        self.lock = threading.Lock()
        self.process: subprocess.Popen[str] | None = None
        self.status = "idle"
        self.command = ""
        self.log_lines: list[str] = []
        self.returncode: int | None = None
        self.started_at: float | None = None
        self.ended_at: float | None = None
        self.manifest = ""

    def snapshot(self) -> dict[str, object]:
        with self.lock:
            elapsed = None
            if self.started_at:
                elapsed = round((self.ended_at or time.time()) - self.started_at, 1)
            return {
                "status": self.status,
                "command": self.command,
                "returncode": self.returncode,
                "elapsed_seconds": elapsed,
                "manifest": self.manifest,
                "log": "\n".join(self.log_lines[-500:]),
            }

    def is_running(self) -> bool:
        with self.lock:
            return self.process is not None and self.process.poll() is None

    def append_log(self, line: str) -> None:
        with self.lock:
            self.log_lines.append(line.rstrip("\n"))
            if len(self.log_lines) > 1000:
                self.log_lines = self.log_lines[-1000:]


JOB = JobState()
MAPPING_JOB = JobState()


INDEX_HTML = r"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>OrthoFinder Runner</title>
  <style>
    :root {
      --bg: #f5f7fb;
      --surface: #ffffff;
      --surface-2: #eef4f8;
      --text: #17212b;
      --muted: #607082;
      --line: #cfd8e3;
      --accent: #18736b;
      --accent-2: #22577a;
      --danger: #a23b3b;
      --warn: #8a5a00;
      --ok: #1e6f43;
      --shadow: 0 8px 24px rgba(20, 34, 48, 0.08);
    }

    * { box-sizing: border-box; }

    body {
      margin: 0;
      min-height: 100vh;
      background: var(--bg);
      color: var(--text);
      font-family: Arial, Helvetica, sans-serif;
      font-size: 14px;
      letter-spacing: 0;
    }

    header {
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 16px;
      padding: 18px 24px;
      border-bottom: 1px solid var(--line);
      background: var(--surface);
    }

    h1 {
      margin: 0;
      font-size: 20px;
      font-weight: 700;
    }

    .status-pill {
      min-width: 112px;
      padding: 8px 12px;
      border: 1px solid var(--line);
      border-radius: 6px;
      color: var(--muted);
      text-align: center;
      background: var(--surface-2);
      text-transform: uppercase;
      font-size: 12px;
      font-weight: 700;
    }

    .status-pill.running { color: #17456b; border-color: #9fc6df; background: #e8f4fb; }
    .status-pill.completed { color: var(--ok); border-color: #a8d4b9; background: #edf8f1; }
    .status-pill.failed { color: var(--danger); border-color: #e1b2b2; background: #fff0f0; }

    .tabs {
      display: flex;
      align-items: center;
      gap: 8px;
      flex-wrap: wrap;
    }

    .tab {
      display: inline-flex;
      align-items: center;
      min-height: 36px;
      padding: 8px 12px;
      border-radius: 6px;
      border: 1px solid var(--line);
      color: var(--text);
      background: var(--surface-2);
      font-weight: 700;
      text-decoration: none;
    }

    .tab.active {
      background: var(--accent);
      border-color: var(--accent);
      color: #ffffff;
    }

    main {
      display: grid;
      grid-template-columns: minmax(320px, 460px) minmax(0, 1fr);
      gap: 18px;
      padding: 18px 24px 24px;
    }

    section {
      background: var(--surface);
      border: 1px solid var(--line);
      border-radius: 8px;
      box-shadow: var(--shadow);
      min-width: 0;
    }

    .panel-header {
      padding: 14px 16px;
      border-bottom: 1px solid var(--line);
      font-weight: 700;
    }

    form, .panel-body {
      padding: 16px;
    }

    label {
      display: block;
      margin-bottom: 6px;
      color: var(--muted);
      font-size: 12px;
      font-weight: 700;
      text-transform: uppercase;
    }

    input, select {
      width: 100%;
      height: 40px;
      border: 1px solid var(--line);
      border-radius: 6px;
      padding: 8px 10px;
      color: var(--text);
      background: #ffffff;
      font: inherit;
    }

    .field {
      margin-bottom: 14px;
    }

    .checkbox-field {
      display: flex;
      align-items: center;
      gap: 10px;
      min-height: 40px;
      margin-bottom: 14px;
    }

    .checkbox-field input {
      width: 18px;
      height: 18px;
      margin: 0;
      accent-color: var(--accent);
    }

    .checkbox-field label {
      margin: 0;
      color: var(--text);
      text-transform: none;
      font-size: 14px;
      font-weight: 700;
    }

    .row {
      display: grid;
      grid-template-columns: 1fr 1fr;
      gap: 12px;
    }

    .actions {
      display: flex;
      flex-wrap: wrap;
      gap: 10px;
      padding-top: 4px;
    }

    button {
      min-height: 40px;
      border: 1px solid var(--line);
      border-radius: 6px;
      padding: 9px 13px;
      font: inherit;
      font-weight: 700;
      cursor: pointer;
      background: var(--surface-2);
      color: var(--text);
    }

    button.primary {
      background: var(--accent);
      border-color: var(--accent);
      color: #ffffff;
    }

    button.secondary {
      background: var(--accent-2);
      border-color: var(--accent-2);
      color: #ffffff;
    }

    button.danger {
      color: #ffffff;
      border-color: var(--danger);
      background: var(--danger);
    }

    button:disabled {
      opacity: 0.5;
      cursor: not-allowed;
    }

    table {
      width: 100%;
      border-collapse: collapse;
      table-layout: fixed;
    }

    th, td {
      padding: 10px 8px;
      border-bottom: 1px solid var(--line);
      text-align: left;
      vertical-align: top;
      word-wrap: break-word;
    }

    th {
      color: var(--muted);
      font-size: 12px;
      text-transform: uppercase;
    }

    .message {
      min-height: 38px;
      margin-top: 14px;
      padding: 10px 12px;
      border-radius: 6px;
      border: 1px solid var(--line);
      background: var(--surface-2);
      color: var(--muted);
      line-height: 1.35;
    }

    .message.error {
      color: var(--danger);
      border-color: #e1b2b2;
      background: #fff0f0;
    }

    .message.ok {
      color: var(--ok);
      border-color: #a8d4b9;
      background: #edf8f1;
    }

    .message.warn {
      color: var(--warn);
      border-color: #dec787;
      background: #fff8e5;
    }

    .console {
      height: 310px;
      overflow: auto;
      margin: 0;
      padding: 12px;
      border-radius: 6px;
      border: 1px solid #101821;
      background: #101821;
      color: #dce7f2;
      font: 13px Consolas, "Courier New", monospace;
      white-space: pre-wrap;
    }

    .command {
      min-height: 54px;
      margin-bottom: 14px;
      padding: 10px 12px;
      border-radius: 6px;
      border: 1px solid var(--line);
      background: #f8fafc;
      font: 13px Consolas, "Courier New", monospace;
      white-space: pre-wrap;
      word-break: break-word;
    }

    .meta {
      display: flex;
      flex-wrap: wrap;
      gap: 10px;
      margin-bottom: 14px;
      color: var(--muted);
    }

    .meta span {
      padding: 6px 8px;
      border-radius: 6px;
      background: var(--surface-2);
      border: 1px solid var(--line);
    }

    @media (max-width: 900px) {
      header {
        align-items: flex-start;
        flex-direction: column;
      }

      main {
        grid-template-columns: 1fr;
        padding: 12px;
      }

      .row {
        grid-template-columns: 1fr;
      }

      .actions button {
        flex: 1 1 150px;
      }
    }
  </style>
</head>
<body>
  <header>
    <h1>OrthoFinder Runner</h1>
    <nav class="tabs" aria-label="Workflow tabs">
      <a class="tab" href="/">Presence/Absence</a>
      <a class="tab" href="/mapping">Known Proteins</a>
      <a class="tab active" href="/orthofinder">OrthoFinder</a>
    </nav>
    <div id="statusPill" class="status-pill">Idle</div>
  </header>

  <main>
    <section>
      <div class="panel-header">Run Setup</div>
      <form id="runForm">
        <div class="field">
          <label for="inputDir">Species FASTA Folder</label>
          <input id="inputDir" name="inputDir" value="__DEFAULT_INPUT__" autocomplete="off">
        </div>
        <div class="field">
          <label for="outputDir">Output Folder</label>
          <input id="outputDir" name="outputDir" value="__DEFAULT_OUTPUT__" autocomplete="off">
        </div>
        <div class="field">
          <label for="orthofinderBin">OrthoFinder Executable</label>
          <input id="orthofinderBin" name="orthofinderBin" value="orthofinder" autocomplete="off">
        </div>
        <div class="row">
          <div class="field">
            <label for="threads">Search Threads</label>
            <input id="threads" name="threads" type="number" min="1" max="128" value="8">
          </div>
          <div class="field">
            <label for="analysisThreads">Analysis Threads</label>
            <input id="analysisThreads" name="analysisThreads" type="number" min="1" max="128" value="8">
          </div>
        </div>
        <div class="field">
          <label for="searchProgram">Search Program</label>
          <select id="searchProgram" name="searchProgram">
            <option value="auto">Auto</option>
            <option value="diamond">DIAMOND</option>
            <option value="mmseqs">MMseqs2</option>
            <option value="blast">BLAST</option>
          </select>
        </div>
        <div class="checkbox-field">
          <input id="orthogroupsOnly" name="orthogroupsOnly" type="checkbox" checked>
          <label for="orthogroupsOnly">Stop after orthogroups</label>
        </div>
        <div class="actions">
          <button type="button" id="validateBtn">Validate</button>
          <button type="button" id="dryRunBtn" class="secondary">Build Command</button>
          <button type="submit" id="runBtn" class="primary">Run</button>
          <button type="button" id="stopBtn" class="danger">Stop</button>
        </div>
        <div id="message" class="message">Ready.</div>
      </form>
    </section>

    <section>
      <div class="panel-header">Inputs</div>
      <div class="panel-body">
        <div class="meta" id="meta"></div>
        <table>
          <thead>
            <tr>
              <th style="width: 24%">Species</th>
              <th style="width: 18%">Sequences</th>
              <th style="width: 18%">Residues</th>
              <th>Warnings</th>
            </tr>
          </thead>
          <tbody id="fastaTable">
            <tr><td colspan="4">No validation run yet.</td></tr>
          </tbody>
        </table>
      </div>
    </section>

    <section style="grid-column: 1 / -1;">
      <div class="panel-header">OrthoFinder Job</div>
      <div class="panel-body">
        <div class="command" id="commandBox">No command built yet.</div>
        <pre class="console" id="logBox"></pre>
      </div>
    </section>
  </main>

  <script>
    const form = document.getElementById("runForm");
    const message = document.getElementById("message");
    const fastaTable = document.getElementById("fastaTable");
    const statusPill = document.getElementById("statusPill");
    const commandBox = document.getElementById("commandBox");
    const logBox = document.getElementById("logBox");
    const meta = document.getElementById("meta");
    const runBtn = document.getElementById("runBtn");
    const stopBtn = document.getElementById("stopBtn");

    function payload() {
      return {
        input_dir: document.getElementById("inputDir").value.trim(),
        output_dir: document.getElementById("outputDir").value.trim(),
        orthofinder_bin: document.getElementById("orthofinderBin").value.trim() || "orthofinder",
        threads: Number(document.getElementById("threads").value || 1),
        analysis_threads: Number(document.getElementById("analysisThreads").value || 1),
        search_program: document.getElementById("searchProgram").value,
        orthogroups_only: document.getElementById("orthogroupsOnly").checked
      };
    }

    function setMessage(text, kind = "") {
      message.className = "message" + (kind ? " " + kind : "");
      message.textContent = text;
    }

    function renderValidation(data) {
      const rows = data.files || [];
      meta.innerHTML = `<span>${rows.length} FASTA file(s)</span><span>${data.total_sequences || 0} sequences</span><span>${data.total_residues || 0} residues</span>`;
      if (!rows.length) {
        fastaTable.innerHTML = '<tr><td colspan="4">No FASTA files found.</td></tr>';
        return;
      }
      fastaTable.innerHTML = rows.map(row => `
        <tr>
          <td>${escapeHtml(row.species_id)}</td>
          <td>${row.sequence_count}</td>
          <td>${row.residue_count}</td>
          <td>${escapeHtml(row.warnings || "")}</td>
        </tr>
      `).join("");
    }

    function escapeHtml(value) {
      return String(value).replace(/[&<>"']/g, char => ({
        "&": "&amp;",
        "<": "&lt;",
        ">": "&gt;",
        '"': "&quot;",
        "'": "&#039;"
      }[char]));
    }

    async function postJson(url, body) {
      const response = await fetch(url, {
        method: "POST",
        headers: {"Content-Type": "application/json"},
        body: JSON.stringify(body)
      });
      const data = await response.json();
      if (!response.ok) {
        throw new Error(data.error || `Request failed with ${response.status}`);
      }
      return data;
    }

    document.getElementById("validateBtn").addEventListener("click", async () => {
      try {
        const data = await postJson("/api/validate", payload());
        renderValidation(data);
        setMessage("Inputs validated.", data.has_warnings ? "warn" : "ok");
      } catch (error) {
        setMessage(error.message, "error");
      }
    });

    document.getElementById("dryRunBtn").addEventListener("click", async () => {
      try {
        const data = await postJson("/api/command", payload());
        renderValidation(data.validation);
        commandBox.textContent = data.command;
        setMessage("Command built.", data.validation.has_warnings ? "warn" : "ok");
      } catch (error) {
        setMessage(error.message, "error");
      }
    });

    form.addEventListener("submit", async event => {
      event.preventDefault();
      try {
        const data = await postJson("/api/run", payload());
        commandBox.textContent = data.command;
        setMessage("OrthoFinder job started.", "ok");
        await refreshStatus();
      } catch (error) {
        setMessage(error.message, "error");
      }
    });

    stopBtn.addEventListener("click", async () => {
      try {
        await postJson("/api/stop", {});
        setMessage("Stop requested.", "warn");
        await refreshStatus();
      } catch (error) {
        setMessage(error.message, "error");
      }
    });

    async function refreshStatus() {
      const response = await fetch("/api/status");
      const data = await response.json();
      statusPill.textContent = data.status || "idle";
      statusPill.className = "status-pill " + (data.status || "idle");
      if (data.command) commandBox.textContent = data.command;
      logBox.textContent = data.log || "";
      if (data.elapsed_seconds !== null) {
        const code = data.returncode === null ? "" : ` | exit ${data.returncode}`;
        statusPill.title = `${data.elapsed_seconds}s${code}`;
      }
      const running = data.status === "running";
      runBtn.disabled = running;
      stopBtn.disabled = !running;
      if (running) {
        logBox.scrollTop = logBox.scrollHeight;
      }
    }

    setInterval(refreshStatus, 2000);
    refreshStatus();
  </script>
</body>
</html>
"""


def json_response(handler: BaseHTTPRequestHandler, status: int, payload: dict[str, object]) -> None:
    body = json.dumps(payload, indent=2).encode("utf-8")
    handler.send_response(status)
    handler.send_header("Content-Type", "application/json; charset=utf-8")
    handler.send_header("Content-Length", str(len(body)))
    handler.end_headers()
    handler.wfile.write(body)


def read_json(handler: BaseHTTPRequestHandler) -> dict[str, object]:
    length = int(handler.headers.get("Content-Length", "0"))
    if length == 0:
        return {}
    raw = handler.rfile.read(length)
    return json.loads(raw.decode("utf-8"))


def analysis_html() -> bytes:
    path = PROJECT_ROOT / "gui" / "orthogroup_analysis.html"
    html = path.read_text(encoding="utf-8")
    html = (
        html.replace("__DEFAULT_OG_TSV__", str(DEFAULT_OG_TSV))
        .replace("__DEFAULT_METADATA__", str(DEFAULT_METADATA))
        .replace("__DEFAULT_ANALYSIS_OUTPUT__", str(DEFAULT_ANALYSIS_OUTPUT))
    )
    return html.encode("utf-8")


def mapping_html() -> bytes:
    path = PROJECT_ROOT / "gui" / "known_protein_mapping.html"
    html = path.read_text(encoding="utf-8")
    html = (
        html.replace("__DEFAULT_OG_TSV__", str(DEFAULT_OG_TSV))
        .replace("__DEFAULT_SPECIES_FASTA_DIR__", str(DEFAULT_SPECIES_FASTA_DIR))
        .replace("__DEFAULT_QUERY_SOURCE__", str(DEFAULT_QUERY_SOURCE))
        .replace("__DEFAULT_MAPPING_OUTPUT__", str(DEFAULT_MAPPING_OUTPUT))
    )
    return html.encode("utf-8")


def validation_payload(input_dir: Path) -> dict[str, object]:
    summaries = validate_fasta_dir(input_dir)
    files = [
        {
            "species_id": summary.species_id,
            "path": str(summary.path),
            "sequence_count": summary.sequence_count,
            "residue_count": summary.residue_count,
            "warnings": "; ".join(summary.warnings),
        }
        for summary in summaries
    ]
    return {
        "files": files,
        "total_sequences": sum(summary.sequence_count for summary in summaries),
        "total_residues": sum(summary.residue_count for summary in summaries),
        "has_warnings": any(summary.warnings for summary in summaries),
    }


def parse_run_payload(data: dict[str, object]) -> dict[str, object]:
    input_dir = Path(str(data.get("input_dir", ""))).expanduser().resolve()
    output_dir = Path(str(data.get("output_dir", ""))).expanduser().resolve()
    orthofinder_bin = str(data.get("orthofinder_bin", "orthofinder") or "orthofinder")
    threads = int(data.get("threads", 8) or 8)
    analysis_threads = int(data.get("analysis_threads", threads) or threads)
    search_program = str(data.get("search_program", "auto") or "auto")
    orthogroups_only = bool(data.get("orthogroups_only", True))
    return {
        "input_dir": input_dir,
        "output_dir": output_dir,
        "orthofinder_bin": orthofinder_bin,
        "threads": threads,
        "analysis_threads": analysis_threads,
        "search_program": search_program,
        "orthogroups_only": orthogroups_only,
    }


def command_payload(data: dict[str, object]) -> tuple[dict[str, object], list[str], Path]:
    parsed = parse_run_payload(data)
    input_dir = parsed["input_dir"]
    output_dir = parsed["output_dir"]
    validation = validation_payload(input_dir)
    command = build_orthofinder_command(**parsed)
    manifest = write_manifest(validate_fasta_dir(input_dir), output_dir)
    return validation, command, manifest


def parse_og_payload(data: dict[str, object]) -> dict[str, object]:
    return {
        "orthogroups_path": Path(str(data.get("orthogroups_path", ""))).expanduser().resolve(),
        "metadata_path": Path(str(data.get("metadata_path", ""))).expanduser().resolve(),
        "output_dir": Path(str(data.get("output_dir", DEFAULT_ANALYSIS_OUTPUT))).expanduser().resolve(),
        "direction": str(data.get("direction", "calcifying") or "calcifying"),
        "min_target_present": int(data.get("min_target_present", 0) or 0),
        "max_background_present": int(data.get("max_background_present", 0) or 0),
    }


def og_filter_payload(data: dict[str, object]) -> tuple[list[dict[str, object]], dict[str, object], dict[str, object]]:
    parsed = parse_og_payload(data)
    rows, summary = filter_orthogroups(
        orthogroups_path=parsed["orthogroups_path"],
        metadata_path=parsed["metadata_path"],
        direction=str(parsed["direction"]),
        min_target_present=int(parsed["min_target_present"]),
        max_background_present=int(parsed["max_background_present"]),
    )
    return rows, summary, parsed


def parse_mapping_payload(data: dict[str, object]) -> dict[str, object]:
    return {
        "orthogroups_path": Path(str(data.get("orthogroups_path", ""))).expanduser().resolve(),
        "species_fasta_dir": Path(str(data.get("species_fasta_dir", ""))).expanduser().resolve(),
        "query_source": Path(str(data.get("query_source", ""))).expanduser().resolve(),
        "output_dir": Path(str(data.get("output_dir", DEFAULT_MAPPING_OUTPUT))).expanduser().resolve(),
        "diamond_bin": str(data.get("diamond_bin", "diamond") or "diamond"),
        "threads": int(data.get("threads", 8) or 8),
        "evalue": str(data.get("evalue", "1e-5") or "1e-5"),
        "max_target_seqs": int(data.get("max_target_seqs", 5) or 5),
        "sensitive": bool(data.get("sensitive", True)),
        "recursive_species": bool(data.get("recursive_species", False)),
    }


def mapping_summary_payload(data: dict[str, object]) -> dict[str, object]:
    parsed = parse_mapping_payload(data)
    return summarize_mapping_inputs(
        orthogroups_tsv=parsed["orthogroups_path"],
        species_fasta_dir=parsed["species_fasta_dir"],
        query_source=parsed["query_source"],
        diamond_bin=str(parsed["diamond_bin"]),
        recursive_species=bool(parsed["recursive_species"]),
    )


def mapping_command(parsed: dict[str, object]) -> list[str]:
    command = [
        sys.executable,
        "scripts/known_protein_mapping.py",
        "--orthogroups",
        str(parsed["orthogroups_path"]),
        "--species-fasta-dir",
        str(parsed["species_fasta_dir"]),
        "--queries",
        str(parsed["query_source"]),
        "--output-dir",
        str(parsed["output_dir"]),
        "--diamond-bin",
        str(parsed["diamond_bin"]),
        "--threads",
        str(parsed["threads"]),
        "--evalue",
        str(parsed["evalue"]),
        "--max-target-seqs",
        str(parsed["max_target_seqs"]),
    ]
    if not parsed["sensitive"]:
        command.append("--fast")
    if parsed["recursive_species"]:
        command.append("--recursive-species")
    return command


def start_job(command: list[str], manifest: Path) -> None:
    with JOB.lock:
        JOB.status = "running"
        JOB.command = command_display(command)
        JOB.log_lines = [f"$ {JOB.command}", f"Input manifest: {manifest}"]
        JOB.returncode = None
        JOB.started_at = time.time()
        JOB.ended_at = None
        JOB.manifest = str(manifest)
        JOB.process = subprocess.Popen(
            command,
            cwd=str(PROJECT_ROOT),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )

    def watch() -> None:
        process = JOB.process
        if process and process.stdout:
            for line in process.stdout:
                JOB.append_log(line)
        returncode = process.wait() if process else 1
        with JOB.lock:
            JOB.returncode = returncode
            JOB.ended_at = time.time()
            JOB.status = "completed" if returncode == 0 else "failed"
            JOB.process = None
        JOB.append_log(f"Process exited with code {returncode}.")

    threading.Thread(target=watch, daemon=True).start()


def start_subprocess_job(job: JobState, command: list[str], initial_lines: list[str]) -> None:
    with job.lock:
        job.status = "running"
        job.command = command_display(command)
        job.log_lines = [f"$ {job.command}", *initial_lines]
        job.returncode = None
        job.started_at = time.time()
        job.ended_at = None
        job.manifest = ""
        job.process = subprocess.Popen(
            command,
            cwd=str(PROJECT_ROOT),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )

    def watch() -> None:
        process = job.process
        if process and process.stdout:
            for line in process.stdout:
                job.append_log(line)
        returncode = process.wait() if process else 1
        with job.lock:
            job.returncode = returncode
            job.ended_at = time.time()
            job.status = "completed" if returncode == 0 else "failed"
            job.process = None
        job.append_log(f"Process exited with code {returncode}.")

    threading.Thread(target=watch, daemon=True).start()


class Handler(BaseHTTPRequestHandler):
    def log_message(self, format: str, *args: object) -> None:
        return

    def do_GET(self) -> None:
        path = urlparse(self.path).path
        if path in {"/", "/analysis"}:
            html = analysis_html()
            self.send_response(HTTPStatus.OK)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.send_header("Content-Length", str(len(html)))
            self.end_headers()
            self.wfile.write(html)
            return
        if path == "/mapping":
            html = mapping_html()
            self.send_response(HTTPStatus.OK)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.send_header("Content-Length", str(len(html)))
            self.end_headers()
            self.wfile.write(html)
            return
        if path == "/orthofinder":
            html = (
                INDEX_HTML.replace("__DEFAULT_INPUT__", str(DEFAULT_INPUT))
                .replace("__DEFAULT_OUTPUT__", str(DEFAULT_OUTPUT))
                .encode("utf-8")
            )
            self.send_response(HTTPStatus.OK)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.send_header("Content-Length", str(len(html)))
            self.end_headers()
            self.wfile.write(html)
            return
        if path == "/api/status":
            json_response(self, HTTPStatus.OK, JOB.snapshot())
            return
        if path == "/api/mapping/status":
            json_response(self, HTTPStatus.OK, MAPPING_JOB.snapshot())
            return
        json_response(self, HTTPStatus.NOT_FOUND, {"error": "Not found"})

    def do_POST(self) -> None:
        path = urlparse(self.path).path
        try:
            data = read_json(self)
            if path == "/api/og/summary":
                parsed = parse_og_payload(data)
                json_response(
                    self,
                    HTTPStatus.OK,
                    summarize(parsed["orthogroups_path"], parsed["metadata_path"]),
                )
                return
            if path == "/api/agent/chat":
                message = str(data.get("message", ""))
                context = data.get("context", {})
                if not isinstance(context, dict):
                    context = {}
                json_response(self, HTTPStatus.OK, agent_chat(message, context))
                return
            if path == "/api/og/filter":
                rows, summary, _parsed = og_filter_payload(data)
                json_response(
                    self,
                    HTTPStatus.OK,
                    {
                        "summary": summary,
                        "rows": rows[:250],
                        "row_limit": 250,
                    },
                )
                return
            if path == "/api/og/export":
                rows, summary, parsed = og_filter_payload(data)
                paths = write_results(rows, summary, parsed["output_dir"])
                json_response(
                    self,
                    HTTPStatus.OK,
                    {
                        "summary": summary,
                        "rows": rows[:250],
                        "row_limit": 250,
                        "paths": paths,
                    },
                )
                return
            if path == "/api/og/report":
                rows, summary, parsed = og_filter_payload(data)
                paths = write_run_report(
                    rows,
                    summary,
                    parsed["output_dir"],
                    parsed["orthogroups_path"],
                    parsed["metadata_path"],
                    str(parsed["direction"]),
                    int(parsed["min_target_present"]),
                    int(parsed["max_background_present"]),
                )
                json_response(
                    self,
                    HTTPStatus.OK,
                    {
                        "summary": summary,
                        "rows": rows[:250],
                        "row_limit": 250,
                        "paths": paths,
                    },
                )
                return
            if path == "/api/mapping/summary":
                json_response(self, HTTPStatus.OK, mapping_summary_payload(data))
                return
            if path == "/api/mapping/command":
                parsed = parse_mapping_payload(data)
                summary = mapping_summary_payload(data)
                json_response(
                    self,
                    HTTPStatus.OK,
                    {
                        "summary": summary,
                        "command": command_display(mapping_command(parsed)),
                    },
                )
                return
            if path == "/api/mapping/run":
                if MAPPING_JOB.is_running():
                    json_response(self, HTTPStatus.CONFLICT, {"error": "A mapping job is already running."})
                    return
                parsed = parse_mapping_payload(data)
                summary = mapping_summary_payload(data)
                if not summary.get("diamond_found"):
                    json_response(
                        self,
                        HTTPStatus.BAD_REQUEST,
                        {"error": "DIAMOND was not found. Install DIAMOND or provide its executable path."},
                    )
                    return
                if summary.get("species_error") or int(summary.get("species_fasta_count", 0)) == 0:
                    json_response(
                        self,
                        HTTPStatus.BAD_REQUEST,
                        {"error": summary.get("species_error") or "No species FASTA files found."},
                    )
                    return
                command = mapping_command(parsed)
                start_subprocess_job(
                    MAPPING_JOB,
                    command,
                    [
                        f"Orthogroups: {summary.get('orthogroup_count', 0)}",
                        f"Queries: {summary.get('query_count', 0)}",
                        f"Species FASTAs: {summary.get('species_fasta_count', 0)}",
                    ],
                )
                json_response(
                    self,
                    HTTPStatus.OK,
                    {"summary": summary, "command": command_display(command)},
                )
                return
            if path == "/api/mapping/stop":
                with MAPPING_JOB.lock:
                    process = MAPPING_JOB.process
                if not process or process.poll() is not None:
                    json_response(self, HTTPStatus.BAD_REQUEST, {"error": "No running mapping job."})
                    return
                process.terminate()
                MAPPING_JOB.append_log("Stop requested by user.")
                json_response(self, HTTPStatus.OK, {"status": "stopping"})
                return
            if path == "/api/validate":
                parsed = parse_run_payload(data)
                json_response(self, HTTPStatus.OK, validation_payload(parsed["input_dir"]))
                return
            if path == "/api/command":
                validation, command, manifest = command_payload(data)
                json_response(
                    self,
                    HTTPStatus.OK,
                    {
                        "validation": validation,
                        "command": command_display(command),
                        "manifest": str(manifest),
                    },
                )
                return
            if path == "/api/run":
                if JOB.is_running():
                    json_response(self, HTTPStatus.CONFLICT, {"error": "A job is already running."})
                    return
                validation, command, manifest = command_payload(data)
                parsed = parse_run_payload(data)
                if not check_orthofinder_available(str(parsed["orthofinder_bin"])):
                    json_response(
                        self,
                        HTTPStatus.BAD_REQUEST,
                        {
                            "error": (
                                "OrthoFinder was not found. Install it or set the executable "
                                "field to the full orthofinder path."
                            )
                        },
                    )
                    return
                start_job(command, manifest)
                json_response(
                    self,
                    HTTPStatus.OK,
                    {
                        "validation": validation,
                        "command": command_display(command),
                        "manifest": str(manifest),
                    },
                )
                return
            if path == "/api/stop":
                with JOB.lock:
                    process = JOB.process
                if not process or process.poll() is not None:
                    json_response(self, HTTPStatus.BAD_REQUEST, {"error": "No running job."})
                    return
                process.terminate()
                JOB.append_log("Stop requested by user.")
                json_response(self, HTTPStatus.OK, {"status": "stopping"})
                return
            json_response(self, HTTPStatus.NOT_FOUND, {"error": "Not found"})
        except Exception as exc:
            json_response(self, HTTPStatus.BAD_REQUEST, {"error": str(exc)})


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8765)
    args = parser.parse_args()

    server = ThreadingHTTPServer((args.host, args.port), Handler)
    print(f"OrthoFinder GUI: http://{args.host}:{args.port}")
    server.serve_forever()


if __name__ == "__main__":
    main()
