#!/usr/bin/env python3
"""Smoke-test the local OrthoFinder GUI HTTP endpoints."""

from __future__ import annotations

import json
import threading
import urllib.request
from http.server import ThreadingHTTPServer

from gui.orthofinder_gui import (
    DEFAULT_ANALYSIS_OUTPUT,
    DEFAULT_INPUT,
    DEFAULT_METADATA,
    DEFAULT_OG_TSV,
    DEFAULT_OUTPUT,
    DEFAULT_MAPPING_OUTPUT,
    DEFAULT_QUERY_SOURCE,
    DEFAULT_SPECIES_FASTA_DIR,
    Handler,
)


def request_json(url: str, payload: dict[str, object] | None = None) -> dict[str, object]:
    data = None
    headers = {}
    if payload is not None:
        data = json.dumps(payload).encode("utf-8")
        headers["Content-Type"] = "application/json"
    request = urllib.request.Request(url, data=data, headers=headers)
    with urllib.request.urlopen(request, timeout=10) as response:
        return json.loads(response.read().decode("utf-8"))


def main() -> None:
    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    host, port = server.server_address
    base_url = f"http://{host}:{port}"
    orthofinder_payload = {
        "input_dir": str(DEFAULT_INPUT),
        "output_dir": str(DEFAULT_OUTPUT),
        "orthofinder_bin": "orthofinder",
        "threads": 2,
        "analysis_threads": 2,
        "search_program": "auto",
        "orthogroups_only": True,
    }
    og_payload = {
        "orthogroups_path": str(DEFAULT_OG_TSV),
        "metadata_path": str(DEFAULT_METADATA),
        "output_dir": str(DEFAULT_ANALYSIS_OUTPUT),
        "direction": "calcifying",
        "min_target_present": 6,
        "max_background_present": 0,
    }
    mapping_payload = {
        "orthogroups_path": str(DEFAULT_OG_TSV),
        "species_fasta_dir": str(DEFAULT_SPECIES_FASTA_DIR),
        "query_source": str(DEFAULT_QUERY_SOURCE),
        "output_dir": str(DEFAULT_MAPPING_OUTPUT),
        "diamond_bin": "diamond",
        "threads": 2,
        "evalue": "1e-5",
        "max_target_seqs": 5,
        "sensitive": True,
        "recursive_species": False,
    }

    status = request_json(f"{base_url}/api/status")
    validation = request_json(f"{base_url}/api/validate", orthofinder_payload)
    command = request_json(f"{base_url}/api/command", orthofinder_payload)
    og_summary = request_json(f"{base_url}/api/og/summary", og_payload)
    og_filter = request_json(f"{base_url}/api/og/filter", og_payload)
    mapping_summary = request_json(f"{base_url}/api/mapping/summary", mapping_payload)
    mapping_command = request_json(f"{base_url}/api/mapping/command", mapping_payload)
    server.shutdown()

    assert status["status"] == "idle"
    assert validation["total_sequences"] == 3
    assert command["command"].endswith(" -og")
    assert og_summary["calcifying_count"] == 7
    assert og_summary["non_calcifying_count"] == 20
    assert og_filter["summary"]["candidate_count"] > 0
    assert mapping_summary["query_file_count"] >= 0
    assert "known_protein_mapping.py" in mapping_command["command"]
    print("GUI smoke test passed.")


if __name__ == "__main__":
    main()
