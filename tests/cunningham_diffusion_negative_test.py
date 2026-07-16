#!/usr/bin/env python3
"""Negative test for Cunningham diffusion closure.

Takes the correct fixture JSON, corrupts the diffusion_coefficient in
fixture A (multiply by 2.0), then runs the closure analysis.
The closure analysis MUST exit non-zero (NOT_CLOSED) because the
mutated Gamma will not match the independent recomputation.

This proves the gate can detect formula violations in the Cunningham
diffusion coefficient.
"""

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description="Cunningham diffusion negative test")
    parser.add_argument("--input", required=True, help="Path to cunningham_diffusion_fixture_results.json")
    args = parser.parse_args()

    json_path = Path(args.input)
    if not json_path.exists():
        print(f"ERROR: JSON not found: {json_path}")
        return 2

    with open(json_path) as f:
        data = json.load(f)

    for fx in data["fixtures"]:
        if fx["label"] == "A_flame_small_soot":
            original = fx["diffusion_coefficient"]
            fx["diffusion_coefficient"] = original * 2.0
            print(f"Mutated A: diffusion_coefficient {original:.6e} -> {fx['diffusion_coefficient']:.6e}")
            break

    script_dir = Path(__file__).parent
    closure_script = script_dir / "cunningham_diffusion_analysis.py"

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".json", dir=json_path.parent, delete=False
    ) as tmp:
        json.dump(data, tmp)
        tmp_path = tmp.name

    try:
        result = subprocess.run(
            [sys.executable, str(closure_script), "--input", tmp_path],
            capture_output=True, text=True
        )
        print(f"CunninghamClosure exit code on mutated JSON: {result.returncode}")
        if result.returncode != 0:
            print("PASS: CunninghamClosure correctly rejected mutated (wrong) data.")
            return 0
        else:
            print("FAIL: CunninghamClosure accepted mutated data (should have failed)!")
            return 1
    finally:
        Path(tmp_path).unlink(missing_ok=True)


if __name__ == "__main__":
    sys.exit(main())
