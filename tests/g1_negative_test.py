#!/usr/bin/env python3
"""G1 negative test: proves G1Closure FAILS when OH stoichiometry is wrong.

Takes the correct oh_fixture_results.json, mutates omega[OH] to simulate
the old per-C bug (1 OH per C2 instead of 2), writes a temp JSON, then
runs g1_closure_analysis.py on it and asserts exit code != 0.

This is the "red" half of the red/green proof: if G1Closure ever passes
on this mutated input, the gate is broken.
"""

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(
        description="G1 negative test (red/green proof)"
    )
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Path to oh_fixture_results.json from the current build tree",
    )
    args = parser.parse_args()
    fixture_path = Path(args.input).resolve()

    if not fixture_path.exists():
        print(f"ERROR: fixture JSON not found: {fixture_path}", file=sys.stderr)
        return 2

    with fixture_path.open("r", encoding="utf-8") as fh:
        data = json.load(fh)

    # Mutate: halve OH consumption (simulate old per-C bug: 1 OH per C2).
    # omega[OH] is negative (consumed); halving its magnitude means less OH
    # consumed -> H and O balances will fail.
    original_oh = data["outputs"]["omega_gas_kg_m3_s"]["OH"]
    data["outputs"]["omega_gas_kg_m3_s"]["OH"] = original_oh * 0.5
    print(f"Mutated omega[OH]: {original_oh:.6e} -> {original_oh * 0.5:.6e}")

    script_dir = Path(__file__).resolve().parent

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".json", delete=False, dir=str(fixture_path.parent)
    ) as tmp:
        json.dump(data, tmp, indent=2)
        tmp_path = tmp.name

    try:
        result = subprocess.run(
            [sys.executable, str(script_dir / "g1_closure_analysis.py"),
             "--input", tmp_path],
            capture_output=True, text=True, timeout=30,
        )
        print(f"G1Closure exit code on mutated JSON: {result.returncode}")
        if result.returncode == 0:
            print("FAIL: G1Closure passed on mutated (wrong) data!", file=sys.stderr)
            print(result.stdout[-500:])
            return 1
        else:
            print("PASS: G1Closure correctly rejected mutated (wrong) data.")
            return 0
    finally:
        Path(tmp_path).unlink(missing_ok=True)


if __name__ == "__main__":
    sys.exit(main())
