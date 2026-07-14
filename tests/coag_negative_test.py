#!/usr/bin/env python3
"""Negative test for coagulation closure.

Takes the correct fixture JSON, injects a fake non-zero S_M10 (mass source)
and a fake non-zero omega_gas into fixture B, then runs the closure analysis.
The closure analysis MUST exit non-zero (NOT_CLOSED).

This proves the gate can detect mass-conservation violations in coagulation.
"""

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description="Coagulation negative test")
    parser.add_argument("--input", required=True, help="Path to coagulation_fixture_results.json")
    args = parser.parse_args()

    json_path = Path(args.input)
    if not json_path.exists():
        print(f"ERROR: JSON not found: {json_path}")
        return 2

    with open(json_path) as f:
        data = json.load(f)

    for fx in data["fixtures"]:
        if fx["label"] == "B_coag_nominal":
            fx["src_coag"][1] = 1.0e-10   # inject fake S_M10 != 0
            fx["omega_gas"][0] = 1.0e-15   # inject fake gas source
            fx["gas_sum"] = sum(fx["omega_gas"])
            fx["soot_mass_rate"] = fx["src_coag"][1] * fx["V0_m3"] * 6.02214076e23 * fx["rho_soot"]
            fx["mass_residual"] = fx["soot_mass_rate"] + fx["gas_sum"]
            print(f"Mutated B: S_M10={fx['src_coag'][1]:.6e}, gas_sum={fx['gas_sum']:.6e}")
            break

    script_dir = Path(__file__).parent
    closure_script = script_dir / "coag_closure_analysis.py"

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
        print(f"CoagClosure exit code on mutated JSON: {result.returncode}")
        if result.returncode != 0:
            print("PASS: CoagClosure correctly rejected mutated (wrong) data.")
            return 0
        else:
            print("FAIL: CoagClosure accepted mutated data (should have failed)!")
            return 1
    finally:
        Path(tmp_path).unlink(missing_ok=True)


if __name__ == "__main__":
    sys.exit(main())
