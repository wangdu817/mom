#!/usr/bin/env python3
"""M1 OH-oxidation Negative Test (red/green proof).

Zeros omega_gas[OH] and omega_gas[CO] for fixture B (simulating no gas
consumption/production), then runs ohox_closure_analysis.py on mutated JSON.
Closure must FAIL.

Exit codes: 0=PASS (closure rejected), 1=FAIL (closure not rejected), 2=ERROR
"""

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="M1 OH-oxidation Negative Test")
    parser.add_argument("--input", type=str, required=True)
    args = parser.parse_args()

    fixture_path = Path(args.input).resolve()
    if not fixture_path.exists():
        print(f"ERROR: JSON not found: {fixture_path}", file=sys.stderr)
        return 2

    script_dir = Path(__file__).resolve().parent
    closure_script = script_dir / "ohox_closure_analysis.py"

    with fixture_path.open("r", encoding="utf-8") as fh:
        data = json.load(fh)

    species = data["fixtures"]["B_oh_ox_above_threshold"]["inputs"]["species_names"]
    iOH = species.index("OH")
    iCO = species.index("CO")

    orig_oh = data["fixtures"]["B_oh_ox_above_threshold"]["outputs"]["omega_gas_kg_m3_s"][iOH]
    orig_co = data["fixtures"]["B_oh_ox_above_threshold"]["outputs"]["omega_gas_kg_m3_s"][iCO]
    data["fixtures"]["B_oh_ox_above_threshold"]["outputs"]["omega_gas_kg_m3_s"][iOH] = 0.0
    data["fixtures"]["B_oh_ox_above_threshold"]["outputs"]["omega_gas_kg_m3_s"][iCO] = 0.0

    print(f"Original omega[OH]={orig_oh}, omega[CO]={orig_co}")
    print(f"Mutated to 0.0 -- simulating no gas consumption/production")
    print(f"Expected: mass closure FAILS")

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".json", dir=fixture_path.parent, delete=False, encoding="utf-8"
    ) as tmp:
        json.dump(data, tmp, indent=2)
        tmp_path = Path(tmp.name)

    try:
        result = subprocess.run(
            [sys.executable, str(closure_script), "--input", str(tmp_path)],
            capture_output=True, text=True,
        )
        print(f"\nClosure exit code on mutated JSON: {result.returncode}")
        if result.returncode != 0:
            print("PASS: OHClosure correctly rejected mutated data.")
            return 0
        else:
            print("FAIL: OHClosure did NOT reject mutated data!")
            return 1
    finally:
        tmp_path.unlink(missing_ok=True)


if __name__ == "__main__":
    sys.exit(main())
