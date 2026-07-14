#!/usr/bin/env python3
"""M1 O2-oxidation Negative Test (red/green proof).

Zeros out omega_gas[O2] and omega_gas[CO] for fixture B (simulating a bug
where oxidation gas consumption is not applied), then runs closure analysis.
Must fail (exit 1).

Exit codes: 0=PASS (correctly rejected), 1=FAIL (not rejected), 2=ERROR
"""

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="M1 O2-oxidation Negative Test")
    parser.add_argument("--input", type=str, required=True)
    args = parser.parse_args()

    fixture_path = Path(args.input).resolve()
    if not fixture_path.exists():
        print(f"ERROR: JSON not found: {fixture_path}", file=sys.stderr)
        return 2

    script_dir = Path(__file__).resolve().parent
    closure_script = script_dir / "o2ox_closure_analysis.py"

    with fixture_path.open("r", encoding="utf-8") as fh:
        data = json.load(fh)

    species = data["fixtures"]["B_ox_above_threshold"]["inputs"]["species_names"]
    iO2 = species.index("O2")
    iCO = species.index("CO")

    # Zero out O2 and CO gas rates for fixture B (simulate no gas consumption)
    omega = data["fixtures"]["B_ox_above_threshold"]["outputs"]["omega_gas_kg_m3_s"]
    orig_o2 = omega[iO2]
    orig_co = omega[iCO]
    omega[iO2] = 0.0
    omega[iCO] = 0.0
    data["fixtures"]["B_ox_above_threshold"]["outputs"]["o2_mass_rate_kg_m3_s"] = 0.0
    data["fixtures"]["B_ox_above_threshold"]["outputs"]["co_mass_rate_kg_m3_s"] = 0.0

    print(f"Original omega[O2]={orig_o2}, omega[CO]={orig_co}")
    print(f"Mutated to 0.0, 0.0 (simulating no gas consumption)")
    print(f"Expected: mass closure FAILS (soot loses mass but no gas gain)")

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
        print(f"\nClosure analysis exit code on mutated JSON: {result.returncode}")
        if result.returncode != 0:
            print("PASS: O2ox closure correctly rejected mutated data.")
            return 0
        else:
            print("FAIL: O2ox closure did NOT reject mutated data!")
            return 1
    finally:
        tmp_path.unlink(missing_ok=True)


if __name__ == "__main__":
    sys.exit(main())
