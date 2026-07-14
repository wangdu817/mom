#!/usr/bin/env python3
"""M1 Condensation Negative Test (red/green proof).

Loads the condensation fixture JSON, zeros out omega_gas[C2H2] for fixture B
(simulating the pre-fix state where PAH gas consumption was gated by
nucleation_model_), then runs cond_closure_analysis.py on the mutated JSON.
The closure analysis must FAIL (exit code 1) on the mutated data.

This proves the closure analysis is not trivially passing -- it genuinely
detects the mass imbalance bug.

Exit codes:
  0 = PASS (closure analysis correctly rejected mutated data)
  1 = FAIL (closure analysis did NOT reject mutated data -- problem!)
  2 = Error
"""

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(
        description="M1 Condensation Negative Test (red/green proof)"
    )
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Path to condensation_fixture_results.json",
    )
    args = parser.parse_args()

    fixture_path = Path(args.input).resolve()
    if not fixture_path.exists():
        print(f"ERROR: fixture JSON not found: {fixture_path}", file=sys.stderr)
        return 2

    script_dir = Path(__file__).resolve().parent
    closure_script = script_dir / "cond_closure_analysis.py"

    # Load original JSON.
    with fixture_path.open("r", encoding="utf-8") as fh:
        data = json.load(fh)

    # Mutate: zero out omega_gas[C2H2] for fixture B (simulate pre-fix bug).
    # Find PAH index.
    species_names = data["fixtures"]["B_cond_only_nominal"]["inputs"]["species_names"]
    pah_idx = species_names.index("C2H2")

    original_omega = data["fixtures"]["B_cond_only_nominal"]["outputs"]["omega_gas_kg_m3_s"][pah_idx]
    data["fixtures"]["B_cond_only_nominal"]["outputs"]["omega_gas_kg_m3_s"][pah_idx] = 0.0
    # Also zero out the derived pah_mass_rate to be consistent.
    data["fixtures"]["B_cond_only_nominal"]["outputs"]["pah_mass_rate_kg_m3_s"] = 0.0
    # Recompute mass_residual (should now be non-zero).
    soot_rate = data["fixtures"]["B_cond_only_nominal"]["outputs"]["soot_mass_rate_kg_m3_s"]
    data["fixtures"]["B_cond_only_nominal"]["outputs"]["mass_residual_kg_m3_s"] = soot_rate

    # Also mutate fixture C.
    data["fixtures"]["C_cond_only_perturbed"]["outputs"]["omega_gas_kg_m3_s"][pah_idx] = 0.0
    data["fixtures"]["C_cond_only_perturbed"]["outputs"]["pah_mass_rate_kg_m3_s"] = 0.0
    soot_rate_c = data["fixtures"]["C_cond_only_perturbed"]["outputs"]["soot_mass_rate_kg_m3_s"]
    data["fixtures"]["C_cond_only_perturbed"]["outputs"]["mass_residual_kg_m3_s"] = soot_rate_c

    print(f"Original omega_gas[C2H2] for B: {original_omega}")
    print(f"Mutated omega_gas[C2H2] for B: 0.0 (simulating pre-fix bug)")
    print(f"Expected: mass closure FAILS (soot mass increases but PAH not consumed)")

    # Write mutated JSON to temp file.
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".json", dir=fixture_path.parent, delete=False, encoding="utf-8"
    ) as tmp:
        json.dump(data, tmp, indent=2)
        tmp_path = Path(tmp.name)

    try:
        # Run closure analysis on mutated JSON.
        result = subprocess.run(
            [sys.executable, str(closure_script), "--input", str(tmp_path)],
            capture_output=True,
            text=True,
        )

        print(f"\nClosure analysis exit code on mutated JSON: {result.returncode}")
        if result.returncode != 0:
            print("PASS: CondClosure correctly rejected mutated (wrong) data.")
            return 0
        else:
            print("FAIL: CondClosure did NOT reject mutated data -- closure check is too weak!")
            print(result.stdout)
            return 1
    finally:
        tmp_path.unlink(missing_ok=True)


if __name__ == "__main__":
    sys.exit(main())
