#!/usr/bin/env python3
"""M1 Surface-growth Negative Test (red/green proof).

Zeros omega_gas[C2H2] and omega_gas[H2] for fixture B (simulating pre-fix
state where gas consumption was missing), runs closure analysis, asserts
exit code != 0.
"""
import argparse, json, subprocess, sys, tempfile
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="M1 SG Negative Test")
    parser.add_argument("--input", type=str, required=True)
    args = parser.parse_args()
    fp = Path(args.input).resolve()
    if not fp.exists(): print(f"ERROR: {fp}", file=sys.stderr); return 2
    sd = Path(__file__).resolve().parent
    cs = sd / "sg_closure_analysis.py"

    with fp.open() as f: data = json.load(f)
    sp = data["fixtures"]["B_sg_above_threshold"]["inputs"]["species_names"]
    iC2H2 = sp.index("C2H2"); iH2 = sp.index("H2")

    # Zero out gas consumption for B
    data["fixtures"]["B_sg_above_threshold"]["outputs"]["omega_gas_kg_m3_s"][iC2H2] = 0.0
    data["fixtures"]["B_sg_above_threshold"]["outputs"]["omega_gas_kg_m3_s"][iH2] = 0.0
    data["fixtures"]["B_sg_above_threshold"]["outputs"]["c2h2_mass_rate_kg_m3_s"] = 0.0
    data["fixtures"]["B_sg_above_threshold"]["outputs"]["h2_mass_rate_kg_m3_s"] = 0.0
    sr = data["fixtures"]["B_sg_above_threshold"]["outputs"]["soot_mass_rate_kg_m3_s"]
    data["fixtures"]["B_sg_above_threshold"]["outputs"]["mass_residual_kg_m3_s"] = sr

    print(f"Zeroed omega_gas[C2H2] and omega_gas[H2] for fixture B")
    print(f"Expected: mass closure FAILS (soot mass increases but no gas consumed)")

    with tempfile.NamedTemporaryFile(mode="w", suffix=".json", dir=fp.parent, delete=False, encoding="utf-8") as tmp:
        json.dump(data, tmp, indent=2); tp = Path(tmp.name)
    try:
        r = subprocess.run([sys.executable, str(cs), "--input", str(tp)], capture_output=True, text=True)
        print(f"\nClosure exit code on mutated JSON: {r.returncode}")
        if r.returncode != 0:
            print("PASS: SGClosure correctly rejected mutated data."); return 0
        else:
            print("FAIL: SGClosure did NOT reject mutated data."); return 1
    finally:
        tp.unlink(missing_ok=True)

if __name__=="__main__": sys.exit(main())
