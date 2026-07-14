#!/usr/bin/env python3
"""M1 Nucleation Negative Test (red/green proof).

Zeros omega_gas[C2H2] for fixtures B and C (simulating a bug where PAH
gas consumption is not computed), then runs nuc_closure_analysis.py.
The closure analysis must FAIL (exit code 1).
"""

import argparse, json, subprocess, sys, tempfile
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="M1 Nucleation Negative Test")
    parser.add_argument("--input", type=str, required=True)
    args = parser.parse_args()
    fp = Path(args.input).resolve()
    if not fp.exists(): print(f"ERROR: {fp}",file=sys.stderr); return 2
    sd = Path(__file__).resolve().parent
    cs = sd / "nuc_closure_analysis.py"

    with fp.open("r",encoding="utf-8") as f: data=json.load(f)
    sn = data["fixtures"]["B_nuc_only_nominal"]["inputs"]["species_names"]
    pi = sn.index("C2H2")
    orig = data["fixtures"]["B_nuc_only_nominal"]["outputs"]["omega_gas_kg_m3_s"][pi]
    for fn in ["B_nuc_only_nominal","C_nuc_only_perturbed"]:
        data["fixtures"][fn]["outputs"]["omega_gas_kg_m3_s"][pi] = 0.0
        data["fixtures"][fn]["outputs"]["pah_mass_rate_kg_m3_s"] = 0.0
        data["fixtures"][fn]["outputs"]["mass_residual_kg_m3_s"] = data["fixtures"][fn]["outputs"]["soot_mass_rate_kg_m3_s"]
    print(f"Original omega_gas[C2H2] for B: {orig}")
    print("Mutated to 0.0 -- expecting closure FAIL")

    with tempfile.NamedTemporaryFile(mode="w",suffix=".json",dir=fp.parent,delete=False,encoding="utf-8") as tmp:
        json.dump(data,tmp,indent=2); tp=Path(tmp.name)
    try:
        r = subprocess.run([sys.executable,str(cs),"--input",str(tp)],capture_output=True,text=True)
        print(f"\nClosure exit code: {r.returncode}")
        if r.returncode!=0: print("PASS: NucClosure rejected mutated data."); return 0
        else: print("FAIL: NucClosure did NOT reject mutated data!"); return 1
    finally:
        tp.unlink(missing_ok=True)

if __name__=="__main__": sys.exit(main())
