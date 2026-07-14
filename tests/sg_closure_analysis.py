#!/usr/bin/env python3
"""M1 Surface-growth-only Closure Analysis.

Read-only analysis of the surface growth fixture JSON. Verifies mass,
element (C/H/O/N) closure for the HACA surface growth pathway.

Reaction: C2H2 + soot-surface -> 2C(soot) + H2
  Gas: C2H2 consumed, H2 produced
  Soot: mass increases (carbon added)
  Mass balance: soot_gain = -(omega_gas[C2H2] + omega_gas[H2])

Also verifies tanh damping behavior: above > near > below threshold.

Exit codes: 0=PASS, 1=FAIL, 2=ERROR
"""

import argparse, json, sys
from pathlib import Path

MW = {"H":1.008,"OH":17.008,"O2":31.999,"H2":2.016,"H2O":18.015,"C2H2":26.038,"CO":28.010,"N2":28.014}
CLOSURE_TOL = 1e-10

def sci(x): return f"{x:.9e}"

def analyze_fixture(name, fx):
    cfg = fx["config"]; inp = fx["inputs"]; out = fx["outputs"]
    species = inp["species_names"]
    omega = out["omega_gas_kg_m3_s"]
    iC2H2 = species.index("C2H2"); iH2 = species.index("H2")

    soot_rate = out["soot_mass_rate_kg_m3_s"]
    c2h2_rate = omega[iC2H2]; h2_rate = omega[iH2]
    mass_residual = soot_rate - (-c2h2_rate - h2_rate)
    mass_imb = mass_residual / soot_rate if abs(soot_rate) > 0 else (0.0 if abs(c2h2_rate+h2_rate)<1e-30 else float("inf"))

    # Element balances [mol/m3/s] -- MW in kg/kmol, omega in kg/m3/s
    # R_C2H2 = -c2h2_rate [kg/m3/s] / (MW_C2H2 [kg/kmol] / 1000) = mol/m3/s
    R_C2H2 = -c2h2_rate * 1000.0 / MW["C2H2"]
    gas_C_loss = R_C2H2 * 2.0
    gas_H_loss_from_c2h2 = R_C2H2 * 2.0
    gas_H_gain_from_h2 = h2_rate * 1000.0 / MW["H2"] * 2.0
    net_gas_H = gas_H_gain_from_h2 - gas_H_loss_from_c2h2  # should be 0

    # Soot gains C mass = 2*WC per C2H2 [kg/m3/s]
    soot_C_mass = R_C2H2 * 2.0 * 12.011 / 1000.0
    soot_C_gain = soot_C_mass / (12.011e-3)  # mol C/m3/s
    C_imb = (soot_C_gain - gas_C_loss) / gas_C_loss if gas_C_loss > 0 else 0.0
    H_imb = net_gas_H / gas_H_loss_from_c2h2 if gas_H_loss_from_c2h2 > 0 else 0.0

    # Other species/processes must be zero
    other_sp = [(s, omega[i]) for i,s in enumerate(species) if s not in ("C2H2","H2") and abs(omega[i])>1e-30]
    other_pr = []
    for pn in ["source_nucleation_mol_m3_s","source_condensation_mol_m3_s","source_oxidation_mol_m3_s","source_coagulation_mol_m3_s"]:
        for i,v in enumerate(out[pn]):
            if abs(v)>1e-30: other_pr.append((pn,i,v))

    return {"name":name,"config":cfg,"soot_rate":soot_rate,"c2h2_rate":c2h2_rate,"h2_rate":h2_rate,
            "mass_residual":mass_residual,"mass_imb":mass_imb,"R_C2H2":R_C2H2,
            "C_imb":C_imb,"H_imb":H_imb,"other_sp":other_sp,"other_pr":other_pr,
            "src_growth":out["source_growth_mol_m3_s"]}

def main():
    parser = argparse.ArgumentParser(description="M1 Surface-growth Closure Analysis")
    parser.add_argument("--input", type=str, required=True)
    args = parser.parse_args()
    path = Path(args.input).resolve()
    if not path.exists(): print(f"ERROR: {path}", file=sys.stderr); return 2
    with path.open() as f: data = json.load(f)
    if "fixtures" not in data: print("ERROR: no fixtures", file=sys.stderr); return 2

    expected = ["A_all_off_control","B_sg_above_threshold","C_sg_below_threshold","D_sg_near_threshold"]
    for n in expected:
        if n not in data["fixtures"]: print(f"ERROR: missing {n}", file=sys.stderr); return 2

    print("="*74); print("M1 SURFACE-GROWTH CLOSURE ANALYSIS"); print(f"Source: {path}"); print("="*74)

    results = {n: analyze_fixture(n, data["fixtures"][n]) for n in expected}

    for r in results.values():
        print(f"\n--- {r['name']} ---")
        c=r["config"]
        print(f"  Config: nuc={c['nucleation_model']} cond={c['condensation_model']} sg={c['surface_growth_model']} ox={c['oxidation_model']} coag={c['coagulation_model']}")
        print(f"  src_growth[1]  = {sci(r['src_growth'][1])}")
        print(f"  soot_rate      = {sci(r['soot_rate'])}")
        print(f"  omega[C2H2]    = {sci(r['c2h2_rate'])}")
        print(f"  omega[H2]      = {sci(r['h2_rate'])}")
        print(f"  mass_residual  = {sci(r['mass_residual'])}")
        print(f"  mass_imb(rel)  = {sci(r['mass_imb'])}")
        print(f"  C_imb(rel)     = {sci(r['C_imb'])}")
        print(f"  H_imb(rel)     = {sci(r['H_imb'])}")
        if r["other_sp"]: print(f"  OTHER SP NONZERO: {r['other_sp']}")
        if r["other_pr"]: print(f"  OTHER PROC NONZERO: {r['other_pr']}")

    print("\n"+"="*74); print("CLOSURE VERDICT"); print("="*74)
    all_ok = True

    rA = results["A_all_off_control"]
    A_ok = abs(rA["soot_rate"])<1e-30 and len(rA["other_sp"])==0 and len(rA["other_pr"])==0
    print(f"  A (all-off): {'PASS' if A_ok else 'FAIL'}")
    if not A_ok: all_ok=False

    rB = results["B_sg_above_threshold"]
    B_ok = (abs(rB["mass_imb"])<CLOSURE_TOL and abs(rB["C_imb"])<CLOSURE_TOL and abs(rB["H_imb"])<CLOSURE_TOL
            and len(rB["other_sp"])==0 and len(rB["other_pr"])==0)
    print(f"  B (above threshold): {'PASS' if B_ok else 'FAIL'}  mass_imb={sci(rB['mass_imb'])}")
    if not B_ok: all_ok=False

    rC = results["C_sg_below_threshold"]
    C_ok = (abs(rC["mass_imb"])<CLOSURE_TOL and rC["src_growth"][1]>=0 and rC["src_growth"][1]<rB["src_growth"][1])
    print(f"  C (below threshold): {'PASS' if C_ok else 'FAIL'}  mass_imb={sci(rC['mass_imb'])}")
    if not C_ok: all_ok=False

    rD = results["D_sg_near_threshold"]
    D_ok = (abs(rD["mass_imb"])<CLOSURE_TOL and rD["src_growth"][1]>rC["src_growth"][1] and rD["src_growth"][1]<rB["src_growth"][1])
    print(f"  D (near threshold): {'PASS' if D_ok else 'FAIL'}  mass_imb={sci(rD['mass_imb'])}")
    if not D_ok: all_ok=False

    print(f"\n  Tanh damping: B({sci(rB['src_growth'][1])}) > D({sci(rD['src_growth'][1])}) > C({sci(rC['src_growth'][1])})")

    print("\n"+"="*74)
    if all_ok:
        print("RESULT: CLOSED / PASS"); print("="*74); return 0
    else:
        print("RESULT: NOT CLOSED"); print("="*74); return 1

if __name__=="__main__": sys.exit(main())
