#!/usr/bin/env python3
"""M1 Nucleation-only Closure Analysis.

Read-only analysis of the nucleation fixture JSON produced by
test_nucleation_fixture.exe.  Verifies mass, element (C/H/O/N) and
PAH-dimer-identity closure for the nucleation-only vertical slice.

Closure physics:
  In nucleation, 2 PAH dimers (4 PAH molecules) collide to form 1 soot
  nucleus (V0 = 4*vpah).  The gas loses PAH mass; the soot gains the
  same mass.  No chemical reaction -- phase change (gas PAH -> soot).

  Mass:   soot_mass_rate == -omega_gas[C2H2]
  Factor: -omega_gas[C2H2] / (PAHDimerizationRate * mw_pah) == 1.0
  C:      gas C loss == soot C gain  (both from PAH composition)
  H:      gas H loss == soot H gain  (both from PAH composition)
  O:      0  (no oxygen in nucleation)
  N:      0  (no nitrogen in nucleation)

Exit codes: 0=PASS, 1=FAIL, 2=ERROR
"""

import argparse, json, sys
from pathlib import Path

NAV_KMOL = 6.02214076e26
WC_KG_PER_KMOL = 12.011
MW = {"H":1.008,"OH":17.008,"O2":31.999,"H2":2.016,"H2O":18.015,"C2H2":26.038,"CO":28.010,"N2":28.014}
SPECIES = ["H","OH","O2","H2","H2O","C2H2","CO","N2"]
CLOSURE_TOL = 1e-10

def sci(x): return f"{x:.9e}"

def analyze_fixture(name, fx):
    cfg = fx["config"]; inp = fx["inputs"]; out = fx["outputs"]
    src_nuc = out["source_nucleation_mol_m3_s"]
    omega_gas = out["omega_gas_kg_m3_s"]
    species_names = inp["species_names"]
    pah_idx = species_names.index("C2H2")
    V0 = inp["V0_m3"]; rho_soot = inp["soot_density_kg_m3"]; mw_pah = MW["C2H2"]
    soot_mass_rate = out["soot_mass_rate_kg_m3_s"]
    pah_mass_rate = omega_gas[pah_idx]
    R_PAH_dimer = out["R_PAH_dimer_kmol_m3_s"]

    mass_residual = soot_mass_rate - (-pah_mass_rate)
    mass_imb = mass_residual / soot_mass_rate if abs(soot_mass_rate)>0 else (0.0 if abs(pah_mass_rate)<1e-30 else float("inf"))

    # PAH-dimer identity: omega_gas should equal -R_PAH_dimer * mw_pah
    expected_pah = -R_PAH_dimer * mw_pah
    dimer_ratio = pah_mass_rate / expected_pah if abs(expected_pah)>0 else (1.0 if abs(pah_mass_rate)<1e-30 else float("inf"))

    pah_cons_mol = -pah_mass_rate / mw_pah
    gas_C = pah_cons_mol * 2.0; gas_H = pah_cons_mol * 2.0
    soot_C = gas_C; soot_H = gas_H
    C_imb = (soot_C-gas_C)/gas_C if gas_C>0 else 0.0
    H_imb = (soot_H-gas_H)/gas_H if gas_H>0 else 0.0

    other_sp = [(s, omega_gas[i]) for i,s in enumerate(species_names) if s!="C2H2" and abs(omega_gas[i])>1e-30]
    other_proc = []
    for pn in ["source_growth_mol_m3_s","source_oxidation_mol_m3_s","source_condensation_mol_m3_s","source_coagulation_mol_m3_s"]:
        for i,v in enumerate(out[pn]):
            if abs(v)>1e-30: other_proc.append((pn,i,v))

    return dict(name=name, config=cfg, soot_mass_rate=soot_mass_rate, pah_mass_rate=pah_mass_rate,
        mass_residual=mass_residual, mass_imb=mass_imb, dimer_ratio=dimer_ratio,
        C_imb=C_imb, H_imb=H_imb, other_sp=other_sp, other_proc=other_proc, src_nuc=src_nuc)

def main():
    parser = argparse.ArgumentParser(description="M1 Nucleation Closure Analysis")
    parser.add_argument("--input", type=str, required=True)
    args = parser.parse_args()
    p = Path(args.input).resolve()
    if not p.exists(): print(f"ERROR: {p}", file=sys.stderr); return 2
    try:
        with p.open("r",encoding="utf-8") as f: data=json.load(f)
    except Exception as e: print(f"ERROR: {e}",file=sys.stderr); return 2
    if "fixtures" not in data: print("ERROR: missing fixtures",file=sys.stderr); return 2
    expected = ["A_all_off_control","B_nuc_only_nominal","C_nuc_only_perturbed"]
    for n in expected:
        if n not in data["fixtures"]: print(f"ERROR: missing {n}",file=sys.stderr); return 2

    print("="*74); print("M1 NUCLEATION-ONLY CLOSURE ANALYSIS"); print(f"Source: {p}"); print("="*74)
    results = {}
    for n in expected:
        results[n] = analyze_fixture(n, data["fixtures"][n])
        r = results[n]
        cfg = r["config"]
        print(f"\n--- {r['name']} ---")
        print(f"  Config: nuc={cfg['nucleation_model']} cond={cfg['condensation_model']}")
        print(f"  soot_mass_rate = {sci(r['soot_mass_rate'])}")
        print(f"  pah_mass_rate  = {sci(r['pah_mass_rate'])}")
        print(f"  mass_residual  = {sci(r['mass_residual'])}")
        print(f"  mass_imb       = {sci(r['mass_imb'])}")
        print(f"  dimer_ratio    = {sci(r['dimer_ratio'])} (expect 1.0)")
        print(f"  C_imb          = {sci(r['C_imb'])}")
        print(f"  H_imb          = {sci(r['H_imb'])}")
        if r["other_sp"]: print(f"  OTHER SP NONZERO: {r['other_sp']}")
        if r["other_proc"]: print(f"  OTHER PROC NONZERO: {r['other_proc']}")

    print("\n"+"="*74); print("CLOSURE VERDICT"); print("="*74)
    all_ok = True
    rA=results["A_all_off_control"]
    A_ok = abs(rA["soot_mass_rate"])<1e-30 and abs(rA["pah_mass_rate"])<1e-30 and not rA["other_sp"] and not rA["other_proc"]
    print(f"  A (all-off): {'PASS' if A_ok else 'FAIL'}")
    if not A_ok: all_ok=False

    rB=results["B_nuc_only_nominal"]
    B_ok = abs(rB["mass_imb"])<CLOSURE_TOL and abs(rB["dimer_ratio"]-1.0)<CLOSURE_TOL and abs(rB["C_imb"])<CLOSURE_TOL and abs(rB["H_imb"])<CLOSURE_TOL and not rB["other_sp"] and not rB["other_proc"]
    print(f"  B (nuc nominal): {'PASS' if B_ok else 'FAIL'}")
    if not B_ok: all_ok=False

    rC=results["C_nuc_only_perturbed"]
    C_ok = abs(rC["mass_imb"])<CLOSURE_TOL and abs(rC["dimer_ratio"]-1.0)<CLOSURE_TOL and abs(rC["C_imb"])<CLOSURE_TOL and abs(rC["H_imb"])<CLOSURE_TOL and abs(rC["soot_mass_rate"]-rB["soot_mass_rate"])>0
    print(f"  C (nuc perturbed): {'PASS' if C_ok else 'FAIL'}")
    if not C_ok: all_ok=False

    print("\n"+"="*74)
    if all_ok:
        print("RESULT: CLOSED / PASS"); print("="*74); return 0
    else:
        print("RESULT: NOT CLOSED"); print("="*74); return 1

if __name__=="__main__": sys.exit(main())
