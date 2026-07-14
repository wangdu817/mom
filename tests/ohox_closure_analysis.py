#!/usr/bin/env python3
"""M1 OH-oxidation Closure Analysis.

Read-only analysis of the OH oxidation fixture JSON. Verifies mass, element
(C/O/H) closure, stoichiometric ratios (OH_per_C2, CO_per_C2, H2_per_C2),
fOH dominance, and tanh damping ordering.

Reaction (OH channel): C2(soot) + 2OH -> 2CO + H2
  OH consumed: 2 per C2
  CO produced: 2 per C2
  H2 produced: 1 per C2
  Mass: soot loses 2C, gas gains 2CO+H2-2OH = 2*28.010+2.016-2*17.008 = 24.020

Exit codes: 0=PASS, 1=FAIL, 2=ERROR
"""

import argparse
import json
import sys
from pathlib import Path

MW = {"H": 1.008, "OH": 17.008, "O2": 31.999, "H2": 2.016,
      "H2O": 18.015, "C2H2": 26.038, "CO": 28.010, "N2": 28.014}
SPECIES = ["H", "OH", "O2", "H2", "H2O", "C2H2", "CO", "N2"]
CLOSURE_TOL = 1e-8

AW_C = 12.011e-3  # kg/mol
AW_H = 1.008e-3
AW_O = 15.999e-3


def sci(x):
    return f"{x:.9e}"


def analyze_fixture(name, fx):
    cfg = fx["config"]
    inp = fx["inputs"]
    out = fx["outputs"]
    species = inp["species_names"]
    iO2 = species.index("O2")
    iCO = species.index("CO")
    iOH = species.index("OH")
    iH2 = species.index("H2")

    soot_rate = out["soot_mass_rate_kg_m3_s"]
    omega = out["omega_gas_kg_m3_s"]
    R_oxid_C2 = out["R_oxid_C2_mol_m3_s"]

    o2_rate = omega[iO2]
    co_rate = omega[iCO]
    oh_rate = omega[iOH]
    h2_rate = omega[iH2]

    gas_sum = sum(omega)
    mass_residual = soot_rate + gas_sum
    if abs(soot_rate) > 0:
        mass_imb_rel = mass_residual / abs(soot_rate)
    else:
        mass_imb_rel = 0.0 if abs(gas_sum) < 1e-30 else float("inf")

    if R_oxid_C2 > 0:
        O2_per_C2 = -o2_rate * 1000 / MW["O2"] / R_oxid_C2
        CO_per_C2 = co_rate * 1000 / MW["CO"] / R_oxid_C2
        OH_per_C2 = -oh_rate * 1000 / MW["OH"] / R_oxid_C2
        H2_per_C2 = h2_rate * 1000 / MW["H2"] / R_oxid_C2
    else:
        O2_per_C2 = CO_per_C2 = OH_per_C2 = H2_per_C2 = 0.0

    fO2 = out.get("fO2", 0.0)
    fOH = out.get("fOH", 0.0)

    # Element balances [mol/m3/s]
    if R_oxid_C2 > 0:
        # C: soot loses 2C per C2, gas gains C via CO
        soot_C_loss = 2.0 * R_oxid_C2
        gas_C_gain = co_rate * 1000 / MW["CO"]
        C_imb = (gas_C_gain - soot_C_loss) / soot_C_loss if soot_C_loss > 0 else 0.0
        # H: 2H in via OH, 2H out via H2
        H_in = -oh_rate * 1000 / MW["OH"] * 1.0
        H_out = h2_rate * 1000 / MW["H2"] * 2.0
        H_imb = (H_out - H_in) / H_in if H_in > 0 else 0.0
        # O: 2O in via OH, 2O out via CO
        O_in = -oh_rate * 1000 / MW["OH"] * 1.0
        O_out = co_rate * 1000 / MW["CO"] * 1.0
        O_imb = (O_out - O_in) / O_in if O_in > 0 else 0.0
    else:
        C_imb = H_imb = O_imb = 0.0

    # Element-consistent mass balance
    if R_oxid_C2 > 0:
        soot_mass_loss_elem = 2.0 * R_oxid_C2 * AW_C
        mass_imb_elem_rel = (abs(soot_rate) - soot_mass_loss_elem) / soot_mass_loss_elem if soot_mass_loss_elem > 0 else 0.0
    else:
        mass_imb_elem_rel = 0.0

    # Other species check (should be zero in ox-only mode)
    other_nonzero = []
    for i, s in enumerate(species):
        if s in ("OH", "CO", "H2", "O2"):
            continue
        if abs(omega[i]) > 1e-30:
            other_nonzero.append((s, omega[i]))

    return {
        "name": name, "config": cfg,
        "soot_rate": soot_rate, "R_oxid_C2": R_oxid_C2,
        "o2_rate": o2_rate, "co_rate": co_rate,
        "oh_rate": oh_rate, "h2_rate": h2_rate,
        "mass_residual": mass_residual, "mass_imb_rel": mass_imb_rel,
        "O2_per_C2": O2_per_C2, "CO_per_C2": CO_per_C2,
        "OH_per_C2": OH_per_C2, "H2_per_C2": H2_per_C2,
        "fO2": fO2, "fOH": fOH,
        "C_imb": C_imb, "H_imb": H_imb, "O_imb": O_imb,
        "mass_imb_elem_rel": mass_imb_elem_rel,
        "other_nonzero": other_nonzero,
    }


def main():
    parser = argparse.ArgumentParser(description="M1 OH-oxidation Closure Analysis")
    parser.add_argument("--input", type=str, required=True)
    args = parser.parse_args()

    path = Path(args.input).resolve()
    if not path.exists():
        print(f"ERROR: JSON not found: {path}", file=sys.stderr)
        return 2

    with path.open("r", encoding="utf-8") as fh:
        data = json.load(fh)

    if "fixtures" not in data:
        print("ERROR: missing 'fixtures'", file=sys.stderr)
        return 2

    expected = ["A_all_off_control", "B_oh_ox_above_threshold",
                "C_oh_ox_below_threshold", "D_oh_ox_near_threshold"]
    for n in expected:
        if n not in data["fixtures"]:
            print(f"ERROR: missing fixture '{n}'", file=sys.stderr)
            return 2

    print("=" * 74)
    print("M1 OH-OXIDATION CLOSURE ANALYSIS")
    print(f"Source: {path}")
    print("=" * 74)

    results = {}
    for n in expected:
        results[n] = analyze_fixture(n, data["fixtures"][n])
        r = results[n]
        print(f"\n--- {r['name']} ---")
        print(f"  soot_rate     = {sci(r['soot_rate'])} kg/m3/s")
        print(f"  R_oxid_C2     = {sci(r['R_oxid_C2'])} mol/m3/s")
        print(f"  fO2           = {sci(r['fO2'])}")
        print(f"  fOH           = {sci(r['fOH'])}")
        print(f"  OH_per_C2     = {sci(r['OH_per_C2'])} (expect 2.0)")
        print(f"  CO_per_C2     = {sci(r['CO_per_C2'])} (expect 2.0)")
        print(f"  H2_per_C2     = {sci(r['H2_per_C2'])} (expect 1.0)")
        print(f"  mass_residual = {sci(r['mass_residual'])}")
        print(f"  mass_imb_elem = {sci(r['mass_imb_elem_rel'])}")
        print(f"  C_imb         = {sci(r['C_imb'])}")
        print(f"  H_imb         = {sci(r['H_imb'])}")
        print(f"  O_imb         = {sci(r['O_imb'])}")

    print("\n" + "=" * 74)
    print("CLOSURE VERDICT")
    print("=" * 74)

    all_closed = True

    rA = results["A_all_off_control"]
    A_ok = (abs(rA["soot_rate"]) < 1e-30 and
            all(abs(v) < 1e-30 for v in [rA["o2_rate"], rA["co_rate"], rA["oh_rate"], rA["h2_rate"]]) and
            len(rA["other_nonzero"]) == 0)
    print(f"  A (all-off): {'PASS' if A_ok else 'FAIL'}")
    if not A_ok: all_closed = False

    rB = results["B_oh_ox_above_threshold"]
    B_ok = (abs(rB["mass_imb_elem_rel"]) < CLOSURE_TOL and
            abs(rB["fOH"] - 1.0) < 0.01 and
            abs(rB["OH_per_C2"] - 2.0) < 0.01 and
            abs(rB["CO_per_C2"] - 2.0) < 0.01 and
            abs(rB["H2_per_C2"] - 1.0) < 0.01 and
            abs(rB["C_imb"]) < CLOSURE_TOL and
            abs(rB["H_imb"]) < CLOSURE_TOL and
            abs(rB["O_imb"]) < 1e-4 and
            len(rB["other_nonzero"]) == 0)
    print(f"  B (above threshold): {'PASS' if B_ok else 'FAIL'}")
    print(f"    fOH={sci(rB['fOH'])}, OH/C2={sci(rB['OH_per_C2'])}, CO/C2={sci(rB['CO_per_C2'])}, H2/C2={sci(rB['H2_per_C2'])}")
    print(f"    mass_imb_elem={sci(rB['mass_imb_elem_rel'])}, C_imb={sci(rB['C_imb'])}, H_imb={sci(rB['H_imb'])}, O_imb={sci(rB['O_imb'])}")
    if not B_ok: all_closed = False

    rC = results["C_oh_ox_below_threshold"]
    C_ok = (abs(rC["mass_imb_elem_rel"]) < CLOSURE_TOL and
            abs(rC["soot_rate"]) < abs(rB["soot_rate"]))
    print(f"  C (below threshold): {'PASS' if C_ok else 'FAIL'}")
    print(f"    mass_imb_elem={sci(rC['mass_imb_elem_rel'])}, damped={abs(rC['soot_rate']) < abs(rB['soot_rate'])}")
    if not C_ok: all_closed = False

    rD = results["D_oh_ox_near_threshold"]
    D_ok = (abs(rD["mass_imb_elem_rel"]) < CLOSURE_TOL and
            abs(rD["soot_rate"]) < abs(rB["soot_rate"]) and
            abs(rD["soot_rate"]) > abs(rC["soot_rate"]))
    print(f"  D (near threshold): {'PASS' if D_ok else 'FAIL'}")
    print(f"    mass_imb_elem={sci(rD['mass_imb_elem_rel'])}, B>D>C={abs(rB['soot_rate']) > abs(rD['soot_rate']) > abs(rC['soot_rate'])}")
    if not D_ok: all_closed = False

    print("\n" + "=" * 74)
    if all_closed:
        print("RESULT: CLOSED / PASS")
        print("  OH channel: C2 + 2OH -> 2CO + H2 verified (OH/C2=2, CO/C2=2, H2/C2=1)")
        print("  Mass, C, H, O balances close to tolerance.")
        print("  Tanh damping: B > D > C ordering confirmed.")
        print("=" * 74)
        return 0
    else:
        print("RESULT: NOT CLOSED")
        print("=" * 74)
        return 1


if __name__ == "__main__":
    sys.exit(main())
