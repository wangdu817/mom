#!/usr/bin/env python3
"""G1 OH Stoichiometry Closure Analysis (Task 1.3).

Read-only analysis of the OH-only oxidation fixture produced by Task 1.2
(mom/build/oh_fixture_results.json). Computes mass, element (C/H/O) and
enthalpy closure for the HMOM soot-model OH oxidation channel
(HMOM.tpp:1103-1129) and determines definitively whether the code applies
per-C or per-C2 stoichiometry.

This script performs NO modification of any source file. Stdlib only.

Exit codes:
  0 = CLOSED / PASS (all closure terms within tolerance)
  1 = NOT CLOSED (one or more closure terms exceed tolerance)
  2 = Error (file not found, JSON parse error, missing fields, etc.)

Author: swarm/hmom
"""

import argparse
import json
import math
import sys
from pathlib import Path

# --------------------------------------------------------------------------
# Physical constants (from HMOM code facts, verified -- do not re-derive)
# --------------------------------------------------------------------------
N_A = 6.02214076e23          # Avogadro number [1/mol]  (MomentMethodBase.hpp:488)
NAV_KMOL = 6.02214076e26     # Avogadro per kmol used by VC2_ (MomentMethodBase.hpp:489)
WC_KG_PER_KMOL = 12.011      # carbon standard atomic weight WC_ (MomentMethodBase.hpp:491)
RHO_SOOT = 1800.0            # soot particle density [kg/m3]
MW_C = 0.012011              # soot carbon molar mass [kg/mol] (= WC_/1000, exact)

# Species molecular weights [kg/mol] (kg/kmol / 1000)
MW = {
    "H": 1.008e-3,
    "OH": 17.008e-3,
    "O2": 31.999e-3,
    "H2": 2.016e-3,
    "H2O": 18.015e-3,
    "C2H2": 26.038e-3,
    "CO": 28.010e-3,
    "N2": 28.014e-3,
}

# Formation enthalpies h_f^0 at 298 K [J/mol]
H_F = {
    "H": 218000.0,
    "OH": 38990.0,
    "O2": 0.0,
    "H2": 0.0,
    "H2O": -241830.0,
    "C2H2": 226730.0,
    "CO": -110530.0,
    "N2": 0.0,
}
H_F_SOOT = 0.0               # graphite reference [J/mol]

SPECIES = ["H", "OH", "O2", "H2", "H2O", "C2H2", "CO", "N2"]


def sci(x: float) -> str:
    """Format a float in scientific notation with 10 significant figures."""
    return f"{x:.9e}"


def main() -> int:
    parser = argparse.ArgumentParser(
        description="G1 OH Stoichiometry Closure Analysis"
    )
    parser.add_argument(
        "--input",
        type=str,
        default=None,
        help="Path to oh_fixture_results.json (default: mom/build/oh_fixture_results.json)",
    )
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    if args.input:
        results_path = Path(args.input).resolve()
    else:
        results_path = (script_dir / ".." / "build" / "oh_fixture_results.json").resolve()

    if not results_path.exists():
        print(f"ERROR: fixture JSON not found: {results_path}", file=sys.stderr)
        return 2

    try:
        with results_path.open("r", encoding="utf-8") as fh:
            data = json.load(fh)
    except (json.JSONDecodeError, OSError) as exc:
        print(f"ERROR: failed to read/parse {results_path}: {exc}", file=sys.stderr)
        return 2

    # Validate required JSON structure.
    required_keys = ["inputs", "outputs"]
    for key in required_keys:
        if key not in data:
            print(f"ERROR: JSON missing required key '{key}'", file=sys.stderr)
            return 2
    output_keys = ["moment_sources_mol_m3_s", "omega_gas_kg_m3_s", "R_oxid_C2_mol_m3_s"]
    for key in output_keys:
        if key not in data["outputs"]:
            print(f"ERROR: JSON missing required outputs.{key}", file=sys.stderr)
            return 2
    input_keys = ["V0_m3", "soot_density_kg_m3"]
    for key in input_keys:
        if key not in data["inputs"]:
            print(f"ERROR: JSON missing required inputs.{key}", file=sys.stderr)
            return 2

    inputs = data["inputs"]
    outputs = data["outputs"]

    V0 = float(inputs["V0_m3"])
    rho_soot = float(inputs["soot_density_kg_m3"])
    S_M10 = float(outputs["moment_sources_mol_m3_s"]["S_M10"])
    omega = {k: float(v) for k, v in outputs["omega_gas_kg_m3_s"].items()}
    R_oxid_C2_json = float(outputs["R_oxid_C2_mol_m3_s"])

    # ------------------------------------------------------------------
    # VC2_: volume of one C2 pair in soot (HMOM.hpp:1066)
    # VC2 = (WC_ / rho_particle_ / Nav_kmol_) * 2.
    # ------------------------------------------------------------------
    VC2 = (WC_KG_PER_KMOL / rho_soot / NAV_KMOL) * 2.0

    # ==================================================================
    # SOOT SIDE (independent of gas stoichiometry)
    # ==================================================================
    # Spec section 8 (spec.md:148): omega_soot = rho_soot * V0 * N_A * S_M10
    omega_soot = rho_soot * V0 * N_A * S_M10                  # [kg/m3/s], < 0
    soot_C_loss = abs(omega_soot) / MW_C                      # [mol C/m3/s]
    # True per-C2 oxidation rate from the SOOT side: -S_M10 * V0 / VC2
    # (HMOM.tpp:1080). This is INVARIANT to the gas-phase OH stoichiometry,
    # unlike R_oxid_C2_json which the fixture derives from omega[OH] (and which
    # therefore scales 1:1 with the applied OH consumption -- buggy 1x, fixed 2x).
    R_oxid_C2_true = -S_M10 * V0 / VC2                        # [mol C2/m3/s]
    C_atoms_per_C2 = soot_C_loss / R_oxid_C2_true            # should be ~2.0

    # ==================================================================
    # GAS SIDE (from omega_gas)
    # ==================================================================
    OH_consumption = -omega["OH"] / MW["OH"]                  # [mol OH/m3/s]
    O2_consumption = -omega["O2"] / MW["O2"]                  # [mol O2/m3/s] ~0
    CO_production = omega["CO"] / MW["CO"]                    # [mol CO/m3/s]
    H2_production = omega["H2"] / MW["H2"]                    # [mol H2/m3/s]

    # ==================================================================
    # RATIOS (the G1 verdict numbers)
    # ==================================================================
    OH_per_C = OH_consumption / soot_C_loss                  # correct 1.0 / buggy 0.5
    CO_per_C = CO_production / soot_C_loss                    # correct 1.0 / buggy 0.5
    H2_per_C = H2_production / soot_C_loss                    # correct 0.5 / buggy 0.25
    OH_per_C2 = OH_consumption / R_oxid_C2_true             # correct 2.0 / buggy 1.0
    CO_per_C2 = CO_production / R_oxid_C2_true               # correct 2.0 / buggy 1.0
    H2_per_C2 = H2_production / R_oxid_C2_true               # correct 1.0 / buggy 0.5

    # ==================================================================
    # MASS BALANCE
    # ==================================================================
    gas_mass_gain = sum(omega[s] for s in SPECIES)           # [kg/m3/s]
    soot_mass_loss = -omega_soot                             # [kg/m3/s], > 0
    # Raw mass balance from the fixture omega sum. NOTE: this carries the
    # fixture MW-table's element inconsistency (CO=28.010 -> O=15.999 while
    # OH=17.008 -> O=16.000), so even a perfectly element-balanced reaction
    # leaves a ~8e-5 residual here. It is a data-table artifact, NOT a
    # stoichiometry error -- see the element-consistent check below.
    mass_imbalance_rel = (gas_mass_gain - soot_mass_loss) / soot_mass_loss

    # ==================================================================
    # ELEMENT BALANCES [mol/m3/s]
    # ==================================================================
    # Carbon: soot loses C, gas gains C as CO (1 C per CO)
    gas_C_gain = CO_production * 1.0
    C_imbalance_rel = (gas_C_gain - soot_C_loss) / soot_C_loss
    # Hydrogen: in via OH (1 H per OH), out via H2 (2 H per H2)
    H_in = OH_consumption * 1.0
    H_out = H2_production * 2.0
    H_imbalance_rel = (H_out - H_in) / H_in
    # Oxygen: in via OH (1 O per OH), out via CO (1 O per CO)
    O_in = OH_consumption * 1.0
    O_out = CO_production * 1.0
    O_imbalance_rel = (O_out - O_in) / O_in

    # ------------------------------------------------------------------
    # ELEMENT-CONSISTENT MASS BALANCE (immune to MW-table rounding)
    # ------------------------------------------------------------------
    # Reconstruct gas/soot mass flux from element moles x standard atomic
    # weights so the balance does not inherit the fixture MW-table's internal
    # O-weight inconsistency. This is the physically authoritative mass check.
    AW_C = 12.011e-3   # [kg/mol] carbon (matches WC_)
    AW_H = 1.008e-3    # [kg/mol] hydrogen
    AW_O = 15.999e-3   # [kg/mol] oxygen
    # Gas element inventory change [mol/m3/s]: C via CO, H via (OH in, H2 out),
    # O via (OH in, CO out).
    gas_mass_gain_elem = (
        CO_production * AW_C
        + (H2_production * 2.0 - OH_consumption * 1.0) * AW_H
        + (CO_production * 1.0 - OH_consumption * 1.0) * AW_O
    )
    soot_mass_loss_elem = soot_C_loss * AW_C
    mass_imbalance_elem_rel = (
        (gas_mass_gain_elem - soot_mass_loss_elem) / soot_mass_loss_elem
    )

    # ==================================================================
    # ENTHALPY CLOSURE [W/m3]
    # ==================================================================
    # Q = heat release = -sum(h_f_i / MW_i * omega_i) - h_f_soot/MW_C * omega_soot
    # Q_HMOM is computed directly from the fixture omega_gas (whatever
    # stoichiometry the code currently applies -- buggy or fixed).
    gas_enthalpy_flux = sum(
        (H_F[s] / MW[s]) * omega[s] for s in SPECIES
    )
    soot_enthalpy_flux = (H_F_SOOT / MW_C) * omega_soot
    Q_HMOM = -gas_enthalpy_flux - soot_enthalpy_flux

    # Cross-check: per-C2 reaction C2 + 2OH -> 2CO + H2 (the physically correct
    # per-C2 stoichiometry).
    dH_C2 = 2.0 * H_F["CO"] + H_F["H2"] - 2.0 * H_F["OH"]     # [J/mol_C2]
    Q_correct_crosscheck = -dH_C2 * R_oxid_C2_true           # [W/m3]

    # Relative enthalpy imbalance of the fixture against the correct per-C2 Q.
    # ~0 when fixed, ~-0.5 when buggy (fixture releases half the heat).
    Q_imbalance_rel = (Q_HMOM - Q_correct_crosscheck) / Q_correct_crosscheck

    # ==================================================================
    # OUTPUT
    # ==================================================================
    print("=" * 74)
    print("G1 OH STOICHIOMETRY CLOSURE ANALYSIS (Task 1.3)")
    print("Source fixture: " + str(results_path))
    print("HMOM OH channel under test: HMOM.tpp:1103-1129")
    print("=" * 74)

    print("\n--- CONSTANTS ---")
    print(f"  VC2 (C2-pair volume)          = {sci(VC2)} m^3")
    print(f"  V0 (nucleus volume)           = {sci(V0)} m^3")
    print(f"  rho_soot                      = {sci(rho_soot)} kg/m^3")
    print(f"  MW_C                          = {sci(MW_C)} kg/mol")
    print(f"  N_A                           = {sci(N_A)} 1/mol")

    print("\n--- SOOT SIDE (independent of gas stoichiometry) ---")
    print(f"  omega_soot  (spec S8)         = {sci(omega_soot)} kg/m^3/s")
    print(f"  soot_C_loss                   = {sci(soot_C_loss)} mol C/m^3/s")
    print(f"  R_oxid_C2 (soot-side, true)   = {sci(R_oxid_C2_true)} mol C2/m^3/s")
    print(f"  R_oxid_C2 (JSON, from omega_OH) = {sci(R_oxid_C2_json)} mol C2/m^3/s")
    print(f"  C_atoms_per_C2 (expect ~2.0)  = {sci(C_atoms_per_C2)}")

    print("\n--- GAS SIDE (from omega_gas) ---")
    print(f"  OH_consumption                = {sci(OH_consumption)} mol OH/m^3/s")
    print(f"  O2_consumption (expect ~0)    = {sci(O2_consumption)} mol O2/m^3/s")
    print(f"  CO_production                 = {sci(CO_production)} mol CO/m^3/s")
    print(f"  H2_production                 = {sci(H2_production)} mol H2/m^3/s")

    print("\n--- RATIOS (G1 VERDICT NUMBERS) ---")
    print(f"  OH_per_C   (correct 1.0 / buggy 0.5)  = {sci(OH_per_C)}")
    print(f"  CO_per_C   (correct 1.0 / buggy 0.5)  = {sci(CO_per_C)}")
    print(f"  H2_per_C   (correct 0.5 / buggy 0.25) = {sci(H2_per_C)}")
    print(f"  OH_per_C2  (correct 2.0 / buggy 1.0)  = {sci(OH_per_C2)}")
    print(f"  CO_per_C2  (correct 2.0 / buggy 1.0)  = {sci(CO_per_C2)}")
    print(f"  H2_per_C2  (correct 1.0 / buggy 0.5)  = {sci(H2_per_C2)}")

    print("\n--- MASS BALANCE ---")
    print(f"  gas_mass_gain                 = {sci(gas_mass_gain)} kg/m^3/s")
    print(f"  soot_mass_loss                = {sci(soot_mass_loss)} kg/m^3/s")
    print(f"  mass_imbalance_rel (raw, MW-table artifact ~-8e-5) = {sci(mass_imbalance_rel)}")
    print(f"  mass_imbalance_elem_rel (element-consistent, authoritative) = {sci(mass_imbalance_elem_rel)}")

    print("\n--- ELEMENT BALANCES ---")
    print(f"  Carbon:   soot_C_loss={sci(soot_C_loss)}  gas_C_gain={sci(gas_C_gain)}")
    print(f"            C_imbalance_rel  (buggy ~-0.5) = {sci(C_imbalance_rel)}")
    print(f"  Hydrogen: H_in={sci(H_in)}  H_out={sci(H_out)}")
    print(f"            H_imbalance_rel  (expect ~0)   = {sci(H_imbalance_rel)}")
    print(f"  Oxygen:   O_in={sci(O_in)}  O_out={sci(O_out)}")
    print(f"            O_imbalance_rel  (expect ~0)   = {sci(O_imbalance_rel)}")

    print("\n--- ENTHALPY CLOSURE ---")
    print(f"  Q_HMOM (from fixture)         = {sci(Q_HMOM)} W/m^3")
    print(f"  Q_correct_crosscheck          = {sci(Q_correct_crosscheck)} W/m^3")
    print(f"  dH_C2 (C2+2OH->2CO+H2)        = {sci(dH_C2)} J/mol_C2")
    print(f"  Q_imbalance_rel (fixed ~0 / buggy ~-0.5) = {sci(Q_imbalance_rel)}")

    # ------------------------------------------------------------------
    # SUMMARY TABLE + VERDICT
    # ------------------------------------------------------------------
    print("\n" + "=" * 74)
    print("SUMMARY TABLE")
    print("=" * 74)
    header = f"{'Quantity':<22}{'Value':<20}{'Correct-per-C2':<18}{'Imbalance(rel)':<16}"
    print(header)
    print("-" * 74)
    rows = [
        ("OH_per_C", OH_per_C, 1.0),
        ("CO_per_C", CO_per_C, 1.0),
        ("H2_per_C", H2_per_C, 0.5),
        ("OH_per_C2", OH_per_C2, 2.0),
        ("CO_per_C2", CO_per_C2, 2.0),
        ("H2_per_C2", H2_per_C2, 1.0),
        ("mass_imbalance(raw)", mass_imbalance_rel, 0.0),
        ("mass_imbalance(elem)", mass_imbalance_elem_rel, 0.0),
        ("C_imbalance", C_imbalance_rel, 0.0),
        ("H_imbalance", H_imbalance_rel, 0.0),
        ("O_imbalance", O_imbalance_rel, 0.0),
    ]
    for name, val, correct in rows:
        print(f"{name:<22}{sci(val):<20}{sci(correct):<18}")

    print("\n" + "=" * 74)
    # Classify the applied stoichiometry from the per-C ratios.
    is_per_c = abs(OH_per_C - 0.5) < 0.05 and abs(CO_per_C - 0.5) < 0.05
    is_per_c2 = abs(OH_per_C - 1.0) < 0.05 and abs(CO_per_C - 1.0) < 0.05
    mass_broken = abs(mass_imbalance_elem_rel + 0.5) < 0.05

    # Closure: all balances must close to < 1e-6 relative and Q must match the
    # per-C2 crosscheck to < 1e-6 relative.
    CLOSURE_TOL = 1e-6
    closure_terms = {
        "mass_imbalance_elem_rel": mass_imbalance_elem_rel,
        "C_imbalance_rel": C_imbalance_rel,
        "H_imbalance_rel": H_imbalance_rel,
        "O_imbalance_rel": O_imbalance_rel,
        "Q_imbalance_rel": Q_imbalance_rel,
    }
    is_closed = all(abs(v) < CLOSURE_TOL for v in closure_terms.values())

    if is_per_c:
        verdict = "PER-C stoichiometry applied at the PER-C2 rate level (R_oxid_C2)"
    elif is_per_c2:
        verdict = "PER-C2 stoichiometry (C2 + 2OH -> 2CO + H2) applied correctly"
    else:
        verdict = (
            f"UNEXPECTED stoichiometry (OH_per_C={OH_per_C:.4f}, "
            f"CO_per_C={CO_per_C:.4f})"
        )
    print("VERDICT: " + verdict)
    print(
        f"  - soot_C_loss = {sci(C_atoms_per_C2)} x R_oxid_C2 "
        "(confirms R_oxid_C2 is genuinely per-C2)"
    )
    print(
        f"  - OH_per_C={OH_per_C:.4f}, CO_per_C={CO_per_C:.4f}, "
        f"H2_per_C={H2_per_C:.4f}"
    )
    print(
        f"  - OH_per_C2={OH_per_C2:.4f}, CO_per_C2={CO_per_C2:.4f}, "
        f"H2_per_C2={H2_per_C2:.4f}"
    )

    print("\n  Closure check (tol < 1e-6 relative):")
    for name, val in closure_terms.items():
        status = "PASS" if abs(val) < CLOSURE_TOL else "FAIL"
        print(f"    {name:<22} = {sci(val):<20} [{status}]")

    print("\n" + "=" * 74)
    if is_closed:
        print(
            "RESULT: CLOSED / PASS -- all mass, element (C/H/O) and enthalpy "
            "balances close to < 1e-6 relative."
        )
        print(
            "  CONCLUSION: HMOM OH channel applies correct per-C2 stoichiometry "
            "(2 OH, 2 CO, 1 H2 per C2). G1 stoichiometry fix (Task 1.4) verified."
        )
        print("=" * 74)
        return 0
    else:
        print(
            "RESULT: NOT CLOSED -- one or more balances exceed 1e-6 relative."
        )
        print(
            f"  Mass {'DOES NOT close (~-50%)' if mass_broken else 'imbalance present'}; "
            "see closure check above."
        )
        print(
            "  CONCLUSION: HMOM.tpp OH channel has a STOICHIOMETRY BUG. Fix (Task 1.4): "
            "multiply OH-channel gas terms by 2 (2 OH, 2 CO, 1 H2 per C2)."
        )
        print("=" * 74)
        return 1


if __name__ == "__main__":
    sys.exit(main())
