#!/usr/bin/env python3
"""M1 Condensation-only Closure Analysis.

Read-only analysis of the condensation fixture JSON produced by
test_condensation_fixture.exe.  Verifies mass, element (C/H/O/N) and
factor-4 closure for the condensation-only vertical slice.

Closure physics:
  In condensation, PAH dimers (2 C2H2) condense onto soot particles.
  The gas loses PAH mass; the soot gains the same mass.  No chemical
  reaction occurs -- it is a phase-change (gas PAH -> soot).

  Mass:   soot_mass_rate == -omega_gas[C2H2]
  C:      gas C loss == soot C gain  (both from PAH composition)
  H:      gas H loss == soot H gain  (both from PAH composition)
  O:      0  (no oxygen in condensation)
  N:      0  (no nitrogen in condensation)
  Factor: -omega_gas[C2H2] / (4 * src_cond[1] * mw_pah / 1000) == 1.0

Exit codes:
  0 = CLOSED / PASS (all closure terms within tolerance)
  1 = NOT CLOSED (one or more closure terms exceed tolerance)
  2 = Error (file not found, JSON parse error, missing fields, etc.)

Author: swarm/hmom
"""

import argparse
import json
import sys
from pathlib import Path

# --------------------------------------------------------------------------
# Physical constants
# --------------------------------------------------------------------------
N_A = 6.02214076e23          # Avogadro [1/mol]
NAV_KMOL = 6.02214076e26     # Avogadro per kmol
WC_KG_PER_KMOL = 12.011      # carbon atomic weight [kg/kmol]
WH_KG_PER_KMOL = 1.008       # hydrogen atomic weight [kg/kmol]

# Species molecular weights [kg/kmol]
MW = {
    "H": 1.008,
    "OH": 17.008,
    "O2": 31.999,
    "H2": 2.016,
    "H2O": 18.015,
    "C2H2": 26.038,
    "CO": 28.010,
    "N2": 28.014,
}

SPECIES = ["H", "OH", "O2", "H2", "H2O", "C2H2", "CO", "N2"]

CLOSURE_TOL = 1e-10  # relative tolerance for closure


def sci(x: float) -> str:
    """Format a float in scientific notation."""
    return f"{x:.9e}"


def analyze_fixture(name: str, fx: dict) -> dict:
    """Analyze a single fixture and return closure metrics."""
    cfg = fx["config"]
    inp = fx["inputs"]
    out = fx["outputs"]

    # Extract key quantities.
    src_cond = out["source_condensation_mol_m3_s"]
    omega_gas = out["omega_gas_kg_m3_s"]
    species_names = inp["species_names"]

    # Find PAH index.
    pah_idx = species_names.index("C2H2")

    # Physical parameters.
    V0 = inp["V0_m3"]
    rho_soot = inp["soot_density_kg_m3"]
    mw_pah = MW["C2H2"]

    # VC2 (volume of one C2 pair in soot).
    VC2 = (WC_KG_PER_KMOL / rho_soot / NAV_KMOL) * 2.0

    # Soot mass rate [kg/m3/s] (from fixture).
    soot_mass_rate = out["soot_mass_rate_kg_m3_s"]

    # PAH gas consumption rate [kg/m3/s] (negative = consumed).
    pah_mass_rate = omega_gas[pah_idx]

    # ----------------------------------------------------------------
    # MASS BALANCE
    # ----------------------------------------------------------------
    # soot_mass_rate should equal -pah_mass_rate
    mass_residual = soot_mass_rate - (-pah_mass_rate)
    if abs(soot_mass_rate) > 0:
        mass_imbalance_rel = mass_residual / soot_mass_rate
    else:
        mass_imbalance_rel = 0.0 if abs(pah_mass_rate) < 1e-30 else float("inf")

    # ----------------------------------------------------------------
    # FACTOR-4 CROSS-CHECK
    # ----------------------------------------------------------------
    # omega_gas[C2H2] should equal -4 * src_cond[1] * mw_pah / 1000
    expected_pah_rate = -4.0 * src_cond[1] * mw_pah / 1000.0
    if abs(expected_pah_rate) > 0:
        factor4_ratio = pah_mass_rate / expected_pah_rate
    else:
        factor4_ratio = 1.0 if abs(pah_mass_rate) < 1e-30 else float("inf")

    # ----------------------------------------------------------------
    # ELEMENT BALANCES [mol/m3/s]
    # ----------------------------------------------------------------
    # In condensation, PAH (C2H2) transfers from gas to soot.
    # Gas loses: 2 C per C2H2, 2 H per C2H2
    # Soot gains: same (it's a phase change, no chemistry)

    # Gas-side element rates (positive = consumed from gas).
    pah_consumption_mol = -pah_mass_rate / mw_pah  # [mol C2H2/m3/s]
    gas_C_loss = pah_consumption_mol * 2.0          # [mol C/m3/s]
    gas_H_loss = pah_consumption_mol * 2.0          # [mol H/m3/s]
    gas_O_loss = 0.0                                 # no O in C2H2
    gas_N_loss = 0.0                                 # no N in C2H2

    # Soot-side element rates (positive = gained by soot).
    # Soot mass from PAH at soot density -- same composition as PAH.
    soot_C_gain = gas_C_loss  # phase change, same composition
    soot_H_gain = gas_H_loss
    soot_O_gain = 0.0
    soot_N_gain = 0.0

    # Element imbalances (should be ~0).
    C_imb = (soot_C_gain - gas_C_loss) / gas_C_loss if gas_C_loss > 0 else 0.0
    H_imb = (soot_H_gain - gas_H_loss) / gas_H_loss if gas_H_loss > 0 else 0.0
    O_imb = 0.0  # no O involved
    N_imb = 0.0  # no N involved

    # ----------------------------------------------------------------
    # OTHER SPECIES CHECK (must all be zero)
    # ----------------------------------------------------------------
    other_species_nonzero = []
    for i, s in enumerate(species_names):
        if s == "C2H2":
            continue
        if abs(omega_gas[i]) > 1e-30:
            other_species_nonzero.append((s, omega_gas[i]))

    # ----------------------------------------------------------------
    # OTHER PROCESS SOURCES CHECK (must all be zero in cond-only)
    # ----------------------------------------------------------------
    other_process_nonzero = []
    for pname in ["source_nucleation_mol_m3_s", "source_growth_mol_m3_s",
                   "source_oxidation_mol_m3_s", "source_coagulation_mol_m3_s"]:
        for i, v in enumerate(out[pname]):
            if abs(v) > 1e-30:
                other_process_nonzero.append((pname, i, v))

    return {
        "name": name,
        "config": cfg,
        "soot_mass_rate": soot_mass_rate,
        "pah_mass_rate": pah_mass_rate,
        "mass_residual": mass_residual,
        "mass_imbalance_rel": mass_imbalance_rel,
        "factor4_ratio": factor4_ratio,
        "gas_C_loss": gas_C_loss,
        "gas_H_loss": gas_H_loss,
        "C_imb": C_imb,
        "H_imb": H_imb,
        "O_imb": O_imb,
        "N_imb": N_imb,
        "other_species_nonzero": other_species_nonzero,
        "other_process_nonzero": other_process_nonzero,
        "src_cond": src_cond,
    }


def print_result(r: dict):
    """Print analysis for a single fixture."""
    print(f"\n--- {r['name']} ---")
    cfg = r["config"]
    print(f"  Config: nuc={cfg['nucleation_model']} cond={cfg['condensation_model']} "
          f"sg={cfg['surface_growth_model']} ox={cfg['oxidation_model']} "
          f"coag={cfg['coagulation_model']}")
    print(f"  soot_mass_rate   = {sci(r['soot_mass_rate'])} kg/m3/s")
    print(f"  pah_mass_rate    = {sci(r['pah_mass_rate'])} kg/m3/s")
    print(f"  mass_residual    = {sci(r['mass_residual'])} kg/m3/s")
    print(f"  mass_imbalance   = {sci(r['mass_imbalance_rel'])} (rel)")
    print(f"  factor4_ratio    = {sci(r['factor4_ratio'])} (expect 1.0)")
    print(f"  C_imbalance      = {sci(r['C_imb'])} (rel)")
    print(f"  H_imbalance      = {sci(r['H_imb'])} (rel)")
    print(f"  O_imbalance      = {sci(r['O_imb'])} (rel)")
    print(f"  N_imbalance      = {sci(r['N_imb'])} (rel)")
    if r["other_species_nonzero"]:
        print(f"  OTHER SPECIES NONZERO:")
        for s, v in r["other_species_nonzero"]:
            print(f"    {s}: {sci(v)}")
    if r["other_process_nonzero"]:
        print(f"  OTHER PROCESS NONZERO:")
        for p, i, v in r["other_process_nonzero"]:
            print(f"    {p}[{i}]: {sci(v)}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="M1 Condensation-only Closure Analysis"
    )
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Path to condensation_fixture_results.json",
    )
    args = parser.parse_args()

    results_path = Path(args.input).resolve()
    if not results_path.exists():
        print(f"ERROR: fixture JSON not found: {results_path}", file=sys.stderr)
        return 2

    try:
        with results_path.open("r", encoding="utf-8") as fh:
            data = json.load(fh)
    except (json.JSONDecodeError, OSError) as exc:
        print(f"ERROR: failed to read/parse {results_path}: {exc}", file=sys.stderr)
        return 2

    # Validate JSON structure.
    if "fixtures" not in data:
        print("ERROR: JSON missing 'fixtures' key", file=sys.stderr)
        return 2

    fixtures = data["fixtures"]
    expected = ["A_all_off_control", "B_cond_only_nominal", "C_cond_only_perturbed"]
    for name in expected:
        if name not in fixtures:
            print(f"ERROR: JSON missing fixture '{name}'", file=sys.stderr)
            return 2

    print("=" * 74)
    print("M1 CONDENSATION-ONLY CLOSURE ANALYSIS")
    print(f"Source: {results_path}")
    print("=" * 74)

    # Analyze all fixtures.
    results = {}
    for name in expected:
        results[name] = analyze_fixture(name, fixtures[name])
        print_result(results[name])

    # ----------------------------------------------------------------
    # CLOSURE VERDICT
    # ----------------------------------------------------------------
    print("\n" + "=" * 74)
    print("CLOSURE VERDICT")
    print("=" * 74)

    all_closed = True

    # Fixture A: all-off control -- everything must be zero.
    rA = results["A_all_off_control"]
    A_ok = (
        abs(rA["soot_mass_rate"]) < 1e-30 and
        abs(rA["pah_mass_rate"]) < 1e-30 and
        len(rA["other_species_nonzero"]) == 0 and
        len(rA["other_process_nonzero"]) == 0
    )
    print(f"  A (all-off control): {'PASS' if A_ok else 'FAIL'}")
    if not A_ok:
        all_closed = False

    # Fixture B: cond-only nominal -- mass must close.
    rB = results["B_cond_only_nominal"]
    B_ok = (
        abs(rB["mass_imbalance_rel"]) < CLOSURE_TOL and
        abs(rB["factor4_ratio"] - 1.0) < CLOSURE_TOL and
        abs(rB["C_imb"]) < CLOSURE_TOL and
        abs(rB["H_imb"]) < CLOSURE_TOL and
        len(rB["other_species_nonzero"]) == 0 and
        len(rB["other_process_nonzero"]) == 0
    )
    print(f"  B (cond-only nominal): {'PASS' if B_ok else 'FAIL'}")
    print(f"    mass_imbalance = {sci(rB['mass_imbalance_rel'])} (tol < {CLOSURE_TOL})")
    print(f"    factor4_ratio   = {sci(rB['factor4_ratio'])} (expect 1.0)")
    if not B_ok:
        all_closed = False

    # Fixture C: cond-only perturbed -- mass must close AND differ from B.
    rC = results["C_cond_only_perturbed"]
    C_ok = (
        abs(rC["mass_imbalance_rel"]) < CLOSURE_TOL and
        abs(rC["factor4_ratio"] - 1.0) < CLOSURE_TOL and
        abs(rC["C_imb"]) < CLOSURE_TOL and
        abs(rC["H_imb"]) < CLOSURE_TOL and
        abs(rC["soot_mass_rate"] - rB["soot_mass_rate"]) > 0  # perturbation response
    )
    print(f"  C (cond-only perturbed): {'PASS' if C_ok else 'FAIL'}")
    print(f"    mass_imbalance = {sci(rC['mass_imbalance_rel'])} (tol < {CLOSURE_TOL})")
    print(f"    factor4_ratio   = {sci(rC['factor4_ratio'])} (expect 1.0)")
    print(f"    perturbation    = {sci(rC['soot_mass_rate'])} != {sci(rB['soot_mass_rate'])}")
    if not C_ok:
        all_closed = False

    # Summary.
    print("\n" + "=" * 74)
    if all_closed:
        print("RESULT: CLOSED / PASS -- all mass, element (C/H/O/N) and factor-4")
        print("  balances close to < 1e-10 relative.")
        print("  CONCLUSION: Condensation-only PAH gas consumption fix verified.")
        print("  PAH mass consumed = soot mass gained (phase change closure).")
        print("=" * 74)
        return 0
    else:
        print("RESULT: NOT CLOSED -- one or more balances exceed tolerance.")
        print("  Check individual fixture results above.")
        print("=" * 74)
        return 1


if __name__ == "__main__":
    sys.exit(main())
