#!/usr/bin/env python3
"""Closure analysis for M1 coagulation fixture.

Coagulation closure contract (review section 5.4):
  - omega_gas     = 0 for ALL species (no gas-phase source)
  - S_M10         = 0 (soot mass conserved)
  - omega_soot    = rho*V0*Nav*S_M10 = 0
  - S_M00         < 0 (number decreases as particles merge)
  - Other process sources (nuc/growth/oxid/cond) all zero
  - mass_residual = soot_mass_rate + gas_sum = 0

Exit codes: 0=CLOSED, 1=NOT_CLOSED, 2=ERROR
"""

import argparse
import json
import sys
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description="Coagulation closure analysis")
    parser.add_argument("--input", required=True, help="Path to coagulation_fixture_results.json")
    args = parser.parse_args()

    json_path = Path(args.input)
    if not json_path.exists():
        print(f"ERROR: JSON file not found: {json_path}")
        return 2

    try:
        with open(json_path) as f:
            data = json.load(f)
    except (json.JSONDecodeError, OSError) as e:
        print(f"ERROR: Cannot parse JSON: {e}")
        return 2

    required_keys = {"fixtures"}
    if not required_keys.issubset(data.keys()):
        print(f"ERROR: Missing required keys. Expected {required_keys}, got {set(data.keys())}")
        return 2

    fixtures = data["fixtures"]
    if not fixtures:
        print("ERROR: No fixtures in JSON")
        return 2

    atol_mass = 1e-18      # kg/m3/s
    atol_moment = 1e-24    # mol/m3/s
    rtol = 1e-6

    all_closed = True
    print("=" * 72)
    print("COAGULATION CLOSURE ANALYSIS")
    print("=" * 72)

    for fx in fixtures:
        label = fx["label"]
        print(f"\n--- {label} ---")
        print(f"  Description: {fx['description']}")

        coag = fx["src_coag"]
        nuc = fx["src_nuc"]
        growth = fx["src_growth"]
        oxid = fx["src_oxid"]
        cond = fx["src_cond"]
        omega_gas = fx["omega_gas"]
        gas_sum = fx["gas_sum"]
        soot_mass_rate = fx["soot_mass_rate"]
        mass_residual = fx["mass_residual"]

        S_M00 = coag[0]
        S_M10 = coag[1]
        S_M01 = coag[2]
        S_N0 = coag[3]

        is_control = (fx["process_switches"]["coagulation"] == 0 and
                      fx["process_switches"]["continuous_coagulation"] == 0 and
                      fx["process_switches"]["nucleation"] == 0 and
                      fx["process_switches"]["condensation"] == 0 and
                      fx["process_switches"]["surface_growth"] == 0 and
                      fx["process_switches"]["oxidation"] == 0)

        checks = []

        if is_control:
            for i in range(4):
                checks.append(("src_all[{}] == 0".format(i), abs(fx["src_all"][i]) < atol_moment, fx["src_all"][i]))
            for k in range(len(omega_gas)):
                checks.append(("omega_gas[{}] == 0".format(k), abs(omega_gas[k]) < atol_mass, omega_gas[k]))
        else:
            checks.append(("S_M10 == 0 (mass conserved)", abs(S_M10) < atol_moment, S_M10))
            checks.append(("S_M00 < 0 (number decreases)", S_M00 < 0., S_M00))
            checks.append(("omega_gas all zero (gas_sum==0)", abs(gas_sum) < atol_mass, gas_sum))
            checks.append(("soot_mass_rate == 0", abs(soot_mass_rate) < atol_mass, soot_mass_rate))
            checks.append(("mass_residual == 0", abs(mass_residual) < atol_mass, mass_residual))

            for i in range(4):
                checks.append(("S_nuc[{}] == 0".format(i), abs(nuc[i]) < atol_moment, nuc[i]))
                checks.append(("S_growth[{}] == 0".format(i), abs(growth[i]) < atol_moment, growth[i]))
                checks.append(("S_oxid[{}] == 0".format(i), abs(oxid[i]) < atol_moment, oxid[i]))
                checks.append(("S_cond[{}] == 0".format(i), abs(cond[i]) < atol_moment, cond[i]))

            sub_ss = fx.get("src_coag_ss", [0,0,0,0])
            sub_sl = fx.get("src_coag_sl", [0,0,0,0])
            sub_ll = fx.get("src_coag_ll", [0,0,0,0])
            sub_cont_ss = fx.get("src_coag_cont_ss", [0,0,0,0])
            sub_cont_sl = fx.get("src_coag_cont_sl", [0,0,0,0])
            sub_cont_ll = fx.get("src_coag_cont_ll", [0,0,0,0])
            for name, val in [("ss_M10", sub_ss[1]), ("sl_M10", sub_sl[1]), ("ll_M10", sub_ll[1]),
                              ("cont_ss_M10", sub_cont_ss[1]), ("cont_sl_M10", sub_cont_sl[1]), ("cont_ll_M10", sub_cont_ll[1])]:
                checks.append(("sub {} S_M10 == 0".format(name), abs(val) < atol_moment, val))

        fx_closed = True
        for name, ok, val in checks:
            status = "PASS" if ok else "FAIL"
            if not ok:
                fx_closed = False
                all_closed = False
            print(f"  [{status}]  {name}  = {val:.6e}")

        print(f"  => {'CLOSED' if fx_closed else 'NOT CLOSED'}")

    # Cross-fixture checks
    print("\n--- Cross-fixture checks ---")
    active = [f for f in fixtures if f["process_switches"]["coagulation"] > 0]

    if len(active) >= 2:
        b = active[0]
        c = active[1] if len(active) > 1 else active[0]

        if abs(b["src_coag"][0]) > 0:
            ratio_c = abs(c["src_coag"][0]) / abs(b["src_coag"][0])
            ok = ratio_c > 1.0
            status = "PASS" if ok else "FAIL"
            if not ok:
                all_closed = False
            print(f"  [{status}]  |C.S_M00| > |B.S_M00|  (ratio={ratio_c:.4e})")

    if len(active) >= 3:
        b = active[0]
        d = active[-1]
        ratio_d = abs(d["src_coag"][0]) / abs(b["src_coag"][0]) if abs(b["src_coag"][0]) > 0 else 0
        ok = ratio_d >= 0.5
        status = "PASS" if ok else "FAIL"
        if not ok:
            all_closed = False
        print(f"  [{status}]  D.S_M00 comparable to B.S_M00  (ratio={ratio_d:.4e})")

    print("\n" + "=" * 72)
    if all_closed:
        print("VERDICT: CLOSED (all coagulation closure checks PASS)")
        return 0
    else:
        print("VERDICT: NOT_CLOSED")
        return 1


if __name__ == "__main__":
    sys.exit(main())
