#!/usr/bin/env python3
"""Closure analysis for M3 Cunningham diffusion fixture.

Cunningham diffusion closure contract:
  - Gamma_model  = model.diffusion_coefficient()
  - Gamma_indep  = max(rho * D_Cunningham, mu / Sc)
  - D_Cunningham = kB*T*Cu / (3*pi*mu*dc_safe)
  - Cu           = 1 + kCunningham * lambda / dc_safe
  - lambda       = (mu/rho) * sqrt(pi * m_gas / (2*kB*T))
  - m_gas        = rho * kB * T / P
  - dc_safe      = max(collision_diameter, 1e-12)
  - rho          = P / (Rgas * T) * MW_mix

The Python script independently recomputes Gamma from first principles
using the same constants and state as the C++ model, then compares.

Exit codes: 0=CLOSED, 1=NOT_CLOSED, 2=ERROR
"""

import argparse
import json
import math
import sys
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description="Cunningham diffusion closure analysis")
    parser.add_argument("--input", required=True, help="Path to cunningham_diffusion_fixture_results.json")
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

    required_keys = {"constants", "fixtures"}
    if not required_keys.issubset(data.keys()):
        print(f"ERROR: Missing required keys. Expected {required_keys}, got {set(data.keys())}")
        return 2

    c = data["constants"]
    kB = c["kB"]
    Rgas = c["Rgas"]
    kCun = c["kCunningham"]
    pi = c["pi"]

    fixtures = data["fixtures"]
    if not fixtures:
        print("ERROR: No fixtures in JSON")
        return 2

    atol = 1e-12
    rtol = 1e-7

    all_closed = True
    print("=" * 72)
    print("CUNNINGHAM DIFFUSION CLOSURE ANALYSIS")
    print("=" * 72)

    indep_results = []

    for fx in fixtures:
        label = fx["label"]
        print(f"\n--- {label} ---")
        print(f"  Description: {fx['description']}")

        T = fx["T_K"]
        P = fx["P_Pa"]
        mu = fx["mu"]
        MW_mix = fx["MW_mix"]
        Sc = fx["schmidt_number"]
        dc_model = fx["collision_diameter"]
        Gamma_model = fx["diffusion_coefficient"]

        rho = P / (Rgas * T) * MW_mix
        dc_safe = max(dc_model, 1e-12)

        m_gas = rho * kB * T / P
        lambda_mfp = (mu / rho) * math.sqrt(pi * m_gas / (2.0 * kB * T))
        Cu = 1.0 + kCun * lambda_mfp / dc_safe
        D_cun = kB * T * Cu / (3.0 * pi * mu * dc_safe)
        Gamma_brownian = rho * D_cun
        Gamma_sc = mu / Sc
        Gamma_indep = max(Gamma_brownian, Gamma_sc)

        indep_results.append({
            "label": label,
            "rho": rho,
            "dc_safe": dc_safe,
            "m_gas": m_gas,
            "lambda": lambda_mfp,
            "Cu": Cu,
            "D_cun": D_cun,
            "Gamma_brownian": Gamma_brownian,
            "Gamma_sc": Gamma_sc,
            "Gamma_indep": Gamma_indep,
            "Gamma_model": Gamma_model,
            "dc_model": dc_model,
            "moment_scale": fx["moment_scale"],
        })

        checks = []

        checks.append(("rho > 0", rho > 0, rho))
        checks.append(("m_gas > 0", m_gas > 0, m_gas))
        checks.append(("lambda > 0", lambda_mfp > 0, lambda_mfp))
        checks.append(("Cu >= 1 (Cunningham >= 1)", Cu >= 1.0, Cu))
        checks.append(("D_Cunningham > 0", D_cun > 0, D_cun))
        checks.append(("Gamma_brownian > 0", Gamma_brownian > 0, Gamma_brownian))
        checks.append(("Gamma_sc > 0", Gamma_sc > 0, Gamma_sc))
        checks.append(("Gamma_indep > 0", Gamma_indep > 0, Gamma_indep))
        checks.append(("Gamma_model > 0", Gamma_model > 0, Gamma_model))

        if abs(Gamma_indep) > 0:
            rel_err = abs(Gamma_model - Gamma_indep) / abs(Gamma_indep)
        else:
            rel_err = abs(Gamma_model - Gamma_indep)
        checks.append(("Gamma_model ~= Gamma_indep (rel_err < rtol)",
                       rel_err < rtol, rel_err))

        if dc_model > 0:
            checks.append(("dc_model > 0 (non-degenerate)", dc_model > 0, dc_model))
            checks.append(("Cu > 1 (non-trivial Cunningham for finite dc)", Cu > 1.0, Cu))
        else:
            checks.append(("dc_model == 0 (degenerate, floor 1e-12)", dc_model == 0, dc_model))
            checks.append(("dc_safe == 1e-12 (floor applied)", abs(dc_safe - 1e-12) < 1e-20, dc_safe))
            checks.append(("Cu >> 1 (tiny dc_safe -> huge Cunningham)", Cu > 1e3, Cu))

        if Gamma_brownian > Gamma_sc:
            checks.append(("Gamma = Gamma_brownian (Brownian dominates)", True, Gamma_brownian))
        else:
            checks.append(("Gamma = Gamma_sc (Schmidt dominates)", True, Gamma_sc))

        fx_closed = True
        for name, ok, val in checks:
            status = "PASS" if ok else "FAIL"
            if not ok:
                fx_closed = False
                all_closed = False
            print(f"  [{status}]  {name}  = {val:.6e}")

        print(f"  rho       = {rho:.6e} kg/m3")
        print(f"  dc_safe   = {dc_safe:.6e} m")
        print(f"  m_gas     = {m_gas:.6e} kg")
        print(f"  lambda    = {lambda_mfp:.6e} m")
        print(f"  Cu        = {Cu:.6e}")
        print(f"  D_cun     = {D_cun:.6e} m2/s")
        print(f"  Gamma_br  = {Gamma_brownian:.6e}")
        print(f"  Gamma_sc  = {Gamma_sc:.6e}")
        print(f"  Gamma_in  = {Gamma_indep:.6e}")
        print(f"  Gamma_mod = {Gamma_model:.6e}")
        print(f"  => {'CLOSED' if fx_closed else 'NOT CLOSED'}")

    print("\n--- Cross-fixture checks ---")

    by_label = {r["label"]: r for r in indep_results}

    if "A_flame_small_soot" in by_label and "B_flame_large_soot" in by_label:
        A = by_label["A_flame_small_soot"]
        B = by_label["B_flame_large_soot"]

        ok = A["dc_safe"] < B["dc_safe"]
        status = "PASS" if ok else "FAIL"
        if not ok: all_closed = False
        print(f"  [{status}]  A.dc < B.dc (small soot < large soot)  A={A['dc_safe']:.4e} B={B['dc_safe']:.4e}")

        ok = A["Cu"] > B["Cu"]
        status = "PASS" if ok else "FAIL"
        if not ok: all_closed = False
        print(f"  [{status}]  A.Cu > B.Cu (stronger Cunningham for small soot)  A={A['Cu']:.4e} B={B['Cu']:.4e}")

        ok = A["Gamma_indep"] > B["Gamma_indep"]
        status = "PASS" if ok else "FAIL"
        if not ok: all_closed = False
        print(f"  [{status}]  A.Gamma > B.Gamma (small soot diffuses faster)  A={A['Gamma_indep']:.4e} B={B['Gamma_indep']:.4e}")

        ok = A["lambda"] > 0 and B["lambda"] > 0
        status = "PASS" if ok else "FAIL"
        if not ok: all_closed = False
        print(f"  [{status}]  lambda > 0 for both A and B  A={A['lambda']:.4e} B={B['lambda']:.4e}")

    if "A_flame_small_soot" in by_label and "C_room_medium_soot" in by_label:
        A = by_label["A_flame_small_soot"]
        C = by_label["C_room_medium_soot"]

        ok = abs(A["lambda"] - C["lambda"]) > 0
        status = "PASS" if ok else "FAIL"
        if not ok: all_closed = False
        print(f"  [{status}]  A.lambda != C.lambda (different T -> different mfp)  A={A['lambda']:.4e} C={C['lambda']:.4e}")

    if "A_flame_small_soot" in by_label and "D_degenerate" in by_label:
        A = by_label["A_flame_small_soot"]
        D = by_label["D_degenerate"]

        ok = D["Cu"] > A["Cu"] * 1e2
        status = "PASS" if ok else "FAIL"
        if not ok: all_closed = False
        print(f"  [{status}]  D.Cu >> A.Cu (tiny dc_safe -> huge Cu)  D={D['Cu']:.4e} A={A['Cu']:.4e}")

        ok = D["Gamma_indep"] > A["Gamma_indep"] * 1e3
        status = "PASS" if ok else "FAIL"
        if not ok: all_closed = False
        print(f"  [{status}]  D.Gamma >> A.Gamma (tiny dc_safe -> huge Gamma)  D={D['Gamma_indep']:.4e} A={A['Gamma_indep']:.4e}")

    print("\n" + "=" * 72)
    if all_closed:
        print("VERDICT: CLOSED (all Cunningham diffusion closure checks PASS)")
        return 0
    else:
        print("VERDICT: NOT_CLOSED")
        return 1


if __name__ == "__main__":
    sys.exit(main())
