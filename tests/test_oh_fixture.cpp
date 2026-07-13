/*-----------------------------------------------------------------------*\
|   MOM Library -- OH-only oxidation fixture harness (task 1.2)             |
|                                                                          |
|   Standalone C++ test harness (no Fluent) that drives                    |
|   HMOM::ComputeSources on a single cell with ONLY the oxidation          |
|   process active AND the O2 channel zeroed out, so that soot oxidation   |
|   proceeds through the OH channel alone.                                 |
|                                                                          |
|   OH-only isolation strategy (see HMOM.tpp:641-642, 1083-1087):          |
|     kox_O2_ = A5 * T^n5 * exp(-E5/T) * conc_O2_ * conc_sootStar          |
|     kox_OH_ = (0.5 / (alpha * chi)) * k6 * conc_OH_                       |
|   The total oxidation rate R_oxid_C2 is split fO2 = kox_O2_/(kox_O2_+    |
|   kox_OH_), fOH = kox_OH_/(kox_O2_+kox_OH_).  Setting the O2 mole         |
|   fraction to a trace (1e-30) drives conc_O2_ -> 0 -> kox_O2_ -> 0 ->     |
|   fO2 -> 0, fOH -> 1.  O2 stays in the thermo (HMOM::SetState looks it   |
|   up by name) with a negligible concentration.                           |
|                                                                          |
|   Representative operating point (sooting laminar flame, post-flame):    |
|     T  = 1800 K, P = 101325 Pa (1 atm), mu = 4.5e-5 kg/(m*s)             |
|     OH-rich, O2-zeroed composition; soot present via scaled moments.     |
|                                                                          |
|   Deliverable: writes mom/build/oh_fixture_results.json (all inputs +     |
|   outputs) for consumption by task 1.3, plus a stdout summary.  No RNG,   |
|   no time-stepping -> deterministic run-to-run.                          |
|                                                                          |
|   Species order (0-based index):                                        |
|     0=H  1=OH  2=O2  3=H2  4=H2O  5=C2H2  6=CO  7=N2                       |
\*-----------------------------------------------------------------------*/

#include "HMOM/HMOM.hpp"

#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <span>
#include <sstream>
#include <string>
#include <vector>

// ============================================================================
// Thermo builder -- mirrors buildSootThermo() in test_verify_all_variants.cpp
// but ADDS CO, which the OH-fixture oxidation channel produces.
//
// HMOM::SetState looks up species BY NAME: H, OH, O2, H2, H2O, C2H2, and the  
// PAH precursor (default "C2H2").  The oxidation product CO is produced by the
// O2 and OH channels, and N2 acts as the inert bulk/closure species.
// ============================================================================
static MOM::BasicThermoData buildOHFixtureThermo()
{
    MOM::BasicThermoData th;
    //            H       OH      O2      H2     H2O     C2H2    CO      N2
    th.names = {"H",    "OH",   "O2",   "H2",  "H2O",  "C2H2", "CO",   "N2"};
    th.mw    = {1.008,  17.008, 31.999, 2.016, 18.015, 26.038, 28.010, 28.014};
    th.nc    = {0,      0,      0,      0,     0,      2,      1,      0};
    th.nh    = {1,      1,      0,      2,     2,      2,      0,      0};
    th.no    = {0,      1,      2,      0,     1,      0,      1,      0};
    th.nn    = {0,      0,      0,      0,     0,      0,      0,      2};
    th.nti   = {0,      0,      0,      0,     0,      0,      0,      0};
    return th;
}

// Mole fractions -> mass fractions helper (mirrors X2Y in the variant test).
static std::vector<double> X2Y(const std::vector<double>& X, const MOM::BasicThermoData& th)
{
    double MW = 0.;
    for (std::size_t k = 0; k < th.names.size(); ++k)
        MW += X[k] * th.mw[k];
    std::vector<double> Y(th.names.size());
    for (std::size_t k = 0; k < th.names.size(); ++k)
        Y[k] = X[k] * th.mw[k] / MW;
    return Y;
}

// Mixture molecular weight from mole fractions [kg/kmol].
static double MixtureMW(const std::vector<double>& X, const MOM::BasicThermoData& th)
{
    double MW = 0.;
    for (std::size_t k = 0; k < th.names.size(); ++k)
        MW += X[k] * th.mw[k];
    return MW;
}

// Format a double with full round-trip precision for a deterministic JSON file.
// (17 significant digits is enough to round-trip any IEEE-754 double.)
static std::string jnum(double v)
{
    std::ostringstream os;
    os << std::setprecision(17) << v;
    return os.str();
}

int main()
{
    std::cout << "=== HMOM OH-only oxidation fixture (task 1.2) ===\n";

    // -- Build the thermo.  It MUST outlive the HMOM model (passed by const ref). --
    const MOM::BasicThermoData thermo = buildOHFixtureThermo();

    std::cout << "Species (" << thermo.NumberOfSpecies() << "): ";
    for (std::size_t k = 0; k < thermo.names.size(); ++k)
        std::cout << thermo.names[k] << (k + 1 < thermo.names.size() ? " " : "\n");

    // -- Construct the HMOM model bound to the thermo. ------------------------
    MOM::HMOM<MOM::BasicThermoData> model(thermo);

    // -- Oxidation-only configuration. ---------------------------------------
    // Only the oxidation process is active; every other process is disabled so
    // the moment sources isolate the OH oxidation contribution.
    MOM::HMOM<MOM::BasicThermoData>::Config cfg;
    cfg.pah_species                  = "C2H2"; // present in the thermo
    cfg.nucleation_model             = 0;
    cfg.surface_growth_model         = 0;
    cfg.condensation_model           = 0;
    cfg.coagulation_model            = 0;
    cfg.continuous_coagulation_model = 0;
    cfg.oxidation_model              = 1;      // <-- ONLY oxidation active
    model.SetupFromConfig(cfg);

    // -- Representative operating point (sooting flame post-flame region). ----
    const double T_K   = 1800.;   // [K]  post-flame zone (soot oxidation dominates)
    const double P_Pa  = 101325.; // [Pa] 1 atm
    const double mu    = 4.5e-5;  // [kg/(m*s)] viscosity

    // Viscosity is not part of Config; set it explicitly BEFORE
    // ComputeSources() (matches the 1.1 harness call order).
    model.SetViscosity(mu);

    // -- Composition (mole fractions): OH-rich, O2-zeroed. --------------------
    //   H     OH    O2(trace)  H2    H2O   C2H2  CO    N2(balance)
    // O2 is set to a trace (1e-30) instead of exactly 0 so that conc_O2_ is
    // negligible (kox_O2_ ~ 0 -> fO2 ~ 0, fOH ~ 1) while avoiding any
    // divide-by-zero edge cases; O2 remains present in the thermo.
    std::vector<double> X = {
        0.001,   // H
        0.020,   // OH
        1.0e-30, // O2 (trace)
        0.005,   // H2
        0.120,   // H2O
        0.002,   // C2H2
        0.030,   // CO
        0.822    // N2 (balance)
    };

    // Normalise so the mole fractions sum to exactly 1.0 (deterministic).
    {
        double sumX = 0.;
        for (double x : X)
            sumX += x;
        for (double& x : X)
            x /= sumX;
    }

    const std::vector<double> Y      = X2Y(X, thermo);
    const double              MW_mix = MixtureMW(X, thermo);

    // -- Call order: SetState -> SetNormalizedMoments -> ComputeSources
    model.SetState(T_K, P_Pa, Y.data());

    // Use the floor initial moments scaled up so soot is present for oxidation.
    // These represent normalised moments [mol/m3]; documented exactly below.
    const auto   ic    = model.initial_moments();
    const double scale = 1.e5;
    const double M00   = ic[0] * scale;
    const double M10   = ic[1] * scale;
    const double M01   = ic[2] * scale;
    const double N0    = ic[3] * scale;
    model.SetNormalizedMoments(M00, M10, M01, N0);

    model.ComputeSources();

    // -- Read back the 4-element moment source span. -------------------------
    const std::span<const double> src = model.sources();

    // -- Read back the gas-phase consumption vector. -------------------------
    const std::span<const double> omega = model.omega_gas();

    // -- Species indices for readback / verification. ------------------------
    const int iH    = thermo.IndexOfSpecies("H");
    const int iOH   = thermo.IndexOfSpecies("OH");
    const int iO2   = thermo.IndexOfSpecies("O2");
    const int iH2   = thermo.IndexOfSpecies("H2");
    const int iH2O  = thermo.IndexOfSpecies("H2O");
    const int iC2H2 = thermo.IndexOfSpecies("C2H2");
    const int iCO   = thermo.IndexOfSpecies("CO");
    const int iN2   = thermo.IndexOfSpecies("N2");

    const double omega_OH = omega[static_cast<std::size_t>(iOH)];
    const double omega_O2 = omega[static_cast<std::size_t>(iO2)];
    const double omega_CO = omega[static_cast<std::size_t>(iCO)];

    // -- Derive R_oxid_C2 and fOH from the PUBLIC omega_gas vector. -----------
    // kox_O2_/kox_OH_ are private with no accessor (HMOM.hpp:1075-1076), so the
    // internal split fractions are not directly observable.  However, with ONLY
    // oxidation active the OH channel is the sole consumer of OH and the O2
    // channel the sole consumer of O2 (HMOM.tpp:1091-1129).
    //
    // O2 channel stoichiometry (HMOM.tpp:1089-1100):
    //   C2 + O2 -> 2CO   (1 O2 consumed per C2)
    //   omega[O2] = -R_O2 * 1 * MW_O2 / 1000
    //   R_O2 = -omega[O2] * 1000 / MW_O2
    //
    // OH channel stoichiometry (HMOM.tpp:1103-1123):
    //   C2 + 2OH -> 2CO + H2   (2 OH consumed per C2)
    //   omega[OH] = -R_OH * 2 * MW_OH / 1000
    //   R_OH = -omega[OH] * 1000 / (2 * MW_OH)
    //
    // R_OH = fOH * R_oxid_C2, R_O2 = fO2 * R_oxid_C2, so
    //   R_oxid_C2 = R_O2 + R_OH   and   fOH = R_OH / R_oxid_C2.
    const double MW_OH = thermo.MolecularWeight(static_cast<unsigned>(iOH));
    const double MW_O2 = thermo.MolecularWeight(static_cast<unsigned>(iO2));
    const double R_OH_channel = -omega_OH * 1000. / (2.0 * MW_OH); // C2 + 2OH -> 2CO + H2
    const double R_O2_channel = -omega_O2 * 1000. / MW_O2;         // C2 + O2  -> 2CO
    const double R_oxid_C2    = R_O2_channel + R_OH_channel;
    const bool   fOH_observable = (R_oxid_C2 > 0.);
    const double fOH = fOH_observable ? (R_OH_channel / R_oxid_C2) : 0.;

    // O2-to-OH ratio for the OH-only verification block.
    const double O2_to_OH_ratio =
        (omega_OH != 0.) ? (omega_O2 / omega_OH) : 0.;

    // Nucleation volume (public accessor) - recorded for task 1.3.
    const double V0_m3 = model.V0();

    // Molar consumption rates for structural cross-checks [mol/m3/s].
    const double OH_consumption_molar = -omega_OH * 1000. / MW_OH;

    // ---------------------------------------------------------------------
    // Structural assertions (M1-0 requirement: fixture must self-verify
    // key invariants, not rely solely on Python post-processing).
    // ---------------------------------------------------------------------
    // Physical constants matching MomentMethodBase.hpp (verified values).
    constexpr double NAV_KMOL_ = 6.02214076e26;  // [1/kmol] MomentMethodBase.hpp:489
    constexpr double WC_KG_KMOL_ = 12.011;        // [kg/kmol] MomentMethodBase.hpp:491
    const double rho_soot_val = model.particle_density();
    const double VC2_check = (WC_KG_KMOL_ / rho_soot_val / NAV_KMOL_) * 2.0;

    // Soot-side R_oxid_C2 (independent of gas-phase stoichiometry).
    const double R_oxid_C2_soot = -src[1] * V0_m3 / VC2_check;  // [mol C2/m3/s]

    // Gas-side per-C2 ratios (must match C2 + 2OH -> 2CO + H2).
    const double MW_CO_v = thermo.MolecularWeight(static_cast<unsigned>(iCO));
    const double MW_H2_v = thermo.MolecularWeight(static_cast<unsigned>(iH2));
    const double OH_per_C2_check =
        (R_oxid_C2_soot > 0.) ? (OH_consumption_molar / R_oxid_C2_soot) : 0.;
    const double CO_per_C2_check =
        (R_oxid_C2_soot > 0.) ? ((omega_CO * 1000. / MW_CO_v) / R_oxid_C2_soot) : 0.;
    const double H2_per_C2_check =
        (R_oxid_C2_soot > 0.) ? ((omega[static_cast<std::size_t>(iH2)] * 1000. / MW_H2_v) / R_oxid_C2_soot) : 0.;

    // Tolerances for structural checks.
    constexpr double STRUCT_TOL = 1e-6;   // relative
    constexpr double TRACE_TOL  = 1e-20;  // absolute (O2 trace)

    int n_fail = 0;
    auto check = [&](const char* name, bool ok, const std::string& detail) {
        std::cout << "  [STRUCT] " << (ok ? "PASS" : "FAIL") << "  " << name;
        if (!ok) { std::cout << "  -- " << detail; ++n_fail; }
        std::cout << "\n";
    };

    std::cout << "\n--- STRUCTURAL ASSERTIONS (C++ self-check) ---\n";
    check("omega_OH < 0 (OH consumed)",      omega_OH < 0.,
          "omega_OH=" + std::to_string(omega_OH));
    check("omega_O2 ~ 0 (OH-only trace)",    std::abs(omega_O2) < TRACE_TOL,
          "omega_O2=" + std::to_string(omega_O2));
    check("omega_CO > 0 (CO produced)",      omega_CO > 0.,
          "omega_CO=" + std::to_string(omega_CO));
    check("S_M10 < 0 (soot mass loss)",      src[1] < 0.,
          "S_M10=" + std::to_string(src[1]));
    check("R_oxid_C2 > 0 (oxidation active)", R_oxid_C2 > 0.,
          "R_oxid_C2=" + std::to_string(R_oxid_C2));
    check("fOH ~ 1.0 (OH-only channel)",     std::abs(fOH - 1.0) < 0.01,
          "fOH=" + std::to_string(fOH));
    check("OH_per_C2 ~ 2.0 (C2+2OH->2CO+H2)", std::abs(OH_per_C2_check - 2.0) < STRUCT_TOL,
          "OH_per_C2=" + std::to_string(OH_per_C2_check));
    check("CO_per_C2 ~ 2.0 (C2+2OH->2CO+H2)", std::abs(CO_per_C2_check - 2.0) < STRUCT_TOL,
          "CO_per_C2=" + std::to_string(CO_per_C2_check));
    check("H2_per_C2 ~ 1.0 (C2+2OH->2CO+H2)", std::abs(H2_per_C2_check - 1.0) < STRUCT_TOL,
          "H2_per_C2=" + std::to_string(H2_per_C2_check));
    // Cross-check: gas-side R_oxid_C2 must match soot-side R_oxid_C2.
    check("R_oxid_C2 gas == soot (cross-check)",
          std::abs(R_oxid_C2 - R_oxid_C2_soot) / std::max(R_oxid_C2_soot, 1e-30) < STRUCT_TOL,
          "gas=" + std::to_string(R_oxid_C2) + " soot=" + std::to_string(R_oxid_C2_soot));

    if (n_fail > 0)
    {
        std::cerr << "ERROR: " << n_fail << " structural assertion(s) failed.\n";
        return 1;
    }
    std::cout << "  All structural assertions passed.\n\n";

    // ---------------------------------------------------------------------
    // Human-readable stdout summary (keep the existing print style).
    // ---------------------------------------------------------------------
    std::cout << std::scientific << std::setprecision(6);

    std::cout << "State: T=" << T_K << " K, P=" << P_Pa << " Pa, mu=" << mu << " kg/m/s\n";
    std::cout << "MW_mix=" << MW_mix << " kg/kmol\n";
    std::cout << "Moments (normalized): M00=" << M00 << " M10=" << M10
              << " M01=" << M01 << " N0=" << N0 << "\n";

    std::cout << "moment_sources[4] (mol/m3/s) =";
    for (std::size_t i = 0; i < src.size(); ++i)
        std::cout << " " << src[i];
    std::cout << "\n";

    std::cout << "omega_gas[" << omega.size() << "] (kg/m3/s) =";
    for (std::size_t i = 0; i < omega.size(); ++i)
        std::cout << " " << thermo.names[i] << "=" << omega[i];
    std::cout << "\n";

    std::cout << "fOH=" << (fOH_observable ? jnum(fOH) : std::string("not-observable"))
              << " (should be ~1.0 for OH-only)\n";
    std::cout << "R_oxid_C2=" << R_oxid_C2 << " mol/m3/s\n";
    std::cout << "OH-only check: omega_O2=" << omega_O2 << " omega_OH=" << omega_OH
              << " ratio(O2/OH)=" << O2_to_OH_ratio << "\n";

    // Structured result line for downstream parsing.
    std::cout << "RESULT: oh_only pipeline ran"
              << " n_moment_sources=" << src.size()
              << " n_omega_gas=" << omega.size()
              << " omega_OH=" << omega_OH
              << " omega_O2=" << omega_O2
              << " omega_CO=" << omega_CO
              << " fOH=" << fOH
              << "\n";

    // ---------------------------------------------------------------------
    // Write the deterministic JSON results file (primary deliverable).
    //
    // Path is relative ("oh_fixture_results.json"): the executable is invoked
    // from the build directory (both by the manual run and by ctest, whose
    // WORKING_DIRECTORY defaults to the build tree), so the file lands at
    // mom/build/oh_fixture_results.json as required.
    // ---------------------------------------------------------------------
    const std::string out_path = "oh_fixture_results.json";
    std::ofstream jf(out_path, std::ios::out | std::ios::trunc);
    if (!jf)
    {
        std::cerr << "ERROR: could not open " << out_path << " for writing\n";
        return 1;
    }

    auto names_json = [&]() {
        std::string s = "[";
        for (std::size_t k = 0; k < thermo.names.size(); ++k)
        {
            s += "\"" + thermo.names[k] + "\"";
            if (k + 1 < thermo.names.size())
                s += ", ";
        }
        s += "]";
        return s;
    };
    auto vec_json = [&](const std::vector<double>& v) {
        std::string s = "[";
        for (std::size_t k = 0; k < v.size(); ++k)
        {
            s += jnum(v[k]);
            if (k + 1 < v.size())
                s += ", ";
        }
        s += "]";
        return s;
    };

    jf << "{\n";
    jf << "  \"task\": \"1.2\",\n";
    jf << "  \"description\": \"OH-only oxidation single-cell fixture\",\n";
    jf << "  \"inputs\": {\n";
    jf << "    \"T_K\": " << jnum(T_K) << ",\n";
    jf << "    \"P_Pa\": " << jnum(P_Pa) << ",\n";
    jf << "    \"mu_kg_ms\": " << jnum(mu) << ",\n";
    jf << "    \"species_names\": " << names_json() << ",\n";
    jf << "    \"mole_fractions\": " << vec_json(X) << ",\n";
    jf << "    \"mass_fractions\": " << vec_json(Y) << ",\n";
    jf << "    \"MW_mix_kg_kmol\": " << jnum(MW_mix) << ",\n";
    jf << "    \"moments_normalized\": {\"M00\": " << jnum(M00) << ", \"M10\": " << jnum(M10)
       << ", \"M01\": " << jnum(M01) << ", \"N0\": " << jnum(N0) << "},\n";
    jf << "    \"soot_density_kg_m3\": " << jnum(model.particle_density()) << ",\n";
    jf << "    \"V0_m3\": " << jnum(V0_m3) << ",\n";
    jf << "    \"config\": {\"oxidation_model\": 1, \"all_others\": 0}\n";
    jf << "  },\n";

    jf << "  \"outputs\": {\n";
    jf << "    \"moment_sources_mol_m3_s\": {\"S_M00\": " << jnum(src[0]) << ", \"S_M10\": "
       << jnum(src[1]) << ", \"S_M01\": " << jnum(src[2]) << ", \"S_N0\": " << jnum(src[3])
       << "},\n";
    jf << "    \"omega_gas_kg_m3_s\": {";
    jf << "\"H\": " << jnum(omega[static_cast<std::size_t>(iH)]) << ", ";
    jf << "\"OH\": " << jnum(omega[static_cast<std::size_t>(iOH)]) << ", ";
    jf << "\"O2\": " << jnum(omega[static_cast<std::size_t>(iO2)]) << ", ";
    jf << "\"H2\": " << jnum(omega[static_cast<std::size_t>(iH2)]) << ", ";
    jf << "\"H2O\": " << jnum(omega[static_cast<std::size_t>(iH2O)]) << ", ";
    jf << "\"C2H2\": " << jnum(omega[static_cast<std::size_t>(iC2H2)]) << ", ";
    jf << "\"CO\": " << jnum(omega[static_cast<std::size_t>(iCO)]) << ", ";
    jf << "\"N2\": " << jnum(omega[static_cast<std::size_t>(iN2)]) << "},\n";
    // kox_O2_ / kox_OH_ are private with no public accessor -> not observable.
    jf << "    \"kox_O2\": null,\n";
    jf << "    \"kox_OH\": null,\n";
    jf << "    \"fOH\": " << (fOH_observable ? jnum(fOH) : std::string("null")) << ",\n";
    jf << "    \"R_oxid_C2_mol_m3_s\": " << jnum(R_oxid_C2) << "\n";
    jf << "  },\n";

    jf << "  \"verification\": {\n";
    jf << "    \"OH_only_check\": \"omega_O2 should be ~0 relative to omega_OH\",\n";
    jf << "    \"omega_O2_kg_m3_s\": " << jnum(omega_O2) << ",\n";
    jf << "    \"omega_OH_kg_m3_s\": " << jnum(omega_OH) << ",\n";
    jf << "    \"O2_to_OH_ratio\": " << jnum(O2_to_OH_ratio) << "\n";
    jf << "  },\n";

    jf << "  \"provenance\": {\n";
    jf << "    \"kox_O2\": \"private member kox_O2_ (HMOM.hpp:1075), no public accessor -> null\",\n";
    jf << "    \"kox_OH\": \"private member kox_OH_ (HMOM.hpp:1076), no public accessor -> null\",\n";
    jf << "    \"fOH\": \"DERIVED from public omega_gas(): fOH = R_OH/(R_O2+R_OH), "
          "R_OH=-omega[OH]*1000/(2*MW_OH) [C2+2OH->2CO+H2], "
          "R_O2=-omega[O2]*1000/MW_O2 [C2+O2->2CO]\",\n";
    jf << "    \"R_oxid_C2\": \"DERIVED from public omega_gas() as R_O2+R_OH; equivalently "
          "-source_oxidation(1)*V0/VC2 but VC2_ is private (HMOM.hpp:1066)\",\n";
    jf << "    \"per_process_source_split\": \"source_oxidation() IS public "
          "(MomentMethodBase::sources_oxidation), but VC2_ scaling is private\"\n";
    jf << "  },\n";

    jf << "  \"reproducibility\": \"deterministic -- no RNG, no time-stepping\"\n";
    jf << "}\n";

    jf.close();
    if (!jf)
    {
        std::cerr << "ERROR: failed while writing/closing " << out_path << "\n";
        return 1;
    }

    std::cout << "Wrote JSON results: " << out_path << "\n";
    std::cout << "=== fixture completed OK ===\n";
    return 0;
}
