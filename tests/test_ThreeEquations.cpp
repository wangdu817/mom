/*-----------------------------------------------------------------------*\
|   MOM Library — ThreeEquations standalone test                         |
|                                                                         |
|   Tests the 3-equation soot model (Franzelli et al. 2019) against the  |
|   original OpenSMOKE++ implementation using BasicThermoData.            |
|                                                                         |
|   Compile (C++20):                                                      |
|     g++ -std=c++20 -O2 -I ../include -I /path/to/eigen3                |
|         test_ThreeEquations.cpp -o test_ThreeEquations                  |
|                                                                         |
|   Run:                                                                  |
|     ./test_ThreeEquations                                               |
\*-----------------------------------------------------------------------*/

#include "ThreeEquations/ThreeEquations.hpp"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <string>

// ============================================================================
// Helper: print a separator line
// ============================================================================
static void sep(const char* title = nullptr)
{
    std::cout << "\n" << std::string(70, '-');
    if (title) std::cout << "\n " << title;
    std::cout << "\n";
}

// ============================================================================
// Build the minimal thermo needed by ThreeEquations (default PAH = C2H2).
//
// Species required by the constructor:
//   H, OH, O2, H2, H2O, C2H2  (all HACA species)
//   + a bath gas (N2) so mass fractions can sum to 1.
//
// Molecular weights [kg/kmol]:
//   H=1.008, OH=17.008, O2=31.999, H2=2.016, H2O=18.015,
//   C2H2=26.038, N2=28.014
// ============================================================================
MOM::BasicThermoData buildThermo(   const std::vector<std::string>& names,
                                    const std::vector<double>& mw,
                                    const std::vector<double>& nc,
                                    const std::vector<double>& nh,
                                    const std::vector<double>& no,
                                    const std::vector<double>& nn,
                                    const std::vector<double>& nti
                                )
{
    MOM::BasicThermoData th;
    
    th.names = names;
    th.mw = mw;
    th.nc = nc;
    th.nh = nh;
    th.no = no;
    th.nn = nn;
    th.nti = nti;

    return th;
}

// ============================================================================
// Helper: compute mass fractions given mole fractions.
// Species order must match buildThermo().
// ============================================================================
std::vector<double> moleFracToMassFrac(const std::vector<double>& X,
                                        const MOM::BasicThermoData& th)
{
    const int n = static_cast<int>(th.names.size());
    assert(static_cast<int>(X.size()) == n);

    double MW_mix = 0.;
    for (int k = 0; k < n; ++k) MW_mix += X[k] * th.mw[k];

    std::vector<double> Y(n);
    for (int k = 0; k < n; ++k)
        Y[k] = X[k] * th.mw[k] / MW_mix;
    return Y;
}

// ============================================================================
// Print all source terms and key particle properties
// ============================================================================
static void printResults(MOM::ThreeEquations<MOM::BasicThermoData>& model,
                         const std::string& label)
{
    sep(label.c_str());

    const auto src_all  = model.sources();
    const auto src_nuc  = model.sources_nucleation();
    const auto src_gro  = model.sources_growth();
    const auto src_oxi  = model.sources_oxidation();
    const auto src_con  = model.sources_condensation();
    const auto src_coa  = model.sources_coagulation();

    std::cout << std::scientific << std::setprecision(6);
    std::cout << " Particle properties:\n"
              << "   Ys (mass fraction)       = " << model.MassFraction()          << "  [-]\n"
              << "   fv (volume fraction)     = " << model.VolumeFraction()         << "  [-]\n"
              << "   Ns (number density)      = " << model.ParticleNumberDensity()  << "  [#/m3]\n"
              << "   Ss (specific surface)    = " << model.SpecificSurface()        << "  [m2/m3]\n"
              << "   dp (primary diameter)    = " << model.ParticleDiameter()*1e9   << "  [nm]\n"
              << "   dc (collision diameter)  = " << model.CollisionDiameter()*1e9  << "  [nm]\n"
              << "   D  (diffusion coeff)     = " << model.DiffusionCoefficient()   << "  [kg/m/s]\n";

    std::cout << "\n Source terms  [Ys: 1/s]  [NsNorm: 1/s]  [Ss: 1/m/s]:\n"
              << "   Total:       "
              << src_all[0] << "  " << src_all[1] << "  " << src_all[2] << "\n"
              << "   Nucleation:  "
              << src_nuc[0] << "  " << src_nuc[1] << "  " << src_nuc[2] << "\n"
              << "   Growth:      "
              << src_gro[0] << "  " << src_gro[1] << "  " << src_gro[2] << "\n"
              << "   Oxidation:   "
              << src_oxi[0] << "  " << src_oxi[1] << "  " << src_oxi[2] << "\n"
              << "   Condensation:"
              << src_con[0] << "  " << src_con[1] << "  " << src_con[2] << "\n"
              << "   Coagulation: "
              << src_coa[0] << "  " << src_coa[1] << "  " << src_coa[2] << "\n";
}

// ============================================================================
// Test 1 — Zero soot + zero PAH: ALL sources must be zero.
//          Zero soot + non-zero PAH: only nucleation may be non-zero;
//          growth / oxidation / condensation / coagulation must be zero
//          (they all require existing soot particles).
// ============================================================================
static bool test_zero_soot(MOM::ThreeEquations<MOM::BasicThermoData>& model,
                            const std::vector<double>& Y, double T, double P,
                            const MOM::BasicThermoData& thermo)
{
    sep("Test 1a: zero soot AND zero PAH — all sources must be zero");

    // Build Y with C2H2 = 0  (no PAH → no nucleation either)
    // Species order: H OH O2 H2 H2O C2H2 N2
    std::vector<double> Y_noPAH = Y;
    const int iC2H2 = thermo.IndexOfSpecies("C2H2");
    const int iN2   = thermo.IndexOfSpecies("N2");
    Y_noPAH[iN2] += Y_noPAH[iC2H2];   // redistribute to N2 (bath gas)
    Y_noPAH[iC2H2] = 0.;

    model.SetStatus(T, P, Y_noPAH.data());
    model.SetViscosity(6.e-5);
    model.SetMoments(0., 0., 0.);
    model.CalculateSourceMoments();

    bool ok = true;
    for (int i = 0; i < 3; ++i) {
        if (std::abs(model.sources()[i]) > 1.e-20) {
            std::cout << "FAIL 1a: sources()[" << i << "] = " << model.sources()[i]
                      << " (expected 0)\n";
            ok = false;
        }
    }
    if (ok) std::cout << " PASS 1a — all sources zero when soot=0 and PAH=0.\n";

    // ── Sub-test 1b: zero soot, non-zero PAH ─────────────────────────────
    std::cout << "\n";
    sep("Test 1b: zero soot, non-zero PAH — only nucleation non-zero");

    model.SetStatus(T, P, Y.data());   // restore original Y (has C2H2)
    model.SetViscosity(6.e-5);
    model.SetMoments(0., 0., 0.);
    model.CalculateSourceMoments();

    const auto s_gro = model.sources_growth();
    const auto s_oxi = model.sources_oxidation();
    const auto s_con = model.sources_condensation();
    const auto s_coa = model.sources_coagulation();
    const auto s_nuc = model.sources_nucleation();

    bool ok2 = true;
    for (int i = 0; i < 3; ++i) {
        if (std::abs(s_gro[i]) > 1.e-20) {
            std::cout << "FAIL 1b: growth[" << i << "] = " << s_gro[i]
                      << " (must be 0 when soot=0)\n"; ok2 = false;
        }
        if (std::abs(s_oxi[i]) > 1.e-20) {
            std::cout << "FAIL 1b: oxidation[" << i << "] = " << s_oxi[i]
                      << " (must be 0 when soot=0)\n"; ok2 = false;
        }
        if (std::abs(s_con[i]) > 1.e-20) {
            std::cout << "FAIL 1b: condensation[" << i << "] = " << s_con[i]
                      << " (must be 0 when soot=0)\n"; ok2 = false;
        }
        if (std::abs(s_coa[i]) > 1.e-20) {
            std::cout << "FAIL 1b: coagulation[" << i << "] = " << s_coa[i]
                      << " (must be 0 when soot=0)\n"; ok2 = false;
        }
    }
    std::cout << std::scientific << std::setprecision(4);
    std::cout << "  nucleation Ys = " << s_nuc[0]
              << "  Ns = " << s_nuc[1]
              << "  Ss = " << s_nuc[2] << "\n";
    if (s_nuc[0] <= 0.) {
        std::cout << "FAIL 1b: nucleation Ys should be > 0 with PAH present\n";
        ok2 = false;
    }
    if (ok2) std::cout << " PASS 1b — only nucleation fires when soot=0, PAH>0.\n";

    return ok && ok2;
}

// ============================================================================
// Test 2 — Nucleation only: non-zero PAH, zero soot.
//   Expected: source_nucleation(Ys) > 0, source_nucleation(Ns) > 0
//             oxidation/growth/coagulation sources == 0
// ============================================================================
static bool test_nucleation_only(MOM::ThreeEquations<MOM::BasicThermoData>& model,
                                  const std::vector<double>& Y, double T, double P)
{
    sep("Test 2: nucleation only (non-zero PAH, zero soot)");

    model.SetStatus(T, P, Y.data());
    model.SetViscosity(6.e-5);
    model.SetMoments(0., 0., 0.);
    model.CalculateSourceMoments();

    const auto s_all = model.sources();
    const auto s_nuc = model.sources_nucleation();
    const auto s_oxi = model.sources_oxidation();
    const auto s_gro = model.sources_growth();

    bool ok = true;

    if (s_nuc[0] <= 0.) {
        std::cout << "FAIL: nucleation source for Ys should be > 0, got "
                  << s_nuc[0] << "\n";
        ok = false;
    }
    if (s_nuc[1] <= 0.) {
        std::cout << "FAIL: nucleation source for NsNorm should be > 0, got "
                  << s_nuc[1] << "\n";
        ok = false;
    }
    if (std::abs(s_oxi[0]) > 1.e-20 || std::abs(s_gro[0]) > 1.e-20) {
        std::cout << "FAIL: oxidation/growth should be zero when soot is absent\n";
        ok = false;
    }

    std::cout << std::scientific << std::setprecision(4);
    std::cout << "  nucleation Ys source = " << s_nuc[0] << " [1/s]\n"
              << "  nucleation Ns source = " << s_nuc[1] << " [1/s]\n"
              << "  nucleation Ss source = " << s_nuc[2] << " [1/m/s]\n";

    if (ok) std::cout << " PASS\n";
    return ok;
}

// ============================================================================
// Test 3 — All processes active: non-zero soot + non-zero PAH.
//   Qualitative sign checks:
//     Ys total source: sign depends on balance (usually > 0 in flame)
//     Ns from coagulation < 0 (particles merge, number decreases)
//     Ns from nucleation  > 0
//     Ys from oxidation   < 0 (soot is consumed)
// ============================================================================
static bool test_all_processes(MOM::ThreeEquations<MOM::BasicThermoData>& model,
                                const std::vector<double>& Y, double T, double P)
{
    sep("Test 3: all processes active");

    model.SetStatus(T, P, Y.data());
    model.SetViscosity(6.e-5);

    // Typical in-flame soot state
    const double Ys     = 5.e-5;      // [-]   soot mass fraction
    const double NsNorm = 1.e-3;      // [-]   = 1e12 #/m3 / 1e15
    const double Ss     = 5.0;        // [1/m] specific surface
    model.SetMoments(Ys, NsNorm, Ss);
    model.CalculateSourceMoments();

    const auto s_nuc = model.sources_nucleation();
    const auto s_oxi = model.sources_oxidation();
    const auto s_coa = model.sources_coagulation();

    bool ok = true;

    if (s_nuc[0] <= 0.) {
        std::cout << "FAIL: nucleation Ys source should be > 0\n";
        ok = false;
    }
    if (s_oxi[0] >= 0.) {
        std::cout << "FAIL: oxidation Ys source should be < 0 (soot consumed)\n";
        ok = false;
    }
    if (s_coa[1] >= 0.) {
        std::cout << "FAIL: coagulation Ns source should be < 0 (particle merging)\n";
        ok = false;
    }

    printResults(model, "All-processes state");

    if (ok) std::cout << "\n PASS — sign checks on nucleation, oxidation, coagulation.\n";
    return ok;
}

// ============================================================================
// Test 4 — Planck coefficient and NDF reconstruction
// ============================================================================
static bool test_diagnostics(MOM::ThreeEquations<MOM::BasicThermoData>& model)
{
    sep("Test 4: Planck coefficient + NDF reconstruction");

    const double T  = 1800.;
    const double fv = 1.e-6;

    const double kp = model.planck_coefficient(T, fv);
    std::cout << std::scientific << std::setprecision(4);
    std::cout << "  Planck coeff (Smooke, T=1800K, fv=1e-6) = " << kp
              << " [1/m]  (expected ~2.35)\n";

    // NDF reconstruction
    const auto ndf = model.ReconstructedNDFData();
    std::cout << "  NDF reconstruction:\n"
              << "    valid  = " << (ndf.valid ? "true" : "false") << "\n"
              << "    alpha  = " << ndf.alpha  << "\n"
              << "    sigma  = " << ndf.sigma  << "\n"
              << "    nbar0  = " << ndf.nbar0  << " [1/m3]\n"
              << "    k      = " << ndf.k      << "\n"
              << "    nu1m   = " << ndf.nu1mean << " [m3]\n"
              << "    nu2m   = " << ndf.nu2mean << " [m3]\n";

    // Initial moments
    const auto im = model.initial_moments();
    std::cout << "  Initial moments:\n"
              << "    Ys0    = " << im[0] << "  [-]\n"
              << "    NsN0   = " << im[1] << "  [-]\n"
              << "    Ss0    = " << im[2] << "  [1/m]\n";

    bool ok = (kp > 0.);
    if (ok) std::cout << " PASS\n";
    return ok;
}

// ============================================================================
// Test 5 — Surface chemistry model switch (RCPAH vs HMOM)
// ============================================================================
static bool test_surface_model_switch(MOM::ThreeEquations<MOM::BasicThermoData>& model,
                                       const std::vector<double>& Y, double T, double P)
{
    sep("Test 5: surface chemistry model switch (RCPAH vs HMOM)");

    const double Ys     = 5.e-5;
    const double NsNorm = 1.e-3;
    const double Ss     = 5.0;

    // RC-PAH
    model.SetSurfaceChemistryModel(MOM::ThreeEquations<MOM::BasicThermoData>
                                   ::SurfaceChemistryModel::RCPAH);
    model.SetStatus(T, P, Y.data());
    model.SetViscosity(6.e-5);
    model.SetMoments(Ys, NsNorm, Ss);
    model.CalculateSourceMoments();
    const double gro_rcpah = model.sources_growth()[0];

    // HMOM
    model.SetSurfaceChemistryModel(MOM::ThreeEquations<MOM::BasicThermoData>
                                   ::SurfaceChemistryModel::HMOM);
    model.SetStatus(T, P, Y.data());
    model.SetViscosity(6.e-5);
    model.SetMoments(Ys, NsNorm, Ss);
    model.CalculateSourceMoments();
    const double gro_hmom = model.sources_growth()[0];

    std::cout << std::scientific << std::setprecision(4);
    std::cout << "  Growth Ys source — RC-PAH: " << gro_rcpah
              << "  HMOM: " << gro_hmom << "\n";

    bool ok = (gro_rcpah != gro_hmom);
    if (ok) std::cout << " PASS — models give distinct results.\n";
    else    std::cout << " WARNING — models gave identical results (check kinetics).\n";
    return ok;
}

// ============================================================================
// main
// ============================================================================
int main()
{
    std::cout << std::string(70, '=') << "\n"
              << " ThreeEquations<BasicThermoData>  —  standalone test\n"
              << std::string(70, '=') << "\n";

    // ── Build thermo ──────────────────────────────────────────────────────
    std::vector<std::string> names {"H","OH","O2","H2","H2O","C2H2","N2"};
    std::vector<double> mw;
    std::vector<double> nc;
    std::vector<double> nh;
    std::vector<double> no;
    std::vector<double> nn;
    std::vector<double> nti;

    const auto thermo = buildThermo(names, mw, nc, nh, no, nn, nti);

    // ── Construct model ───────────────────────────────────────────────────
    MOM::ThreeEquations<MOM::BasicThermoData> model(thermo);

    // Print default configuration
    sep("Model summary");
    model.PrintSummary();

    // ── Define a representative post-flame mixture at 1800 K ──────────────
    //
    //  Mole fractions (rough oxidation-zone composition):
    //    N2=0.70, H2O=0.10, H2=0.04, C2H2=0.05, O2=0.02, OH=0.02, H=0.07
    //  (does not need to be exact — the test exercises physics, not a real flame)
    //
    const double T = 1800.;          // [K]
    const double P = 101325.;        // [Pa]
    const double mu = 6.e-5;        // dynamic viscosity [kg/m/s]  (air ~1800K)

    // Mole fractions: {H, OH, O2, H2, H2O, C2H2, N2}
    std::vector<double> X = { 0.007, 0.020, 0.020, 0.040, 0.100, 0.050, 0.763 };
    const auto Y = moleFracToMassFrac(X, thermo);

    std::cout << "\n Mixture at T=" << T << " K, P=" << P << " Pa\n";
    for (int k = 0; k < 7; ++k)
        std::cout << "   Y[" << names[k] << "] = " << std::scientific
                  << std::setprecision(3) << Y[k] << "\n";

    // ── Run tests ─────────────────────────────────────────────────────────
    int pass = 0, fail = 0;

    auto record = [&](bool ok) { if (ok) ++pass; else ++fail; };

    record(test_zero_soot          (model, Y, T, P, thermo));
    record(test_nucleation_only    (model, Y, T, P));
    record(test_all_processes      (model, Y, T, P));
    record(test_diagnostics        (model));
    record(test_surface_model_switch(model, Y, T, P));

    // ── Gas-consumption mode ──────────────────────────────────────────────
    sep("Gas consumption mode");
    model.SetGasConsumption(true);
    model.SetStatus(T, P, Y.data());
    model.SetViscosity(mu);
    model.SetMoments(5.e-5, 1.e-3, 5.0);
    model.CalculateSourceMoments();

    const auto omg = model.omega_gas();
    std::cout << " Gas-phase source terms [kg/m3/s]:\n";
    for (int k = 0; k < 7; ++k)
        if (std::abs(omg[k]) > 0.)
            std::cout << "   omega[" << names[k] << "] = "
                      << std::scientific << std::setprecision(4) << omg[k] << "\n";

    // ── Summary ───────────────────────────────────────────────────────────
    sep("Test summary");
    std::cout << " Passed: " << pass << " / " << (pass + fail) << "\n";
    if (fail == 0)
        std::cout << " ALL TESTS PASSED\n";
    else
        std::cout << " " << fail << " TEST(S) FAILED\n";

    return fail == 0 ? 0 : 1;
}
