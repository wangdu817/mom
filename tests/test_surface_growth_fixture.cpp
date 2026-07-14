/*-----------------------------------------------------------------------*\
|   MOM Library -- Surface-growth-only fixture harness (M1 vertical slice)  |
|                                                                          |
|   Standalone C++ test harness (no Fluent) that drives                    |
|   HMOM::ComputeSources on a single cell with ONLY surface growth         |
|   active (nucleation OFF, condensation OFF, oxidation OFF,               |
|   coagulation OFF), isolating the HACA surface growth pathway for        |
|   mass/element closure verification.                                     |
|                                                                          |
|   Four fixture variants:                                                 |
|     A. All-off control: every soot process disabled.                     |
|     B. SG-only, radicals ABOVE threshold: H and OH mass fractions        |
|        well above tanh damping thresholds (YH>2e-9, YOH>2e-8).           |
|        Full-strength surface growth.                                     |
|     C. SG-only, radicals BELOW threshold: H and OH mass fractions        |
|        below damping thresholds. Sources heavily damped.                 |
|     D. SG-only, radicals NEAR threshold: H and OH at threshold.          |
|        Partial damping (smooth tanh transition).                         |
|                                                                          |
|   Reaction: C2H2 + soot-surface -> 2C(soot) + H2                        |
|   Gas: C2H2 consumed, H2 produced.                                      |
|   Soot: mass increases (carbon added).                                   |
|   Mass balance: soot_gain = -omega_gas[C2H2] + omega_gas[H2]            |
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
// Thermo builder -- same species set as condensation/nucleation fixtures.
// ============================================================================
static MOM::BasicThermoData buildSGFixtureThermo()
{
    MOM::BasicThermoData th;
    th.names = {"H",    "OH",   "O2",   "H2",  "H2O",  "C2H2", "CO",   "N2"};
    th.mw    = {1.008,  17.008, 31.999, 2.016, 18.015, 26.038, 28.010, 28.014};
    th.nc    = {0,      0,      0,      0,     0,      2,      1,      0};
    th.nh    = {1,      1,      0,      2,     2,      2,      0,      0};
    th.no    = {0,      1,      2,      0,     1,      0,      1,      0};
    th.nn    = {0,      0,      0,      0,     0,      0,      0,      2};
    th.nti   = {0,      0,      0,      0,     0,      0,      0,      0};
    return th;
}

static std::vector<double> X2Y(const std::vector<double>& X, const MOM::BasicThermoData& th)
{
    double MW = 0.;
    for (std::size_t k = 0; k < th.names.size(); ++k) MW += X[k] * th.mw[k];
    std::vector<double> Y(th.names.size());
    for (std::size_t k = 0; k < th.names.size(); ++k) Y[k] = X[k] * th.mw[k] / MW;
    return Y;
}

static double MixtureMW(const std::vector<double>& X, const MOM::BasicThermoData& th)
{
    double MW = 0.;
    for (std::size_t k = 0; k < th.names.size(); ++k) MW += X[k] * th.mw[k];
    return MW;
}

static std::string jnum(double v) { std::ostringstream os; os << std::setprecision(17) << v; return os.str(); }

constexpr double NAV_MOL_  = 6.02214076e23;
constexpr double NAV_KMOL_ = 6.02214076e26;
constexpr double WC_KG_KMOL_ = 12.011;

// ============================================================================
// Per-fixture result structure.
// ============================================================================
struct FixtureResult
{
    std::string label, description;
    int nucleation_model, condensation_model, surface_growth_model;
    int oxidation_model, coagulation_model, continuous_coagulation_model;
    double T_K, P_Pa, mu;
    std::vector<double> X, Y;
    double MW_mix, M00, M10, M01, N0, rho_soot, V0_m3;
    double src_nuc[4], src_growth[4], src_oxid[4], src_cond[4], src_coag[4], src_all[4];
    std::vector<double> omega_gas;
    double soot_mass_rate, c2h2_mass_rate, h2_mass_rate, mass_residual;
    int n_struct_fail;
};

static FixtureResult run_fixture(
    const std::string& label, const std::string& description,
    const MOM::BasicThermoData& thermo,
    int nuc_model, int cond_model, int sg_model, int ox_model, int coag_model, int cont_coag_model,
    double T_K, double P_Pa, double mu,
    const std::vector<double>& X_in,
    double moment_scale)
{
    FixtureResult r;
    r.label = label; r.description = description;
    r.nucleation_model = nuc_model; r.condensation_model = cond_model;
    r.surface_growth_model = sg_model; r.oxidation_model = ox_model;
    r.coagulation_model = coag_model; r.continuous_coagulation_model = cont_coag_model;
    r.T_K = T_K; r.P_Pa = P_Pa; r.mu = mu;

    r.X = X_in;
    { double s=0; for(double x:r.X) s+=x; for(double& x:r.X) x/=s; }
    r.Y = X2Y(r.X, thermo);
    r.MW_mix = MixtureMW(r.X, thermo);

    MOM::HMOM<MOM::BasicThermoData> model(thermo);
    MOM::HMOM<MOM::BasicThermoData>::Config cfg;
    cfg.pah_species = "C2H2";
    cfg.nucleation_model = nuc_model; cfg.condensation_model = cond_model;
    cfg.surface_growth_model = sg_model; cfg.oxidation_model = ox_model;
    cfg.coagulation_model = coag_model; cfg.continuous_coagulation_model = cont_coag_model;
    model.SetupFromConfig(cfg);
    model.SetViscosity(mu);
    model.SetState(T_K, P_Pa, r.Y.data());

    const auto ic = model.initial_moments();
    r.M00 = ic[0] * moment_scale; r.M10 = ic[1] * moment_scale;
    r.M01 = ic[2] * moment_scale; r.N0  = ic[3] * moment_scale;
    model.SetNormalizedMoments(r.M00, r.M10, r.M01, r.N0);
    model.ComputeSources();

    const auto src_all = model.sources();
    const auto omega   = model.omega_gas();
    const auto src_nuc = model.sources_nucleation_impl();
    const auto src_gr  = model.sources_growth_impl();
    const auto src_ox  = model.sources_oxidation_impl();
    const auto src_co  = model.sources_condensation_impl();
    const auto src_cg  = model.sources_coagulation_impl();

    for (int i = 0; i < 4; ++i) {
        r.src_nuc[i]    = src_nuc[static_cast<std::size_t>(i)];
        r.src_growth[i] = src_gr[static_cast<std::size_t>(i)];
        r.src_oxid[i]   = src_ox[static_cast<std::size_t>(i)];
        r.src_cond[i]   = src_co[static_cast<std::size_t>(i)];
        r.src_coag[i]   = src_cg[static_cast<std::size_t>(i)];
        r.src_all[i]    = src_all[static_cast<std::size_t>(i)];
    }
    r.omega_gas.assign(omega.begin(), omega.end());
    r.rho_soot = model.particle_density();
    r.V0_m3    = model.V0();

    const int iC2H2 = 5, iH2 = 3; // species indices
    r.soot_mass_rate  = r.src_growth[1] * r.V0_m3 * NAV_MOL_ * r.rho_soot;
    r.c2h2_mass_rate  = r.omega_gas[static_cast<std::size_t>(iC2H2)];
    r.h2_mass_rate    = r.omega_gas[static_cast<std::size_t>(iH2)];
    // Mass balance: soot_gain = -(omega_gas[C2H2] + omega_gas[H2])
    // C2H2 consumed (negative), H2 produced (positive). Net gas loss = soot gain.
    // soot_mass_rate should equal -(c2h2 + h2) = -c2h2 - h2
    r.mass_residual = r.soot_mass_rate - (-r.c2h2_mass_rate - r.h2_mass_rate);

    r.n_struct_fail = 0;
    auto check = [&](const char* name, bool ok, const std::string& detail) {
        std::cout << "  [STRUCT] " << (ok ? "PASS" : "FAIL") << "  " << name;
        if (!ok) { std::cout << "  -- " << detail; ++r.n_struct_fail; }
        std::cout << "\n";
    };

    std::cout << "\n--- " << label << ": structural assertions ---\n";

    if (sg_model == 0 && nuc_model == 0 && cond_model == 0 && ox_model == 0 && coag_model == 0) {
        for (int i = 0; i < 4; ++i) {
            std::string nm = "src_all[" + std::to_string(i) + "] == 0 (all-off)";
            check(nm.c_str(), r.src_all[i] == 0., "val=" + std::to_string(r.src_all[i]));
        }
        check("omega_gas all zero (all-off)",
              std::all_of(r.omega_gas.begin(), r.omega_gas.end(), [](double v){ return v == 0.; }),
              "nonzero entries present");
    } else {
        check("src_nuc[1] == 0 (nucleation off)",
              r.src_nuc[1] == 0., "val=" + std::to_string(r.src_nuc[1]));
        check("src_cond[1] == 0 (condensation off)",
              r.src_cond[1] == 0., "val=" + std::to_string(r.src_cond[1]));
        check("src_oxid[1] == 0 (oxidation off)",
              r.src_oxid[1] == 0., "val=" + std::to_string(r.src_oxid[1]));
        check("src_coag[1] == 0 (coagulation off)",
              r.src_coag[1] == 0., "val=" + std::to_string(r.src_coag[1]));
    }

    return r;
}

// ============================================================================
// JSON serialization.
// ============================================================================
static std::string names_json(const MOM::BasicThermoData& th) {
    std::string s = "[";
    for (std::size_t k = 0; k < th.names.size(); ++k) { s += "\"" + th.names[k] + "\""; if (k+1<th.names.size()) s+=", "; }
    s += "]"; return s;
}
static std::string vec_json(const std::vector<double>& v) {
    std::string s = "[";
    for (std::size_t k = 0; k < v.size(); ++k) { s += jnum(v[k]); if (k+1<v.size()) s+=", "; }
    s += "]"; return s;
}
static std::string src4_json(const double s[4]) {
    return "[" + jnum(s[0]) + ", " + jnum(s[1]) + ", " + jnum(s[2]) + ", " + jnum(s[3]) + "]";
}

static void write_fixture_json(std::ofstream& jf, const FixtureResult& r, const MOM::BasicThermoData& thermo, bool is_last)
{
    const int iC2H2 = 5, iH2 = 3;
    jf << "  \"" << r.label << "\": {\n";
    jf << "    \"description\": \"" << r.description << "\",\n";
    jf << "    \"config\": {\n";
    jf << "      \"nucleation_model\": " << r.nucleation_model << ",\n";
    jf << "      \"condensation_model\": " << r.condensation_model << ",\n";
    jf << "      \"surface_growth_model\": " << r.surface_growth_model << ",\n";
    jf << "      \"oxidation_model\": " << r.oxidation_model << ",\n";
    jf << "      \"coagulation_model\": " << r.coagulation_model << ",\n";
    jf << "      \"continuous_coagulation_model\": " << r.continuous_coagulation_model << "\n";
    jf << "    },\n";
    jf << "    \"inputs\": {\n";
    jf << "      \"T_K\": " << jnum(r.T_K) << ",\n";
    jf << "      \"P_Pa\": " << jnum(r.P_Pa) << ",\n";
    jf << "      \"mu_kg_ms\": " << jnum(r.mu) << ",\n";
    jf << "      \"species_names\": " << names_json(thermo) << ",\n";
    jf << "      \"mole_fractions\": " << vec_json(r.X) << ",\n";
    jf << "      \"mass_fractions\": " << vec_json(r.Y) << ",\n";
    jf << "      \"MW_mix_kg_kmol\": " << jnum(r.MW_mix) << ",\n";
    jf << "      \"moments_normalized\": {\"M00\": " << jnum(r.M00) << ", \"M10\": " << jnum(r.M10)
       << ", \"M01\": " << jnum(r.M01) << ", \"N0\": " << jnum(r.N0) << "},\n";
    jf << "      \"soot_density_kg_m3\": " << jnum(r.rho_soot) << ",\n";
    jf << "      \"V0_m3\": " << jnum(r.V0_m3) << "\n";
    jf << "    },\n";
    jf << "    \"outputs\": {\n";
    jf << "      \"source_nucleation_mol_m3_s\": " << src4_json(r.src_nuc) << ",\n";
    jf << "      \"source_growth_mol_m3_s\": " << src4_json(r.src_growth) << ",\n";
    jf << "      \"source_oxidation_mol_m3_s\": " << src4_json(r.src_oxid) << ",\n";
    jf << "      \"source_condensation_mol_m3_s\": " << src4_json(r.src_cond) << ",\n";
    jf << "      \"source_coagulation_mol_m3_s\": " << src4_json(r.src_coag) << ",\n";
    jf << "      \"source_all_mol_m3_s\": " << src4_json(r.src_all) << ",\n";
    jf << "      \"omega_gas_kg_m3_s\": " << vec_json(r.omega_gas) << ",\n";
    jf << "      \"soot_mass_rate_kg_m3_s\": " << jnum(r.soot_mass_rate) << ",\n";
    jf << "      \"c2h2_mass_rate_kg_m3_s\": " << jnum(r.c2h2_mass_rate) << ",\n";
    jf << "      \"h2_mass_rate_kg_m3_s\": " << jnum(r.h2_mass_rate) << ",\n";
    jf << "      \"mass_residual_kg_m3_s\": " << jnum(r.mass_residual) << "\n";
    jf << "    },\n";
    jf << "    \"verification\": {\n";
    jf << "      \"n_struct_fail\": " << r.n_struct_fail << ",\n";
    jf << "      \"mass_closure_formula\": \"soot_mass_rate - (-(omega_gas[C2H2]) - omega_gas[H2])\",\n";
    jf << "      \"mass_closure_residual\": " << jnum(r.mass_residual) << ",\n";
    jf << "      \"c2h2_species_index\": " << iC2H2 << ",\n";
    jf << "      \"h2_species_index\": " << iH2 << ",\n";
    jf << "      \"reaction\": \"C2H2 + soot-surface -> 2C(soot) + H2\"\n";
    jf << "    }\n";
    jf << "  }" << (is_last ? "" : ",") << "\n";
}

// ============================================================================
int main()
{
    std::cout << "=== HMOM Surface-growth-only fixture (M1 vertical slice) ===\n";
    const MOM::BasicThermoData thermo = buildSGFixtureThermo();

    const double T_K = 1800., P_Pa = 101325., mu = 4.5e-5, moment_scale = 1.e5;

    // tanh thresholds: YH=2e-9, YOH=2e-8 (mass fractions)
    // Need mass fractions, not mole fractions. MW_mix ~ 26.6 kg/kmol.
    // Y_H = X_H * MW_H / MW_mix. For Y_H = 2e-9: X_H = 2e-9 * 26.6 / 1.008 ~ 5.3e-8
    // For Y_OH = 2e-8: X_OH = 2e-8 * 26.6 / 17.008 ~ 3.1e-8

    // A. All-off control
    std::vector<double> X_A = {
        0.001, 0.001, 1e-30, 0.005, 0.120, 0.020, 0.030, 0.823
    };
    // B. Above threshold: H=0.001 (Y_H ~ 3.8e-5 >> 2e-9), OH=0.001 (Y_OH ~ 6.4e-4 >> 2e-8)
    std::vector<double> X_B = X_A; // same composition, radicals well above threshold
    // C. Below threshold: H and OH at trace levels (Y_H < 2e-9, Y_OH < 2e-8)
    std::vector<double> X_C = {
        1e-10, 1e-10, 1e-30, 0.005, 0.120, 0.020, 0.030, 0.825
    };
    // D. Near threshold: H ~ 5.3e-8 mol frac (Y_H ~ 2e-9), OH ~ 3.1e-8 (Y_OH ~ 2e-8)
    std::vector<double> X_D = {
        5.3e-8, 3.1e-8, 1e-30, 0.005, 0.120, 0.020, 0.030, 0.825
    };

    auto fa = run_fixture("A_all_off_control", "All soot processes OFF.", thermo,
        0,0,0,0,0,0, T_K,P_Pa,mu, X_A, moment_scale);
    auto fb = run_fixture("B_sg_above_threshold", "SG=1, radicals above tanh threshold.", thermo,
        0,0,1,0,0,0, T_K,P_Pa,mu, X_B, moment_scale);
    auto fc = run_fixture("C_sg_below_threshold", "SG=1, radicals below tanh threshold (damped).", thermo,
        0,0,1,0,0,0, T_K,P_Pa,mu, X_C, moment_scale);
    auto fd = run_fixture("D_sg_near_threshold", "SG=1, radicals near tanh threshold (partial).", thermo,
        0,0,1,0,0,0, T_K,P_Pa,mu, X_D, moment_scale);

    // Cross-fixture assertions
    int n_cross = 0;
    auto cross = [&](const char* name, bool ok, const std::string& detail) {
        std::cout << "  [CROSS] " << (ok ? "PASS" : "FAIL") << "  " << name;
        if (!ok) { std::cout << "  -- " << detail; ++n_cross; }
        std::cout << "\n";
    };
    std::cout << "\n--- Cross-fixture assertions ---\n";
    cross("B.src_growth[1] > 0 (above threshold, active)",
          fb.src_growth[1] > 0., "val=" + std::to_string(fb.src_growth[1]));
    cross("B.src_growth[1] > C.src_growth[1] (above > below)",
          fb.src_growth[1] > fc.src_growth[1],
          "B=" + std::to_string(fb.src_growth[1]) + " C=" + std::to_string(fc.src_growth[1]));
    cross("C.src_growth[1] >= 0 (below threshold, damped but non-negative)",
          fc.src_growth[1] >= 0., "val=" + std::to_string(fc.src_growth[1]));
    cross("D.src_growth[1] > C.src_growth[1] (near > below)",
          fd.src_growth[1] > fc.src_growth[1],
          "D=" + std::to_string(fd.src_growth[1]) + " C=" + std::to_string(fc.src_growth[1]));
    cross("B.src_growth[1] > D.src_growth[1] (above > near)",
          fb.src_growth[1] > fd.src_growth[1],
          "B=" + std::to_string(fb.src_growth[1]) + " D=" + std::to_string(fd.src_growth[1]));

    const int total_fail = fa.n_struct_fail + fb.n_struct_fail + fc.n_struct_fail + fd.n_struct_fail + n_cross;
    std::cout << "\n--- Summary ---\n";
    std::cout << "  Total: " << total_fail << " failures\n";

    std::cout << std::scientific << std::setprecision(6);
    for (const auto* r : {&fa, &fb, &fc, &fd}) {
        std::cout << "\n  " << r->label << ":\n";
        std::cout << "    src_growth[1] = " << r->src_growth[1] << "\n";
        std::cout << "    omega[C2H2]   = " << r->c2h2_mass_rate << "\n";
        std::cout << "    omega[H2]     = " << r->h2_mass_rate << "\n";
        std::cout << "    soot_rate     = " << r->soot_mass_rate << "\n";
        std::cout << "    mass_residual = " << r->mass_residual << "\n";
    }

    // Write JSON
    const std::string out_path = "surface_growth_fixture_results.json";
    std::ofstream jf(out_path, std::ios::out | std::ios::trunc);
    if (!jf) { std::cerr << "ERROR: could not open " << out_path << "\n"; return 1; }
    jf << "{\n";
    jf << "  \"task\": \"M1-surface-growth\",\n";
    jf << "  \"description\": \"Surface-growth-only single-cell fixtures (4 variants)\",\n";
    jf << "  \"species_names\": " << names_json(thermo) << ",\n";
    jf << "  \"fixtures\": {\n";
    write_fixture_json(jf, fa, thermo, false);
    write_fixture_json(jf, fb, thermo, false);
    write_fixture_json(jf, fc, thermo, false);
    write_fixture_json(jf, fd, thermo, true);
    jf << "  },\n";
    jf << "  \"closure\": {\n";
    jf << "    \"mass_formula\": \"soot_mass_rate = -omega_gas[C2H2] + omega_gas[H2]\",\n";
    jf << "    \"reaction\": \"C2H2 + soot-surface -> 2C(soot) + H2\",\n";
    jf << "    \"tanh_thresholds\": {\"YH\": 2e-9, \"YOH\": 2e-8}\n";
    jf << "  },\n";
    jf << "  \"total_struct_failures\": " << total_fail << "\n";
    jf << "}\n";
    jf.close();
    if (!jf) { std::cerr << "ERROR: failed writing " << out_path << "\n"; return 1; }
    std::cout << "\nWrote JSON results: " << out_path << "\n";

    if (total_fail > 0) { std::cerr << "ERROR: " << total_fail << " assertion(s) failed.\n"; return 1; }
    std::cout << "=== surface growth fixture completed OK ===\n";
    return 0;
}
