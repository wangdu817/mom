/*-----------------------------------------------------------------------*\
|   MOM Library -- O2 oxidation-only fixture harness (M1 vertical slice)    |
|                                                                          |
|   Standalone C++ test harness (no Fluent) that drives                    |
|   HMOM::ComputeSources on a single cell with ONLY oxidation              |
|   active (nucleation OFF, condensation OFF, surface growth OFF,          |
|   coagulation OFF), isolating the O2 oxidation pathway for               |
|   mass/element closure verification.                                     |
|                                                                          |
|   Four fixture variants:                                                 |
|     A. All-off control: every soot process disabled.                     |
|     B. O2-ox only, radicals ABOVE threshold: H and OH mass fractions     |
|        well above tanh damping thresholds. High O2, low OH -> fO2 ~ 1.   |
|     C. O2-ox only, radicals BELOW threshold: H and OH below damping.     |
|        Sources heavily damped by tanh.                                   |
|     D. O2-ox only, radicals NEAR threshold: partial damping.             |
|                                                                          |
|   Reaction (O2 channel): C2(soot) + O2 -> 2CO                            |
|   Gas: O2 consumed, CO produced.                                         |
|   Soot: mass decreases (carbon removed).                                 |
|   Mass balance: -soot_mass_rate = omega_gas[O2] + omega_gas[CO]          |
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

static MOM::BasicThermoData buildOxFixtureThermo()
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
    double soot_mass_rate;
    double o2_mass_rate, co_mass_rate, oh_mass_rate, h2_mass_rate;
    double mass_residual;
    double R_oxid_C2;
    double fO2, fOH;
    double O2_per_C2, CO_per_C2, OH_per_C2, H2_per_C2;
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

    const int iO2 = 2, iCO = 6, iOH = 1, iH2 = 3;

    // Soot mass rate [kg/m3/s] (negative = soot loses mass)
    r.soot_mass_rate = r.src_oxid[1] * r.V0_m3 * NAV_MOL_ * r.rho_soot;

    // Gas species rates
    r.o2_mass_rate  = r.omega_gas[static_cast<std::size_t>(iO2)];
    r.co_mass_rate  = r.omega_gas[static_cast<std::size_t>(iCO)];
    r.oh_mass_rate  = r.omega_gas[static_cast<std::size_t>(iOH)];
    r.h2_mass_rate  = r.omega_gas[static_cast<std::size_t>(iH2)];

    // Mass balance: soot_mass_rate + sum(omega_gas) = 0
    // For O2-ox: soot_mass_rate (neg) + o2 (neg) + co (pos) + oh (neg) + h2 (pos) = 0
    double gas_sum = 0.;
    for (double v : r.omega_gas) gas_sum += v;
    r.mass_residual = r.soot_mass_rate + gas_sum;

    // R_oxid_C2 from soot side [mol C2/m3/s]
    const double VC2 = (WC_KG_KMOL_ / r.rho_soot / NAV_KMOL_) * 2.0;
    r.R_oxid_C2 = -r.src_oxid[1] * r.V0_m3 / VC2;

    // fO2, fOH from gas side
    if (r.R_oxid_C2 > 0.) {
        // O2 consumed: -omega[O2] / MW_O2 * 1000 = fO2 * R_oxid_C2 [mol/m3/s]
        r.O2_per_C2 = -r.o2_mass_rate * 1000. / thermo.mw[iO2] / r.R_oxid_C2;
        r.CO_per_C2 = r.co_mass_rate * 1000. / thermo.mw[iCO] / r.R_oxid_C2;
        r.OH_per_C2 = -r.oh_mass_rate * 1000. / thermo.mw[iOH] / r.R_oxid_C2;
        r.H2_per_C2 = r.h2_mass_rate * 1000. / thermo.mw[iH2] / r.R_oxid_C2;
        r.fO2 = r.O2_per_C2;       // 1 O2 per C2 via O2 channel
        r.fOH = r.OH_per_C2 / 2.;  // 2 OH per C2 via OH channel
    } else {
        r.O2_per_C2 = 0.; r.CO_per_C2 = 0.; r.OH_per_C2 = 0.; r.H2_per_C2 = 0.;
        r.fO2 = 0.; r.fOH = 0.;
    }

    r.n_struct_fail = 0;
    auto check = [&](const char* name, bool ok, const std::string& detail) {
        std::cout << "  [STRUCT] " << (ok ? "PASS" : "FAIL") << "  " << name;
        if (!ok) { std::cout << "  -- " << detail; ++r.n_struct_fail; }
        std::cout << "\n";
    };

    std::cout << "\n--- " << label << ": structural assertions ---\n";

    if (ox_model == 0 && nuc_model == 0 && cond_model == 0 && sg_model == 0 && coag_model == 0) {
        for (int i = 0; i < 4; ++i) {
            std::string nm = "src_all[" + std::to_string(i) + "] == 0 (all-off)";
            check(nm.c_str(), r.src_all[i] == 0., "val=" + std::to_string(r.src_all[i]));
        }
        check("omega_gas all zero (all-off)",
              std::all_of(r.omega_gas.begin(), r.omega_gas.end(), [](double v){ return v == 0.; }),
              "nonzero entries present");
    } else {
        check("src_oxid[1] < 0 (oxidation active, soot mass loss)",
              r.src_oxid[1] < 0., "val=" + std::to_string(r.src_oxid[1]));
        check("src_nuc[1] == 0 (nucleation off)",
              r.src_nuc[1] == 0., "val=" + std::to_string(r.src_nuc[1]));
        check("src_cond[1] == 0 (condensation off)",
              r.src_cond[1] == 0., "val=" + std::to_string(r.src_cond[1]));
        check("src_growth[1] == 0 (growth off)",
              r.src_growth[1] == 0., "val=" + std::to_string(r.src_growth[1]));
        check("src_coag[1] == 0 (coagulation off)",
              r.src_coag[1] == 0., "val=" + std::to_string(r.src_coag[1]));
        check("R_oxid_C2 > 0 (positive C2 removal rate)",
              r.R_oxid_C2 > 0., "val=" + std::to_string(r.R_oxid_C2));
    }

    return r;
}

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
    jf << "      \"o2_mass_rate_kg_m3_s\": " << jnum(r.o2_mass_rate) << ",\n";
    jf << "      \"co_mass_rate_kg_m3_s\": " << jnum(r.co_mass_rate) << ",\n";
    jf << "      \"oh_mass_rate_kg_m3_s\": " << jnum(r.oh_mass_rate) << ",\n";
    jf << "      \"h2_mass_rate_kg_m3_s\": " << jnum(r.h2_mass_rate) << ",\n";
    jf << "      \"mass_residual_kg_m3_s\": " << jnum(r.mass_residual) << ",\n";
    jf << "      \"R_oxid_C2_mol_m3_s\": " << jnum(r.R_oxid_C2) << ",\n";
    jf << "      \"fO2\": " << jnum(r.fO2) << ",\n";
    jf << "      \"fOH\": " << jnum(r.fOH) << ",\n";
    jf << "      \"O2_per_C2\": " << jnum(r.O2_per_C2) << ",\n";
    jf << "      \"CO_per_C2\": " << jnum(r.CO_per_C2) << ",\n";
    jf << "      \"OH_per_C2\": " << jnum(r.OH_per_C2) << ",\n";
    jf << "      \"H2_per_C2\": " << jnum(r.H2_per_C2) << "\n";
    jf << "    },\n";
    jf << "    \"verification\": {\n";
    jf << "      \"n_struct_fail\": " << r.n_struct_fail << ",\n";
    jf << "      \"mass_closure_formula\": \"soot_mass_rate + sum(omega_gas) = 0\",\n";
    jf << "      \"mass_closure_residual\": " << jnum(r.mass_residual) << ",\n";
    jf << "      \"reaction_O2\": \"C2(soot) + O2 -> 2CO\",\n";
    jf << "      \"reaction_OH\": \"C2(soot) + 2OH -> 2CO + H2\"\n";
    jf << "    }\n";
    jf << "  }" << (is_last ? "" : ",") << "\n";
}

int main()
{
    std::cout << "=== HMOM O2-oxidation-only fixture (M1 vertical slice) ===\n";
    const MOM::BasicThermoData thermo = buildOxFixtureThermo();

    const double T_K = 1800., P_Pa = 101325., mu = 4.5e-5, moment_scale = 1.e5;

    // Species: H, OH, O2, H2, H2O, C2H2, CO, N2
    // Need significant O2 for O2-oxidation, and radicals (H, OH) above tanh threshold.
    // OH kept small to get fO2 ~ 1.0 (O2-dominant).

    // A. All-off control
    std::vector<double> X_A = {
        0.001, 0.001, 0.02, 0.005, 0.120, 0.020, 0.030, 0.803
    };
    // B. Above threshold: H=0.001, OH=1e-5 (small but Y_OH ~ 6.4e-6 >> 2e-8), O2=0.02
    std::vector<double> X_B = {
        0.001, 1e-5, 0.02, 0.005, 0.120, 0.020, 0.030, 0.804
    };
    // C. Below threshold: H and OH trace (damped)
    std::vector<double> X_C = {
        1e-10, 1e-10, 0.02, 0.005, 0.120, 0.020, 0.030, 0.805
    };
    // D. Near threshold: H ~ 5.3e-8, OH ~ 3.1e-8
    std::vector<double> X_D = {
        5.3e-8, 3.1e-8, 0.02, 0.005, 0.120, 0.020, 0.030, 0.805
    };

    auto fa = run_fixture("A_all_off_control", "All soot processes OFF.", thermo,
        0,0,0,0,0,0, T_K,P_Pa,mu, X_A, moment_scale);
    auto fb = run_fixture("B_ox_above_threshold", "OX=1, radicals above threshold, O2-dominant.", thermo,
        0,0,0,1,0,0, T_K,P_Pa,mu, X_B, moment_scale);
    auto fc = run_fixture("C_ox_below_threshold", "OX=1, radicals below threshold (damped).", thermo,
        0,0,0,1,0,0, T_K,P_Pa,mu, X_C, moment_scale);
    auto fd = run_fixture("D_ox_near_threshold", "OX=1, radicals near threshold (partial).", thermo,
        0,0,0,1,0,0, T_K,P_Pa,mu, X_D, moment_scale);

    int n_cross = 0;
    auto cross = [&](const char* name, bool ok, const std::string& detail) {
        std::cout << "  [CROSS] " << (ok ? "PASS" : "FAIL") << "  " << name;
        if (!ok) { std::cout << "  -- " << detail; ++n_cross; }
        std::cout << "\n";
    };
    std::cout << "\n--- Cross-fixture assertions ---\n";
    cross("B.src_oxid[1] < 0 (above threshold, active oxidation)",
          fb.src_oxid[1] < 0., "val=" + std::to_string(fb.src_oxid[1]));
    cross("|B.src_oxid[1]| > |C.src_oxid[1]| (above > below)",
          std::abs(fb.src_oxid[1]) > std::abs(fc.src_oxid[1]),
          "B=" + std::to_string(fb.src_oxid[1]) + " C=" + std::to_string(fc.src_oxid[1]));
    cross("|B.src_oxid[1]| > |D.src_oxid[1]| (above > near)",
          std::abs(fb.src_oxid[1]) > std::abs(fd.src_oxid[1]),
          "B=" + std::to_string(fb.src_oxid[1]) + " D=" + std::to_string(fd.src_oxid[1]));
    cross("B.fO2 > 0.5 (O2-dominant channel)",
          fb.fO2 > 0.5, "fO2=" + std::to_string(fb.fO2));

    const int total_fail = fa.n_struct_fail + fb.n_struct_fail + fc.n_struct_fail + fd.n_struct_fail + n_cross;
    std::cout << "\n--- Summary ---\n";
    std::cout << "  Total: " << total_fail << " failures\n";

    std::cout << std::scientific << std::setprecision(6);
    for (const auto* r : {&fa, &fb, &fc, &fd}) {
        std::cout << "\n  " << r->label << ":\n";
        std::cout << "    src_oxid[1]  = " << r->src_oxid[1] << "\n";
        std::cout << "    R_oxid_C2    = " << r->R_oxid_C2 << "\n";
        std::cout << "    fO2          = " << r->fO2 << "\n";
        std::cout << "    fOH          = " << r->fOH << "\n";
        std::cout << "    O2_per_C2    = " << r->O2_per_C2 << "\n";
        std::cout << "    CO_per_C2    = " << r->CO_per_C2 << "\n";
        std::cout << "    omega[O2]    = " << r->o2_mass_rate << "\n";
        std::cout << "    omega[CO]    = " << r->co_mass_rate << "\n";
        std::cout << "    soot_rate    = " << r->soot_mass_rate << "\n";
        std::cout << "    mass_residual= " << r->mass_residual << "\n";
    }

    const std::string out_path = "o2_oxidation_fixture_results.json";
    std::ofstream jf(out_path, std::ios::out | std::ios::trunc);
    if (!jf) { std::cerr << "ERROR: could not open " << out_path << "\n"; return 1; }
    jf << "{\n";
    jf << "  \"task\": \"M1-o2-oxidation\",\n";
    jf << "  \"description\": \"O2-oxidation-only single-cell fixtures (4 variants)\",\n";
    jf << "  \"species_names\": " << names_json(thermo) << ",\n";
    jf << "  \"fixtures\": {\n";
    write_fixture_json(jf, fa, thermo, false);
    write_fixture_json(jf, fb, thermo, false);
    write_fixture_json(jf, fc, thermo, false);
    write_fixture_json(jf, fd, thermo, true);
    jf << "  },\n";
    jf << "  \"closure\": {\n";
    jf << "    \"mass_formula\": \"soot_mass_rate + sum(omega_gas) = 0\",\n";
    jf << "    \"reaction_O2\": \"C2(soot) + O2 -> 2CO\",\n";
    jf << "    \"reaction_OH\": \"C2(soot) + 2OH -> 2CO + H2\",\n";
    jf << "    \"tanh_thresholds\": {\"YH\": 2e-9, \"YOH\": 2e-8}\n";
    jf << "  },\n";
    jf << "  \"total_struct_failures\": " << total_fail << "\n";
    jf << "}\n";
    jf.close();
    if (!jf) { std::cerr << "ERROR: failed writing " << out_path << "\n"; return 1; }
    std::cout << "\nWrote JSON results: " << out_path << "\n";

    if (total_fail > 0) { std::cerr << "ERROR: " << total_fail << " assertion(s) failed.\n"; return 1; }
    std::cout << "=== O2 oxidation fixture completed OK ===\n";
    return 0;
}
