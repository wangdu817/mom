/*-----------------------------------------------------------------------*\
|   MOM Library -- Coagulation-only fixture harness (M1 vertical slice)    |
|                                                                          |
|   Drives HMOM::ComputeSources with ONLY coagulation active.              |
|   Closure contract (review 5.4):                                         |
|     - omega_gas     = 0 for ALL species (no gas-phase source)            |
|     - S_M10         = 0 (soot mass conserved)                            |
|     - omega_soot    = rho*V0*Nav*S_M10 = 0                              |
|     - Q_HMOM        = 0 (no heat release)                                |
|     - S_M00 < 0     (number decreases as particles merge)                |
|     - S_M01, S_N0   may change (morphology / primary-particle count)     |
|                                                                          |
|   Four fixture variants:                                                 |
|     A. All-off control.                                                  |
|     B. Coagulation (discrete only), nominal soot loading.                |
|     C. Coagulation (discrete only), high soot loading (more collisions). |
|     D. Coagulation (discrete + continuous), nominal soot loading.        |
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

static MOM::BasicThermoData buildCoagFixtureThermo()
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
    double src_coag_ss[4], src_coag_sl[4], src_coag_ll[4];
    double src_coag_cont_ss[4], src_coag_cont_sl[4], src_coag_cont_ll[4];
    std::vector<double> omega_gas;
    double soot_mass_rate;
    double mass_residual;
    double gas_sum;
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

    const auto cd = model.coagulation_detail();
    for (int i = 0; i < 4; ++i) {
        r.src_coag_ss[i]      = cd.disc_ss[static_cast<std::size_t>(i)];
        r.src_coag_sl[i]      = cd.disc_sl[static_cast<std::size_t>(i)];
        r.src_coag_ll[i]      = cd.disc_ll[static_cast<std::size_t>(i)];
        r.src_coag_cont_ss[i] = cd.cont_ss[static_cast<std::size_t>(i)];
        r.src_coag_cont_sl[i] = cd.cont_sl[static_cast<std::size_t>(i)];
        r.src_coag_cont_ll[i] = cd.cont_ll[static_cast<std::size_t>(i)];
    }

    r.omega_gas.assign(omega.begin(), omega.end());
    r.rho_soot = model.particle_density();
    r.V0_m3    = model.V0();

    r.soot_mass_rate = r.src_coag[1] * r.V0_m3 * NAV_MOL_ * r.rho_soot;

    r.gas_sum = 0.;
    for (double v : r.omega_gas) r.gas_sum += v;
    r.mass_residual = r.soot_mass_rate + r.gas_sum;

    r.n_struct_fail = 0;
    auto check = [&](const char* name, bool ok, const std::string& detail) {
        std::cout << "  [STRUCT] " << (ok ? "PASS" : "FAIL") << "  " << name;
        if (!ok) { std::cout << "  -- " << detail; ++r.n_struct_fail; }
        std::cout << "\n";
    };

    std::cout << "\n--- " << label << ": structural assertions ---\n";

    if (coag_model == 0 && cont_coag_model == 0 && nuc_model == 0 && cond_model == 0 && sg_model == 0 && ox_model == 0) {
        for (int i = 0; i < 4; ++i) {
            check(("src_all[" + std::to_string(i) + "] == 0").c_str(),
                  r.src_all[i] == 0., jnum(r.src_all[i]));
        }
        for (std::size_t k = 0; k < r.omega_gas.size(); ++k) {
            check(("omega_gas[" + std::to_string(k) + "] == 0").c_str(),
                  r.omega_gas[k] == 0., jnum(r.omega_gas[k]));
        }
    } else {
        check("S_M10 (coag mass source) == 0",
              std::abs(r.src_coag[1]) < 1e-30,
              "S_M10=" + jnum(r.src_coag[1]));

        check("S_M00 (number source) < 0",
              r.src_coag[0] < 0.,
              "S_M00=" + jnum(r.src_coag[0]));

        check("omega_gas all zero",
              std::abs(r.gas_sum) < 1e-30,
              "gas_sum=" + jnum(r.gas_sum));

        check("soot_mass_rate == 0",
              std::abs(r.soot_mass_rate) < 1e-30,
              "soot_mass_rate=" + jnum(r.soot_mass_rate));

        check("mass_residual == 0",
              std::abs(r.mass_residual) < 1e-30,
              "mass_residual=" + jnum(r.mass_residual));

        check("S_nuc all zero",
              r.src_nuc[0] == 0. && r.src_nuc[1] == 0. && r.src_nuc[2] == 0. && r.src_nuc[3] == 0.,
              "");

        check("S_growth all zero",
              r.src_growth[0] == 0. && r.src_growth[1] == 0. && r.src_growth[2] == 0. && r.src_growth[3] == 0.,
              "");

        check("S_oxid all zero",
              r.src_oxid[0] == 0. && r.src_oxid[1] == 0. && r.src_oxid[2] == 0. && r.src_oxid[3] == 0.,
              "");

        check("S_cond all zero",
              r.src_cond[0] == 0. && r.src_cond[1] == 0. && r.src_cond[2] == 0. && r.src_cond[3] == 0.,
              "");

        check("S_M10 sub-components all zero (ss/sl/ll/cont_ss/cont_sl/cont_ll)",
              std::abs(r.src_coag_ss[1]) < 1e-30 &&
              std::abs(r.src_coag_sl[1]) < 1e-30 &&
              std::abs(r.src_coag_ll[1]) < 1e-30 &&
              std::abs(r.src_coag_cont_ss[1]) < 1e-30 &&
              std::abs(r.src_coag_cont_sl[1]) < 1e-30 &&
              std::abs(r.src_coag_cont_ll[1]) < 1e-30,
              "");
    }

    return r;
}

static void write_json(const std::vector<FixtureResult>& results, const std::string& path)
{
    std::ofstream f(path);
    f << std::setprecision(17);
    f << "{\n  \"fixtures\": [\n";
    for (std::size_t i = 0; i < results.size(); ++i) {
        const auto& r = results[i];
        f << "    {\n";
        f << "      \"label\": \"" << r.label << "\",\n";
        f << "      \"description\": \"" << r.description << "\",\n";
        f << "      \"process_switches\": {"
          << "\"nucleation\":" << r.nucleation_model
          << ",\"condensation\":" << r.condensation_model
          << ",\"surface_growth\":" << r.surface_growth_model
          << ",\"oxidation\":" << r.oxidation_model
          << ",\"coagulation\":" << r.coagulation_model
          << ",\"continuous_coagulation\":" << r.continuous_coagulation_model << "},\n";
        f << "      \"T_K\": " << jnum(r.T_K) << ",\n";
        f << "      \"P_Pa\": " << jnum(r.P_Pa) << ",\n";
        f << "      \"mu\": " << jnum(r.mu) << ",\n";
        f << "      \"MW_mix\": " << jnum(r.MW_mix) << ",\n";
        f << "      \"Y\": ["; for (std::size_t k=0;k<r.Y.size();++k){if(k)f<<", ";f<<jnum(r.Y[k]);} f << "],\n";
        f << "      \"moments\": {"
          << "\"M00\":" << jnum(r.M00) << ",\"M10\":" << jnum(r.M10)
          << ",\"M01\":" << jnum(r.M01) << ",\"N0\":" << jnum(r.N0) << "},\n";
        f << "      \"rho_soot\": " << jnum(r.rho_soot) << ",\n";
        f << "      \"V0_m3\": " << jnum(r.V0_m3) << ",\n";
        f << "      \"src_nuc\": ["; for(int k=0;k<4;++k){if(k)f<<", ";f<<jnum(r.src_nuc[k]);} f << "],\n";
        f << "      \"src_growth\": ["; for(int k=0;k<4;++k){if(k)f<<", ";f<<jnum(r.src_growth[k]);} f << "],\n";
        f << "      \"src_oxid\": ["; for(int k=0;k<4;++k){if(k)f<<", ";f<<jnum(r.src_oxid[k]);} f << "],\n";
        f << "      \"src_cond\": ["; for(int k=0;k<4;++k){if(k)f<<", ";f<<jnum(r.src_cond[k]);} f << "],\n";
        f << "      \"src_coag\": ["; for(int k=0;k<4;++k){if(k)f<<", ";f<<jnum(r.src_coag[k]);} f << "],\n";
        f << "      \"src_all\": ["; for(int k=0;k<4;++k){if(k)f<<", ";f<<jnum(r.src_all[k]);} f << "],\n";
        f << "      \"src_coag_ss\": ["; for(int k=0;k<4;++k){if(k)f<<", ";f<<jnum(r.src_coag_ss[k]);} f << "],\n";
        f << "      \"src_coag_sl\": ["; for(int k=0;k<4;++k){if(k)f<<", ";f<<jnum(r.src_coag_sl[k]);} f << "],\n";
        f << "      \"src_coag_ll\": ["; for(int k=0;k<4;++k){if(k)f<<", ";f<<jnum(r.src_coag_ll[k]);} f << "],\n";
        f << "      \"src_coag_cont_ss\": ["; for(int k=0;k<4;++k){if(k)f<<", ";f<<jnum(r.src_coag_cont_ss[k]);} f << "],\n";
        f << "      \"src_coag_cont_sl\": ["; for(int k=0;k<4;++k){if(k)f<<", ";f<<jnum(r.src_coag_cont_sl[k]);} f << "],\n";
        f << "      \"src_coag_cont_ll\": ["; for(int k=0;k<4;++k){if(k)f<<", ";f<<jnum(r.src_coag_cont_ll[k]);} f << "],\n";
        f << "      \"omega_gas\": ["; for(std::size_t k=0;k<r.omega_gas.size();++k){if(k)f<<", ";f<<jnum(r.omega_gas[k]);} f << "],\n";
        f << "      \"soot_mass_rate\": " << jnum(r.soot_mass_rate) << ",\n";
        f << "      \"gas_sum\": " << jnum(r.gas_sum) << ",\n";
        f << "      \"mass_residual\": " << jnum(r.mass_residual) << ",\n";
        f << "      \"n_struct_fail\": " << r.n_struct_fail << "\n";
        f << "    }";
        if (i + 1 < results.size()) f << ",";
        f << "\n";
    }
    f << "  ]\n}\n";
}

int main()
{
    const auto thermo = buildCoagFixtureThermo();

    constexpr double T  = 1800.;
    constexpr double P  = 101325.;
    constexpr double mu = 4.5e-5;

    std::vector<double> X = {0.01, 0.001, 0.02, 0.01, 0.05, 0.02, 0.01, 0.70};

    std::vector<FixtureResult> results;

    std::cout << "=== M1 Coagulation Fixture ===\n";

    results.push_back(run_fixture(
        "A_all_off", "All processes off (control)",
        thermo, 0,0,0,0,0,0, T,P,mu, X, 1e5));

    results.push_back(run_fixture(
        "B_coag_nominal", "Coagulation discrete only, nominal soot",
        thermo, 0,0,0,0,1,0, T,P,mu, X, 1e5));

    results.push_back(run_fixture(
        "C_coag_high_soot", "Coagulation discrete only, high soot loading",
        thermo, 0,0,0,0,1,0, T,P,mu, X, 1e7));

    results.push_back(run_fixture(
        "D_coag_continuous", "Coagulation discrete + continuous, nominal soot",
        thermo, 0,0,0,0,1,1, T,P,mu, X, 1e5));

    std::cout << "\n=== Cross-fixture assertions ===\n";
    int n_cross_fail = 0;
    auto xcheck = [&](const char* name, bool ok, const std::string& detail) {
        std::cout << "  [CROSS] " << (ok ? "PASS" : "FAIL") << "  " << name;
        if (!ok) { std::cout << "  -- " << detail; ++n_cross_fail; }
        std::cout << "\n";
    };

    const auto& A = results[0];
    const auto& B = results[1];
    const auto& C = results[2];
    const auto& D = results[3];

    xcheck("A: all sources zero",
           A.src_all[0]==0. && A.src_all[1]==0. && A.src_all[2]==0. && A.src_all[3]==0., "");

    xcheck("B: S_M00 < 0 (number decreases)",
           B.src_coag[0] < 0., "S_M00=" + jnum(B.src_coag[0]));

    xcheck("B: S_M10 == 0 (mass conserved)",
           std::abs(B.src_coag[1]) < 1e-30, "S_M10=" + jnum(B.src_coag[1]));

    xcheck("B: gas_sum == 0 (no gas source)",
           std::abs(B.gas_sum) < 1e-30, "gas_sum=" + jnum(B.gas_sum));

    xcheck("C: |S_M00| > |B.S_M00| (more soot -> more coagulation)",
           std::abs(C.src_coag[0]) > std::abs(B.src_coag[0]),
           "|C|=" + jnum(std::abs(C.src_coag[0])) + " |B|=" + jnum(std::abs(B.src_coag[0])));

    xcheck("C: S_M10 == 0 (mass conserved)",
           std::abs(C.src_coag[1]) < 1e-30, "S_M10=" + jnum(C.src_coag[1]));

    xcheck("D: S_M10 == 0 (mass conserved, discrete+continuous)",
           std::abs(D.src_coag[1]) < 1e-30, "S_M10=" + jnum(D.src_coag[1]));

    xcheck("D: S_M00 < 0 (number decreases)",
           D.src_coag[0] < 0., "S_M00=" + jnum(D.src_coag[0]));

    xcheck("D: gas_sum == 0 (no gas source, discrete+continuous)",
           std::abs(D.gas_sum) < 1e-30, "gas_sum=" + jnum(D.gas_sum));

    xcheck("D: |S_M00| >= |B.S_M00| (continuous adds to discrete)",
           std::abs(D.src_coag[0]) >= std::abs(B.src_coag[0]) * 0.99,
           "|D|=" + jnum(std::abs(D.src_coag[0])) + " |B|=" + jnum(std::abs(B.src_coag[0])));

    int total_fail = 0;
    for (const auto& r : results) total_fail += r.n_struct_fail;
    total_fail += n_cross_fail;

    std::cout << "\n=== Summary ===\n";
    std::cout << "  Structural failures: " << total_fail << "\n";
    for (const auto& r : results) {
        std::cout << "  " << r.label
                  << "  S_M00=" << jnum(r.src_coag[0])
                  << "  S_M10=" << jnum(r.src_coag[1])
                  << "  S_M01=" << jnum(r.src_coag[2])
                  << "  S_N0="  << jnum(r.src_coag[3])
                  << "  gas_sum=" << jnum(r.gas_sum)
                  << "  mass_res=" << jnum(r.mass_residual)
                  << "  fails=" << r.n_struct_fail << "\n";
    }

    write_json(results, "coagulation_fixture_results.json");
    std::cout << "\n  JSON written: coagulation_fixture_results.json\n";

    if (total_fail > 0) {
        std::cerr << "FAILED with " << total_fail << " structural/cross errors\n";
        return 1;
    }
    std::cout << "\nALL PASS\n";
    return 0;
}
