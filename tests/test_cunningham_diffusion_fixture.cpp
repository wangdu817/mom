/*-----------------------------------------------------------------------*\
|   MOM Library -- Cunningham diffusion fixture harness (M3 vertical slice) |
|                                                                          |
|   Drives HMOM::diffusion_coefficient() with all soot processes off.      |
|   Closure contract:                                                      |
|     - Gamma = max(rho * D_Cunningham, mu / Sc)                           |
|     - D_Cunningham = kB*T*Cu / (3*pi*mu*dc_safe)                        |
|     - Cu = 1 + kCunningham * lambda / dc_safe                           |
|     - lambda = (mu/rho) * sqrt(pi * m_gas / (2*kB*T))                  |
|     - m_gas = rho * kB * T / P                                          |
|     - dc_safe = max(collision_diameter, 1e-12)                          |
|     - rho = P / (Rgas * T) * MW_mix                                    |
|                                                                          |
|   Four fixture variants:                                                 |
|     A. Flame (1800 K), small soot -- high Kn, strong Cunningham          |
|     B. Flame (1800 K), large soot -- low Kn, weak Cunningham             |
|     C. Room temp (300 K), medium soot -- different gas properties        |
|     D. Degenerate (zero moments) -- dc=0, diameter floor 1e-12           |
\*-----------------------------------------------------------------------*/

#include "HMOM/HMOM.hpp"

#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

static MOM::BasicThermoData buildCunninghamFixtureThermo()
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

constexpr double kB_         = 1.380649e-23;
constexpr double Rgas_       = 8314.46261815324;
constexpr double kCunningham_= 2.154;
constexpr double PI_         = 3.14159265358979323846;
constexpr double NAV_MOL_    = 6.02214076e23;
constexpr double WC_KG_KMOL_ = 12.011;

struct CunninghamResult
{
    std::string label, description;
    double T_K, P_Pa, mu;
    std::vector<double> X, Y;
    double MW_mix, moment_scale, mass_factor;
    double M00, M10, M01, N0;
    double rho_soot, V0_m3;
    double collision_diameter;
    double diffusion_coefficient;
    double schmidt_number;
    double particle_diameter;
    double volume_fraction;
    int n_struct_fail;
};

static CunninghamResult run_cunningham_fixture(
    const std::string& label, const std::string& description,
    const MOM::BasicThermoData& thermo,
    double T_K, double P_Pa, double mu,
    const std::vector<double>& X_in,
    double moment_scale,
    double mass_factor = 1.0)
{
    CunninghamResult r;
    r.label = label; r.description = description;
    r.T_K = T_K; r.P_Pa = P_Pa; r.mu = mu;
    r.moment_scale = moment_scale;
    r.mass_factor = mass_factor;

    r.X = X_in;
    { double s=0; for(double x:r.X) s+=x; for(double& x:r.X) x/=s; }
    r.Y = X2Y(r.X, thermo);
    r.MW_mix = MixtureMW(r.X, thermo);

    MOM::HMOM<MOM::BasicThermoData> model(thermo);
    MOM::HMOM<MOM::BasicThermoData>::Config cfg;
    cfg.pah_species = "C2H2";
    cfg.nucleation_model = 0; cfg.condensation_model = 0;
    cfg.surface_growth_model = 0; cfg.oxidation_model = 0;
    cfg.coagulation_model = 0; cfg.continuous_coagulation_model = 0;
    model.SetupFromConfig(cfg);
    model.SetViscosity(mu);
    model.SetState(T_K, P_Pa, r.Y.data());

    const auto ic = model.initial_moments();
    r.M00 = ic[0] * moment_scale;
    r.M10 = ic[1] * moment_scale * mass_factor;
    r.M01 = ic[2] * moment_scale * mass_factor;
    r.N0  = ic[3] * moment_scale;
    model.SetNormalizedMoments(r.M00, r.M10, r.M01, r.N0);

    r.collision_diameter    = model.collision_diameter();
    r.diffusion_coefficient = model.diffusion_coefficient();
    r.schmidt_number        = model.schmidt_number();
    r.particle_diameter     = model.particle_diameter();
    r.volume_fraction       = model.volume_fraction();
    r.rho_soot              = model.particle_density();
    r.V0_m3                 = model.V0();

    r.n_struct_fail = 0;
    auto check = [&](const char* name, bool ok, const std::string& detail) {
        std::cout << "  [STRUCT] " << (ok ? "PASS" : "FAIL") << "  " << name;
        if (!ok) { std::cout << "  -- " << detail; ++r.n_struct_fail; }
        std::cout << "\n";
    };

    std::cout << "\n--- " << label << ": structural assertions ---\n";

    check("diffusion_coefficient > 0",
          r.diffusion_coefficient > 0., jnum(r.diffusion_coefficient));
    check("schmidt_number > 0",
          r.schmidt_number > 0., jnum(r.schmidt_number));

    if (moment_scale > 0.) {
        check("collision_diameter > 0 (non-degenerate)",
              r.collision_diameter > 0., jnum(r.collision_diameter));
        check("particle_diameter > 0 (non-degenerate)",
              r.particle_diameter > 0., jnum(r.particle_diameter));
    } else {
        check("collision_diameter == 0 (degenerate, floor applies)",
              r.collision_diameter == 0., jnum(r.collision_diameter));
    }

    check("volume_fraction >= 0",
          r.volume_fraction >= 0., jnum(r.volume_fraction));

    return r;
}

static void write_json(const std::vector<CunninghamResult>& results, const std::string& path)
{
    std::ofstream f(path);
    f << std::setprecision(17);
    f << "{\n";
    f << "  \"constants\": {\n";
    f << "    \"kB\": " << jnum(kB_) << ",\n";
    f << "    \"Rgas\": " << jnum(Rgas_) << ",\n";
    f << "    \"kCunningham\": " << jnum(kCunningham_) << ",\n";
    f << "    \"pi\": " << jnum(PI_) << ",\n";
    f << "    \"Nav_mol\": " << jnum(NAV_MOL_) << ",\n";
    f << "    \"WC_kg_kmol\": " << jnum(WC_KG_KMOL_) << "\n";
    f << "  },\n";
    f << "  \"fixtures\": [\n";
    for (std::size_t i = 0; i < results.size(); ++i) {
        const auto& r = results[i];
        f << "    {\n";
        f << "      \"label\": \"" << r.label << "\",\n";
        f << "      \"description\": \"" << r.description << "\",\n";
        f << "      \"T_K\": " << jnum(r.T_K) << ",\n";
        f << "      \"P_Pa\": " << jnum(r.P_Pa) << ",\n";
        f << "      \"mu\": " << jnum(r.mu) << ",\n";
        f << "      \"MW_mix\": " << jnum(r.MW_mix) << ",\n";
        f << "      \"moment_scale\": " << jnum(r.moment_scale) << ",\n";
        f << "      \"mass_factor\": " << jnum(r.mass_factor) << ",\n";
        f << "      \"Y\": ["; for (std::size_t k=0;k<r.Y.size();++k){if(k)f<<", ";f<<jnum(r.Y[k]);} f << "],\n";
        f << "      \"moments\": {"
          << "\"M00\":" << jnum(r.M00) << ",\"M10\":" << jnum(r.M10)
          << ",\"M01\":" << jnum(r.M01) << ",\"N0\":" << jnum(r.N0) << "},\n";
        f << "      \"rho_soot\": " << jnum(r.rho_soot) << ",\n";
        f << "      \"V0_m3\": " << jnum(r.V0_m3) << ",\n";
        f << "      \"collision_diameter\": " << jnum(r.collision_diameter) << ",\n";
        f << "      \"diffusion_coefficient\": " << jnum(r.diffusion_coefficient) << ",\n";
        f << "      \"schmidt_number\": " << jnum(r.schmidt_number) << ",\n";
        f << "      \"particle_diameter\": " << jnum(r.particle_diameter) << ",\n";
        f << "      \"volume_fraction\": " << jnum(r.volume_fraction) << ",\n";
        f << "      \"n_struct_fail\": " << r.n_struct_fail << "\n";
        f << "    }";
        if (i + 1 < results.size()) f << ",";
        f << "\n";
    }
    f << "  ]\n}\n";
}

int main()
{
    const auto thermo = buildCunninghamFixtureThermo();

    std::vector<double> X = {0.01, 0.001, 0.02, 0.01, 0.05, 0.02, 0.01, 0.70};

    std::vector<CunninghamResult> results;

    std::cout << "=== M3 Cunningham Diffusion Fixture ===\n";

    results.push_back(run_cunningham_fixture(
        "A_flame_small_soot", "Flame 1800 K, small soot (high Kn, strong Cunningham)",
        thermo, 1800., 101325., 4.5e-5, X, 1e5, 1.0));

    results.push_back(run_cunningham_fixture(
        "B_flame_large_soot", "Flame 1800 K, large soot (low Kn, weak Cunningham)",
        thermo, 1800., 101325., 4.5e-5, X, 1e5, 1000.0));

    results.push_back(run_cunningham_fixture(
        "C_room_medium_soot", "Room 300 K, medium soot (different gas properties)",
        thermo, 300., 101325., 1.8e-5, X, 1e7));

    results.push_back(run_cunningham_fixture(
        "D_degenerate", "Flame 1800 K, zero moments (dc=0, diameter floor 1e-12)",
        thermo, 1800., 101325., 4.5e-5, X, 0.));

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

    xcheck("A: collision_diameter < B: collision_diameter (small soot < large soot)",
           A.collision_diameter < B.collision_diameter,
           "A=" + jnum(A.collision_diameter) + " B=" + jnum(B.collision_diameter));

    xcheck("A: diffusion > B: diffusion (small soot diffuses faster)",
           A.diffusion_coefficient > B.diffusion_coefficient,
           "A=" + jnum(A.diffusion_coefficient) + " B=" + jnum(B.diffusion_coefficient));

    xcheck("D: collision_diameter == 0 (degenerate)",
           D.collision_diameter == 0.,
           "D=" + jnum(D.collision_diameter));

    xcheck("D: diffusion > 0 (floor applies, Gamma computed with dc_safe=1e-12)",
           D.diffusion_coefficient > 0.,
           "D=" + jnum(D.diffusion_coefficient));

    xcheck("D: diffusion >> A: diffusion (tiny dc_safe -> huge D)",
           D.diffusion_coefficient > A.diffusion_coefficient * 1e2,
           "D=" + jnum(D.diffusion_coefficient) + " A=" + jnum(A.diffusion_coefficient));

    xcheck("C: diffusion > 0 (room temp still produces valid Gamma)",
           C.diffusion_coefficient > 0.,
           "C=" + jnum(C.diffusion_coefficient));

    int total_fail = 0;
    for (const auto& r : results) total_fail += r.n_struct_fail;
    total_fail += n_cross_fail;

    std::cout << "\n=== Summary ===\n";
    std::cout << "  Structural failures: " << total_fail << "\n";
    for (const auto& r : results) {
        std::cout << "  " << r.label
                  << "  dc=" << jnum(r.collision_diameter)
                  << "  Gamma=" << jnum(r.diffusion_coefficient)
                  << "  Sc=" << jnum(r.schmidt_number)
                  << "  dp=" << jnum(r.particle_diameter)
                  << "  fv=" << jnum(r.volume_fraction)
                  << "  fails=" << r.n_struct_fail << "\n";
    }

    write_json(results, "cunningham_diffusion_fixture_results.json");
    std::cout << "\n  JSON written: cunningham_diffusion_fixture_results.json\n";

    if (total_fail > 0) {
        std::cerr << "FAILED with " << total_fail << " structural/cross errors\n";
        return 1;
    }
    std::cout << "\nALL PASS\n";
    return 0;
}
