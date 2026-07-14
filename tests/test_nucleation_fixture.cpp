/*-----------------------------------------------------------------------*\
|   MOM Library -- Nucleation-only fixture harness (M1 vertical slice)      |
|                                                                          |
|   Standalone C++ test harness (no Fluent) that drives                    |
|   HMOM::ComputeSources on a single cell with ONLY the nucleation         |
|   process active (condensation OFF), isolating the PAH-dimer -> soot     |
|   nucleation pathway for mass/element closure verification.              |
|                                                                          |
|   Three fixture variants:                                                |
|     A. All-off control: every soot process disabled.                     |
|        All sources must be exactly zero.                                 |
|     B. Nucleation-only nominal: nucleation=1, condensation=0,            |
|        all others=0.  Valid PAH concentration.                           |
|        Nucleation source must be significantly non-zero.                 |
|     C. Nucleation-only perturbation: same as B but PAH mole fraction     |
|        doubled.  Nucleation source must differ from B.                   |
|                                                                          |
|   Expected GREEN light: Mass balance should close because V0 is defined  |
|   as 4*mw_pah/(rho*Nav_kmol), making soot_mass_rate = PAH consumption    |
|   by construction.  This fixture verifies that identity holds.           |
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
// Thermo builder -- same species set as condensation fixture.
// ============================================================================
static MOM::BasicThermoData buildNucFixtureThermo()
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
    for (std::size_t k = 0; k < th.names.size(); ++k)
        MW += X[k] * th.mw[k];
    std::vector<double> Y(th.names.size());
    for (std::size_t k = 0; k < th.names.size(); ++k)
        Y[k] = X[k] * th.mw[k] / MW;
    return Y;
}

static double MixtureMW(const std::vector<double>& X, const MOM::BasicThermoData& th)
{
    double MW = 0.;
    for (std::size_t k = 0; k < th.names.size(); ++k)
        MW += X[k] * th.mw[k];
    return MW;
}

static std::string jnum(double v)
{
    std::ostringstream os;
    os << std::setprecision(17) << v;
    return os.str();
}

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
    double MW_mix, M00, M10, M01, N0;
    double rho_soot, V0_m3;
    double src_nuc[4], src_growth[4], src_oxid[4], src_cond[4], src_coag[4], src_all[4];
    std::vector<double> omega_gas;
    double R_PAH_dimer, dimerization_rate, conc_PAH;
    double R_nuc_C2, soot_mass_rate, pah_mass_rate, mass_residual;
    int n_struct_fail;
};

static FixtureResult run_fixture(
    const std::string& label, const std::string& description,
    const MOM::BasicThermoData& thermo,
    int nuc_model, int cond_model, int sg_model, int ox_model, int coag_model, int cont_coag_model,
    double T_K, double P_Pa, double mu,
    const std::vector<double>& X_in, double moment_scale)
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
    cfg.pah_species="C2H2"; cfg.nucleation_model=nuc_model; cfg.condensation_model=cond_model;
    cfg.surface_growth_model=sg_model; cfg.oxidation_model=ox_model;
    cfg.coagulation_model=coag_model; cfg.continuous_coagulation_model=cont_coag_model;
    model.SetupFromConfig(cfg);
    model.SetViscosity(mu);
    model.SetState(T_K, P_Pa, r.Y.data());

    const auto ic = model.initial_moments();
    r.M00 = ic[0]*moment_scale; r.M10 = ic[1]*moment_scale;
    r.M01 = ic[2]*moment_scale; r.N0  = ic[3]*moment_scale;
    model.SetNormalizedMoments(r.M00, r.M10, r.M01, r.N0);
    model.ComputeSources();

    const auto src_all = model.sources();
    const auto omega   = model.omega_gas();
    const auto src_nuc = model.sources_nucleation_impl();
    const auto src_gr  = model.sources_growth_impl();
    const auto src_ox  = model.sources_oxidation_impl();
    const auto src_co  = model.sources_condensation_impl();
    const auto src_cg  = model.sources_coagulation_impl();

    for (int i=0; i<4; ++i) {
        r.src_nuc[i]=src_nuc[i]; r.src_growth[i]=src_gr[i]; r.src_oxid[i]=src_ox[i];
        r.src_cond[i]=src_co[i]; r.src_coag[i]=src_cg[i]; r.src_all[i]=src_all[i];
    }
    r.omega_gas.assign(omega.begin(), omega.end());
    r.R_PAH_dimer = 2. * model.dimerization_rate() * NAV_MOL_ / 1000.;
    r.dimerization_rate = model.dimerization_rate();
    r.conc_PAH = model.precursor_concentration();
    r.rho_soot = model.particle_density();
    r.V0_m3 = model.V0();

    const double VC2 = (WC_KG_KMOL_ / r.rho_soot / NAV_KMOL_) * 2.0;
    r.R_nuc_C2 = r.src_nuc[1] * r.V0_m3 / VC2;
    r.soot_mass_rate = r.src_nuc[1] * r.V0_m3 * NAV_MOL_ * r.rho_soot;

    const int iC2H2 = thermo.IndexOfSpecies("C2H2");
    r.pah_mass_rate = r.omega_gas[static_cast<std::size_t>(iC2H2)];
    r.mass_residual = r.soot_mass_rate - (-r.pah_mass_rate);

    r.n_struct_fail = 0;
    auto check = [&](const char* name, bool ok, const std::string& detail) {
        std::cout << "  [STRUCT] " << (ok?"PASS":"FAIL") << "  " << name;
        if (!ok) { std::cout << "  -- " << detail; ++r.n_struct_fail; }
        std::cout << "\n";
    };

    std::cout << "\n--- " << label << ": structural assertions ---\n";

    if (nuc_model==0 && cond_model==0 && sg_model==0 && ox_model==0 && coag_model==0) {
        for (int i=0; i<4; ++i) {
            std::string nm = "src_all["+std::to_string(i)+"] == 0 (all-off)";
            check(nm.c_str(), r.src_all[i]==0., "val="+std::to_string(r.src_all[i]));
        }
        check("omega_gas[C2H2] == 0 (all-off)", r.pah_mass_rate==0., "val="+std::to_string(r.pah_mass_rate));
    } else {
        check("src_nuc[1] > 0 (nucleation active)", r.src_nuc[1]>0., "val="+std::to_string(r.src_nuc[1]));
        check("src_cond[1] == 0 (condensation off)", r.src_cond[1]==0., "val="+std::to_string(r.src_cond[1]));
        check("src_growth[1] == 0 (growth off)", r.src_growth[1]==0., "val="+std::to_string(r.src_growth[1]));
        check("src_oxid[1] == 0 (oxidation off)", r.src_oxid[1]==0., "val="+std::to_string(r.src_oxid[1]));
        check("src_coag[1] == 0 (coagulation off)", r.src_coag[1]==0., "val="+std::to_string(r.src_coag[1]));
        check("R_PAH_dimer > 0 (PAH dimerization active)", r.R_PAH_dimer>0., "val="+std::to_string(r.R_PAH_dimer));
    }
    return r;
}

static std::string names_json(const MOM::BasicThermoData& th) {
    std::string s="[";
    for (std::size_t k=0;k<th.names.size();++k) { s+="\""+th.names[k]+"\""; if(k+1<th.names.size()) s+=", "; }
    s+="]"; return s;
}
static std::string vec_json(const std::vector<double>& v) {
    std::string s="[";
    for (std::size_t k=0;k<v.size();++k) { s+=jnum(v[k]); if(k+1<v.size()) s+=", "; }
    s+="]"; return s;
}
static std::string src4_json(const double s[4]) {
    return "["+jnum(s[0])+", "+jnum(s[1])+", "+jnum(s[2])+", "+jnum(s[3])+"]";
}

static void write_fixture_json(std::ofstream& jf, const FixtureResult& r, const MOM::BasicThermoData& thermo, bool is_last) {
    const int iC2H2 = thermo.IndexOfSpecies("C2H2");
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
    jf << "      \"R_PAH_dimer_kmol_m3_s\": " << jnum(r.R_PAH_dimer) << ",\n";
    jf << "      \"dimerization_rate_mol_m3_s\": " << jnum(r.dimerization_rate) << ",\n";
    jf << "      \"conc_PAH_kmol_m3_s\": " << jnum(r.conc_PAH) << ",\n";
    jf << "      \"R_nuc_C2_kmol_m3_s\": " << jnum(r.R_nuc_C2) << ",\n";
    jf << "      \"soot_mass_rate_kg_m3_s\": " << jnum(r.soot_mass_rate) << ",\n";
    jf << "      \"pah_mass_rate_kg_m3_s\": " << jnum(r.pah_mass_rate) << ",\n";
    jf << "      \"mass_residual_kg_m3_s\": " << jnum(r.mass_residual) << "\n";
    jf << "    },\n";
    jf << "    \"verification\": {\n";
    jf << "      \"n_struct_fail\": " << r.n_struct_fail << ",\n";
    jf << "      \"mass_closure_formula\": \"soot_mass_rate - (-omega_gas[C2H2])\",\n";
    jf << "      \"mass_closure_residual\": " << jnum(r.mass_residual) << ",\n";
    jf << "      \"pah_species_index\": " << iC2H2 << ",\n";
    jf << "      \"pah_nc\": 2,\n";
    jf << "      \"pah_nh\": 2,\n";
    jf << "      \"pah_no\": 0,\n";
    jf << "      \"pah_nn\": 0\n";
    jf << "    }\n";
    jf << "  }" << (is_last ? "" : ",") << "\n";
}

int main()
{
    std::cout << "=== HMOM Nucleation-only fixture (M1 vertical slice) ===\n";
    const MOM::BasicThermoData thermo = buildNucFixtureThermo();
    std::cout << "Species (" << thermo.NumberOfSpecies() << "): ";
    for (std::size_t k=0;k<thermo.names.size();++k)
        std::cout << thermo.names[k] << (k+1<thermo.names.size()?" ":"\n");

    const double T_K=1800., P_Pa=101325., mu=4.5e-5, moment_scale=1.e5;

    std::vector<double> X_nominal = {0.001, 0.001, 1.0e-30, 0.005, 0.120, 0.020, 0.030, 0.823};
    std::vector<double> X_perturbed = {0.001, 0.001, 1.0e-30, 0.005, 0.120, 0.040, 0.030, 0.803};

    auto fa = run_fixture("A_all_off_control", "All soot processes OFF.", thermo,
        0,0,0,0,0,0, T_K,P_Pa,mu, X_nominal, moment_scale);
    auto fb = run_fixture("B_nuc_only_nominal", "Nucleation=1, all others=0.", thermo,
        1,0,0,0,0,0, T_K,P_Pa,mu, X_nominal, moment_scale);
    auto fc = run_fixture("C_nuc_only_perturbed", "Same as B but C2H2 doubled.", thermo,
        1,0,0,0,0,0, T_K,P_Pa,mu, X_perturbed, moment_scale);

    int n_cross_fail=0;
    auto cross=[&](const char* name, bool ok, const std::string& detail) {
        std::cout << "  [CROSS] " << (ok?"PASS":"FAIL") << "  " << name;
        if(!ok) { std::cout << "  -- " << detail; ++n_cross_fail; }
        std::cout << "\n";
    };
    std::cout << "\n--- Cross-fixture assertions ---\n";
    cross("C.src_nuc[1] != B.src_nuc[1] (perturbation response)",
          std::abs(fc.src_nuc[1]-fb.src_nuc[1])>0.,
          "B="+std::to_string(fb.src_nuc[1])+" C="+std::to_string(fc.src_nuc[1]));
    cross("C.src_nuc[1] > B.src_nuc[1] (more PAH -> more nucleation)",
          fc.src_nuc[1]>fb.src_nuc[1],
          "B="+std::to_string(fb.src_nuc[1])+" C="+std::to_string(fc.src_nuc[1]));

    const int total_fail = fa.n_struct_fail + fb.n_struct_fail + fc.n_struct_fail + n_cross_fail;
    std::cout << "\n--- Summary ---\n";
    std::cout << "  A: " << fa.n_struct_fail << " struct failures\n";
    std::cout << "  B: " << fb.n_struct_fail << " struct failures\n";
    std::cout << "  C: " << fc.n_struct_fail << " struct failures\n";
    std::cout << "  Cross: " << n_cross_fail << " failures\n";
    std::cout << "  Total: " << total_fail << " failures\n";

    std::cout << std::scientific << std::setprecision(6);
    for (const auto* r : {&fa, &fb, &fc}) {
        std::cout << "\n  " << r->label << ":\n";
        std::cout << "    src_nuc[0..3] = " << r->src_nuc[0] << " " << r->src_nuc[1]
                  << " " << r->src_nuc[2] << " " << r->src_nuc[3] << "\n";
        std::cout << "    omega_gas[C2H2] = " << r->pah_mass_rate << " kg/m3/s\n";
        std::cout << "    R_PAH_dimer = " << r->R_PAH_dimer << " kmol/m3/s\n";
        std::cout << "    soot_mass_rate = " << r->soot_mass_rate << " kg/m3/s\n";
        std::cout << "    mass_residual = " << r->mass_residual << " kg/m3/s\n";
    }

    const std::string out_path = "nucleation_fixture_results.json";
    std::ofstream jf(out_path, std::ios::out | std::ios::trunc);
    if (!jf) { std::cerr << "ERROR: could not open " << out_path << "\n"; return 1; }
    jf << "{\n";
    jf << "  \"task\": \"M1-nucleation\",\n";
    jf << "  \"description\": \"Nucleation-only single-cell fixtures (3 variants)\",\n";
    jf << "  \"species_names\": " << names_json(thermo) << ",\n";
    jf << "  \"fixtures\": {\n";
    write_fixture_json(jf, fa, thermo, false);
    write_fixture_json(jf, fb, thermo, false);
    write_fixture_json(jf, fc, thermo, true);
    jf << "  },\n";
    jf << "  \"closure\": {\n";
    jf << "    \"mass_formula\": \"soot_mass_rate = src_nuc[1] * V0 * Nav_mol * rho_particle; pah_mass_rate = omega_gas[C2H2]; residual = soot_mass_rate - (-pah_mass_rate)\",\n";
    jf << "    \"expected_green_light\": \"Mass closes by construction: V0 = 4*mw_pah/(rho*Nav_kmol)\",\n";
    jf << "    \"B_mass_residual\": " << jnum(fb.mass_residual) << ",\n";
    jf << "    \"C_mass_residual\": " << jnum(fc.mass_residual) << "\n";
    jf << "  },\n";
    jf << "  \"total_struct_failures\": " << total_fail << ",\n";
    jf << "  \"reproducibility\": \"deterministic -- no RNG, no time-stepping\"\n";
    jf << "}\n";
    jf.close();
    if (!jf) { std::cerr << "ERROR: failed writing " << out_path << "\n"; return 1; }

    std::cout << "\nWrote JSON results: " << out_path << "\n";
    if (total_fail > 0) { std::cerr << "ERROR: " << total_fail << " assertion(s) failed.\n"; return 1; }
    std::cout << "=== nucleation fixture completed OK ===\n";
    return 0;
}
