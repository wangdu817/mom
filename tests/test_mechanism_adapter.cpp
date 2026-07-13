/*-----------------------------------------------------------------------*\
|   MOM Library -- MechanismAdapter unit tests (task 3.4)                  |
|                                                                          |
|   Standalone C++ test (plain main(), manual checks with std::exit(1) on  |
|   failure, return 0 on success) exercising MOM::MechanismAdapter's three  |
|   pure static conversions build_thermo / build_config / build_descriptor  |
|   plus the spec §5.1 guard rails, and verifying HMOM accepts the built    |
|   Thermo + Config.                                                        |
|                                                                          |
|   Species order (0-based): H OH O2 H2 H2O C2H2 CO A4 N2 (N2 last, bulk). |
\*-----------------------------------------------------------------------*/

#include "MOM/MechanismAdapter.hpp"

#include "HMOM/HMOM.hpp"
#include "MOM/ManifestReader.hpp"
#include "MOM/ThermoProxy.hpp"
#include "test_manifest_fixture.hpp"

#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <string>

using mom_test::make_test_manifest;

namespace
{

int g_checks = 0;

void check(bool cond, const std::string& what)
{
    ++g_checks;
    if (!cond)
    {
        std::cerr << "FAIL: " << what << "\n";
        std::exit(1);
    }
}

bool almost_equal(double a, double b, double tol = 1e-9)
{
    return std::fabs(a - b) <= tol * (1.0 + std::fabs(a) + std::fabs(b));
}

} // namespace

int main()
{
    std::cout << "=== MechanismAdapter unit tests (task 3.4) ===\n";

    const MOM::Manifest manifest = make_test_manifest();

    // -- 1. build_thermo: count, names, MW, element counts (double), nti=0. --
    {
        const MOM::BasicThermoData thermo = MOM::MechanismAdapter::build_thermo(manifest);

        check(thermo.NumberOfSpecies() == 9u, "thermo: species count");
        check(thermo.names.size() == 9, "thermo: names size");
        check(thermo.names[0] == "H", "thermo: name[0] H");
        check(thermo.names[5] == "C2H2", "thermo: name[5] C2H2");
        check(thermo.names[7] == "A4", "thermo: name[7] A4");
        check(thermo.names[8] == "N2", "thermo: name[8] N2");

        // 4. MW units: BasicThermoData.mw is kg/kmol, direct copy from manifest.
        check(almost_equal(thermo.MolecularWeight(5), 26.038), "thermo: C2H2 mw kg/kmol");
        check(almost_equal(thermo.MolecularWeight(7), 202.256), "thermo: A4 mw kg/kmol");
        check(almost_equal(thermo.MolecularWeight(8), 28.014), "thermo: N2 mw kg/kmol");

        // Element counts widened to double.
        check(almost_equal(thermo.NumberOfCarbonAtoms(7), 16.0), "thermo: A4 nc");
        check(almost_equal(thermo.NumberOfHydrogenAtoms(7), 10.0), "thermo: A4 nh");
        check(almost_equal(thermo.NumberOfCarbonAtoms(5), 2.0), "thermo: C2H2 nc");
        check(almost_equal(thermo.NumberOfOxygenAtoms(2), 2.0), "thermo: O2 no");
        check(almost_equal(thermo.NumberOfNitrogenAtoms(8), 2.0), "thermo: N2 nn");

        // nti = 0 for every species (no titanium in this mechanism).
        for (unsigned i = 0; i < thermo.NumberOfSpecies(); ++i)
            check(almost_equal(thermo.NumberOfTitaniumAtoms(i), 0.0), "thermo: nti==0");
        std::cout << "  [ok] 1/4. build_thermo\n";
    }

    // -- 2. build_config: pah_species, thermophoretic, sticking, gas_consump. -
    {
        const MOM::HMOM<MOM::BasicThermoData>::Config cfg =
            MOM::MechanismAdapter::build_config(manifest);

        check(cfg.pah_species == "A4", "config: pah_species == A4");
        check(cfg.thermophoretic_model == 0, "config: thermophoretic_model == 0");
        check(cfg.sticking_model == "constant", "config: sticking_model == constant");
        check(cfg.gas_consumption == true, "config: gas_consumption == true");
        check(cfg.nucleation_model == 1, "config: nucleation_model");
        check(cfg.oxidation_model == 1, "config: oxidation_model");
        check(cfg.fractal_diameter_model == 1, "config: fractal_diameter_model");
        check(cfg.collision_diameter_model == 2, "config: collision_diameter_model");
        check(almost_equal(cfg.soot_density_kg_m3, 1800.0), "config: soot_density");
        check(almost_equal(cfg.sticking_coeff_constant, 0.002), "config: sticking_coeff");
        check(almost_equal(cfg.schmidt_number, 50.0), "config: schmidt_number");
        std::cout << "  [ok] 2. build_config\n";
    }

    // -- 3. build_descriptor: indices, operating_pressure, h_f0 values. ------
    {
        const MOM::MechanismDescriptor desc =
            MOM::MechanismAdapter::build_descriptor(manifest);

        check(desc.pah_index == 7, "descriptor: pah_index");
        check(desc.bulk_index == 8, "descriptor: bulk_index");
        check(desc.c2h2_index == 5, "descriptor: c2h2_index");
        check(desc.h2_index == 3, "descriptor: h2_index");
        check(desc.o2_index == 2, "descriptor: o2_index");
        check(desc.oh_index == 1, "descriptor: oh_index");
        check(desc.co_index == 6, "descriptor: co_index");
        check(almost_equal(desc.operating_pressure_pa, 101325.0), "descriptor: operating_pressure");

        check(desc.h_f0_j_mol.size() == 9, "descriptor: h_f0_j_mol size");
        check(desc.h_f0_j_kg.size() == 9, "descriptor: h_f0_j_kg size");
        check(almost_equal(desc.h_f0_j_mol[6], -110529.0), "descriptor: CO h_f0_j_mol");
        check(almost_equal(desc.h_f0_j_mol[0], 217997.0), "descriptor: H h_f0_j_mol");
        check(almost_equal(desc.h_f0_j_kg[6], -110529.0 * 1000.0 / 28.010, 1e-6),
              "descriptor: CO h_f0_j_kg");

        // NASA-7 / temperature bounds carried through in mechanism order.
        check(desc.nasa7_low.size() == 9, "descriptor: nasa7_low size");
        check(almost_equal(desc.nasa7_low[0][5], -1000.0), "descriptor: nasa7_low coeff");
        check(almost_equal(desc.temp_common[0], 1000.0), "descriptor: temp_common");
        std::cout << "  [ok] 3. build_descriptor\n";
    }

    // -- 5. h_f0 sign: CO/H2O negative, H2/O2/N2 zero. -----------------------
    {
        const MOM::MechanismDescriptor desc =
            MOM::MechanismAdapter::build_descriptor(manifest);
        check(desc.h_f0_j_mol[6] < 0.0, "sign: CO h_f0 negative");   // CO
        check(desc.h_f0_j_mol[4] < 0.0, "sign: H2O h_f0 negative");  // H2O
        check(almost_equal(desc.h_f0_j_mol[3], 0.0), "sign: H2 h_f0 zero");   // H2
        check(almost_equal(desc.h_f0_j_mol[2], 0.0), "sign: O2 h_f0 zero");   // O2
        check(almost_equal(desc.h_f0_j_mol[8], 0.0), "sign: N2 h_f0 zero");   // N2
        std::cout << "  [ok] 5. h_f0 signs\n";
    }

    // -- 6. surface_density_correction: 0.0 -> false, 1.0 -> true. -----------
    {
        MOM::Manifest m0 = make_test_manifest();
        m0.surface_density_correction = 0.0;
        const auto cfg0 = MOM::MechanismAdapter::build_config(m0);
        check(cfg0.surface_density_correction == false, "sdc: 0.0 -> false");

        MOM::Manifest m1 = make_test_manifest();
        m1.surface_density_correction = 1.0;
        const auto cfg1 = MOM::MechanismAdapter::build_config(m1);
        check(cfg1.surface_density_correction == true, "sdc: 1.0 -> true");
        std::cout << "  [ok] 6. surface_density_correction bool conversion\n";
    }

    // -- 7. Guard rail: sticking_model == "golaut" throws. -------------------
    {
        MOM::Manifest m = make_test_manifest();
        m.sticking_model = "golaut";
        bool threw = false;
        try
        {
            MOM::MechanismAdapter::build_config(m);
        }
        catch (const MOM::ManifestParseError&)
        {
            threw = true;
        }
        check(threw, "guard: golaut sticking throws");
        std::cout << "  [ok] 7. golaut sticking guard rail\n";
    }

    // -- 8. Guard rail: thermophoretic_model != 0 throws. --------------------
    {
        MOM::Manifest m = make_test_manifest();
        m.thermophoretic_model = 1;
        bool threw = false;
        try
        {
            MOM::MechanismAdapter::build_config(m);
        }
        catch (const MOM::ManifestParseError&)
        {
            threw = true;
        }
        check(threw, "guard: thermophoretic!=0 throws");
        std::cout << "  [ok] 8. thermophoretic guard rail\n";
    }

    // -- 9. HMOM accepts the constructed Thermo + Config (no throw). ---------
    {
        const MOM::BasicThermoData thermo = MOM::MechanismAdapter::build_thermo(manifest);
        const MOM::HMOM<MOM::BasicThermoData>::Config cfg =
            MOM::MechanismAdapter::build_config(manifest);

        bool ok = true;
        try
        {
            MOM::HMOM<MOM::BasicThermoData> model(thermo);
            model.SetupFromConfig(cfg);
            check(model.omega_gas().size() == 9, "hmom: omega_gas sized to species");
            check(model.precursor_species() == "A4", "hmom: precursor is A4");
        }
        catch (const std::exception& e)
        {
            ok = false;
            std::cerr << "  unexpected exception: " << e.what() << "\n";
        }
        check(ok, "hmom: SetupFromConfig no throw");
        std::cout << "  [ok] 9. HMOM accepts thermo + config\n";
    }

    std::cout << "=== MechanismAdapter: all checks passed (" << g_checks << ") ===\n";
    return 0;
}
