/*-----------------------------------------------------------------------*\
|   MOM Library -- shared test manifest fixture (task 3.4)                  |
|                                                                          |
|   Header-only helper shared by test_manifest_reader / test_mechanism_    |
|   adapter / test_mom_context.  Provides make_test_manifest(): a fully-    |
|   populated, VALID 9-species Manifest (H OH O2 H2 H2O C2H2 CO A4 N2,      |
|   N2 last / bulk).  Everything is inline; no separate translation unit.   |
\*-----------------------------------------------------------------------*/

#pragma once

#include "MOM/ManifestReader.hpp"

#include <array>
#include <string>

namespace mom_test
{

// Fully-populated, valid 9-species manifest used across the task-3.4 tests.
inline MOM::Manifest make_test_manifest()
{
    struct SpeciesSeed
    {
        const char* name;
        double mw_kg_kmol;
        int nc, nh, no, nn;
        double h_f0_j_mol;
    };

    const SpeciesSeed seeds[] = {
        {"H",    1.008,   0, 1, 0, 0, 217997.0},
        {"OH",   17.008,  0, 1, 1, 0, 37300.0},
        {"O2",   31.999,  0, 0, 2, 0, 0.0},
        {"H2",   2.016,   0, 2, 0, 0, 0.0},
        {"H2O",  18.015,  0, 2, 1, 0, -241825.0},
        {"C2H2", 26.038,  2, 2, 0, 0, 228199.0},
        {"CO",   28.010,  1, 0, 1, 0, -110529.0},
        {"A4",   202.256, 16, 10, 0, 0, 200000.0},
        {"N2",   28.014,  0, 0, 0, 2, 0.0},
    };

    const std::array<double, 7> nasa7_low  = {3.5, 0.0, 0.0, 0.0, 0.0, -1000.0, 3.5};
    const std::array<double, 7> nasa7_high = {3.5, 0.0, 0.0, 0.0, 0.0, -1000.0, 3.5};

    MOM::Manifest m;
    m.magic             = "HMOM-MANIFEST";
    m.schema_version    = 1;
    m.generator_version = "test";
    m.generated_at      = "2026-01-01T00:00:00Z";

    m.config_file   = "test.toml";
    m.config_sha256 = std::string(64, '0');
    m.thermo_file   = "thermo.dat";
    m.thermo_sha256 = std::string(64, '0');

    m.species_count = 9;
    for (const auto& s : seeds)
    {
        MOM::ManifestSpecies sp;
        sp.name        = s.name;
        sp.mw_kg_kmol  = s.mw_kg_kmol;
        sp.mw_kg_mol   = s.mw_kg_kmol / 1000.0;
        sp.nc          = s.nc;
        sp.nh          = s.nh;
        sp.no          = s.no;
        sp.nn          = s.nn;
        sp.h_f0_j_mol  = s.h_f0_j_mol;
        sp.h_f0_j_kg   = s.h_f0_j_mol * 1000.0 / s.mw_kg_kmol;
        sp.nasa7_low   = nasa7_low;
        sp.nasa7_high  = nasa7_high;
        sp.temp_low    = 200.0;
        sp.temp_common = 1000.0;
        sp.temp_high   = 5000.0;
        m.species.push_back(sp);
    }

    m.pah_index  = 7; // A4
    m.bulk_index = 8; // N2
    m.c2h2_index = 5;
    m.h2_index   = 3;
    m.o2_index   = 2;
    m.oh_index   = 1;
    m.co_index   = 6;

    m.nucleation_model             = 1;
    m.condensation_model           = 1;
    m.surface_growth_model         = 1;
    m.oxidation_model              = 1;
    m.coagulation_model            = 1;
    m.continuous_coagulation_model = 1;
    m.thermophoretic_model         = 0;
    m.fractal_diameter_model       = 1;
    m.collision_diameter_model     = 2;
    m.soot_density_kg_m3           = 1800.0;
    m.surface_density_per_m2       = 1.7e19;
    m.surface_density_correction   = 0.0;
    m.sticking_model               = "constant";
    m.sticking_coeff_constant      = 0.002;
    m.schmidt_number               = 50.0;
    m.gas_consumption              = true;

    m.udm_mode  = "release";
    m.udm_count = 18;
    m.uds_count = 4;

    m.operating_pressure_pa = 101325.0;

    m.radiation_zones   = {1};
    m.radiation_enabled = true;
    m.radiation_model   = "Smooke";

    m.reaction_heat_enabled              = true;
    m.soot_formation_enthalpy_j_kg       = 0.0;
    m.require_sensible_enthalpy_equation = true;

    m.manifest_sha256 = "";

    return m;
}

} // namespace mom_test
