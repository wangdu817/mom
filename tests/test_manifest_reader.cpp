/*-----------------------------------------------------------------------*\
|   MOM Library -- ManifestReader unit tests (task 3.4)                    |
|                                                                          |
|   Standalone C++ test (plain main(), manual checks with std::exit(1) on  |
|   failure, return 0 on success) exercising MOM::ManifestReader::parse /   |
|   load / validate against a valid 9-species manifest and against a set    |
|   of deliberately-broken variants that must throw ManifestParseError.     |
|                                                                          |
|   Species order (0-based): H OH O2 H2 H2O C2H2 CO A4 N2 (N2 last, bulk). |
\*-----------------------------------------------------------------------*/

#include "MOM/ManifestReader.hpp"
#include "test_manifest_fixture.hpp"

#include <nlohmann/json.hpp>

#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using mom_test::make_test_manifest;

// ============================================================================
// Minimal check helpers (no external test framework, per existing pattern).
// ============================================================================
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

// ----------------------------------------------------------------------------
// to_json(): serialise a Manifest into the JSON layout ManifestReader expects.
// ----------------------------------------------------------------------------
nlohmann::json to_json(const MOM::Manifest& m)
{
    nlohmann::json j;
    j["magic"]             = m.magic;
    j["schema_version"]    = m.schema_version;
    j["generator_version"] = m.generator_version;
    j["generated_at"]      = m.generated_at;

    j["sources"] = {
        {"config_file", m.config_file},
        {"config_sha256", m.config_sha256},
        {"thermo_file", m.thermo_file},
        {"thermo_sha256", m.thermo_sha256},
    };

    j["species_count"] = m.species_count;

    nlohmann::json species = nlohmann::json::array();
    for (const auto& sp : m.species)
    {
        nlohmann::json s;
        s["name"]        = sp.name;
        s["mw_kg_kmol"]  = sp.mw_kg_kmol;
        s["mw_kg_mol"]   = sp.mw_kg_mol;
        s["nc"]          = sp.nc;
        s["nh"]          = sp.nh;
        s["no"]          = sp.no;
        s["nn"]          = sp.nn;
        s["h_f0_j_kg"]   = sp.h_f0_j_kg;
        s["h_f0_j_mol"]  = sp.h_f0_j_mol;
        s["nasa7_low"]   = sp.nasa7_low;
        s["nasa7_high"]  = sp.nasa7_high;
        s["temp_low"]    = sp.temp_low;
        s["temp_common"] = sp.temp_common;
        s["temp_high"]   = sp.temp_high;
        species.push_back(s);
    }
    j["species"] = species;

    j["indices"] = {
        {"pah_index", m.pah_index},
        {"bulk_index", m.bulk_index},
        {"source_species", {
            {"C2H2", m.c2h2_index},
            {"H2", m.h2_index},
            {"O2", m.o2_index},
            {"OH", m.oh_index},
            {"CO", m.co_index},
        }},
    };

    j["hmom_options"] = {
        {"nucleation_model", m.nucleation_model},
        {"condensation_model", m.condensation_model},
        {"surface_growth_model", m.surface_growth_model},
        {"oxidation_model", m.oxidation_model},
        {"coagulation_model", m.coagulation_model},
        {"continuous_coagulation_model", m.continuous_coagulation_model},
        {"thermophoretic_model", m.thermophoretic_model},
        {"fractal_diameter_model", m.fractal_diameter_model},
        {"collision_diameter_model", m.collision_diameter_model},
        {"soot_density_kg_m3", m.soot_density_kg_m3},
        {"surface_density_per_m2", m.surface_density_per_m2},
        {"surface_density_correction", m.surface_density_correction},
        {"sticking_model", m.sticking_model},
        {"sticking_coeff_constant", m.sticking_coeff_constant},
        {"schmidt_number", m.schmidt_number},
        {"gas_consumption", m.gas_consumption},
    };

    j["fluent_options"] = {
        {"udm_mode", m.udm_mode},
        {"udm_count", m.udm_count},
        {"uds_count", m.uds_count},
    };

    j["operating_pressure_pa"] = m.operating_pressure_pa;

    j["radiation"] = {
        {"zones", m.radiation_zones},
        {"enabled", m.radiation_enabled},
        {"model", m.radiation_model},
    };

    j["energy_options"] = {
        {"reaction_heat_enabled", m.reaction_heat_enabled},
        {"soot_formation_enthalpy_j_kg", m.soot_formation_enthalpy_j_kg},
        {"require_sensible_enthalpy_equation", m.require_sensible_enthalpy_equation},
    };

    j["manifest_sha256"] = m.manifest_sha256;
    return j;
}

// Returns true if parsing @p j threw ManifestParseError.
bool parse_throws(const nlohmann::json& j)
{
    try
    {
        MOM::ManifestReader::parse(j.dump());
        return false;
    }
    catch (const MOM::ManifestParseError&)
    {
        return true;
    }
}

} // namespace

int main()
{
    std::cout << "=== ManifestReader unit tests (task 3.4) ===\n";

    // -- 1. Valid manifest parse: every field round-trips. -------------------
    {
        const MOM::Manifest ref = make_test_manifest();
        const MOM::Manifest m   = MOM::ManifestReader::parse(to_json(ref).dump());

        check(m.magic == "HMOM-MANIFEST", "valid: magic");
        check(m.schema_version == 1, "valid: schema_version");
        check(m.generator_version == "test", "valid: generator_version");
        check(m.generated_at == "2026-01-01T00:00:00Z", "valid: generated_at");
        check(m.config_file == "test.toml", "valid: config_file");
        check(m.config_sha256 == std::string(64, '0'), "valid: config_sha256");
        check(m.thermo_file == "thermo.dat", "valid: thermo_file");

        check(m.species_count == 9, "valid: species_count");
        check(m.species.size() == 9, "valid: species array size");
        check(m.species[0].name == "H", "valid: first species H");
        check(m.species[8].name == "N2", "valid: last species N2");
        check(almost_equal(m.species[5].mw_kg_kmol, 26.038), "valid: C2H2 mw");
        check(almost_equal(m.species[5].mw_kg_mol, 0.026038), "valid: C2H2 mw_kg_mol");
        check(m.species[7].nc == 16 && m.species[7].nh == 10, "valid: A4 element counts");
        check(m.species[8].nn == 2, "valid: N2 nitrogen count");
        check(almost_equal(m.species[6].h_f0_j_mol, -110529.0), "valid: CO h_f0_j_mol");
        check(almost_equal(m.species[6].h_f0_j_kg, -110529.0 * 1000.0 / 28.010, 1e-6),
              "valid: CO h_f0_j_kg derived");
        check(almost_equal(m.species[0].nasa7_low[5], -1000.0), "valid: nasa7_low coeff");
        check(almost_equal(m.species[0].temp_common, 1000.0), "valid: temp_common");

        check(m.pah_index == 7, "valid: pah_index");
        check(m.bulk_index == 8, "valid: bulk_index");
        check(m.c2h2_index == 5, "valid: c2h2_index");
        check(m.h2_index == 3, "valid: h2_index");
        check(m.o2_index == 2, "valid: o2_index");
        check(m.oh_index == 1, "valid: oh_index");
        check(m.co_index == 6, "valid: co_index");

        check(m.thermophoretic_model == 0, "valid: thermophoretic_model");
        check(m.sticking_model == "constant", "valid: sticking_model");
        check(m.gas_consumption == true, "valid: gas_consumption");
        check(almost_equal(m.soot_density_kg_m3, 1800.0), "valid: soot_density");

        check(m.udm_mode == "release", "valid: udm_mode");
        check(m.udm_count == 18, "valid: udm_count");
        check(m.uds_count == 4, "valid: uds_count");
        check(almost_equal(m.operating_pressure_pa, 101325.0), "valid: operating_pressure_pa");
        check(m.radiation_zones.size() == 1 && m.radiation_zones[0] == 1, "valid: radiation_zones");
        check(m.radiation_enabled == true, "valid: radiation_enabled");
        check(m.reaction_heat_enabled == true, "valid: reaction_heat_enabled");
        std::cout << "  [ok] 1. valid manifest parse\n";
    }

    // -- 2. Invalid magic -> throws. -----------------------------------------
    {
        nlohmann::json j = to_json(make_test_manifest());
        j["magic"]       = "NOT-HMOM";
        check(parse_throws(j), "invalid magic throws");
        std::cout << "  [ok] 2. invalid magic throws\n";
    }

    // -- 3. Invalid schema_version -> throws. --------------------------------
    {
        nlohmann::json j     = to_json(make_test_manifest());
        j["schema_version"]  = 99;
        check(parse_throws(j), "invalid schema_version throws");
        std::cout << "  [ok] 3. invalid schema_version throws\n";
    }

    // -- 4. Species count mismatch -> throws. --------------------------------
    {
        nlohmann::json j   = to_json(make_test_manifest());
        j["species_count"] = 8; // array still has 9
        check(parse_throws(j), "species count mismatch throws");
        std::cout << "  [ok] 4. species count mismatch throws\n";
    }

    // -- 5. Bulk index not last -> throws. -----------------------------------
    {
        nlohmann::json j            = to_json(make_test_manifest());
        j["indices"]["bulk_index"]  = 4; // not species_count - 1
        check(parse_throws(j), "bulk index not last throws");
        std::cout << "  [ok] 5. bulk index not last throws\n";
    }

    // -- 6. MW <= 0 -> throws. -----------------------------------------------
    {
        nlohmann::json j             = to_json(make_test_manifest());
        j["species"][0]["mw_kg_kmol"] = 0.0;
        check(parse_throws(j), "non-positive mw throws");
        std::cout << "  [ok] 6. mw <= 0 throws\n";
    }

    // -- 7. Invalid udm_count for mode -> throws. ----------------------------
    {
        nlohmann::json j                   = to_json(make_test_manifest());
        j["fluent_options"]["udm_count"]   = 38; // 38 is the debug count, mode is release
        check(parse_throws(j), "invalid udm_count for mode throws");
        std::cout << "  [ok] 7. invalid udm_count for mode throws\n";
    }

    // -- 8. Invalid uds_count -> throws. -------------------------------------
    {
        nlohmann::json j                 = to_json(make_test_manifest());
        j["fluent_options"]["uds_count"] = 3; // must be 4
        check(parse_throws(j), "invalid uds_count throws");
        std::cout << "  [ok] 8. invalid uds_count throws\n";
    }

    // -- 9. File not found -> throws. ----------------------------------------
    {
        bool threw = false;
        try
        {
            MOM::ManifestReader::load("this_manifest_does_not_exist_3p4.json");
        }
        catch (const MOM::ManifestParseError&)
        {
            threw = true;
        }
        check(threw, "load of missing file throws");
        std::cout << "  [ok] 9. file not found throws\n";
    }

    // -- 10. N2 not at bulk_index -> throws. ---------------------------------
    {
        nlohmann::json j            = to_json(make_test_manifest());
        j["species"][8]["name"]     = "AR"; // bulk species is no longer N2
        check(parse_throws(j), "N2 not at bulk_index throws");
        std::cout << "  [ok] 10. N2 not at bulk_index throws\n";
    }

    // -- 11. load() from a real temp file succeeds (round-trip via disk). ----
    {
        const MOM::Manifest ref = make_test_manifest();
        const std::string path  = "test_manifest_reader_tmp.json";
        {
            std::ofstream out(path, std::ios::out | std::ios::trunc);
            check(static_cast<bool>(out), "temp manifest file open");
            out << to_json(ref).dump();
        }
        const MOM::Manifest m = MOM::ManifestReader::load(path);
        check(m.species_count == 9, "load: species_count");
        check(m.species[8].name == "N2", "load: N2 last");
        std::remove(path.c_str());
        std::cout << "  [ok] 11. load() from temp file\n";
    }

    std::cout << "=== ManifestReader: all checks passed (" << g_checks << ") ===\n";
    return 0;
}
