/*-----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Author: Alberto Cuoci <alberto.cuoci@polimi.it>                       |
|   CRECK Modeling Lab <https://www.creckmodeling.polimi.it>              |
|   Department of Chemistry, Materials, and Chemical Engineering          |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of the OpenSMOKEpp library.                         |
|                                                                         |
|   Copyright (C) 2026 Alberto Cuoci.                                     |
|                                                                         |
|   OpenSMOKEpp is free software: you can redistribute it and/or modify   |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   OpenSMOKEpp is distributed in the hope that it will be useful,        |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with OpenSMOKEpp. If not, see <https://www.gnu.org/licenses/>.  |
|                                                                         |
\*-----------------------------------------------------------------------*/

#include "MOM/ManifestReader.hpp"

#include <nlohmann/json.hpp>

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

namespace MOM
{

Manifest ManifestReader::load(const std::string& path)
{
    std::ifstream file(path, std::ios::binary);
    if (!file.is_open())
        throw ManifestParseError("Cannot open manifest file: " + path);

    std::ostringstream buffer;
    buffer << file.rdbuf();
    return parse(buffer.str());
}

Manifest ManifestReader::parse(const std::string& json_text)
{
    nlohmann::json j;
    try
    {
        j = nlohmann::json::parse(json_text);
    }
    catch (const nlohmann::json::exception& e)
    {
        throw ManifestParseError(std::string("JSON parse error: ") + e.what());
    }

    Manifest m;
    try
    {
        m.magic = j.at("magic").get<std::string>();
        m.schema_version = j.at("schema_version").get<int>();
        m.generator_version = j.at("generator_version").get<std::string>();
        m.generated_at = j.at("generated_at").get<std::string>();

        const auto& sources = j.at("sources");
        m.config_file = sources.at("config_file").get<std::string>();
        m.config_sha256 = sources.at("config_sha256").get<std::string>();
        m.thermo_file = sources.at("thermo_file").get<std::string>();
        m.thermo_sha256 = sources.at("thermo_sha256").get<std::string>();

        m.species_count = j.at("species_count").get<int>();

        const auto& species_array = j.at("species");
        m.species.reserve(species_array.size());
        for (const auto& s : species_array)
        {
            ManifestSpecies sp;
            sp.name = s.at("name").get<std::string>();
            sp.mw_kg_kmol = s.at("mw_kg_kmol").get<double>();
            sp.mw_kg_mol = s.at("mw_kg_mol").get<double>();
            sp.nc = s.at("nc").get<int>();
            sp.nh = s.at("nh").get<int>();
            sp.no = s.at("no").get<int>();
            sp.nn = s.at("nn").get<int>();
            sp.h_f0_j_kg = s.at("h_f0_j_kg").get<double>();
            sp.h_f0_j_mol = s.at("h_f0_j_mol").get<double>();

            const auto& low = s.at("nasa7_low");
            const auto& high = s.at("nasa7_high");
            if (low.size() != 7)
                throw ManifestParseError("Species '" + sp.name +
                    "' nasa7_low must have 7 coefficients, got " +
                    std::to_string(low.size()));
            if (high.size() != 7)
                throw ManifestParseError("Species '" + sp.name +
                    "' nasa7_high must have 7 coefficients, got " +
                    std::to_string(high.size()));
            for (std::size_t k = 0; k < 7; ++k)
            {
                sp.nasa7_low[k] = low[k].get<double>();
                sp.nasa7_high[k] = high[k].get<double>();
            }

            sp.temp_low = s.at("temp_low").get<double>();
            sp.temp_common = s.at("temp_common").get<double>();
            sp.temp_high = s.at("temp_high").get<double>();

            m.species.push_back(std::move(sp));
        }

        const auto& indices = j.at("indices");
        m.pah_index = indices.at("pah_index").get<int>();
        m.bulk_index = indices.at("bulk_index").get<int>();
        const auto& src_species = indices.at("source_species");
        m.c2h2_index = src_species.at("C2H2").get<int>();
        m.h2_index = src_species.at("H2").get<int>();
        m.o2_index = src_species.at("O2").get<int>();
        m.oh_index = src_species.at("OH").get<int>();
        m.co_index = src_species.at("CO").get<int>();

        const auto& hmom = j.at("hmom_options");
        m.nucleation_model = hmom.at("nucleation_model").get<int>();
        m.condensation_model = hmom.at("condensation_model").get<int>();
        m.surface_growth_model = hmom.at("surface_growth_model").get<int>();
        m.oxidation_model = hmom.at("oxidation_model").get<int>();
        m.coagulation_model = hmom.at("coagulation_model").get<int>();
        m.continuous_coagulation_model = hmom.at("continuous_coagulation_model").get<int>();
        m.thermophoretic_model = hmom.at("thermophoretic_model").get<int>();
        m.fractal_diameter_model = hmom.at("fractal_diameter_model").get<int>();
        m.collision_diameter_model = hmom.at("collision_diameter_model").get<int>();
        m.soot_density_kg_m3 = hmom.at("soot_density_kg_m3").get<double>();
        m.surface_density_per_m2 = hmom.at("surface_density_per_m2").get<double>();
        m.surface_density_correction = hmom.at("surface_density_correction").get<double>();
        m.sticking_model = hmom.at("sticking_model").get<std::string>();
        m.sticking_coeff_constant = hmom.at("sticking_coeff_constant").get<double>();
        m.schmidt_number = hmom.at("schmidt_number").get<double>();
        m.gas_consumption = hmom.at("gas_consumption").get<bool>();

        const auto& fluent = j.at("fluent_options");
        m.udm_mode = fluent.at("udm_mode").get<std::string>();
        m.udm_count = fluent.at("udm_count").get<int>();
        m.uds_count = fluent.at("uds_count").get<int>();

        m.operating_pressure_pa = j.at("operating_pressure_pa").get<double>();

        const auto& radiation = j.at("radiation");
        m.radiation_zones = radiation.at("zones").get<std::vector<int>>();
        m.radiation_enabled = radiation.at("enabled").get<bool>();
        m.radiation_model = radiation.at("model").get<std::string>();

        const auto& energy = j.at("energy_options");
        m.reaction_heat_enabled = energy.at("reaction_heat_enabled").get<bool>();
        m.soot_formation_enthalpy_j_kg = energy.at("soot_formation_enthalpy_j_kg").get<double>();
        m.require_sensible_enthalpy_equation =
            energy.at("require_sensible_enthalpy_equation").get<bool>();

        m.manifest_sha256 = j.at("manifest_sha256").get<std::string>();
    }
    catch (const ManifestParseError&)
    {
        throw;
    }
    catch (const nlohmann::json::exception& e)
    {
        throw ManifestParseError(std::string("JSON field error: ") + e.what());
    }

    validate(m);
    return m;
}

void ManifestReader::validate(const Manifest& m)
{
    if (m.magic != kMagic)
        throw ManifestParseError("Invalid magic: expected '" + std::string(kMagic) +
            "', got '" + m.magic + "'");

    if (m.schema_version != kSchemaVersion)
        throw ManifestParseError("Unsupported schema_version: expected " +
            std::to_string(kSchemaVersion) + ", got " + std::to_string(m.schema_version));

    if (m.species_count <= 0)
        throw ManifestParseError("species_count must be positive, got " +
            std::to_string(m.species_count));

    if (static_cast<std::size_t>(m.species_count) != m.species.size())
        throw ManifestParseError("species_count (" + std::to_string(m.species_count) +
            ") does not match species array size (" + std::to_string(m.species.size()) + ")");

    if (m.bulk_index != m.species_count - 1)
        throw ManifestParseError("bulk_index must be species_count - 1 (" +
            std::to_string(m.species_count - 1) + "), got " + std::to_string(m.bulk_index));

    struct NamedIndex { const char* name; int value; };
    const NamedIndex source_indices[] = {
        {"c2h2_index", m.c2h2_index},
        {"h2_index", m.h2_index},
        {"o2_index", m.o2_index},
        {"oh_index", m.oh_index},
        {"co_index", m.co_index},
        {"pah_index", m.pah_index},
    };

    for (const auto& idx : source_indices)
    {
        if (idx.value < 0 || idx.value >= m.species_count)
            throw ManifestParseError(std::string(idx.name) + " out of range [0, " +
                std::to_string(m.species_count) + "), got " + std::to_string(idx.value));
        if (idx.value == m.bulk_index)
            throw ManifestParseError(std::string(idx.name) + " (" + std::to_string(idx.value) +
                ") must be distinct from bulk_index (" + std::to_string(m.bulk_index) + ")");
    }

    const std::size_t n = sizeof(source_indices) / sizeof(source_indices[0]);
    for (std::size_t a = 0; a < n; ++a)
        for (std::size_t b = a + 1; b < n; ++b)
            if (source_indices[a].value == source_indices[b].value)
                throw ManifestParseError(std::string(source_indices[a].name) + " and " +
                    source_indices[b].name + " share the same index (" +
                    std::to_string(source_indices[a].value) + ")");

    if (m.species[static_cast<std::size_t>(m.bulk_index)].name != "N2")
        throw ManifestParseError("bulk species (index " + std::to_string(m.bulk_index) +
            ") must be 'N2', got '" +
            m.species[static_cast<std::size_t>(m.bulk_index)].name + "'");

    for (std::size_t i = 0; i < m.species.size(); ++i)
    {
        const ManifestSpecies& sp = m.species[i];
        if (!(sp.mw_kg_kmol > 0.0))
            throw ManifestParseError("Species '" + sp.name + "' has non-positive mw_kg_kmol: " +
                std::to_string(sp.mw_kg_kmol));
        if (!(sp.mw_kg_mol > 0.0))
            throw ManifestParseError("Species '" + sp.name + "' has non-positive mw_kg_mol: " +
                std::to_string(sp.mw_kg_mol));
        if (!std::isfinite(sp.h_f0_j_kg))
            throw ManifestParseError("Species '" + sp.name + "' has non-finite h_f0_j_kg");
        if (!std::isfinite(sp.h_f0_j_mol))
            throw ManifestParseError("Species '" + sp.name + "' has non-finite h_f0_j_mol");

        for (std::size_t k = 0; k < 7; ++k)
        {
            if (!std::isfinite(sp.nasa7_low[k]))
                throw ManifestParseError("Species '" + sp.name +
                    "' has non-finite nasa7_low[" + std::to_string(k) + "]");
            if (!std::isfinite(sp.nasa7_high[k]))
                throw ManifestParseError("Species '" + sp.name +
                    "' has non-finite nasa7_high[" + std::to_string(k) + "]");
        }
    }

    if (m.uds_count != kUdsCount)
        throw ManifestParseError("uds_count must be " + std::to_string(kUdsCount) +
            ", got " + std::to_string(m.uds_count));

    if (m.udm_mode == "release")
    {
        if (m.udm_count != kUdmReleaseCount)
            throw ManifestParseError("udm_count for release mode must be " +
                std::to_string(kUdmReleaseCount) + ", got " + std::to_string(m.udm_count));
    }
    else if (m.udm_mode == "debug")
    {
        if (m.udm_count != kUdmDebugCount)
            throw ManifestParseError("udm_count for debug mode must be " +
                std::to_string(kUdmDebugCount) + ", got " + std::to_string(m.udm_count));
    }
    else
    {
        throw ManifestParseError("Invalid udm_mode: expected 'release' or 'debug', got '" +
            m.udm_mode + "'");
    }

    if (!(m.operating_pressure_pa > 0.0))
        throw ManifestParseError("operating_pressure_pa must be positive, got " +
            std::to_string(m.operating_pressure_pa));
}

} // namespace MOM
