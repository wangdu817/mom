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

#pragma once

#include <array>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace MOM
{

/**
 * @struct ManifestSpecies
 * @brief Per-species record parsed from the HMOM runtime manifest.
 *
 * Mirrors one entry of the manifest ``species`` array emitted by
 * ``tools/manifest_writer.py``. All indices are 0-based and molecular weights
 * are provided in both [kg/kmol] and [kg/mol]; formation enthalpies in both
 * [J/kg] and [J/mol].
 */
struct ManifestSpecies
{
    std::string name;         //!< Species name (mechanism order).
    double mw_kg_kmol;        //!< Molecular weight [kg/kmol].
    double mw_kg_mol;         //!< Molecular weight [kg/mol].
    int nc;                   //!< Carbon atom count per molecule.
    int nh;                   //!< Hydrogen atom count per molecule.
    int no;                   //!< Oxygen atom count per molecule.
    int nn;                   //!< Nitrogen atom count per molecule.
    double h_f0_j_kg;         //!< Standard formation enthalpy [J/kg].
    double h_f0_j_mol;        //!< Standard formation enthalpy [J/mol].
    std::array<double, 7> nasa7_low;  //!< NASA-7 low-temperature coefficients.
    std::array<double, 7> nasa7_high; //!< NASA-7 high-temperature coefficients.
    double temp_low;          //!< Lower temperature bound [K].
    double temp_common;       //!< Common (switch) temperature [K].
    double temp_high;         //!< Upper temperature bound [K].
};

/**
 * @struct Manifest
 * @brief Fully-parsed HMOM runtime manifest consumed by the Fluent UDF bridge.
 *
 * Corresponds to the JSON document emitted by ``tools/manifest_writer.py``
 * (see design-final sections 7-8). The @c manifest_sha256 is computed over the
 * manifest with that field set to the empty string and provides the
 * determinism / integrity guarantee.
 */
struct Manifest
{
    // Header / provenance.
    std::string magic;              //!< Must equal ManifestReader::kMagic.
    int schema_version;             //!< Must equal ManifestReader::kSchemaVersion.
    std::string generator_version;  //!< Writer version string.
    std::string generated_at;       //!< ISO 8601 UTC provenance timestamp.

    // Sources.
    std::string config_file;        //!< Relative path to the source config.
    std::string config_sha256;      //!< SHA-256 of the source config file.
    std::string thermo_file;        //!< Relative path to the thermo file.
    std::string thermo_sha256;      //!< SHA-256 of the thermo file.

    // Species table.
    int species_count;                     //!< Number of species entries.
    std::vector<ManifestSpecies> species;  //!< Species records (mechanism order).

    // Index map.
    int pah_index;   //!< 0-based index of the PAH species.
    int bulk_index;  //!< 0-based index of the bulk (N2) species.
    int c2h2_index;  //!< 0-based index of C2H2 source species.
    int h2_index;    //!< 0-based index of H2 source species.
    int o2_index;    //!< 0-based index of O2 source species.
    int oh_index;    //!< 0-based index of OH source species.
    int co_index;    //!< 0-based index of CO source species.

    // HMOM options (all 16).
    int nucleation_model;
    int condensation_model;
    int surface_growth_model;
    int oxidation_model;
    int coagulation_model;
    int continuous_coagulation_model;
    int thermophoretic_model;
    int fractal_diameter_model;
    int collision_diameter_model;
    double soot_density_kg_m3;
    double surface_density_per_m2;
    double surface_density_correction;
    std::string sticking_model;
    double sticking_coeff_constant;
    double schmidt_number;
    bool gas_consumption;

    // Fluent options.
    std::string udm_mode;  //!< "release" or "debug".
    int udm_count;         //!< Derived UDM slot count (18 release / 38 debug).
    int uds_count;         //!< UDS count (4 for HMOM NEq=4).

    // Operating conditions.
    double operating_pressure_pa;  //!< Fixed operating pressure [Pa].

    // Radiation.
    std::vector<int> radiation_zones;  //!< Radiation zone IDs.
    bool radiation_enabled;            //!< Radiation coupling enabled flag.
    std::string radiation_model;       //!< Radiation model name.

    // Energy.
    bool reaction_heat_enabled;             //!< Reaction-heat coupling flag.
    double soot_formation_enthalpy_j_kg;    //!< Soot formation enthalpy [J/kg].
    bool require_sensible_enthalpy_equation;//!< Requires sensible-enthalpy eqn.

    // Self-verification.
    std::string manifest_sha256;  //!< SHA-256 over the manifest (hash field empty).
};

/**
 * @class ManifestParseError
 * @brief Thrown when a manifest cannot be parsed or fails validation.
 */
class ManifestParseError : public std::runtime_error
{
public:
    explicit ManifestParseError(const std::string& message)
        : std::runtime_error(message)
    {
    }
};

/**
 * @class ManifestReader
 * @brief Loads, parses, and validates an HMOM runtime manifest.
 *
 * Declaration only; implementations are provided in task 3.1. The reader
 * verifies the magic string, schema version, self SHA-256, and structural
 * invariants (species count, N2 last, distinct source indices, UDM count
 * consistency with @c udm_mode).
 */
class ManifestReader
{
public:
    /**
     * @brief Load and parse a manifest from a file path.
     * @param path Filesystem path to the manifest JSON.
     * @return The fully-parsed and validated Manifest.
     * @throws ManifestParseError on I/O, parse, or validation failure.
     */
    static Manifest load(const std::string& path);

    /**
     * @brief Parse and validate a manifest from an in-memory JSON string.
     * @param json_text UTF-8 JSON document.
     * @return The fully-parsed and validated Manifest.
     * @throws ManifestParseError on parse or validation failure.
     */
    static Manifest parse(const std::string& json_text);

    /**
     * @brief Validate a parsed manifest's structural invariants.
     * @param m Manifest to validate.
     * @throws ManifestParseError if any invariant is violated.
     */
    static void validate(const Manifest& m);

    static constexpr const char* kMagic = "HMOM-MANIFEST";  //!< Manifest magic string.
    static constexpr int kSchemaVersion = 1;                //!< Supported schema version.
    static constexpr int kUdsCount = 4;                     //!< Fixed UDS count.
    static constexpr int kUdmReleaseCount = 18;             //!< Release-mode UDM count.
    static constexpr int kUdmDebugCount = 38;               //!< Debug-mode UDM count.
};

} // namespace MOM
