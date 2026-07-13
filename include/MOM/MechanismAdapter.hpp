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
#include <string>
#include <vector>

#include "HMOM/HMOM.hpp"
#include "MOM/ManifestReader.hpp"
#include "MOM/ThermoProxy.hpp"

namespace MOM
{

/**
 * @struct MechanismDescriptor
 * @brief Mechanism-specific data that BasicThermoData does not store.
 *
 * BasicThermoData only carries the fields required by the ThermoMap concept
 * (names, molecular weights, atom counts). The HMOM soot model itself needs no
 * more than that, but downstream energy coupling (Q_HMOM, Phase 8) requires the
 * formation enthalpies and NASA-7 polynomial coefficients, together with the
 * index map of the key soot precursor / product species and the fixed operating
 * pressure. This descriptor collects exactly those extra quantities, in
 * mechanism (0-based) order, straight from the parsed Manifest.
 */
struct MechanismDescriptor
{
    // -- Enthalpy data (not in BasicThermoData) ------------------------------
    std::vector<double> h_f0_j_kg;   //!< Standard formation enthalpy [J/kg], mechanism order.
    std::vector<double> h_f0_j_mol;  //!< Standard formation enthalpy [J/mol], mechanism order.

    // -- NASA-7 coefficients (not in BasicThermoData) ------------------------
    std::vector<std::array<double, 7>> nasa7_low;   //!< Low-T NASA-7 coefficients.
    std::vector<std::array<double, 7>> nasa7_high;  //!< High-T NASA-7 coefficients.
    std::vector<double> temp_low;                   //!< Lower temperature bound [K].
    std::vector<double> temp_common;                //!< Common (switch) temperature [K].
    std::vector<double> temp_high;                  //!< Upper temperature bound [K].

    // -- Index map -----------------------------------------------------------
    int pah_index = -1;   //!< 0-based index of the PAH species.
    int bulk_index = -1;  //!< 0-based index of the bulk (N2) species.
    int c2h2_index = -1;  //!< 0-based index of C2H2 source species.
    int h2_index = -1;    //!< 0-based index of H2 source species.
    int o2_index = -1;    //!< 0-based index of O2 source species.
    int oh_index = -1;    //!< 0-based index of OH source species.
    int co_index = -1;    //!< 0-based index of CO source species.

    // -- Operating conditions ------------------------------------------------
    double operating_pressure_pa = 0.0;  //!< Fixed operating pressure [Pa].
};

/**
 * @class MechanismAdapter
 * @brief Bridges a parsed Manifest to the MOM library's runtime data structures.
 *
 * Provides pure static conversions from a validated @c Manifest (task 3.1) to
 * the objects the HMOM soot model consumes:
 *  - @ref build_thermo    → a @c BasicThermoData satisfying the ThermoMap concept;
 *  - @ref build_config    → an @c HMOM<BasicThermoData>::Config;
 *  - @ref build_descriptor→ a @ref MechanismDescriptor with enthalpy / NASA-7 data.
 *
 * All conversions preserve mechanism (0-based) species ordering. Unit handling
 * follows the manifest contract: molecular weights are carried in [kg/kmol] and
 * copied directly; integer atom counts are widened to @c double; the
 * @c surface_density_correction flag is converted from the manifest's @c double
 * (1.0 / 0.0) to the Config @c bool. Spec §5.1 guard rails are enforced: the
 * "golaut" sticking model is rejected and a non-zero thermophoretic model is
 * rejected (v1 must use 0).
 */
class MechanismAdapter
{
public:
    /**
     * @brief Build a BasicThermoData from the manifest species array.
     *
     * Copies names and molecular weights [kg/kmol] verbatim, widens the integer
     * carbon/hydrogen/oxygen/nitrogen atom counts to @c double, and sets the
     * titanium count to 0.0 for every species (no Ti in this mechanism).
     *
     * @param manifest Validated manifest.
     * @return A BasicThermoData in mechanism order.
     */
    static BasicThermoData build_thermo(const Manifest& manifest);

    /**
     * @brief Build an HMOM<BasicThermoData>::Config from the manifest options.
     *
     * Maps every HMOM option present in the manifest onto the Config struct and
     * leaves all Config fields not represented in the manifest at their library
     * defaults.
     *
     * @param manifest Validated manifest.
     * @return A populated Config.
     * @throws ManifestParseError if @c sticking_model is "golaut" (spec §5.1) or
     *         if @c thermophoretic_model is non-zero (spec §5.1, v1 must be 0).
     */
    static HMOM<BasicThermoData>::Config build_config(const Manifest& manifest);

    /**
     * @brief Build a MechanismDescriptor from the manifest.
     *
     * Collects the formation enthalpies, NASA-7 coefficients, temperature
     * bounds, index map, and operating pressure in mechanism order.
     *
     * @param manifest Validated manifest.
     * @return A populated MechanismDescriptor.
     */
    static MechanismDescriptor build_descriptor(const Manifest& manifest);
};

} // namespace MOM
