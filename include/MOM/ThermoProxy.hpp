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

#include <concepts>
#include <string_view>
#include <span>
#include <string>
#include <vector>
#include <functional>

namespace MOM
{

// ============================================================================
// ThermoMap concept
// ============================================================================
//
// Minimal thermodynamics interface required by any moment method.
// This concept fully decouples the MOM library core from OpenSMOKE++.
// Any thermodynamics backend (OpenSMOKE++, Cantera, custom) can be used
// by providing an adapter that satisfies this concept.
//
// All indices are 0-based. IndexOfSpecies returns -1 if absent.
// Molecular weights are in [kg/kmol].
// Atom counts return 0.0 for species that do not contain that element.
// ============================================================================

template <typename T>
concept ThermoMap = requires(const T& t, unsigned i, std::string_view name) {
    // Total number of species in the kinetic mechanism
    { t.NumberOfSpecies() } -> std::convertible_to<unsigned>;

    // 0-based species index; returns -1 if species is absent in the mechanism
    { t.IndexOfSpecies(name) } -> std::convertible_to<int>;

    // Molecular weight of species i [kg/kmol]
    { t.MolecularWeight(i) } -> std::same_as<double>;

    // Atom counts per molecule of species i
    { t.NumberOfCarbonAtoms(i) } -> std::same_as<double>;
    { t.NumberOfHydrogenAtoms(i) } -> std::same_as<double>;
    { t.NumberOfOxygenAtoms(i) } -> std::same_as<double>;
    { t.NumberOfNitrogenAtoms(i) } -> std::same_as<double>;
    { t.NumberOfTitaniumAtoms(i) } -> std::same_as<double>; // returns 0.0 for non-Ti species
};

// ============================================================================
// BasicThermoData
// ============================================================================
//
// A lightweight, self-contained struct satisfying ThermoMap.
// Useful for unit tests, standalone executables, and CI.
// Populate by hand or parse from a minimal input file.
// ============================================================================

struct BasicThermoData
{
    std::vector<std::string> names;
    std::vector<double> mw;
    std::vector<double> nc;
    std::vector<double> nh;
    std::vector<double> no;
    std::vector<double> nn;
    std::vector<double> nti;

    [[nodiscard]] unsigned NumberOfSpecies() const noexcept
    {
        return static_cast<unsigned>(names.size());
    }

    [[nodiscard]] int IndexOfSpecies(std::string_view name) const noexcept
    {
        for (unsigned i = 0; i < names.size(); ++i)
            if (names[i] == name)
                return static_cast<int>(i);
        return -1;
    }

    [[nodiscard]] double MolecularWeight(unsigned i) const noexcept { return mw[i]; }

    [[nodiscard]] double NumberOfCarbonAtoms(unsigned i) const noexcept { return nc[i]; }

    [[nodiscard]] double NumberOfHydrogenAtoms(unsigned i) const noexcept { return nh[i]; }

    [[nodiscard]] double NumberOfOxygenAtoms(unsigned i) const noexcept { return no[i]; }

    [[nodiscard]] double NumberOfNitrogenAtoms(unsigned i) const noexcept { return nn[i]; }

    [[nodiscard]] double NumberOfTitaniumAtoms(unsigned i) const noexcept { return nti[i]; }
};

static_assert(ThermoMap<BasicThermoData>, "BasicThermoData must satisfy ThermoMap");

// ============================================================================
// Additional thermo adapters
// ============================================================================
//
// Wraps OpenSMOKE::ThermodynamicsMap_CHEMKIN to satisfy the ThermoMap concept.
// This adapter lives in the compatibility layer; the MOM library core never
// includes OpenSMOKE++ headers directly.
//
// USAGE:
//   #include "OpenSMOKE++ headers..."
//   #include "MOM/ThermoProxy.hpp"
//
//   OpenSMOKE::ThermodynamicsMap_CHEMKIN osThermo(...);
//   MOM::OpenSMOKEThermo thermo(osThermo);
//   MOM::HMOM model(thermo);
//
// The implementation of each method is in ThermoProxy.cpp
// ============================================================================

} // namespace MOM
