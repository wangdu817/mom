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

#include <cstdint>
#include <string_view>

namespace MOM
{

// ============================================================================
// Process model flags
// ============================================================================
//
// Strongly-typed enum classes for all physical sub-processes.
// Using enum class (rather than raw int) gives type safety at call sites
// and clear compiler diagnostics when an invalid value is passed.
//
// All concrete moment method classes accept both the enum class and a plain
// int (via converting constructors) for backward compatibility with existing
// CFD code that uses integer flags.
//
// Convention:
//   Off     = 0  (process disabled)
//   Standard = 1 (primary model variant, method-specific)
//   Extended = 2 (alternative / extended model variant, method-specific)
// ============================================================================

enum class NucleationModel : int
{
    Off      = 0,
    Standard = 1, //!< Default nucleation model for each method
    Extended = 2  //!< Alternative nucleation model (e.g. BrookesMoss-Hall)
};

enum class CoagulationModel : int
{
    Off      = 0,
    Standard = 1 //!< Free-molecular + continuum coagulation
};

enum class SurfaceGrowthModel : int
{
    Off      = 0,
    Standard = 1 //!< HACA-based surface growth
};

enum class OxidationModel : int
{
    Off      = 0,
    Standard = 1, //!< OH + O2 oxidation
    Extended = 2  //!< Alternative oxidation model (e.g. BrookesMoss-Hall)
};

enum class CondensationModel : int
{
    Off      = 0,
    Standard = 1 //!< PAH condensation on soot
};

enum class SinteringModel : int
{
    Off      = 0,
    Standard = 1 //!< Viscous flow sintering (TiO2)
};

enum class ThermophoreticModel : int
{
    Off      = 0,
    Standard = 1 //!< Thermophoretic drift in diffusion coefficient
};

// ============================================================================
// Planck mean absorption coefficient models for radiative heat transfer
// ============================================================================

enum class PlanckCoeffModel : int
{
    None   = 0, //!< No radiative contribution from particles
    Smooke = 1, //!< Smooke et al. correlation (default for soot)
    Kent   = 2, //!< Kent & Honnery correlation
    Sazhin = 3  //!< Sazhin et al. correlation
};

// ============================================================================
// String-to-enum parsers
// ============================================================================
//
// Used by SetupFromDictionary() implementations. Returns the Off variant
// if the label is not recognised (never throws).
// ============================================================================

[[nodiscard]] constexpr PlanckCoeffModel PlanckCoeffModelFromString(std::string_view s) noexcept
{
    if (s == "Smooke" || s == "smooke" || s == "SMOOKE")
        return PlanckCoeffModel::Smooke;
    if (s == "Kent" || s == "kent" || s == "KENT")
        return PlanckCoeffModel::Kent;
    if (s == "Sazhin" || s == "sazhin" || s == "SAZHIN")
        return PlanckCoeffModel::Sazhin;
    if (s == "None" || s == "none" || s == "NONE")
        return PlanckCoeffModel::None;
    return PlanckCoeffModel::None;
}

[[nodiscard]] constexpr NucleationModel NucleationModelFromString(std::string_view s) noexcept
{
    if (s == "Standard" || s == "standard" || s == "1")
        return NucleationModel::Standard;
    if (s == "Extended" || s == "extended" || s == "2")
        return NucleationModel::Extended;
    return NucleationModel::Off;
}

[[nodiscard]] constexpr OxidationModel OxidationModelFromString(std::string_view s) noexcept
{
    if (s == "Standard" || s == "standard" || s == "1")
        return OxidationModel::Standard;
    if (s == "Extended" || s == "extended" || s == "2")
        return OxidationModel::Extended;
    return OxidationModel::Off;
}

} // namespace MOM
