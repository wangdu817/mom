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

// ============================================================================
// MomVariantList.hpp — authoritative registry of all MOM variants
// ============================================================================
//
// THIS IS THE ONLY FILE THAT NEEDS TO BE EDITED WHEN ADDING A NEW VARIANT.
//
// To register a new variant "FooMOM":
//
//   Step 1.  Create include/FooMOM/FooMOM.hpp (CRTP subclass of MomentMethodBase)
//            and add the required static member:
//
//              static constexpr std::array<std::string_view, N> variant_labels
//                  { "FooMOM", "foomom" };
//
//   Step 2.  Add one #include line below (under "Registered variant headers").
//
//   Step 3.  Add FooMOM to the AllVariants type alias below.
//
// That is all. The AnyMomentMethod variant type, the MakeAnyMomentMethod
// factory, and the MomentMethod concept checks are all derived automatically
// from AllVariants.
// ============================================================================

// ── Infrastructure ────────────────────────────────────────────────────────────
#include "MomentMethodConcept.hpp" // MomentMethod concept (used in ConceptCheck)

// ── Registered variant headers ────────────────────────────────────────────────
// Add one #include per new variant here.
#include "HMOM/HMOM.hpp"
#include "BrookesMoss/BrookesMoss.hpp"
#include "ThreeEquations/ThreeEquations.hpp"
#include "TiO2/TiO2.hpp"

#include <variant>

namespace MOM
{

// ============================================================================
// detail::TypeList<Variants...>
// ============================================================================
//
// Compile-time list of template-template parameters (uninstantiated variant
// class templates).  Provides:
//
//   TypeList<Vs...>::AsVariant<Thermo>
//       → std::variant<Vs<Thermo>...>
//
//   TypeList<Vs...>::ConceptCheck<TestThermo>
//       A struct whose instantiation static_asserts that every Vs<TestThermo>
//       satisfies the MomentMethod concept.  Trigger by explicit instantiation
//       in MOM.hpp:
//
//           template struct MOM::AllVariants::ConceptCheck<MOM::BasicThermoData>;
//
// ============================================================================

namespace detail
{

template <template <typename> class... Variants> struct TypeList
{
    // ── Type-list operations ─────────────────────────────────────────────────

    /// Number of registered variants.
    static constexpr std::size_t size = sizeof...(Variants);

    /// Instantiate as std::variant<Variants<Thermo>...>.
    template <typename Thermo> using AsVariant = std::variant<Variants<Thermo>...>;

    // ── Compile-time concept gate ─────────────────────────────────────────────
    //
    // Instantiate ConceptCheck<BasicThermoData> (via explicit template
    // instantiation in MOM.hpp) to fire the static_assert for every registered
    // variant in a single statement.  The fold expression fires individually
    // for each Vs, so the compiler error points at the offending variant.

    template <typename TestThermo> struct ConceptCheck
    {
        static_assert((MomentMethod<Variants<TestThermo>> && ...),
                      "[MOM] One or more registered variants in AllVariants do not satisfy "
                      "the MomentMethod concept. Check: SetMoments(span), "
                      "CalculateSourceMoments() noexcept, DiffusionCoefficient() const, "
                      "CollisionDiameter(), SpecificSurface(), PrintSummary().");
    };
};

} // namespace detail

// ============================================================================
// AllVariants — THE ONLY LINE TO CHANGE WHEN ADDING A VARIANT
// ============================================================================
//
// The canonical registry.  Every type listed here:
//   • appears as an alternative in AnyMomentMethod<Thermo>
//   • is tried by MakeAnyMomentMethod via its variant_labels static member
//   • is concept-checked by the explicit instantiation in MOM.hpp
//
// Order is preserved in the variant (first alternative is the default-active
// state after std::variant default construction, which is not meaningful here
// since the factory always constructs explicitly).

using AllVariants = detail::TypeList<HMOM, BrookesMoss, ThreeEquations, TiO2>;

} // namespace MOM
