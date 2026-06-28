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

// MomVariantList.hpp is the single authoritative registry of all variants.
// It transitively provides: all variant headers, MomentMethodConcept.hpp,
// detail::TypeList, AllVariants, and <variant>.
#include "MomVariantList.hpp"
#include <functional>
#include <span>
#include <string_view>
#include <stdexcept>
#include <string>

namespace MOM {

// ============================================================================
// AnyMomentMethod<Thermo> — runtime-selectable moment method
// ============================================================================
//
// A std::variant over all four concrete moment method types. Enables
// runtime selection of the method from an input file or command-line flag
// without duplicating CFD solver code.
//
// Dispatch is via std::visit, which generates a jump-table (one indirect
// branch) — comparable to a vtable call but without aliasing penalties that
// prevent auto-vectorisation of the source computation loop.
//
// USAGE:
//
//   // Construction (method chosen at runtime):
//   MOM::AnyMomentMethod<MyThermo> model =
//       MOM::MakeAnyMomentMethod<MyThermo>(thermo, "HMOM");
//
//   // Preferred per-cell call — one std::visit for the full update:
//   MOM::ComputeCell(model, T, P, Y, mu, cell_moments);
//   auto src = MOM::GetSources(model);       // std::span<const double>
//   auto n   = MOM::GetNEquations(model);    // unsigned
//
//   // Individual setters remain available when only part of the state changes:
//   MOM::SetState(model, T, P, Y);
//   MOM::SetMoments(model, cell_moments);
//   MOM::Compute(model);
//
// ============================================================================

// ── AnyMomentMethod<Thermo> — derived automatically from AllVariants ──────────
//
// Expands to std::variant<HMOM<Thermo>, BrookesMoss<Thermo>, ...> with the
// exact set of types registered in MomVariantList.hpp::AllVariants.
// No manual edit required here when adding a new variant.

template <ThermoMap Thermo>
using AnyMomentMethod = typename AllVariants::template AsVariant<Thermo>;

// ============================================================================
// detail::FactoryHelper — compile-time recursive label dispatcher
// ============================================================================
//
// Iterates the registered variants at compile time, checking each type's
// variant_labels member against the runtime label string.  Falls through
// recursively until a match is found or the list is exhausted.
//
// Derived from AllVariants in MakeAnyMomentMethod, so it automatically
// covers every registered variant without manual if-chains.
// ============================================================================

namespace detail {

// Base case: empty list — no variant matched the label.
template <template<typename> class... Vs>
struct FactoryHelper
{
    template <typename Thermo>
    [[noreturn]] static AnyMomentMethod<Thermo>
    make(const Thermo&, std::string_view label)
    {
        throw std::invalid_argument(
            "MOM::MakeAnyMomentMethod: unknown method label '" +
            std::string(label) +
            "'. See MomVariantList.hpp for registered variants and their labels.");
    }
};

// Recursive case: try First's labels, then recurse into Rest...
template <template<typename> class First, template<typename> class... Rest>
struct FactoryHelper<First, Rest...>
{
    template <typename Thermo>
    static AnyMomentMethod<Thermo>
    make(const Thermo& thermo, std::string_view label)
    {
        for (std::string_view lbl : First<Thermo>::variant_labels)
            if (lbl == label)
                return AnyMomentMethod<Thermo>{ std::in_place_type<First<Thermo>>, thermo };
        return FactoryHelper<Rest...>::template make<Thermo>(thermo, label);
    }
};

// Unpack TypeList<Vs...> → FactoryHelper<Vs...>::make
// This allows MakeAnyMomentMethod to delegate to the correct FactoryHelper
// specialisation without naming the individual variant types.
template <template<typename> class... Vs, typename Thermo>
inline AnyMomentMethod<Thermo>
make_from_type_list(TypeList<Vs...>, const Thermo& thermo, std::string_view label)
{
    return FactoryHelper<Vs...>::template make<Thermo>(thermo, label);
}

} // namespace detail

// ============================================================================
// Factory function
// ============================================================================
//
// Constructs an AnyMomentMethod holding the concrete type whose variant_labels
// contains @p label (exact, case-sensitive match).
//
// @param thermo  Thermodynamics object (must outlive the returned variant)
// @param label   Any label registered in a variant's variant_labels member.
//                See MomVariantList.hpp for all registered variants.
//
// @throws std::invalid_argument if label is not recognised.
// ============================================================================

template <ThermoMap Thermo>
[[nodiscard]] AnyMomentMethod<Thermo>
MakeAnyMomentMethod(const Thermo& thermo, std::string_view label);

// ============================================================================
// Generic dispatch helpers (free functions, not member functions)
// ============================================================================
//
// These reproduce the full MomentMethod concept interface as free functions
// operating on AnyMomentMethod. Each is a one-liner std::visit wrapper.
// They are the "uniform call site" for CFD solver code that works with any
// method without template parameters.
//
// All helpers have O(1) dispatch cost (one indirect branch via jump table).
// ============================================================================

/// Returns the number of transported equations for the active method.
template <ThermoMap Thermo>
[[nodiscard]] inline unsigned GetNEquations(const AnyMomentMethod<Thermo>& m) noexcept {
    return std::visit([](const auto& mm) -> unsigned { return mm.n_equations; }, m);
}

/// Injects thermodynamic state into the active method.
template <ThermoMap Thermo>
inline void SetState(AnyMomentMethod<Thermo>& m,
                     double T, double P_Pa, const double* Y) noexcept {
    std::visit([&](auto& mm) { mm.SetStatus(T, P_Pa, Y); }, m);
}

/// Sets the mixture dynamic viscosity [kg/m/s].
template <ThermoMap Thermo>
inline void SetViscosity(AnyMomentMethod<Thermo>& m, double mu) noexcept {
    std::visit([mu](auto& mm) { mm.SetViscosity(mu); }, m);
}

/// Sets current moment values from a contiguous span.
template <ThermoMap Thermo>
inline void SetMoments(AnyMomentMethod<Thermo>& m,
                       std::span<const double> moments) noexcept {
    std::visit([moments](auto& mm) { mm.SetMoments(moments); }, m);
}

/// Computes all source terms and gas-phase consumption.
template <ThermoMap Thermo>
inline void Compute(AnyMomentMethod<Thermo>& m) {
    std::visit([](auto& mm) { mm.CalculateSourceMoments(); }, m);
}

// ============================================================================
// ComputeCell — single-call per-cell entry point (runtime path)
// ============================================================================
//
// Collapses the four-call sequence (SetState/SetMoments/SetViscosity/Compute)
// into a single std::visit, reducing jump-table dispatches from four to one.
//
// For the compile-time path (type known statically), use the MomentMethod-
// constrained overload in MomentMethodConcept.hpp instead.
// ============================================================================

template <ThermoMap Thermo>
inline void ComputeCell(AnyMomentMethod<Thermo>& m,
                        double                   T,
                        double                   P_Pa,
                        const double*            Y,
                        double                   mu,
                        std::span<const double>  moments) noexcept
{
    std::visit([T, P_Pa, Y, mu, moments](auto& mm) noexcept {
        mm.SetStatus(T, P_Pa, Y);
        mm.SetMoments(moments);
        mm.SetViscosity(mu);
        mm.CalculateSourceMoments();
    }, m);
}

/// Write header line
template <ThermoMap Thermo>
inline void WriteHeaderLine(AnyMomentMethod<Thermo>& m, MOM::OutputFileColumns& fOutput, const unsigned int precision) {
    std::visit([&fOutput, precision](auto& mm) { mm.WriteHeaderLine(fOutput, precision); }, m);
}

/// Write header line
template <ThermoMap Thermo>
inline void WriteOutputLine(AnyMomentMethod<Thermo>& m, 
								MOM::OutputFileColumns& fOutput, const double T, const double P_Pa, const double* Y, const double mu, const double* M) {
    std::visit([&fOutput, T, P_Pa, Y, mu, M](auto& mm) { mm.WriteOutputLine(fOutput, T, P_Pa, Y, mu, M); }, m);
}

#if defined(MOM_USE_DICTIONARY)
/// Setup from a dictionary
template <ThermoMap Thermo, typename Dictionary>
[[nodiscard]] inline std::expected<void, std::string>
SetupFromDictionary(AnyMomentMethod<Thermo>& m, Dictionary& dict)
{
    return std::visit([&dict](auto& mm) -> std::expected<void, std::string> {
        return mm.SetupFromDictionary(dict);
    }, m);
}
#endif

template <ThermoMap Thermo>
[[nodiscard]] inline int GetThermophoreticModel(const AnyMomentMethod<Thermo>& m) noexcept {
    return std::visit([](const auto& mm) { return mm.thermophoretic_model(); }, m);
}

template <ThermoMap Thermo>
[[nodiscard]] inline bool GetGasConsumption(const AnyMomentMethod<Thermo>& m) noexcept {
    return std::visit([](const auto& mm) { return mm.GasConsumption(); }, m);
}

/// Returns a span over the total source vector (all processes summed).
template <ThermoMap Thermo>
[[nodiscard]] inline std::span<const double>
GetSources(const AnyMomentMethod<Thermo>& m) noexcept {
    return std::visit([](const auto& mm) { return mm.sources(); }, m);
}

/// Returns a span over the gas-phase source terms [kg/m3/s].
template <ThermoMap Thermo>
[[nodiscard]] inline std::span<const double>
GetOmegaGas(const AnyMomentMethod<Thermo>& m) noexcept {
    return std::visit([](const auto& mm) { return mm.omega_gas(); }, m);
}

/// Returns soot/particle volume fraction [-].
template <ThermoMap Thermo>
[[nodiscard]] inline double GetVolumeFraction(const AnyMomentMethod<Thermo>& m) noexcept {
    return std::visit([](const auto& mm) { return mm.VolumeFraction(); }, m);
}

/// Returns mean primary particle diameter [m].
template <ThermoMap Thermo>
[[nodiscard]] inline double GetParticleDiameter(const AnyMomentMethod<Thermo>& m) noexcept {
    return std::visit([](const auto& mm) { return mm.ParticleDiameter(); }, m);
}

/// Returns particle number density [#/m3].
template <ThermoMap Thermo>
[[nodiscard]] inline double GetParticleNumberDensity(const AnyMomentMethod<Thermo>& m) noexcept {
    return std::visit([](const auto& mm) { return mm.ParticleNumberDensity(); }, m);
}

/// Returns particle mass fraction [-].
template <ThermoMap Thermo>
[[nodiscard]] inline double GetMassFraction(const AnyMomentMethod<Thermo>& m) noexcept {
    return std::visit([](const auto& mm) { return mm.MassFraction(); }, m);
}

/// Returns the Schmidt number for the particle phase [-].
template <ThermoMap Thermo>
[[nodiscard]] inline double GetSchmidtNumber(const AnyMomentMethod<Thermo>& m) noexcept {
    return std::visit([](const auto& mm) { return mm.schmidt_number(); }, m);
}

/// Returns the effective diffusion coefficient [kg/m/s].
template <ThermoMap Thermo>
[[nodiscard]] inline double GetDiffusionCoefficient(const AnyMomentMethod<Thermo>& m) noexcept {
    return std::visit([](const auto& mm) { return mm.DiffusionCoefficient(); }, m);
}

/// Returns true if the Planck absorption coefficient should be included.
template <ThermoMap Thermo>
[[nodiscard]] inline bool GetRadiativeHeatTransfer(const AnyMomentMethod<Thermo>& m) noexcept {
    return std::visit([](const auto& mm) { return mm.radiative_heat_transfer(); }, m);
}

/// Returns the Planck absorption coefficient [1/m].
template <ThermoMap Thermo>
[[nodiscard]] inline double GetPlanckCoefficient(const AnyMomentMethod<Thermo>& m,
                                                  double T, double fv) noexcept {
    return std::visit([T, fv](const auto& mm) {
        return mm.planck_coefficient(T, fv);
    }, m);
}

/// Returns the initial moment values for solver initialisation.
template <ThermoMap Thermo>
[[nodiscard]] inline std::span<const double>
GetInitialMoments(const AnyMomentMethod<Thermo>& m) noexcept {
    return std::visit([](const auto& mm) { return mm.initial_moments(); }, m);
}

/// Returns true if a gas-closure dummy species is active.
template <ThermoMap Thermo>
[[nodiscard]] inline bool GetClosureDummySpeciesIsActive(
    const AnyMomentMethod<Thermo>& m) noexcept
{
    return std::visit([](const auto& mm) {
        return mm.ClosureDummySpeciesIsActive();
    }, m);
}

/// Returns 0-based index of the gas-closure dummy species (-1 if inactive).
template <ThermoMap Thermo>
[[nodiscard]] inline int GetClosureDummyIndex(
    const AnyMomentMethod<Thermo>& m) noexcept
{
    return std::visit([](const auto& mm) { return mm.closure_dummy_index(); }, m);
}

/// Returns 0-based precursor species index.
template <ThermoMap Thermo>
[[nodiscard]] inline int GetPrecursorIndex(const AnyMomentMethod<Thermo>& m) noexcept {
    return std::visit([](const auto& mm) { return mm.precursor_index(); }, m);
}

/// Prints the model configuration summary.
template <ThermoMap Thermo>
inline void PrintSummary(const AnyMomentMethod<Thermo>& m) {
    std::visit([](const auto& mm) { mm.PrintSummary(); }, m);
}

// ============================================================================
// Compile-time variant for performance-critical inner loops
// ============================================================================
//
// When the method IS known at compile time, use this alias directly in a
// templated cell-loop function. This allows full inlining and vectorisation
// of CalculateSourceMoments().
//
// Pattern:
//
//   using ParticleModel = MOM::ThreeEquations<MyThermo>;  // ← one line to switch
//   static_assert(MOM::MomentMethod<ParticleModel>);
//
//   template <MOM::MomentMethod M>
//   void CellLoop(M& model, /* ... */) {
//       for (auto& cell : cells) {
//           model.SetStatus(cell.T, cell.P, cell.Y.data());
//           model.SetMoments(cell.moments);
//           model.CalculateSourceMoments();
//           // access model.sources() — zero-copy span
//       }
//   }
//
// ============================================================================

} // namespace MOM

// AnyMomentMethod<Thermo> is a thin std::variant alias — its only template body
// is the MakeAnyMomentMethod factory (AnyMomentMethod.tpp).
// That factory is trivially inlined and does not benefit from pre-compilation,
// so we always include it (no extern template needed here).
#include "AnyMomentMethod.tpp"
