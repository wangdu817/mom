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

/**
 * @file Sources.hpp
 * @brief Source term accessors and process activation queries for AnyMomentMethod.
 *
 * Provides:
 * - The total moment source vector and total gas-phase consumption vector.
 * - Zero-copy per-process source spans (nucleation, growth, coagulation,
 *   condensation, sintering).  For models that do not implement a given
 *   process, the CRTP base class returns a static zero span of size
 *   `n_equations` — the caller never needs to check for model type.
 * - Process activation flag queries — return the integer model flag
 *   (0 = off, >0 = active variant index) for each physical process.
 *
 * Included automatically by `MOM/MOM.hpp`.
 *
 * @par Distinguishing capability from activation
 * Two separate questions can be asked per process:
 *
 * - **Type capability** (compile time): `MOM::ModelsOxidation<M>` — does this
 *   model *type* support the process at all?  Answered by the capability
 *   concepts in `MomentMethodConcept.hpp`.
 *
 * - **Instance activation** (runtime): `MOM::GetOxidationModel(m) > 0` — is
 *   the process currently enabled in this model *instance*?  Answered by the
 *   functions in this file.
 *
 * Always prefer the activation check in CFD code; the capability check is
 * primarily used by `MomentMethodReporter` for zero-fallback column tagging.
 *
 * @par Operator-splitting source functions
 * The oxidation-specific splitting functions (`GetOxidationSources`,
 * `GetSourcesWithoutOxidation`, `GetOxidationRateCoefficients`, etc.) are
 * in `Splitting.hpp`, not here, because they carry additional semantics
 * (stiffness removal, analytical sub-step) beyond simple source retrieval.
 */

#include "AnyMomentMethod.hpp"

namespace MOM
{

/**
 * @name Total source term accessors
 * @{
 */

/**
 * @brief Returns a zero-copy span over the total source vector [model-specific units].
 * @return Span of size `n_equations`; valid until next `CalculateSourceMoments()`.
 */
template <ThermoMap Thermo>
[[nodiscard]] inline std::span<const double> GetSources(const AnyMomentMethod<Thermo>& m) noexcept
{
    return std::visit([](const auto& mm) { return mm.sources(); }, m);
}

/** @brief Returns a zero-copy span over gas-phase consumption terms [kg/m³/s]. */
template <ThermoMap Thermo>
[[nodiscard]] inline std::span<const double> GetOmegaGas(const AnyMomentMethod<Thermo>& m) noexcept
{
    return std::visit([](const auto& mm) { return mm.omega_gas(); }, m);
}

/** @} */

/**
 * @name Per-process source term accessors
 *
 * Zero-copy spans into the model's internal per-process source storage.
 * For models that do not implement a given process, the CRTP base-class
 * fallback returns a span over a `static constexpr` zero array — no runtime
 * branch, no allocation, span size equals `n_equations`.
 *
 * Always check `GetXxxModel() > 0` before consuming these values if you need
 * to distinguish "process is off" from "process is on but rate happens to be zero".
 * @{
 */

/** @brief Zero-copy span over nucleation source terms. Zero span if not modelled. */
template <ThermoMap Thermo>
[[nodiscard]] inline std::span<const double>
GetNucleationSources(const AnyMomentMethod<Thermo>& m) noexcept
{
    return std::visit([](const auto& mm) { return mm.sources_nucleation(); }, m);
}

/** @brief Zero-copy span over surface-growth source terms. Zero span if not modelled. */
template <ThermoMap Thermo>
[[nodiscard]] inline std::span<const double>
GetGrowthSources(const AnyMomentMethod<Thermo>& m) noexcept
{
    return std::visit([](const auto& mm) { return mm.sources_growth(); }, m);
}

/** @brief Zero-copy span over coagulation source terms. Zero span if not modelled. */
template <ThermoMap Thermo>
[[nodiscard]] inline std::span<const double>
GetCoagulationSources(const AnyMomentMethod<Thermo>& m) noexcept
{
    return std::visit([](const auto& mm) { return mm.sources_coagulation(); }, m);
}

/** @brief Zero-copy span over condensation source terms. Zero span if not modelled (e.g. BrookesMoss). */
template <ThermoMap Thermo>
[[nodiscard]] inline std::span<const double>
GetCondensationSources(const AnyMomentMethod<Thermo>& m) noexcept
{
    return std::visit([](const auto& mm) { return mm.sources_condensation(); }, m);
}

/** @brief Zero-copy span over sintering source terms. Zero span if not modelled (MetalOxide only). */
template <ThermoMap Thermo>
[[nodiscard]] inline std::span<const double>
GetSinteringSources(const AnyMomentMethod<Thermo>& m) noexcept
{
    return std::visit([](const auto& mm) { return mm.sources_sintering(); }, m);
}

/** @} */

/**
 * @name Process activation flag queries
 *
 * Return the integer model flag for each physical process (0 = off, >0 = active
 * variant index).  These route through `model_X()` CRTP dispatchers in
 * `MomentMethodBase`, which return 0 for models that do not implement that
 * process (e.g. `GetSinteringModel` always returns 0 for all soot models).
 *
 * @par Interpretation of return values
 * - `0` — process is off (either not supported by the model type, or disabled
 *          at runtime via configuration).
 * - `1` — primary/standard model variant.
 * - `2` — alternative/extended model variant (e.g. BrookesMoss-Hall oxidation).
 *
 * @par Example
 * @code
 *   if (MOM::GetOxidationModel(mom_) > 0) {
 *       // oxidation is active — use operator splitting (see Splitting.hpp)
 *       MOM::GetSourcesWithoutOxidation(mom_, src_buf);
 *   }
 * @endcode
 * @{
 */

/** @brief Returns the nucleation model flag (0 = off, >0 = active variant). */
template <ThermoMap Thermo>
[[nodiscard]] inline int GetNucleationModel(const AnyMomentMethod<Thermo>& m) noexcept
{
    return std::visit([](const auto& mm) { return mm.model_nucleation(); }, m);
}

/** @brief Returns the surface-growth model flag (0 = off, >0 = active variant). */
template <ThermoMap Thermo>
[[nodiscard]] inline int GetGrowthModel(const AnyMomentMethod<Thermo>& m) noexcept
{
    return std::visit([](const auto& mm) { return mm.model_growth(); }, m);
}

/** @brief Returns the coagulation model flag (0 = off, >0 = active variant). */
template <ThermoMap Thermo>
[[nodiscard]] inline int GetCoagulationModel(const AnyMomentMethod<Thermo>& m) noexcept
{
    return std::visit([](const auto& mm) { return mm.model_coagulation(); }, m);
}

/** @brief Returns the condensation model flag (0 = off, >0 = active variant). */
template <ThermoMap Thermo>
[[nodiscard]] inline int GetCondensationModel(const AnyMomentMethod<Thermo>& m) noexcept
{
    return std::visit([](const auto& mm) { return mm.model_condensation(); }, m);
}

/** @brief Returns the oxidation model flag (0 = off, >0 = active variant). */
template <ThermoMap Thermo>
[[nodiscard]] inline int GetOxidationModel(const AnyMomentMethod<Thermo>& m) noexcept
{
    return std::visit([](const auto& mm) { return mm.model_oxidation(); }, m);
}

/** @brief Returns the sintering model flag (0 = off, >0 = active; MetalOxide only). */
template <ThermoMap Thermo>
[[nodiscard]] inline int GetSinteringModel(const AnyMomentMethod<Thermo>& m) noexcept
{
    return std::visit([](const auto& mm) { return mm.model_sintering(); }, m);
}

/** @} */

} // namespace MOM
