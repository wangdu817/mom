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
 * @file Splitting.hpp
 * @brief Operator-splitting support functions for AnyMomentMethod.
 *
 * Provides the complete set of functions required to implement Lie–Trotter
 * (or Strang) operator splitting that removes the stiff oxidation eigenvalue
 * from the moment ODE system.
 *
 * Included automatically by `MOM/MOM.hpp`.
 *
 * @par Motivation
 * Oxidation is a first-order sink on all moment variables.  Its linearised
 * rate coefficient κ_i = |source_oxidation[i]| / |M_i| can reach 10⁵–10⁶ s⁻¹
 * in a flame — orders of magnitude larger than the CFD time step allows in an
 * implicit BDF solver.  The Patankar limiter bounds values but not eigenvalues,
 * so the BDF solver still brackets the C⁰ kink and produces tiny sub-steps.
 * Operator splitting removes the eigenvalue entirely:
 *
 * **Step A** — ODE sub-step with non-oxidation sources only:
 * @code
 *   double src_buf[MAX_MOM];
 *   MOM::Compute(mom);
 *   MOM::GetSourcesWithoutOxidation(mom, {src_buf, n_mom});
 *   // pass src_buf to the stiff ODE solver — no large eigenvalue
 * @endcode
 *
 * **Step B** — analytical oxidation sub-step after the ODE completes:
 * @code
 *   double kappa[MAX_MOM];
 *   MOM::GetOxidationRateCoefficients(mom, {M, n_mom}, {kappa, n_mom});
 *   for (int i = 0; i < n_mom; ++i)
 *       M[i] *= std::exp(-kappa[i] * dt);
 * @endcode
 *
 * **Step C** — gas-phase correction with saturation factor φ:
 * @code
 *   // Inside Equations(), pass the non-oxidation gas sources only:
 *   double gas_buf[MAX_SPECIES];
 *   MOM::GetOmegaGasWithoutOxidation(mom, {gas_buf, n_species});
 *
 *   // After Step B, apply the integrated oxidation contribution:
 *   const double x   = kappa[MASS_IDX] * dt;   // use the mass-fraction moment
 *   const double phi = (x > 1e-8) ? (1. - std::exp(-x)) / x : 1.;
 *   auto ox_gas = MOM::GetOmegaGasOxidation(mom);
 *   for (int k = 0; k < n_species; ++k)
 *       Y[k] += ox_gas[k] / rho * dt * phi;
 * @endcode
 *
 * @par When to use this header
 * Include `MOM/MOM.hpp` (which pulls in `Splitting.hpp` automatically) when
 * compiling any translation unit that uses operator splitting.  If only dispatch
 * or simple source retrieval is needed, those translation units compile faster
 * including only `Dispatch.hpp` or `Sources.hpp` directly.
 *
 * @par Functions provided
 * - `GetOxidationSources`             — zero-copy span of oxidation-only moment sources
 * - `GetSourcesWithoutOxidation`      — total minus oxidation, written into caller buffer
 * - `GetOxidationRateCoefficients`    — first-order rate coefficients κ_i [1/s]
 * - `GetOmegaGasOxidation`            — zero-copy span of oxidation-only gas-phase sources
 * - `GetOmegaGasWithoutOxidation`     — total gas minus oxidation, written into caller buffer
 */

#include <algorithm>   // std::min
#include <cmath>       // std::abs, std::exp (std::exp used in doc; not called here)

#include "AnyMomentMethod.hpp"

namespace MOM
{

/**
 * @name Operator-splitting — moment-space functions
 * @{
 */

/**
 * @brief Returns a zero-copy span over the **oxidation-only** moment source vector.
 *
 * Points directly into the model's internal `source_oxidation_` storage —
 * zero-overhead.  All entries are ≤ 0 (oxidation is always a sink for soot).
 *
 * @pre  Compute() must have been called at the current state.
 * @return Span of size `n_equations`; valid until next Compute().
 */
template <ThermoMap Thermo>
[[nodiscard]] inline std::span<const double>
GetOxidationSources(const AnyMomentMethod<Thermo>& m) noexcept
{
    return std::visit([](const auto& mm) { return mm.sources_oxidation(); }, m);
}

/**
 * @brief Writes `source_all[i] - source_oxidation[i]` into @p out for each moment.
 *
 * For operator splitting of the stiff oxidation terms: pass these reduced sources
 * to the stiff ODE solver instead of the full `GetSources()` vector, then apply
 * the oxidation sub-step analytically with `GetOxidationRateCoefficients()`.
 *
 * @pre  Compute() must have been called at the current state.
 * @param[in]  m    Model variant.
 * @param[out] out  Caller-allocated buffer; size must be >= n_equations.
 */
template <ThermoMap Thermo>
inline void GetSourcesWithoutOxidation(const AnyMomentMethod<Thermo>& m,
                                        std::span<double> out) noexcept
{
    std::visit(
        [&out](const auto& mm)
        {
            const auto total = mm.sources();           // source_all_ — zero-copy
            const auto ox    = mm.sources_oxidation(); // source_oxidation_ — zero-copy
            const std::size_t N = std::min(total.size(), out.size());
            for (std::size_t i = 0; i < N; ++i)
                out[i] = total[i] - ox[i];
        },
        m);
}

/**
 * @brief Computes the effective first-order oxidation rate coefficient [1/s] per moment.
 *
 * Linearises the (generally nonlinear) oxidation depletion as a first-order decay:
 *
 *   κ_i = max(−source_oxidation[i], 0) / max(|M_i|, ε)
 *
 * The exact analytical solution for the oxidation sub-step then reads:
 *
 *   M_i(t + Δt) = M_i(t) · exp(−κ_i · Δt)
 *
 * This is unconditionally stable regardless of how large κ_i is — it removes the
 * stiff oxidation eigenvalue from the ODE system entirely.
 *
 * @note  κ_i is evaluated at the **current** stored state (the last Compute() call).
 *        For Lie–Trotter splitting, call this after the ODE step has finished.
 *        For symmetric Strang splitting, evaluate at t + Δt/2.
 *
 * @pre   Compute() must have been called at the current state.
 * @param[in]  m                Model variant.
 * @param[in]  current_moments  Transported moment values M_i at current time.
 *                              Size must be >= n_equations.
 * @param[out] kappa_out        Output buffer for κ_i [1/s].  Size >= n_equations.
 */
template <ThermoMap Thermo>
inline void GetOxidationRateCoefficients(const AnyMomentMethod<Thermo>& m,
                                          std::span<const double> current_moments,
                                          std::span<double>       kappa_out) noexcept
{
    std::visit(
        [&current_moments, &kappa_out](const auto& mm)
        {
            const auto ox = mm.sources_oxidation();
            const std::size_t N =
                std::min({ox.size(), current_moments.size(), kappa_out.size()});
            constexpr double eps = 1.e-300;
            for (std::size_t i = 0; i < N; ++i)
            {
                // Oxidation is a sink: source_ox[i] ≤ 0.  Rate coeff is non-negative.
                const double neg_src = -ox[i];
                const double M       = std::max(std::abs(current_moments[i]), eps);
                kappa_out[i] = (neg_src > 0.) ? neg_src / M : 0.;
            }
        },
        m);
}

/** @} */

/**
 * @name Operator-splitting — gas-phase functions
 * @{
 */

/**
 * @brief Returns a zero-copy span over the **oxidation-only** gas-phase source vector [kg/m³/s].
 *
 * Points directly into internal storage — zero-overhead.
 * Returns an empty span for models without oxidation gas coupling (e.g. MetalOxide/TiO₂).
 *
 * @pre Compute() must have been called at the current state.
 */
template <ThermoMap Thermo>
[[nodiscard]] inline std::span<const double>
GetOmegaGasOxidation(const AnyMomentMethod<Thermo>& m) noexcept
{
    return std::visit([](const auto& mm) { return mm.omega_gas_oxidation(); }, m);
}

/**
 * @brief Writes `omega_gas[k] − omega_gas_oxidation[k]` into @p out for each species.
 *
 * For operator splitting: pass these reduced gas sources inside `Equations()`,
 * then apply the oxidation contribution analytically after the ODE step completes.
 *
 * @par Post-step gas-phase correction with saturation factor φ
 * After the outer CFD timestep @p dt has completed and the soot moments have been
 * updated with the exponential decay (Step B above), apply the integrated oxidation
 * effect on gas species:
 *
 * @code
 *   // kappa[MASS_IDX]: from GetOxidationRateCoefficients() for the mass-fraction moment
 *   //   MASS_IDX = 0 for ThreeEquations / BrookesMoss; = 1 for HMOM
 *   const double x   = kappa[MASS_IDX] * dt;
 *   const double phi = (x > 1e-8) ? (1. - std::exp(-x)) / x : 1.;
 *   auto ox_gas = MOM::GetOmegaGasOxidation(mom);
 *   for (int k = 0; k < n_species; ++k)
 *       Y[k] += ox_gas[k] / rho * dt * phi;
 * @endcode
 *
 * The factor φ ∈ (0, 1] collapses to 1 for small `κ·Δt` (first-order approximation)
 * and falls toward `1/(κ·Δt)` for fast oxidation, bounding the integrated gas
 * consumption to at most the available reactant mass.
 *
 * @pre  Compute() must have been called at the current state.
 * @param[in]  m    Model variant.
 * @param[out] out  Caller-allocated buffer; size must be >= n_species.
 */
template <ThermoMap Thermo>
inline void GetOmegaGasWithoutOxidation(const AnyMomentMethod<Thermo>& m,
                                         std::span<double> out) noexcept
{
    std::visit(
        [&out](const auto& mm)
        {
            const auto total = mm.omega_gas();           // full omega_gas_   — zero-copy
            const auto ox    = mm.omega_gas_oxidation(); // oxidation-only    — zero-copy
            const std::size_t N = std::min(total.size(), out.size());
            for (std::size_t k = 0; k < N; ++k)
                out[k] = total[k] - (k < ox.size() ? ox[k] : 0.);
        },
        m);
}

/** @} */

} // namespace MOM
