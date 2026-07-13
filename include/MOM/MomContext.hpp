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
#include <memory>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "HMOM/HMOM.hpp"
#include "MOM/ManifestReader.hpp"
#include "MOM/MechanismAdapter.hpp"
#include "MOM/ThermoProxy.hpp"

namespace MOM
{

/**
 * @struct MomContext
 * @brief Rank-local owner of the Thermo + HMOM model + scratch buffers.
 *
 * A @c MomContext bundles everything a single MPI rank needs to evaluate HMOM
 * soot source terms cell-by-cell:
 *  - the @ref BasicThermoData built from the runtime manifest (task 3.2);
 *  - the configured @c HMOM<BasicThermoData> model;
 *  - the @ref MechanismDescriptor metadata (enthalpy / NASA-7 / index map);
 *  - pre-allocated scratch buffers to avoid per-cell heap allocation.
 *
 * @par Lifetime (design-final §5.1)
 * @c thermo_ is declared before @c model_ so that it is constructed first and
 * destroyed last: the HMOM constructor stores a @c const @c Thermo& that must
 * outlive the model (HMOM.hpp:288). @c model_ is a @c std::unique_ptr because
 * the HMOM constructor requires the already-constructed @c thermo_.
 *
 * @par Error handling
 * The constructor never throws or crashes. If manifest loading, mechanism
 * adaptation, or HMOM setup fails, the exception is caught, @c valid_ is set to
 * @c false, and the message is stored in @c error_message_. The caller must
 * check @ref is_valid() before using the context; @ref compute() is a no-op on
 * an invalid context.
 *
 * @par Threading
 * Not copyable, not movable, and not thread-safe. Each MPI rank constructs its
 * own @c MomContext from the manifest file and drives it single-threaded. No
 * static or global state is used.
 */
struct MomContext final
{
    // -- Construction --------------------------------------------------------

    /**
     * @brief Construct a context by loading a manifest from a file.
     *
     * Loads and validates the manifest, builds the thermo/config/descriptor via
     * @ref MechanismAdapter, constructs and configures the HMOM model, and
     * pre-allocates scratch buffers. On any failure the context is left in the
     * invalid state (see @ref is_valid, @ref error_message).
     *
     * @param manifest_path Filesystem path to the manifest JSON.
     */
    explicit MomContext(const std::string& manifest_path);

    /**
     * @brief Construct a context directly from an already-parsed manifest.
     *
     * Identical to the file-path constructor but skips the I/O step. Intended
     * for testing without a manifest file on disk.
     *
     * @param manifest Fully-parsed and validated manifest.
     */
    explicit MomContext(const Manifest& manifest);

    // -- Deleted copy/move (lifetime safety, design-final P0-03) -------------
    MomContext(const MomContext&)            = delete;
    MomContext& operator=(const MomContext&) = delete;
    MomContext(MomContext&&)                 = delete;
    MomContext& operator=(MomContext&&)      = delete;

    // -- Per-cell computation -------------------------------------------------

    /**
     * @brief Evaluate HMOM source terms for a single cell.
     *
     * Wraps the noexcept @ref ComputeCell free function. On an invalid context
     * this is a no-op (the model is not touched). After a successful call the
     * results are read via the accessor methods below.
     *
     * @param T       Gas temperature [K].
     * @param P_Pa    Gas pressure [Pa].
     * @param Y       Species mass fractions (size = @ref n_species).
     * @param mu      Mixture dynamic viscosity [kg/m/s].
     * @param moments Current moment values [M00, M10, M01, N0].
     */
    void compute(double T, double P_Pa, std::span<const double> Y, double mu,
                 std::span<const double> moments) noexcept;

    // -- Accessors (delegate to model_) --------------------------------------

    /** @brief Source moments [S_M00, S_M10, S_M01, S_N0] [mol/m³/s], size 4. */
    [[nodiscard]] std::span<const double> sources() const noexcept;

    /** @brief Gas-phase source terms [kg/m³/s], size = @ref n_species. */
    [[nodiscard]] std::span<const double> omega_gas() const noexcept;

    /** @brief Primary particle (nucleation-mode) volume [m³]. */
    [[nodiscard]] double V0() const noexcept;

    /** @brief Soot material density [kg/m³]. */
    [[nodiscard]] double particle_density() const noexcept;

    /** @brief Effective soot diffusion coefficient Γ_mom [kg/m/s]. */
    [[nodiscard]] double diffusion_coefficient() const noexcept;

    /** @brief Planck mean absorption coefficient [1/m] for radiation coupling.
     *  Dispatches to the model's planck_coefficient (default: Smooke). */
    [[nodiscard]] double planck_coefficient(double T, double fv) const noexcept;

    /** @brief Soot volume fraction [-] (model's M10-based volume fraction). */
    [[nodiscard]] double volume_fraction() const noexcept;

    /** @brief Primary particle diameter [m] (model's particle_diameter). */
    [[nodiscard]] double particle_diameter() const noexcept;

    /** @brief Collision diameter [m] (model's collision_diameter). */
    [[nodiscard]] double collision_diameter() const noexcept;

    /** @brief Initial moments [M00, M10, M01, N0] [mol/m^3] for fresh init. */
    [[nodiscard]] std::span<const double> initial_moments() const noexcept;

    /** @brief Mechanism metadata (enthalpy / NASA-7 / index map). */
    [[nodiscard]] const MechanismDescriptor& descriptor() const noexcept { return descriptor_; }

    /** @brief Radiation zone IDs from the manifest. */
    [[nodiscard]] std::span<const int> zone_ids() const noexcept { return zone_ids_; }

    /** @brief Manifest SHA-256 (integrity / determinism guarantee). */
    [[nodiscard]] std::string_view manifest_hash() const noexcept { return manifest_hash_; }

    /** @brief True if construction succeeded and the context is usable. */
    [[nodiscard]] bool is_valid() const noexcept { return valid_; }

    /** @brief Error description if invalid; empty if valid. */
    [[nodiscard]] std::string_view error_message() const noexcept { return error_message_; }

    /** @brief Number of species in the mechanism. */
    [[nodiscard]] unsigned n_species() const noexcept;

    /** @brief Species names in mechanism order. */
    [[nodiscard]] std::span<const std::string> species_names() const noexcept { return thermo_.names; }

    // -- Owned data (declaration order matters for lifetime!) ----------------

    BasicThermoData thermo_;                        //!< Declared FIRST — outlives model_.
    MechanismDescriptor descriptor_;                //!< Mechanism metadata.
    std::unique_ptr<HMOM<BasicThermoData>> model_;  //!< Declared AFTER thermo_.

    // -- Scratch buffers (pre-allocated, avoid per-cell allocation) ----------

    std::vector<double> y_buffer_;         //!< size = n_species (mass fractions).
    std::array<double, 4> moments_buffer_; //!< M00, M10, M01, N0.

    // -- Metadata ------------------------------------------------------------

    std::vector<int> zone_ids_;    //!< Radiation zone IDs.
    std::string manifest_hash_;    //!< Manifest SHA-256.
    std::string error_message_;    //!< Error description (empty if valid_).
    bool valid_ = false;           //!< True if construction succeeded.

private:
    /** @brief Shared build path used by both public constructors. */
    void build_from_manifest(const Manifest& manifest);
};

} // namespace MOM
