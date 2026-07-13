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

#include "MOM/MomContext.hpp"

#include <exception>

#include "MOM/MomentMethodConcept.hpp"

namespace MOM
{

MomContext::MomContext(const std::string& manifest_path)
{
    try
    {
        const Manifest manifest = ManifestReader::load(manifest_path);
        build_from_manifest(manifest);
    }
    catch (const std::exception& e)
    {
        valid_         = false;
        error_message_ = e.what();
    }
    catch (...)
    {
        valid_         = false;
        error_message_ = "unknown error while constructing MomContext";
    }
}

MomContext::MomContext(const Manifest& manifest)
{
    try
    {
        build_from_manifest(manifest);
    }
    catch (const std::exception& e)
    {
        valid_         = false;
        error_message_ = e.what();
    }
    catch (...)
    {
        valid_         = false;
        error_message_ = "unknown error while constructing MomContext";
    }
}

void MomContext::build_from_manifest(const Manifest& manifest)
{
    // Build the thermo map first: model_ stores a const reference to it.
    thermo_ = MechanismAdapter::build_thermo(manifest);

    // Mechanism metadata (enthalpy / NASA-7 / index map / operating pressure).
    descriptor_ = MechanismAdapter::build_descriptor(manifest);

    // HMOM configuration (validated by the adapter against spec §5.1 guards).
    const HMOM<BasicThermoData>::Config config = MechanismAdapter::build_config(manifest);

    // Construct the model bound to the already-constructed thermo_, then apply
    // the configuration. thermo_ must be alive before this point (lifetime).
    model_ = std::make_unique<HMOM<BasicThermoData>>(thermo_);
    model_->SetupFromConfig(config);

    // Pre-allocate scratch buffers to avoid per-cell heap allocation.
    y_buffer_.assign(thermo_.NumberOfSpecies(), 0.0);
    moments_buffer_.fill(0.0);

    // Metadata carried through from the manifest.
    zone_ids_      = manifest.radiation_zones;
    manifest_hash_ = manifest.manifest_sha256;

    error_message_.clear();
    valid_ = true;
}

void MomContext::compute(double T, double P_Pa, std::span<const double> Y, double mu,
                         std::span<const double> moments) noexcept
{
    if (!valid_)
        return;

    MOM::ComputeCell(*model_, T, P_Pa, Y.data(), mu, moments);
}

std::span<const double> MomContext::sources() const noexcept
{
    return model_->sources();
}

std::span<const double> MomContext::omega_gas() const noexcept
{
    return model_->omega_gas();
}

double MomContext::V0() const noexcept
{
    return model_->V0();
}

double MomContext::particle_density() const noexcept
{
    return model_->particle_density();
}

double MomContext::diffusion_coefficient() const noexcept
{
    return model_->diffusion_coefficient();
}

double MomContext::planck_coefficient(double T, double fv) const noexcept
{
    if (!model_) return 0.0;
    return model_->planck_coefficient(T, fv);
}

double MomContext::volume_fraction() const noexcept
{
    if (!model_) return 0.0;
    return model_->volume_fraction();
}

double MomContext::particle_diameter() const noexcept
{
    if (!model_) return 0.0;
    return model_->particle_diameter();
}

double MomContext::collision_diameter() const noexcept
{
    if (!model_) return 0.0;
    return model_->collision_diameter();
}

std::span<const double> MomContext::initial_moments() const noexcept
{
    if (!model_) return {};
    return model_->initial_moments();
}

unsigned MomContext::n_species() const noexcept
{
    return thermo_.NumberOfSpecies();
}

} // namespace MOM
