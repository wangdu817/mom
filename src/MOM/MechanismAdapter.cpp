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

#include "MOM/MechanismAdapter.hpp"

#include <cstddef>

namespace MOM
{

BasicThermoData MechanismAdapter::build_thermo(const Manifest& manifest)
{
    BasicThermoData thermo;

    const std::size_t n = manifest.species.size();
    thermo.names.reserve(n);
    thermo.mw.reserve(n);
    thermo.nc.reserve(n);
    thermo.nh.reserve(n);
    thermo.no.reserve(n);
    thermo.nn.reserve(n);
    thermo.nti.reserve(n);

    for (const ManifestSpecies& sp : manifest.species)
    {
        thermo.names.push_back(sp.name);
        // Manifest mw_kg_kmol and BasicThermoData.mw are both [kg/kmol]: direct copy.
        thermo.mw.push_back(sp.mw_kg_kmol);
        // Integer atom counts widen to double.
        thermo.nc.push_back(static_cast<double>(sp.nc));
        thermo.nh.push_back(static_cast<double>(sp.nh));
        thermo.no.push_back(static_cast<double>(sp.no));
        thermo.nn.push_back(static_cast<double>(sp.nn));
        // No titanium in this mechanism.
        thermo.nti.push_back(0.0);
    }

    return thermo;
}

HMOM<BasicThermoData>::Config MechanismAdapter::build_config(const Manifest& manifest)
{
    HMOM<BasicThermoData>::Config cfg;

    // ---- Guard rails (spec §5.1) ------------------------------------------
    if (manifest.sticking_model == "golaut")
        throw ManifestParseError(
            "sticking_model 'golaut' is not permitted (spec §5.1)");

    if (manifest.thermophoretic_model != 0)
        throw ManifestParseError(
            "thermophoretic_model must be 0 for v1 (spec §5.1), got " +
            std::to_string(manifest.thermophoretic_model));

    // ---- Activation / PAH setup -------------------------------------------
    // pah_index is validated in-range by ManifestReader::validate.
    cfg.pah_species = manifest.species[static_cast<std::size_t>(manifest.pah_index)].name;

    // ---- Geometry models --------------------------------------------------
    cfg.fractal_diameter_model = manifest.fractal_diameter_model;
    cfg.collision_diameter_model = manifest.collision_diameter_model;

    // ---- Soot/particle properties -----------------------------------------
    cfg.soot_density_kg_m3 = manifest.soot_density_kg_m3;
    cfg.surface_density_per_m2 = manifest.surface_density_per_m2;
    // Manifest carries a double (1.0 / 0.0); Config wants a bool.
    cfg.surface_density_correction = (manifest.surface_density_correction != 0.0);

    // ---- Process model selection ------------------------------------------
    cfg.nucleation_model = manifest.nucleation_model;
    cfg.condensation_model = manifest.condensation_model;
    cfg.surface_growth_model = manifest.surface_growth_model;
    cfg.oxidation_model = manifest.oxidation_model;
    cfg.coagulation_model = manifest.coagulation_model;
    cfg.continuous_coagulation_model = manifest.continuous_coagulation_model;
    cfg.thermophoretic_model = manifest.thermophoretic_model;

    // ---- Sticking coefficient ---------------------------------------------
    cfg.sticking_model = manifest.sticking_model;
    cfg.sticking_coeff_constant = manifest.sticking_coeff_constant;

    // ---- Gas consumption --------------------------------------------------
    cfg.gas_consumption = manifest.gas_consumption;

    // ---- Transport --------------------------------------------------------
    cfg.schmidt_number = manifest.schmidt_number;

    // All remaining Config fields (is_active, simplified_pah_mass,
    // surf_dens_a1/a2/b1/b2, gas_closure_dummy_species, radiative_heat_transfer,
    // planck_coefficient, HACA kinetics, debug_mode) are not represented in the
    // manifest and keep their library defaults.

    return cfg;
}

MechanismDescriptor MechanismAdapter::build_descriptor(const Manifest& manifest)
{
    MechanismDescriptor desc;

    const std::size_t n = manifest.species.size();
    desc.h_f0_j_kg.reserve(n);
    desc.h_f0_j_mol.reserve(n);
    desc.nasa7_low.reserve(n);
    desc.nasa7_high.reserve(n);
    desc.temp_low.reserve(n);
    desc.temp_common.reserve(n);
    desc.temp_high.reserve(n);

    for (const ManifestSpecies& sp : manifest.species)
    {
        desc.h_f0_j_kg.push_back(sp.h_f0_j_kg);
        desc.h_f0_j_mol.push_back(sp.h_f0_j_mol);
        desc.nasa7_low.push_back(sp.nasa7_low);
        desc.nasa7_high.push_back(sp.nasa7_high);
        desc.temp_low.push_back(sp.temp_low);
        desc.temp_common.push_back(sp.temp_common);
        desc.temp_high.push_back(sp.temp_high);
    }

    desc.pah_index = manifest.pah_index;
    desc.bulk_index = manifest.bulk_index;
    desc.c2h2_index = manifest.c2h2_index;
    desc.h2_index = manifest.h2_index;
    desc.o2_index = manifest.o2_index;
    desc.oh_index = manifest.oh_index;
    desc.co_index = manifest.co_index;

    desc.operating_pressure_pa = manifest.operating_pressure_pa;

    return desc;
}

} // namespace MOM
