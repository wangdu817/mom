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

// Included at the bottom of AnyMomentMethod.hpp — not compiled standalone.

#pragma once

#include <stdexcept>
#include <string>

namespace MOM {

template <ThermoMap Thermo>
AnyMomentMethod<Thermo> MakeAnyMomentMethod(const Thermo& thermo,
                                             std::string_view label)
{
    if (label == "HMOM"           || label == "hmom")
        return AnyMomentMethod<Thermo>{ std::in_place_type<HMOM<Thermo>>, thermo };

    if (label == "BrookesMoss"    || label == "brookesmoss" || label == "BM")
        return AnyMomentMethod<Thermo>{ std::in_place_type<BrookesMoss<Thermo>>, thermo };

    if (label == "ThreeEquations" || label == "threeequations" || label == "3Eq")
        return AnyMomentMethod<Thermo>{ std::in_place_type<ThreeEquations<Thermo>>, thermo };

    if (label == "TiO2"           || label == "tio2")
        return AnyMomentMethod<Thermo>{ std::in_place_type<TiO2<Thermo>>, thermo };

    throw std::invalid_argument(
        std::string("MOM::MakeAnyMomentMethod: unknown method label '")
        + std::string(label) + "'. "
        "Valid options: HMOM, BrookesMoss, ThreeEquations, TiO2."
    );
}

} // namespace MOM
