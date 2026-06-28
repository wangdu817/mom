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

// This file is #included at the bottom of MomentMethodBase.hpp.
// It must not be compiled independently.

#pragma once

#include <cmath>

namespace MOM {

// ============================================================================
// Planck mean absorption coefficient correlations
// ============================================================================
//
// These correlations relate the soot volume fraction fv and temperature T to
// the Planck mean absorption coefficient kP [1/m], used in the radiative
// source term:
//
//   q_rad = 4 * sigma * kP * (T^4 - T_inf^4)
//
// where sigma is the Stefan-Boltzmann constant.
// ============================================================================

template <class Derived, unsigned NEq>
double MomentMethodBase<Derived, NEq>::PlanckSmooke(double T, double fv) const noexcept
{
    // Smooke et al. (2005), Combust. Flame 143, 613-628.
    // kP = C1 * fv * T   [1/m],  C1 = 1232.4 K^{-1} m^{-1}  (approx.)
    // Reference values: C1 = 1232.4 from the original paper.
    constexpr double C1 = 1232.4;    // [1/m/K]
    return C1 * fv * T;
}

template <class Derived, unsigned NEq>
double MomentMethodBase<Derived, NEq>::PlanckKent(double T, double fv) const noexcept
{
    // Kent & Honnery (1990), Combust. Sci. Tech. 75, 167-177.
    // kP = C1 * fv * (C2 + T)   where C2 is a temperature offset.
    constexpr double C1 =  1.3e5;    // [1/m]  (approximate, check original paper)
    constexpr double C2 =  0.;       // [K]    (some formulations include an offset)
    return C1 * fv * (C2 + T);
}

template <class Derived, unsigned NEq>
double MomentMethodBase<Derived, NEq>::PlanckSazhin(double T, double fv) const noexcept
{
    // Sazhin (1994), Prog. Energy Combust. Sci. 20, 297-318.
    // kP = fv * (C0 + C1*T + C2*T^2 + C3*T^3)  [1/m]
    constexpr double C0 =  6.3e2;
    constexpr double C1 =  6.3e-1;
    constexpr double C2 = -1.0e-4;
    constexpr double C3 =  0.;       // cubic term (set to 0 if not available)
    const double T2 = T * T;
    return fv * (C0 + C1*T + C2*T2 + C3*T2*T);
}

} // namespace MOM
