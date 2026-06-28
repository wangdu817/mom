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

/*-----------------------------------------------------------------------*\
|   MOM Library — Method of Moments for particle population dynamics      |
|   CRECK Modeling Group, Politecnico di Milano                           |
|                                                                         |
|   Master include. In most cases, including only this header is enough.  |
|                                                                         |
|   USAGE (compile-time method selection — zero overhead):                |
|                                                                         |
|     #include "MOM/MOM.hpp"                                              |
|                                                                         |
|     using ParticleModel = MOM::ThreeEquations<MOM::Thermo>;             |
|     // ↑ change this one line to switch to HMOM, BrookesMoss, or TiO2  |
|     static_assert(MOM::MomentMethod<ParticleModel>);                   |
|                                                                         |
|   USAGE (runtime method selection — one indirect branch per call):      |
|                                                                         |
|     auto model = MOM::MakeAnyMomentMethod(thermo, "HMOM");             |
|     MOM::SetState(model, T, P, Y);                                      |
|     MOM::SetMoments(model, moments_span);                               |
|     MOM::Compute(model);                                                |
|     auto src = MOM::GetSources(model);                                  |
|                                                                         |
\*-----------------------------------------------------------------------*/

#pragma once

// ── Core infrastructure ───────────────────────────────────────────────────────
#include "ThermoProxy.hpp"           // ThermoMap concept + Thermo adapter
#include "ProcessFlags.hpp"          // Shared process model enum classes
#include "MomentMethodBase.hpp"      // CRTP base: shared state, Planck, zero sources
#include "MomentMethodConcept.hpp"   // MomentMethod C++20 concept (the contract)

// ── Concrete method implementations ───────────────────────────────────────────
#include "HMOM.hpp"                  // 4-moment Hybrid MOM (Mueller et al. 2009)
#include "BrookesMoss.hpp"           // 2-equation Brookes-Moss soot model
#include "ThreeEquations.hpp"        // 3-equation soot (Franzelli et al. 2019)
#include "TiO2.hpp"                  // 3-equation TiO2 nanoparticle model

// ── Runtime selection ─────────────────────────────────────────────────────────
#include "AnyMomentMethod.hpp"       // std::variant wrapper + dispatch helpers

// ============================================================================
// Compile-time concept satisfaction checks
// ============================================================================
//
// These static_asserts verify that every concrete class fully satisfies the
// MomentMethod concept when instantiated with BasicThermoData. 
// They fire at #include time, providing immediate
// feedback if a refactoring accidentally removes a required method.
//
// The checks use BasicThermoData (the standalone test thermo) so that they
// compile in any environment.
// ============================================================================

static_assert(MOM::MomentMethod<MOM::HMOM<MOM::BasicThermoData>>,
    "[MOM] HMOM does not satisfy the MomentMethod concept. "
    "Check: SetMoments(span), CalculateSourceMoments(), CollisionDiameter(), "
    "SpecificSurface(), DiffusionCoefficient(), PrintSummary().");

static_assert(MOM::MomentMethod<MOM::BrookesMoss<MOM::BasicThermoData>>,
    "[MOM] BrookesMoss does not satisfy the MomentMethod concept.");

static_assert(MOM::MomentMethod<MOM::ThreeEquations<MOM::BasicThermoData>>,
    "[MOM] ThreeEquations does not satisfy the MomentMethod concept.");

static_assert(MOM::MomentMethod<MOM::TiO2<MOM::BasicThermoData>>,
    "[MOM] TiO2 does not satisfy the MomentMethod concept.");
