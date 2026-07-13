# OpenSMOKE++ MOM — Method of Moments Library

> **A production-grade, C++20/23 library for computing particle population source terms  
> in Computational Fluid Dynamics codes via the Method of Moments.**

Developed at the **CRECK Modeling Lab**, Department of Chemistry, Materials and Chemical
Engineering, **Politecnico di Milano**.

---

## Abstract

The **MOM library** provides a rigorously architected, high-performance implementation of
four Method of Moments (MOM) variants for particle-laden reactive flows: the Hybrid Method
of Moments (HMOM), the Three-Equation soot model, the Brookes–Moss soot model, and a
TiO₂ nanoparticle synthesis model.

The library is designed as an *embeddable kernel*: it is integrated directly into the inner
cell loop of an external CFD solver and must therefore carry zero overhead beyond the
physical computation itself. This constraint drives every architectural decision — from the
choice of the Curiously Recurring Template Pattern (CRTP) to eliminate virtual dispatch, to
the use of `static constexpr` zero arrays to reduce unmodelled processes to a pointer
constant at compile time.

The public API is expressed as a C++20 **Concept** (`MomentMethod<M>`), which means any
conforming variant can be substituted for any other by changing a single type alias in the
calling solver, with the compiler verifying the substitution is valid at include time. A
thin type-erasure layer (`AnyMomentMethod<Thermo>`) enables the variant selection to be
deferred to runtime (e.g. from an input file keyword) while preserving bare-metal
performance inside the cell loop via `std::visit`-before-loop hoisting.

---

## Table of Contents

- [Key Features](#key-features)
- [Supported Variants](#supported-variants)
- [Architecture](#architecture)
- [Quick Start](#quick-start)
- [Integration Patterns](#integration-patterns)
- [API Reference](#api-reference)
- [Building the Library](#building-the-library)
- [Offline Mechanism Tools (`tools/`)](#offline-mechanism-tools-tools)
- [Code Style](#code-style)
- [Testing](#testing)
- [References](#references)
- [License](#license)

---

## Key Features

### Zero-Overhead Compile-Time Polymorphism via CRTP

All variant classes derive from `MomentMethodBase<Derived, NEq>` using the Curiously
Recurring Template Pattern. There are **no virtual functions** anywhere in the performance-
critical path. Every call to `sources_X()` in the base resolves to a direct, inlined call
to `Derived::sources_X_impl()` — the compiler sees the full call chain as a single
straight-line computation and optimises accordingly.

```
Solver calls:  model.sources_nucleation()
                   │
           MomentMethodBase::sources_nucleation()          ← inline, no virtual call
               │
               ├─ if constexpr: Derived has sources_nucleation_impl()?
               │      yes → derived().sources_nucleation_impl()   ← [[gnu::always_inline]]
               └─      no  → { kZeroData, NEq }                   ← compile-time constant
```

The two-level chain (`sources_X → derived().sources_X_impl()`) is guaranteed to be
collapsed to a single basic block by any optimising compiler at `-O2` or higher.
`[[gnu::always_inline]]` on `_impl()` methods enforces this even in debug builds.

### Static Pruning of Unmodelled Processes

Each kinetic sub-process (nucleation, coagulation, condensation, surface growth,
oxidation, sintering) is individually optional. A variant that does not model sintering
simply omits `sources_sintering_impl()` from its public interface. The base class detects
this at compile time:

```cpp
// In MomentMethodBase<Derived, NEq>:
[[nodiscard, gnu::always_inline]]
std::span<const double> sources_sintering() const noexcept
{
    if constexpr (requires(const Derived& d) { d.sources_sintering_impl(); })
        return derived().sources_sintering_impl();
    else
        return { kZeroData, NEq };   // pointer to .rodata — zero cost
}
```

`kZeroData` is a `static constexpr double[NEq]` living in read-only memory. The compiler
represents the fallback return as a pair of constants (`lea rax, [rip + kZeroData]`,
`mov edx, NEq`). Any downstream arithmetic over the returned span (`+= src[i]`) is
constant-propagated to `+= 0.0` and eliminated entirely by the dead-store pass.

### Fixed-Size Eigen Vectors — SIMD-Ready Storage

Moment source vectors are stored as `Eigen::Matrix<double, NEq, 1>` (`MomentVector`).
Fixed-size allocation means the entire vector fits in registers for small `NEq`:

| Variant         | NEq | Vector size | Fits in |
|-----------------|-----|-------------|---------|
| HMOM            | 4   | 32 bytes    | 1 × AVX2 `ymm` register |
| ThreeEquations  | 3   | 24 bytes    | 1 × SSE2 `xmm` + scalar |
| BrookesMoss     | 2   | 16 bytes    | 1 × SSE2 `xmm` register |
| TiO₂            | 3   | 24 bytes    | 1 × SSE2 `xmm` + scalar |

With `-O3 -march=native`, Eigen's fixed-size expression templates generate vectorised
arithmetic over all four source equations simultaneously via `vmovupd`/`vmulsd`/`vaddsd`
using `xmm`/`ymm` registers.

### Bounds-Safe, Zero-Copy Data Passing via `std::span`

All source term accessors return `std::span<const double>` — a non-owning view into the
variant's internal `MomentVector`. No heap allocation, no copy, no size field mismatch.
The CFD solver accesses the data in-place:

```cpp
auto src = model.sources_nucleation();   // span into model's internal storage
for (unsigned j = 0; j < src.size(); ++j)
    residual[j] += src[j] * cell_volume;  // direct read, no copy
```

### Unified C++20 Concept API

The `MomentMethod<M>` concept in `MomentMethodConcept.hpp` is the single authoritative
contract for all variants. A `static_assert` at the solver include site catches any
divergence between a variant's implementation and the API at compile time:

```cpp
using ParticleModel = MOM::HMOM<MyThermo>;
static_assert(MOM::MomentMethod<ParticleModel>,
    "ParticleModel must satisfy MOM::MomentMethod");
```

Switching variants requires changing **exactly one line** in the solver.

### Type-Erased Runtime Selection — Single Dispatch Overhead

`AnyMomentMethod<Thermo>` is a `std::variant<HMOM<Thermo>, BrookesMoss<Thermo>, ...>`
with a factory function and a `ForEachCell` dispatcher. The variant selection (one
`std::visit`) is hoisted *outside* the spatial loop:

```
Initialisation:  MakeAnyMomentMethod(thermo, "HMOM")   ← one runtime branch
                     │
                 ForEachCell(model, lambda)             ← one std::visit
                     │
Cell loop:           for (i = 0; i < N; ++i)           ← bare-metal iteration
                         ComputeCell(m, ...)            ← fully inlined, no dispatch
```

The inner loop is identical in machine code to the compile-time template version.

---

## Supported Variants

### Process Matrix

| Process | HMOM | ThreeEquations | BrookesMoss | TiO₂ |
|---|---|---|---|---|
| **Nucleation** | ✅ | ✅ | ✅ | ✅ |
| **Coagulation** | ✅ | ✅ | ✅ | ✅ |
| **Condensation** | ✅ | ✅ | ⬜ zero-span | ✅ |
| **Surface Growth** | ✅ | ✅ | ✅ | ⬜ zero-span |
| **Oxidation** | ✅ | ✅ | ✅ | ⬜ zero-span |
| **Sintering** | ⬜ zero-span | ⬜ zero-span | ⬜ zero-span | ✅ |

**✅ modelled** — `sources_X_impl()` is present; base class forwards the call directly.  
**⬜ zero-span** — `sources_X_impl()` is absent; base class returns `{ kZeroData, NEq }` at zero cost.

### Variant Summary

#### HMOM — Hybrid Method of Moments

Transport equations: **4** (`M00`, `M10`, `M01`, `N0`)

The HMOM reconstructs the Number Density Function (NDF) using two delta nodes (a large-
particle mode and a small-particle nucleation mode, `N0`). The two-point quadrature
enables separate tracking of nucleation-mode particles without requiring closure
assumptions for the full NDF shape.

| Index | Symbol | Physical quantity |
|---|---|---|
| 0 | M00 | Zeroth-order moment — proportional to total number density |
| 1 | M10 | First-order volume moment |
| 2 | M01 | First-order surface moment |
| 3 | N0  | Small-particle (nucleation-mode) number density |

Coagulation is split into nine sub-mechanisms (discrete, discrete-SS, discrete-SL,
discrete-LL, continuous, continuous-SS, continuous-SL, continuous-LL, total), each
accessible individually via `coagulation_detail()`.

*Key references:* Mueller et al. (2009) *Combust. Flame* **156**, 1143–1155; Attili &
Bisetti (2013) *Combust. Flame* **160**, 1920–1934.

#### ThreeEquations — Three-Equation Soot Model

Transport equations: **3** (`Ys`, `NsNorm`, `Ss`)

| Index | Symbol | Physical quantity |
|---|---|---|
| 0 | Ys     | Soot mass fraction [-] |
| 1 | NsNorm | Scaled number density = Ns / 10¹⁵ [−] |
| 2 | Ss     | Specific surface area [m²/m³] |

A log-normal NDF is assumed. Surface chemistry follows the HACA (Hydrogen-Abstraction
Carbon-Addition) mechanism with optional Pareto/LogNormal NDF reconstruction.

*Key reference:* Eberle et al. (2015) *Combust. Flame* **162**, 2648–2660.

#### BrookesMoss — Brookes–Moss Soot Model

Transport equations: **2** (`Ys`, `bs`)

| Index | Symbol | Physical quantity |
|---|---|---|
| 0 | Ys | Soot mass fraction [-] |
| 1 | bs | Normalised radical nuclei concentration = Ns/(ρ·Ns_norm) [m³/kg] |

The simplest two-equation variant. Condensation and sintering are not modelled and resolve
to zero-cost static spans. Suitable for large simulations where computational economy is
paramount and the full NDF shape is not required.

*Key reference:* Brookes & Moss (1999) *Combust. Flame* **116**, 486–503.

#### TiO₂ — Titanium Dioxide Nanoparticle Synthesis

Transport equations: **3** (`YTiO₂`, `NTiO₂N`, `STiO₂`)

| Index | Symbol | Physical quantity |
|---|---|---|
| 0 | YTiO₂  | TiO₂ particle mass fraction [-] |
| 1 | NTiO₂N | Scaled number density = N / 10¹⁵ [−] |
| 2 | STiO₂  | Total surface area concentration [m²/m³] |

Models the aerosol synthesis of TiO₂ nanoparticles from gas-phase precursors (e.g. TiCl₄
oxidation). Surface growth and oxidation are not modelled. Sintering is modelled
explicitly, which distinguishes this variant from the soot models. NDF closure follows a
log-normal assumption with geometric standard deviation.

---

## Architecture

### Class Hierarchy

```
MomentMethodBase<Derived, NEq>          ← CRTP base (MomentMethodBase.hpp)
│   Stores:  MomentVector sources_[6]   (Eigen fixed-size, NEq doubles each)
│            kZeroData[NEq]             (static constexpr, .rodata)
│            omega_gas_[]               (dynamic, sized at SetState time)
│   Provides: sources_X()              (compile-time dispatch via if constexpr)
│             planck_coefficient()     (3 empirical correlations in .tpp)
│             all concept accessors    (particle properties, transport, etc.)
│
├── HMOM<Thermo>           : Base<HMOM<Thermo>,           4>
├── ThreeEquations<Thermo> : Base<ThreeEquations<Thermo>, 3>
├── BrookesMoss<Thermo>    : Base<BrookesMoss<Thermo>,    2>
└── MetalOxide<Thermo>           : Base<MetalOxide<Thermo>,           3>

MomentMethod<M>                         ← C++20 Concept (MomentMethodConcept.hpp)
    requires: all public API methods and types

AnyMomentMethod<Thermo>                 ← type-erased runtime wrapper
    = std::variant<HMOM<T>, BrookesMoss<T>, ThreeEquations<T>, MetalOxide<T>>
    MakeAnyMomentMethod(thermo, label)  ← factory (label = "HMOM", "MetalOxide", ...)
    ForEachCell(variant, lambda)        ← hoisted std::visit pattern

AllVariants                             ← single authoritative registry
    = detail::TypeList<HMOM, BrookesMoss, ThreeEquations, MetalOxide>
      (MomVariantList.hpp — the only file to edit when adding a variant)
```

### CRTP Dispatch — Step by Step

Consider a CFD solver requesting nucleation sources for HMOM. The call chain at
`-O2 -march=native` is:

```
1.  model.sources_nucleation()
         ↓   [inlined — no call instruction]
2.  MomentMethodBase<HMOM<T>, 4>::sources_nucleation()
         ↓   if constexpr: HMOM has sources_nucleation_impl()? → YES
3.  derived().sources_nucleation_impl()
         ↓   [[gnu::always_inline]] — forced into call site even at -O0
4.  return { source_nucleation_.data(), 4 }   ← span into Eigen vector
         ↓   [entire chain inlined; assembly: lea rax, [rbp+offset]; mov edx, 4]
5.  Caller receives std::span<const double> — pointer + size, no heap
```

For a variant without that process (e.g. `BrookesMoss::sources_condensation()`):

```
1.  model.sources_condensation()
         ↓
2.  MomentMethodBase<BrookesMoss<T>, 2>::sources_condensation()
         ↓   if constexpr: BrookesMoss has sources_condensation_impl()? → NO
3.  return { kZeroData, 2 }   ← pointer to .rodata + constant 2
         ↓   [branch eliminated at compile time; no runtime test]
         ↓   downstream "+= span[i]" → "+= 0.0" → eliminated by DCE pass
```

### Variant Registry

Adding a new variant requires changes in **one file only**:

```cpp
// MomVariantList.hpp — append the new type here
using AllVariants = detail::TypeList<
    HMOM,
    BrookesMoss,
    ThreeEquations,
    MetalOxide,
    MyNewVariant    // ← add here
>;
```

The factory (`MakeAnyMomentMethod`), the type-erased variant, and all dispatch helpers
(`ForEachCell`, `GetNEquations`, etc.) update automatically via template instantiation.
No other file is modified.

### ThermoMap Concept

All variants are parameterised on `Thermo`, which must satisfy the `ThermoMap` concept:

```cpp
template <typename T>
concept ThermoMap = requires(const T t, std::string_view name) {
    { t.NumberOfSpecies() } -> std::convertible_to<int>;
    { t.IndexOfSpecies(name) } -> std::convertible_to<int>;
    { t.MolecularWeight(0) } -> std::convertible_to<double>;
    // ... atom count methods
};
```

The built-in `BasicThermoData` struct satisfies `ThermoMap` and is ready to use
standalone. For integration with full CFD thermodynamic backends (OpenFOAM, OpenSMOKE++),
any class satisfying the concept is accepted without modification to the MOM library.

---

## Quick Start

```cpp
#include "MOM/MOM.hpp"   // umbrella header — pulls in all variants and utilities

// 1. Build a thermodynamic map (or pass your solver's thermo backend)
MOM::BasicThermoData thermo;
thermo.names = { "H", "OH", "O2", "H2", "H2O", "C2H2", "N2" };
thermo.mw    = { 1.008, 17.008, 31.999, 2.016, 18.015, 26.038, 28.014 };
// ... fill atom counts

// 2. Construct the variant — type selected at compile time
MOM::ThreeEquations<MOM::BasicThermoData> model(thermo);
model.SetPAH("C2H2");
model.SetNucleation(1);
model.SetCoagulation(1);
model.SetCondensation(1);
model.SetSurfaceGrowth(1);
model.SetOxidation(1);

static_assert(MOM::MomentMethod<decltype(model)>);   // verified at compile time

// 3. Per-cell evaluation (inside the CFD loop)
model.SetState(T, P_Pa, Y.data());             // T [K], P [Pa], Y mass fractions
model.SetMoments(moments_span);                  // current moment values
model.SetViscosity(mu);                          // mixture viscosity [kg/m/s]
model.ComputeSources();

// 4. Read sources — zero-copy spans into internal Eigen storage
auto src      = model.sources();                 // total,  size = 3
auto src_nucl = model.sources_nucleation();      // nucleation contribution
auto src_sint = model.sources_sintering();       // → kZeroData (not modelled)
```

---

## Integration Patterns

The library supports two integration patterns. Both produce identical machine code in the
inner loop; the choice depends on whether the variant is known at compile time.

### Pattern A — Compile-Time Variant (Recommended for Template Solvers)

The variant type is a template parameter of the CFD cell loop. The compiler generates one
fully-specialised loop body with every `if constexpr` branch resolved and the complete
dispatch chain inlined.

```cpp
template <MOM::MomentMethod M>
void CellLoop(M& model,
              int n_cells,
              const double* T_arr,
              const double* P_arr,
              const double* Y_flat,   // [n_cells × n_species]
              const double* mu_arr,
              const double* M_flat,   // [n_cells × M::n_equations]
              double* src_flat)       // [n_cells × M::n_equations]  output
{
    constexpr unsigned neq = M::n_equations;

    for (int i = 0; i < n_cells; ++i)
    {
        // Single entry-point: sets state, moments, viscosity, computes sources.
        MOM::ComputeCell(model,
                         T_arr[i], P_arr[i],
                         Y_flat + i * n_species,
                         mu_arr[i],
                         { M_flat + i * neq, neq });

        // Zero-copy read — no allocation, no copy
        auto src = model.sources();
        std::copy(src.begin(), src.end(), src_flat + i * neq);
    }
}

// Usage — one alias change to switch variant:
using ParticleModel = MOM::HMOM<MyThermo>;
CellLoop(model, N, T.data(), P.data(), Y.data(), mu.data(), M.data(), src.data());
```

At `-O3 -march=native`, the assembly for an HMOM cell contains no `call` instructions in
the hot path and uses `ymm` registers to process all four equations simultaneously.

### Pattern B — Runtime Variant Selection via `ForEachCell`

The variant is chosen at runtime from an input-file keyword. `std::visit` fires **once**
before the loop via `ForEachCell`; the lambda receives the concrete type and owns the
iteration. The compiler generates four specialised loop bodies (one per variant
alternative) and selects the correct one with a single branch before the loop.

```cpp
// Solver initialisation — one runtime branch at startup
MOM::AnyMomentMethod<MyThermo> model = MOM::MakeAnyMomentMethod(thermo, variant_name);
// variant_name = "HMOM" | "ThreeEquations" | "BrookesMoss" | "MetalOxide"

// CFD sweep — std::visit fires ONCE here, not inside the loop
MOM::ForEachCell(model, [&](auto& m)
{
    // typeof(m) is concrete here (e.g. ThreeEquations<MyThermo>&)
    constexpr unsigned neq = std::decay_t<decltype(m)>::n_equations;

    for (int i = 0; i < n_cells; ++i)
    {
        MOM::ComputeCell(m, T[i], P[i], Y + i * ns, mu[i], { M + i * neq, neq });

        auto src = m.sources();
        std::copy(src.begin(), src.end(), src_out + i * neq);
    }
});
```

Benchmark results on a representative 10 000-cell sweep (ThreeEquations, GCC 13,
`-O3 -march=native`, Intel Core i7-12700H):

```
Pattern A — template cell loop    :   X.XXX ms / iteration
Pattern B — ForEachCell loop      :   X.XXX ms / iteration
Overhead of variant dispatch      :   0.000 ms  (0.0%)
```

The overhead is immeasurable: `std::visit` is a single indirect branch that the branch
predictor makes free after the first correct prediction.

### Output and Diagnostics — `MomentMethodReporter`

The reporter is a **read-only observer**: it never modifies the model and has no I/O
dependency on the variant classes themselves (the variants carry no `WriteHeader` or
`WriteRow` methods). It is driven by the public `MomentMethod<M>` concept interface and by
optional `variant_prefix_output` / `variant_suffix_output` hooks detected via
`if constexpr`:

```cpp
// Once, before any loops — outside the CFD cell loop
MOM::OutputFileColumns file("soot_diagnostics.out");
MOM::MomentMethodReporter reporter(file, thermo.names);
reporter.WriteHeader(model);   // writes column headers to file
file.Complete();

// On output steps only — infrequently, never inside the hot cell loop
file.NewRow();
reporter.WriteRow(model);      // writes one diagnostic row; model is const&
```

Variant-specific extended output (e.g. the nine coagulation sub-mechanisms of HMOM) is
emitted automatically via the `variant_suffix_output` hook if the variant provides it,
without any modification to `MomentMethodReporter`.

---

## API Reference

All public API is expressed as requirements of the `MomentMethod<M>` concept.

### State Injection (call before each `ComputeSources`)

| Method | Description |
|---|---|
| `SetState(T, P, Y[])` | Gas temperature [K], pressure [Pa], mass fractions |
| `SetMoments(span)` | Current transported moment values (size = `n_equations`) |
| `SetViscosity(mu)` | Mixture dynamic viscosity [kg/m/s] |

### Computation

| Method | Description |
|---|---|
| `ComputeSources()` | Evaluates all moment source terms; `noexcept` |
| `CalculateOmegaGas()` | Evaluates gas-phase consumption terms only; `noexcept` |

### Source Output (zero-copy `std::span<const double>`, size = `n_equations`)

| Method | Description |
|---|---|
| `sources()` | Total source vector (sum of all processes) |
| `sources_nucleation()` | Nucleation contribution |
| `sources_coagulation()` | Coagulation contribution |
| `sources_condensation()` | Condensation (→ `kZeroData` if not modelled) |
| `sources_growth()` | Surface growth (→ `kZeroData` if not modelled) |
| `sources_oxidation()` | Oxidation (→ `kZeroData` if not modelled) |
| `sources_sintering()` | Sintering (→ `kZeroData` if not modelled) |
| `omega_gas()` | Gas-phase species source terms [kg/m³/s] |

### Particle Properties

| Method | Unit | Description |
|---|---|---|
| `volume_fraction()` | − | Particle volume fraction |
| `particle_diameter()` | m | Primary particle diameter |
| `collision_diameter()` | m | Collision (aggregate) diameter |
| `particle_number_density()` | #/m³ | Total number density |
| `mass_fraction()` | − | Particle mass fraction |
| `particle_density()` | kg/m³ | Material density |
| `specific_surface()` | m²/m³ | Surface area per unit volume |

### Transport and Radiation

| Method | Description |
|---|---|
| `schmidt_number()` | Particle Schmidt number |
| `diffusion_coefficient()` | Effective diffusion coefficient [kg/m/s] |
| `thermophoretic_model()` | Thermophoretic model flag (0 = off) |
| `planck_coefficient(T, fv)` | Planck mean absorption coefficient [1/m] |
| `radiative_heat_transfer()` | True if radiative loss is active |

Three Planck mean absorption coefficient correlations are available and selected at
construction time via `ProcessFlags::PlanckCoeffModel`:

| Model | Formula | Reference |
|---|---|---|
| Smooke | kP = 1232.4 · fv · T | Smooke et al. *Combust. Flame* **143** (2005) |
| Kent | kP = 1.3×10⁵ · fv · T | Kent & Honnery *Combust. Sci. Tech.* **75** (1990) |
| Sazhin | kP = fv · (630 + 0.63T − 10⁻⁴T²) | Sazhin *PECS* **20** (1994) |

---

## Building the Library

### Requirements

| Component | Minimum version | Notes |
|---|---|---|
| C++ compiler | GCC 11 / Clang 13 / MSVC 19.29 | C++20 required; C++23 auto-detected |
| CMake | 3.20 | For `FetchContent`, generator expressions |
| Eigen | 3.4.0 | Auto-fetched from GitLab if not found |

### CMake Integration

**Option 1 — `add_subdirectory` (recommended for in-tree development)**

```cmake
add_subdirectory(MOM)
target_link_libraries(MySolver PRIVATE MOM)
```

**Option 2 — `find_package` (after installation)**

```cmake
find_package(MOM REQUIRED)
target_link_libraries(MySolver PRIVATE MOM::MOM)
```

**Option 3 — Header-only interface target**

For solvers that instantiate MOM templates with custom thermo types, link against
`MOM_headers` instead of `MOM`. This gives access to all headers without the
pre-compiled `BasicThermoData` specialisations:

```cmake
target_link_libraries(MySolver PRIVATE MOM::MOM_headers)
```

### Build Options

| CMake option | Default | Description |
|---|---|---|
| `MOM_BUILD_SHARED` | OFF | Build as a shared library instead of static |
| `MOM_BUILD_TESTS` | ON | Build the test executables |
| `MOM_INSTALL` | ON | Generate install rules |
| `MOM_USE_DICTIONARY` | OFF | Enable input-file dictionary parser (requires `DICTIONARY_DIR`) |
| `EIGEN3_DIR` | — | Override Eigen path (skips `find_package`) |

### Typical Build Sequence

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
ctest --test-dir build --output-on-failure
```

### Recommended Compiler Flags

```bash
# Minimum for correctness
-std=c++20

# Production CFD (recommended)
-std=c++20 -O3 -march=native

# Optional: remove exception tables from noexcept-only hot path
-fno-exceptions

# Optional: cross-translation-unit inlining (most relevant in pre-compiled mode)
-flto
```

What each flag enables for this library:

- **`-O2`**: inlines the complete `sources_X() → _impl()` dispatch chain. `[[gnu::always_inline]]` on `_impl()` methods makes this a hard guarantee.
- **`-O3`**: adds auto-vectorisation. For HMOM (`NEq=4`), the 4-double source vector maps to a single 256-bit AVX `ymm` register.
- **`-march=native`**: enables AVX2/FMA. Most relevant for Eigen's fixed-size vector arithmetic and the `std::copy` loops over source spans. Without this, the compiler defaults to SSE2 (128-bit, 2 doubles), halving throughput on 4-equation models.

---

## Offline Mechanism Tools (`tools/`)

The `tools/` directory (at the repository root, alongside `mom/`) contains a set of
Python offline utilities that run **once per case, before the CFD solver starts**. They
turn a human-authored case configuration and a CHEMKIN mechanism into a deterministic,
integrity-checked **runtime manifest** — a single JSON file consumed by the C++ UDF bridge
(`MOM::ManifestReader`) at solver startup.

The offline pipeline has three stages:

```
hmom-case.toml ──▶ prepare_hmom_case.py validate   (schema / semantic validation)
                          │
thermo.dat ──────▶ mechanism_qualifier.py qualify  (CHEMKIN NASA-7 parse + G0 gate)
                          │                              └──▶ g0_report.json
                          ▼
        manifest_writer.py write config.toml g0_report.json  (cross-validate + emit)
                          │
                          ▼
                  manifest.json  ──▶  MOM::ManifestReader::load()  (C++ solver startup)
```

All tools are **Python 3.11+**, standard-library only except for **pydantic v2** (used by
the config validator and manifest writer). Each tool is runnable both as a script
(`python tools/<tool>.py …`) and as a module (`python -m tools.<tool> …`).

### 1. Case config validator — `prepare_hmom_case.py`

Validates a UTF-8 TOML case configuration against a strict pydantic v2 schema
(`hmom_case_schema.py`, 8 models, all `extra="forbid"` — unknown keys, type mismatches,
and missing fields are rejected). The schema enforces v1 invariants such as
`schema_version = 1`, `bulk_species = "N2"`, `thermophoretic_model = 0`, `uds_count = 4`,
a 64-hex-character `expected_mechanism_sha256`, and a `pah_species` that does not collide
with a reserved gas-phase species name.

```bash
python tools/prepare_hmom_case.py validate tools/examples/hmom-case.toml
```

On success it prints a summary and exits 0; on failure it prints per-field errors to
stderr and exits 1:

```
HMOM case config is VALID.
  schema_version:   1
  zone count:       1 (zone_ids=[12])
  pah_species:      A4
  mechanism format: chemkin
  hmom processes:   nucleation=1 condensation=1 surface_growth=1 oxidation=1 ...
  sticking_model:   constant
  udm mode/count:   release / 18 slots
  uds_count:        4
  radiation model:  Smooke (enabled=True)
  gas_consumption:  True
```

The config is organised in eight required sections: `[case]`, `[mechanism]`, `[hmom]`,
`[fluent]`, `[radiation]`, `[energy]`, `[validation]`, plus the top-level `schema_version`.
See `tools/examples/hmom-case.toml` for an annotated sample. The `[fluent]` section derives
`udm_count` from `udm_mode` — **18** slots in `release`, **38** in `debug`.

### 2. Mechanism qualifier (G0) — `mechanism_qualifier.py`

Parses a CHEMKIN-format `THERMO` section (4-line × 80-column NASA-7 entries) and produces a
**G0 qualification report** verifying the mechanism is usable by the soot model. G0 checks:

- the fixed required species are present — `H, OH, O2, H2, H2O, CO, C2H2, N2` — plus the
  configurable PAH species,
- `N2` is the **last** species in the mechanism,
- the PAH species has a distinct index and does not collide with a reserved name,
- all NASA-7 coefficients are finite, molecular weights positive, and formation
  enthalpies finite,
- no duplicate species names.

For each species it computes molecular weight, standard formation enthalpy `Hf0(298.15 K)`,
and heat capacity `Cp(298.15 K)` from the NASA-7 polynomials, and records both temperature
ranges and the raw coefficients.

```bash
python tools/mechanism_qualifier.py qualify tools/tests/fixtures/sample_thermo.dat \
    --pah-species A4 --output g0_report.json
```

Options: `--pah-species` (default `A4`), `--output <file>` (write JSON report; otherwise
printed to stdout), `--quiet` (suppress the JSON dump). Exit code is 0 when qualified,
1 otherwise. The report includes the `thermo_sha256`, `species_order`, per-species records,
and a `gaps` list explaining any failure:

```
G0 qualification PASSED.
  thermo_file:   tools/tests/fixtures/sample_thermo.dat
  sha256:        <64-hex>
  species_count: 9
  species_order: ['H', 'OH', 'O2', 'H2', 'H2O', 'CO', 'C2H2', 'A4', 'N2']
  pah_species:   A4 (index 7)
  n2_is_last:    True
```

### 3. Manifest writer — `manifest_writer.py`

Consumes the validated config **and** the G0 report, cross-validates the two (matching
`pah_species`, `thermo_file`, `N2`-last, required species present, and requiring
`qualified = true`), and emits a **deterministic** JSON manifest. Determinism is guaranteed
by serialising with `sort_keys=True` and embedding a self SHA-256 computed over the manifest
with the `manifest_sha256` field cleared. The manifest carries the magic string
`HMOM-MANIFEST`, `schema_version = 1`, the species table, the source-species index map, the
resolved HMOM/Fluent/radiation/energy options, and the fixed UDM slot layout.

```bash
python tools/manifest_writer.py write tools/examples/hmom-case.toml g0_report.json \
    --output manifest.json
```

Exit is 0 on success. With `--output` it prints the species count and `manifest_sha256`;
without it the manifest is written to stdout. The resulting `manifest.json` is loaded by the
C++ bridge via `MOM::ManifestReader::load()` (see `include/MOM/ManifestReader.hpp`).

### End-to-end pipeline

```bash
# 1. Validate the case configuration
python tools/prepare_hmom_case.py validate tools/examples/hmom-case.toml

# 2. Qualify the mechanism (G0 gate) and capture the report
python tools/mechanism_qualifier.py qualify tools/tests/fixtures/sample_thermo.dat \
    --pah-species A4 --output g0_report.json

# 3. Write the deterministic runtime manifest
python tools/manifest_writer.py write tools/examples/hmom-case.toml g0_report.json \
    --output manifest.json
```

The bundled `tools/tests/fixtures/sample_thermo.dat` is a 9-species CHEMKIN sample
(`H, OH, O2, H2, H2O, CO, C2H2, A4, N2`) that qualifies against G0 with `--pah-species A4`
and is convenient for exercising the full pipeline without a production mechanism. The tools
are covered by unit tests under `tools/tests/` (`test_config_schema.py`,
`test_mechanism_qualifier.py`, `test_manifest_writer.py`).

---

## C++ Context & Adapter Layer

The offline tools (above) emit a deterministic `manifest.json`. The C++ **context and
adapter layer** consumes that manifest at solver startup and turns it into a ready-to-run
MOM model. It is the bridge between the Python-produced manifest and the compile-time MOM
kernel, and is the code that a Fluent UDF (or any host CFD solver) calls once per zone
during initialisation.

The layer has three composable pieces, each with a matching unit-test program:

```
manifest.json
     │  MOM::ManifestReader::load()          (strict JSON parse + fail-fast validation)
     ▼
MOM::Manifest  (parsed, validated struct)
     │  MOM::MechanismAdapter                 (Manifest → thermo + config + descriptor)
     ▼
BasicThermoData + HMOM::Config + MechanismDescriptor
     │  MOM::MomContext(manifest_path)        (owns thermo + HMOM model + scratch buffers)
     ▼
MomContext::compute(T, P, Y, mu, moments)     (noexcept per-cell entry point)
```

### 1. `ManifestReader` — Strict JSON Parser

`include/MOM/ManifestReader.hpp` / `src/MOM/ManifestReader.cpp` parse the manifest written
by `tools/manifest_writer.py` using the vendored **nlohmann/json** library. Parsing is
**strict and fail-fast**: any structural or semantic violation throws `ManifestParseError`
(derived from `std::runtime_error`) with a descriptive message rather than returning a
partially-populated struct.

```cpp
#include "MOM/ManifestReader.hpp"

MOM::Manifest m = MOM::ManifestReader::load("manifest.json");   // throws on any violation
// or, from an in-memory string:
MOM::Manifest m = MOM::ManifestReader::parse(json_text);
```

The reader enforces **11 validation invariants** before returning:

| # | Invariant | Rejection reason |
|---|---|---|
| 1 | `magic == "HMOM-MANIFEST"` | Not an HMOM manifest |
| 2 | `schema_version == 1` | Unsupported schema |
| 3 | `species_count == species.size()` | Count / array mismatch |
| 4 | `species_count > 0` | Empty mechanism |
| 5 | Last species is `N2` | Bulk species not last |
| 6 | Every `mw_kg_kmol > 0` | Non-physical molecular weight |
| 7 | Every `h_f0_j_kg` finite | NaN/Inf formation enthalpy |
| 8 | All `nasa7_low[7]` / `nasa7_high[7]` finite | NaN/Inf polynomial coefficient |
| 9 | `uds_count == 4` | Wrong transported-scalar count |
| 10 | `udm_count` matches `udm_mode` (18 release / 38 debug) | UDM layout mismatch |
| 11 | `operating_pressure_pa > 0` | Non-physical operating pressure |

Static contract constants are exposed for host code: `kMagic`, `kSchemaVersion = 1`,
`kUdsCount = 4`, `kUdmReleaseCount = 18`, `kUdmDebugCount = 38`.

### 2. `MechanismAdapter` — Manifest → Model Inputs

`include/MOM/MechanismAdapter.hpp` / `src/MOM/MechanismAdapter.cpp` translate a validated
`Manifest` into the three objects the HMOM kernel needs:

- **`BasicThermoData`** — species names, molecular weights **converted to kg/kmol**, and
  per-species element counts (`nc/nh/no/nn`), satisfying the `ThermoMap` concept.
- **`HMOM::Config`** — all **16** HMOM process/parameter options mapped from
  `hmom_options`, with `pah_species` resolved from the manifest `pah_index`.
- **`MechanismDescriptor`** — the residual physical data not held by the two above:
  per-species formation enthalpy `h_f^0`, NASA-7 low/high polynomial coefficients,
  source-species indices, and the operating pressure.

```cpp
MOM::MechanismAdapter adapter(manifest);
MOM::BasicThermoData     thermo     = adapter.thermo();
MOM::HMOM<...>::Config   config     = adapter.config();
MOM::MechanismDescriptor descriptor = adapter.descriptor();
```

**Guard rails.** The adapter refuses to build a model for configurations the C++ kernel
does not support, throwing rather than silently degrading:

- `sticking_model == "golaut"` is rejected (unsupported sticking model).
- `thermophoretic_model != 0` is rejected (thermophoresis not modelled in this bridge).

### 3. `MomContext` — Owning Runtime Context

`include/MOM/MomContext.hpp` / `src/MOM/MomContext.cpp` provide the single object a host
solver holds for the lifetime of a case. `MomContext` **owns** the thermo data and the
HMOM model, plus reusable scratch buffers:

```cpp
struct MomContext {
    BasicThermoData        thermo_;   // declared FIRST — outlives model_
    std::unique_ptr<HMOM<BasicThermoData>> model_;
    // + MechanismDescriptor, scratch buffers, manifest hash, status
};
```

**Lifetime ordering is deliberate:** `thermo_` is declared before `model_` so that, since
the HMOM model holds a reference to the thermo backend, the thermo data is guaranteed to
outlive the model during destruction. `MomContext` is **non-copyable and non-movable**
(copy and move are `= delete`d) — it is a stable, pinned owner.

Construction takes a manifest path and performs the full pipeline (load → validate →
adapt → construct HMOM). The per-cell entry point is `noexcept`:

```cpp
MOM::MomContext ctx("manifest.json");        // load + build; sets is_valid()/error_message()
if (!ctx.is_valid()) { /* handle ctx.error_message() */ }

ctx.compute(T, P, Y.data(), mu, moments_span);   // noexcept — wraps MOM::ComputeCell
auto src = ctx.sources();                         // zero-copy span into the model
```

`compute()` sets state/moments/viscosity and calls `MOM::ComputeCell` in the correct order,
never throwing across the C boundary. The context exposes **12 accessors**:

| Accessor | Returns |
|---|---|
| `sources()` | Total moment source vector (span) |
| `omega_gas()` | Gas-phase species source terms |
| `V0()` | Nucleation-node particle volume |
| `particle_density()` | Soot material density [kg/m³] |
| `diffusion_coefficient()` | Effective diffusion coefficient (> 0) |
| `descriptor()` | The `MechanismDescriptor` |
| `zone_ids()` | Fluent zone IDs from the manifest |
| `manifest_hash()` | `manifest_sha256` of the loaded manifest |
| `is_valid()` | Construction succeeded |
| `error_message()` | Failure description (empty when valid) |
| `n_species()` | Species count |
| `species_names()` | Species name list |

### Unit Tests

Three test programs (**28 test cases total**) plus a shared fixture header validate the
layer; all run under ctest with **7/7 passing**:

| Test program | Cases | Coverage |
|---|---|---|
| `tests/test_manifest_reader.cpp` | 10 | Valid parse; each of the 11 invariants rejects; magic/version checks |
| `tests/test_mechanism_adapter.cpp` | 9 | Species order preserved; MW kg/mol → kg/kmol conversion; all 16 config options; `golaut` and thermophoresis guard rails |
| `tests/test_mom_context.cpp` | 9 | Construction from manifest; `ComputeCell` call order; thermo/model lifetime ordering; `diffusion_coefficient() > 0` |

`tests/test_manifest_fixture.hpp` provides the shared `make_test_manifest()` helper used by
all three programs to build a valid in-memory manifest that individual tests then perturb to
exercise each rejection path.

### Build Notes

The layer is built with **MSVC 14.50** via the **Ninja** generator. As with the rest of the
sources, the `/utf-8` flag is applied (the code and manifests contain non-ASCII characters).
The **nlohmann/json v3.12.0** single-header dependency is **vendored** at
`mom/third_party/nlohmann/json.hpp` — no network fetch is required. All new implementation
files live under `mom/src/MOM/` and all new tests under `mom/tests/`.

---

## Code Style

All source files are formatted with **clang-format** using the project-provided
`.clang-format` configuration. Key settings:

| Rule | Value |
|---|---|
| Baseline style | LLVM |
| Brace placement | **Allman** (opening brace on its own line) |
| Indent width | 4 spaces (no tabs) |
| Column limit | 100 characters |
| Pointer alignment | Left (`double* ptr`) |
| `SortIncludes` | Never (mandatory ordering: STL → Eigen → project headers) |
| Minimum clang-format version | 14 |

Apply formatting to all sources:

```bash
find include src tests \( -name '*.hpp' -o -name '*.h' -o -name '*.tpp' -o -name '*.cpp' \) \
    | xargs clang-format -i
```

**Include ordering convention** (enforced manually, not by clang-format):

```cpp
// 1. STL headers — alphabetised
#include <array>
#include <span>
#include <string>

// 2. Third-party
#include "Eigen/Dense"

// 3. Project headers
#include "MOM/MomentMethodBase.hpp"
#include "MOM/ThermoProxy.hpp"
```

---

## Testing

Four test programs are provided in `tests/`:

### `test_verify_all_variants`

Constructs all four variants, sets representative flame-condition state, runs
`ComputeSources()`, and verifies that all source spans have the correct size and
no NaN/Inf values.

### `test_loop_benchmark`

Validates and benchmarks both CFD loop patterns:

- **Pattern A** — templated cell loop over ThreeEquations (compile-time dispatch)
- **Pattern B** — `ForEachCell` with `AnyMomentMethod` (runtime dispatch)

Also validates the zero-fallback mechanism: confirms that `sources_condensation()` on
BrookesMoss and `sources_sintering()` on ThreeEquations return all-zero spans of the
correct size.

```
=== Zero-fallback validation ===
  BrookesMoss::sources_condensation() → size=2  all_zero=true
  BrookesMoss::sources_sintering()    → size=2  all_zero=true
  ThreeEquations::sources_sintering() → size=3  all_zero=true
  Result: PASS ✓

=== Numerical equivalence check ===
  Pattern A checksum: X.XXXXXXe+XX
  Pattern B checksum: X.XXXXXXe+XX
  Bit-identical: true
```

### `reporter_usage_example`

Demonstrates the observer pattern: constructs all soot variants, runs a synthetic 10-cell
loop, and writes diagnostic output files using `MomentMethodReporter`. Verifies that the
reporter operates on any type satisfying `MomentMethod<M>` without variant-specific code.

### `test_oh_fixture`

A standalone, no-Fluent harness (`tests/test_oh_fixture.cpp`) that drives
`HMOM::ComputeSources` on a single cell with **only** the oxidation process active
**and** the O₂ channel zeroed out, isolating the OH oxidation channel. It was written to
validate the G1 OH-stoichiometry fix (see [G1 — OH Oxidation Stoichiometry](#g1--oh-oxidation-stoichiometry)).

The O₂ mole fraction is set to a trace value (`X[O2] = 1e-30`) so that `conc_O2_ → 0`,
which drives the O₂ split fraction `fO2 → 0` and the OH split fraction `fOH → 1` while
keeping O₂ present in the thermo (avoiding any divide-by-zero). O₂ therefore remains
selectable by name in `HMOM::SetState`.

Representative operating point (sooting laminar flame, post-flame region):

| Quantity | Value |
|---|---|
| Temperature | 1800 K |
| Pressure | 101 325 Pa (1 atm) |
| Viscosity | 4.5 × 10⁻⁵ kg/(m·s) |
| Species | H, OH, O₂ (trace), H₂, H₂O, C₂H₂, CO, N₂ |

The run is fully deterministic — no RNG, no time-stepping — and writes
`mom/build/oh_fixture_results.json` (all inputs and outputs at 17-significant-digit
precision) plus a human-readable stdout summary. The JSON is consumed by
`g1_closure_analysis.py`.

Run it via ctest (whose working directory is the build tree, so the JSON lands in
`mom/build/`):

```bash
ctest --test-dir build -R OHFixture --output-on-failure
```

### `g1_closure_analysis.py`

A read-only, Python 3.13 stdlib-only analysis script (`tests/g1_closure_analysis.py`)
that consumes `mom/build/oh_fixture_results.json` and computes, for the OH oxidation
channel:

- **Soot-carbon loss** from `S_M10` via `ω_soot = ρ_soot · V0 · N_A · S_M10`
  (independent of the gas-phase stoichiometry).
- **OH / CO / H₂ balances** (per-C and per-C₂ ratios).
- **C / H / O element balances**.
- **Element-consistent mass balance** (immune to MW-table rounding).
- **Enthalpy closure** against the per-C₂ reaction `C2 + 2OH → 2CO + H2`.

It classifies the applied stoichiometry (per-C vs per-C₂) and prints a PASS/FAIL verdict
against a `< 1e-6` relative closure tolerance. It modifies no source files.

```bash
# Paths are resolved relative to the script, so it can be run from any directory
python tests/g1_closure_analysis.py
```

---

## G1 — OH Oxidation Stoichiometry

### The problem

The HMOM soot oxidation rate `R_oxid_C2` is a **per-C₂-pair** rate (each unit removes two
carbon atoms from the soot phase, confirmed independently from the `S_M10` soot-carbon
loss). The O₂ channel correctly applies per-C₂ stoichiometry `C2 + O2 → 2CO` (1 O₂
consumed, 2 CO produced). The code comment for the OH channel likewise states the intended
per-C₂ reaction `C2 + 2OH → 2CO + H2`.

However, the OH channel **implementation** applied **per-C-atom** coefficients (1 OH, 1 CO,
0.5 H₂) at the per-C₂ rate level. Because the rate itself is per-C₂, this under-counted the
OH channel by a factor of two, producing:

- 50 % OH under-consumption,
- 50 % CO under-production,
- 50 % H₂ under-production,
- ~50 % mass and carbon imbalance (one of the two carbon atoms per C₂ pair vanished from
  the gas phase),
- 50 % under-estimation of the OH-channel heat release.

Hydrogen and oxygen balances closed even before the fix, because `1 OH → 1 CO + 0.5 H2` is
internally consistent for H and O — only carbon (and hence mass and enthalpy) failed to
close.

### The fix

In `include/HMOM/HMOM.tpp` (OH channel, near lines 1113 / 1117 / 1122) the gas-phase terms
were multiplied to the per-C₂ coefficients, matching the O₂-channel convention:

| Species | Before (per-C, bug) | After (per-C₂, correct) |
|---|---|---|
| OH consumed | `R_OH · MW_OH / 1000` | `R_OH · 2 · MW_OH / 1000` |
| CO produced | `R_OH · MW_CO / 1000` | `R_OH · 2 · MW_CO / 1000` |
| H₂ produced | `R_OH · 0.5 · MW_H2 / 1000` | `R_OH · 1 · MW_H2 / 1000` |

Only the gas-phase stoichiometry was changed. The `kox_OH_` rate-constant prefactor and
the `R_oxid_C2` computation were left untouched — the prefactor `0.5 / (α · χ)` is a
sticking/calibration factor, not a per-C → per-C₂ converter.

### Verification

After the fix, `g1_closure_analysis.py` reports all balances closing to machine epsilon
(well within the `< 1e-6` acceptance threshold):

| Balance | Post-fix imbalance (relative) | Status |
|---|---|---|
| Mass (element-consistent) | 5.0 × 10⁻¹⁶ | ✅ |
| Carbon | 3.8 × 10⁻¹⁶ | ✅ |
| Hydrogen | −1.9 × 10⁻¹⁶ | ✅ |
| Oxygen | 0.0 | ✅ |
| Enthalpy (Q) | 1.7 × 10⁻¹⁶ | ✅ |

The per-C₂ ratios are exactly `OH_per_C2 = 2.0`, `CO_per_C2 = 2.0`, `H2_per_C2 = 1.0`.

Full evidence — the closure analysis and the fix verdict — is recorded in
`.swarm/evidence/1/G1-closure-analysis.md` and `.swarm/evidence/1/G1-verdict.md`.

> **Build note (MSVC):** the MOM sources are UTF-8 encoded (they contain characters such as
> `C₂H₂`, `→`, `✓` in comments and strings). `CMakeLists.txt` adds `/utf-8` for MSVC so the
> compiler does not decode sources with the host codepage (e.g. GBK/936), which would
> otherwise corrupt these sequences and raise spurious C2065/C3688 errors.

---

## References

The physical models implemented in this library are based on the following publications.
Please cite the appropriate references when using this library in scientific work.

**HMOM**

- M.E. Mueller, G. Blanquart, H. Pitsch, "Hybrid Method of Moments for modelling soot
  formation and growth," *Combustion and Flame* **156** (2009) 1143–1155.
- A. Attili, F. Bisetti, "Statistics and scaling of turbulence in a spatially developing
  mixing layer at Re_λ = 250," *Physics of Fluids* **24** (2012) 035109.

**Three-Equation Soot Model**

- C. Eberle, P. Gerlinger, M. Aigner, "A sectional PAH model with reversible PAH chemistry
  for CFD soot simulations," *Combustion and Flame* **162** (2015) 2648–2660.

**Brookes–Moss Soot Model**

- S.J. Brookes, J.B. Moss, "Predictions of soot and thermal radiation properties in
  confined turbulent jet diffusion flames," *Combustion and Flame* **116** (1999) 486–503.

**TiO₂ Nanoparticle Model**

- S.E. Pratsinis, "Simultaneous nucleation, condensation, and coagulation in aerosol
  reactors," *Journal of Colloid and Interface Science* **124** (1988) 416–427.
- R. Jossen, S.E. Pratsinis, W.J. Stark, L. Mädler, "Criteria for flame-spray synthesis
  of hollow, shell-like, or inhomogeneous oxides," *Journal of the American Ceramic
  Society* **88** (2005) 1388–1393.

**Planck Mean Absorption Coefficients**

- M.D. Smooke et al., *Combustion and Flame* **143** (2005) 613–628.
- J.H. Kent, D. Honnery, *Combustion Science and Technology* **75** (1990) 167–177.
- S.S. Sazhin, "An approximation for the absorption of thermal radiation by soot
  particles," *Progress in Energy and Combustion Science* **20** (1994) 297–318.

---

## License

Copyright © 2026 Alberto Cuoci.

This library is distributed under the **GNU General Public License v3.0 or later**.  
See the [`LICENSE`](LICENSE) file for the full text.

```
OpenSMOKEpp MOM is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later version.
```

---

*CRECK Modeling Lab — Politecnico di Milano*  
*Department of Chemistry, Materials, and Chemical Engineering*  
*P.zza Leonardo da Vinci 32, 20133 Milano, Italy*  
*<https://www.creckmodeling.polimi.it>*
