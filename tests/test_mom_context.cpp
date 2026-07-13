/*-----------------------------------------------------------------------*\
|   MOM Library -- MomContext unit tests (task 3.4)                        |
|                                                                          |
|   Standalone C++ test (plain main(), manual checks with std::exit(1) on  |
|   failure, return 0 on success) exercising MOM::MomContext: construction  |
|   from a valid manifest and from an invalid path, a finite compute() on   |
|   one cell, no-op compute() on an invalid context, particle-property      |
|   accessors, species names, zone IDs, and manifest hash.                  |
|                                                                          |
|   Species order (0-based): H OH O2 H2 H2O C2H2 CO A4 N2 (N2 last, bulk). |
\*-----------------------------------------------------------------------*/

#include "MOM/MomContext.hpp"

#include "MOM/ManifestReader.hpp"
#include "test_manifest_fixture.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <span>
#include <string>
#include <type_traits>
#include <vector>

using mom_test::make_test_manifest;

namespace
{

int g_checks = 0;

void check(bool cond, const std::string& what)
{
    ++g_checks;
    if (!cond)
    {
        std::cerr << "FAIL: " << what << "\n";
        std::exit(1);
    }
}

bool almost_equal(double a, double b, double tol = 1e-9)
{
    return std::fabs(a - b) <= tol * (1.0 + std::fabs(a) + std::fabs(b));
}

} // namespace

// 9. Copy/move must not be available (design-final P0-03). Verified at compile
// time via static_assert rather than attempting a failing compilation.
static_assert(!std::is_copy_constructible_v<MOM::MomContext>,
              "MomContext must not be copy-constructible");
static_assert(!std::is_copy_assignable_v<MOM::MomContext>,
              "MomContext must not be copy-assignable");
static_assert(!std::is_move_constructible_v<MOM::MomContext>,
              "MomContext must not be move-constructible");

int main()
{
    std::cout << "=== MomContext unit tests (task 3.4) ===\n";

    const MOM::Manifest manifest = make_test_manifest();

    // -- 1. Construct from a valid Manifest: is_valid, n_species. ------------
    {
        MOM::MomContext ctx(manifest);
        check(ctx.is_valid(), "valid: is_valid() true");
        check(ctx.error_message().empty(), "valid: error_message empty");
        check(ctx.n_species() == 9u, "valid: n_species == 9");
        std::cout << "  [ok] 1. construct from valid manifest\n";
    }

    // -- 2. Construct from invalid path: not valid, message non-empty. -------
    {
        MOM::MomContext ctx("this_manifest_does_not_exist_momctx.json");
        check(!ctx.is_valid(), "invalid: is_valid() false");
        check(!ctx.error_message().empty(), "invalid: error_message non-empty");
        std::cout << "  [ok] 2. construct from invalid path\n";
    }

    // -- 3. compute(): sources (4) and omega_gas (9) all finite. -------------
    {
        MOM::MomContext ctx(manifest);
        check(ctx.is_valid(), "compute: context valid");

        // Representative sooting-flame post-flame state.
        const double T    = 1800.0;    // [K]
        const double P    = 101325.0;  // [Pa]
        const double mu   = 4.5e-5;    // [kg/m/s]

        // Mass fractions (size = n_species), roughly normalised. Order:
        //   H OH O2 H2 H2O C2H2 CO A4 N2
        std::vector<double> Y = {
            1.0e-4, // H
            1.0e-3, // OH
            5.0e-3, // O2
            5.0e-4, // H2
            0.10,   // H2O
            2.0e-3, // C2H2
            2.0e-2, // CO
            1.0e-5, // A4 (PAH precursor)
            0.86,   // N2 (balance)
        };
        double sumY = 0.0;
        for (double y : Y)
            sumY += y;
        for (double& y : Y)
            y /= sumY;

        // Moments [M00, M10, M01, N0] — use the model's own initial (floor)
        // moments scaled up so soot is present. model_ is public on the struct.
        const std::span<const double> ic = ctx.model_->initial_moments();
        const double scale = 1.0e5;
        std::array<double, 4> moments = {
            ic[0] * scale, ic[1] * scale, ic[2] * scale, ic[3] * scale};

        ctx.compute(T, P, std::span<const double>(Y), mu,
                    std::span<const double>(moments));

        const std::span<const double> src   = ctx.sources();
        const std::span<const double> omega = ctx.omega_gas();

        check(src.size() == 4, "compute: sources size == 4");
        check(omega.size() == 9, "compute: omega_gas size == 9");
        for (std::size_t i = 0; i < src.size(); ++i)
            check(std::isfinite(src[i]), "compute: source finite");
        for (std::size_t i = 0; i < omega.size(); ++i)
            check(std::isfinite(omega[i]), "compute: omega finite");
        std::cout << "  [ok] 3. compute() produces finite sources + omega\n";
    }

    // -- 4. compute() on an invalid context is a no-op (no crash). -----------
    {
        MOM::MomContext ctx("still_no_such_manifest.json");
        check(!ctx.is_valid(), "noop: invalid context");
        std::vector<double> Y(9, 0.0);
        std::array<double, 4> moments = {0.0, 0.0, 0.0, 0.0};
        // Must not crash or throw (compute is noexcept and guarded).
        ctx.compute(1800.0, 101325.0, std::span<const double>(Y), 4.5e-5,
                    std::span<const double>(moments));
        std::cout << "  [ok] 4. compute() no-op on invalid context\n";
    }

    // -- 5. V0() > 0, particle_density() == 1800, diffusion_coefficient() > 0. -
    // diffusion_coefficient() depends on the gas state (mu, rho, T, P) that is
    // established by compute()/SetState, so drive one cell first.
    {
        MOM::MomContext ctx(manifest);
        check(ctx.is_valid(), "props: context valid");

        std::vector<double> Y(9, 0.0);
        Y[8] = 1.0; // pure N2 bulk is enough to establish a finite gas state
        std::array<double, 4> moments = {0.0, 0.0, 0.0, 0.0};
        ctx.compute(1800.0, 101325.0, std::span<const double>(Y), 4.5e-5,
                    std::span<const double>(moments));

        check(ctx.V0() > 0.0, "props: V0() > 0");
        check(almost_equal(ctx.particle_density(), 1800.0), "props: particle_density == 1800");
        check(ctx.diffusion_coefficient() > 0.0, "props: diffusion_coefficient > 0");
        std::cout << "  [ok] 5. particle-property accessors\n";
    }

    // -- 6. species_names(): 9 names, first H, last N2. ----------------------
    {
        MOM::MomContext ctx(manifest);
        const std::span<const std::string> names = ctx.species_names();
        check(names.size() == 9, "names: size == 9");
        check(names.front() == "H", "names: first H");
        check(names.back() == "N2", "names: last N2");
        std::cout << "  [ok] 6. species_names()\n";
    }

    // -- 7. zone_ids() matches manifest radiation zones. ---------------------
    {
        MOM::MomContext ctx(manifest);
        const std::span<const int> zones = ctx.zone_ids();
        check(zones.size() == manifest.radiation_zones.size(), "zones: size matches");
        check(zones.size() == 1 && zones[0] == 1, "zones: {1}");
        std::cout << "  [ok] 7. zone_ids()\n";
    }

    // -- 8. manifest_hash() matches manifest_sha256. -------------------------
    {
        MOM::Manifest m  = make_test_manifest();
        m.manifest_sha256 = std::string(64, 'a');
        MOM::MomContext ctx(m);
        check(ctx.manifest_hash() == std::string(64, 'a'), "hash: matches manifest_sha256");

        MOM::MomContext ctx_empty(manifest); // empty hash in default fixture
        check(ctx_empty.manifest_hash().empty(), "hash: empty when manifest hash empty");
        std::cout << "  [ok] 8. manifest_hash()\n";
    }

    std::cout << "=== MomContext: all checks passed (" << g_checks << ") ===\n";
    return 0;
}
