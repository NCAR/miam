// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <miam/util/constants.hpp>
#include <miam/math/condensation_rate.hpp>

#include <cmath>
#include <numbers>

#include <gtest/gtest.h>

using namespace miam;

namespace
{
  // Physical parameters for test cases
  constexpr double D_g = 1.3e-5;       // Gas diffusion coefficient [m² s⁻¹] (typical for SO2 in air)
  constexpr double alpha = 0.1;        // Mass accommodation coefficient
  constexpr double gas_molecular_weight = 0.064064;  // Molecular weight of SO2 [kg mol⁻¹]
  constexpr double T = 298.15;         // Temperature [K]
  constexpr double r_eff = 1.0e-6;     // Effective radius [m] (1 μm)
  constexpr double N = 1.0e9;          // Number concentration [# m⁻³]

  // Helper: compute mean molecular speed
  double mean_molecular_speed(double T, double Mw)
  {
    return std::sqrt(8.0 * miam::GAS_CONSTANT * T / (std::numbers::pi * Mw));
  }

  // Helper: compute mean free path
  double mean_free_path(double D_g, double c_bar)
  {
    return 3.0 * D_g / c_bar;
  }

  // Helper: compute Fuchs-Sutugin correction factor (Fuchs1971, Zaveri2008)
  double fuchs_sutugin(double Kn, double alpha)
  {
    double denom = Kn * Kn + (1.0 + 0.283 * alpha) * Kn + 0.75 * alpha;
    return 0.75 * alpha * (1.0 + Kn) / denom;
  }
}  // namespace

TEST(CondensationRate, ComputeValueMatchesHandCalculation)
{
  auto provider = MakeCondensationRateProvider(D_g, alpha, gas_molecular_weight);

  double c_bar = mean_molecular_speed(T, gas_molecular_weight);
  double lambda = mean_free_path(D_g, c_bar);
  double Kn = lambda / r_eff;
  double f = fuchs_sutugin(Kn, alpha);
  double expected_k_cond = 4.0 * std::numbers::pi * r_eff * N * D_g * f;

  double k_cond = provider.ComputeValue(r_eff, N, T);
  EXPECT_DOUBLE_EQ(k_cond, expected_k_cond);
}

TEST(CondensationRate, ContinuumLimit)
{
  // For very large particles (Kn → 0), f(Kn) → 1
  // k_cond → 4π · r_eff · N · D_g
  double large_r = 1.0e-3;  // 1 mm radius — continuum regime
  auto provider = MakeCondensationRateProvider(D_g, alpha, gas_molecular_weight);

  double k_cond = provider.ComputeValue(large_r, N, T);
  double continuum_limit = 4.0 * std::numbers::pi * large_r * N * D_g;

  // f(Kn) should be very close to 1 for large particles
  EXPECT_NEAR(k_cond / continuum_limit, 1.0, 5.0e-3);
}

TEST(CondensationRate, FullAccommodation)
{
  // With α = 1, the correction factor should be closer to 1
  double alpha_full = 1.0;
  auto provider_full = MakeCondensationRateProvider(D_g, alpha_full, gas_molecular_weight);
  auto provider_low = MakeCondensationRateProvider(D_g, 0.01, gas_molecular_weight);

  double k_full = provider_full.ComputeValue(r_eff, N, T);
  double k_low = provider_low.ComputeValue(r_eff, N, T);

  // Higher accommodation coefficient should give higher condensation rate
  EXPECT_GT(k_full, k_low);
}

TEST(CondensationRate, LinearInN)
{
  auto provider = MakeCondensationRateProvider(D_g, alpha, gas_molecular_weight);

  double k1 = provider.ComputeValue(r_eff, N, T);
  double k2 = provider.ComputeValue(r_eff, 2.0 * N, T);

  // k_cond is linear in N, so doubling N should double k_cond
  EXPECT_DOUBLE_EQ(k2, 2.0 * k1);
}

TEST(CondensationRate, TemperatureDependence)
{
  auto provider = MakeCondensationRateProvider(D_g, alpha, gas_molecular_weight);

  double k_cold = provider.ComputeValue(r_eff, N, 250.0);
  double k_warm = provider.ComputeValue(r_eff, N, 350.0);

  // Higher temperature → larger mean molecular speed → smaller mean free path
  // → smaller Kn → f(Kn) closer to 1 → generally different rate
  // Both should be positive
  EXPECT_GT(k_cold, 0.0);
  EXPECT_GT(k_warm, 0.0);
}

TEST(CondensationRate, ComputeValueAndDerivativesConsistency)
{
  auto provider = MakeCondensationRateProvider(D_g, alpha, gas_molecular_weight);

  double k_cond_val = provider.ComputeValue(r_eff, N, T);

  double k_cond, dk_dr, dk_dN;
  provider.ComputeValueAndDerivatives(r_eff, N, T, k_cond, dk_dr, dk_dN);

  // k_cond from both methods should match exactly
  EXPECT_DOUBLE_EQ(k_cond, k_cond_val);
}

TEST(CondensationRate, DerivativesDkDrFiniteDifference)
{
  auto provider = MakeCondensationRateProvider(D_g, alpha, gas_molecular_weight);

  double k_cond, dk_dr, dk_dN;
  provider.ComputeValueAndDerivatives(r_eff, N, T, k_cond, dk_dr, dk_dN);

  // Finite difference approximation for dk/dr
  double h = r_eff * 1.0e-6;
  double k_plus = provider.ComputeValue(r_eff + h, N, T);
  double k_minus = provider.ComputeValue(r_eff - h, N, T);
  double dk_dr_fd = (k_plus - k_minus) / (2.0 * h);

  EXPECT_NEAR(dk_dr, dk_dr_fd, std::abs(dk_dr_fd) * 1.0e-5);
}

TEST(CondensationRate, DerivativesDkDnFiniteDifference)
{
  auto provider = MakeCondensationRateProvider(D_g, alpha, gas_molecular_weight);

  double k_cond, dk_dr, dk_dN;
  provider.ComputeValueAndDerivatives(r_eff, N, T, k_cond, dk_dr, dk_dN);

  // dk/dN = k/N (linear), so dk/dN should match k_cond/N exactly
  EXPECT_DOUBLE_EQ(dk_dN, k_cond / N);

  // Also verify with finite difference
  double h = N * 1.0e-6;
  double k_plus = provider.ComputeValue(r_eff, N + h, T);
  double k_minus = provider.ComputeValue(r_eff, N - h, T);
  double dk_dN_fd = (k_plus - k_minus) / (2.0 * h);

  EXPECT_NEAR(dk_dN, dk_dN_fd, std::abs(dk_dN_fd) * 1.0e-5);
}

TEST(CondensationRate, DerivativesDkDrAtDifferentRegimes)
{
  auto provider = MakeCondensationRateProvider(D_g, alpha, gas_molecular_weight);

  // Small particle (large Kn — free molecular regime)
  {
    double r_small = 1.0e-8;  // 10 nm
    double k_cond, dk_dr, dk_dN;
    provider.ComputeValueAndDerivatives(r_small, N, T, k_cond, dk_dr, dk_dN);

    // Use h = r * 1e-4 here: the standard 1e-6 produces step sizes near floating-point
    // precision for a 10 nm radius and causes catastrophic cancellation.
    double h = r_small * 1.0e-4;
    double k_plus = provider.ComputeValue(r_small + h, N, T);
    double k_minus = provider.ComputeValue(r_small - h, N, T);
    double dk_dr_fd = (k_plus - k_minus) / (2.0 * h);

    EXPECT_NEAR(dk_dr, dk_dr_fd, std::abs(dk_dr_fd) * 1.0e-4);
  }

  // Large particle (small Kn — continuum regime)
  {
    double r_large = 1.0e-4;  // 100 μm
    double k_cond, dk_dr, dk_dN;
    provider.ComputeValueAndDerivatives(r_large, N, T, k_cond, dk_dr, dk_dN);

    double h = r_large * 1.0e-6;
    double k_plus = provider.ComputeValue(r_large + h, N, T);
    double k_minus = provider.ComputeValue(r_large - h, N, T);
    double dk_dr_fd = (k_plus - k_minus) / (2.0 * h);

    EXPECT_NEAR(dk_dr, dk_dr_fd, std::abs(dk_dr_fd) * 1.0e-4);
  }
}

TEST(CondensationRate, DerivativesContinuumLimit)
{
  // For very large radius (Kn → 0, f → 1):
  // k_cond ≈ 4π · r · N · D_g
  // dk/dr ≈ 4π · N · D_g
  double r_large = 1.0e-3;
  auto provider = MakeCondensationRateProvider(D_g, alpha, gas_molecular_weight);

  double k_cond, dk_dr, dk_dN;
  provider.ComputeValueAndDerivatives(r_large, N, T, k_cond, dk_dr, dk_dN);

  double expected_dk_dr = 4.0 * std::numbers::pi * N * D_g;
  EXPECT_NEAR(dk_dr / expected_dk_dr, 1.0, 1.0e-5);

  // Also verify against finite difference
  double h = r_large * 1.0e-6;
  double k_plus = provider.ComputeValue(r_large + h, N, T);
  double k_minus = provider.ComputeValue(r_large - h, N, T);
  double dk_dr_fd = (k_plus - k_minus) / (2.0 * h);
  EXPECT_NEAR(dk_dr, dk_dr_fd, std::abs(dk_dr_fd) * 1.0e-5);
}

TEST(CondensationRate, DerivativesFreeMolecularFiniteDifference)
{
  // Deep free-molecular regime (Kn ≈ 124 at r = 1 nm): verify dk/dr against finite
  // difference. At this radius the standard h = r*1e-6 causes catastrophic cancellation;
  // h = r*1e-4 gives ~1e-8 truncation error which is well within tolerance.
  double r_small = 1.0e-9;  // 1 nm
  auto provider = MakeCondensationRateProvider(D_g, alpha, gas_molecular_weight);

  double k_cond, dk_dr, dk_dN;
  provider.ComputeValueAndDerivatives(r_small, N, T, k_cond, dk_dr, dk_dN);

  double h = r_small * 1.0e-4;
  double k_plus = provider.ComputeValue(r_small + h, N, T);
  double k_minus = provider.ComputeValue(r_small - h, N, T);
  double dk_dr_fd = (k_plus - k_minus) / (2.0 * h);

  EXPECT_NEAR(dk_dr, dk_dr_fd, std::abs(dk_dr_fd) * 1.0e-4);
}

TEST(CondensationRate, DifferentGasSpecies)
{
  // Compare two gases with different molecular weights but same D_g and alpha
  double molecular_weight_light = 0.018015;  // Water [kg mol⁻¹]
  double molecular_weight_heavy = 0.098079;  // H2SO4 [kg mol⁻¹]

  auto provider_light = MakeCondensationRateProvider(D_g, alpha, molecular_weight_light);
  auto provider_heavy = MakeCondensationRateProvider(D_g, alpha, molecular_weight_heavy);

  double k_light = provider_light.ComputeValue(r_eff, N, T);
  double k_heavy = provider_heavy.ComputeValue(r_eff, N, T);

  // Both should be positive and different
  EXPECT_GT(k_light, 0.0);
  EXPECT_GT(k_heavy, 0.0);
  EXPECT_NE(k_light, k_heavy);
}

TEST(CondensationRate, PositiveKcondForAllInputs)
{
  auto provider = MakeCondensationRateProvider(D_g, alpha, gas_molecular_weight);

  // Test a range of realistic physical conditions
  for (double r : { 1.0e-9, 1.0e-8, 1.0e-7, 1.0e-6, 1.0e-5, 1.0e-4 })
  {
    for (double n : { 1.0e6, 1.0e9, 1.0e12 })
    {
      for (double t : { 200.0, 298.15, 400.0 })
      {
        EXPECT_GT(provider.ComputeValue(r, n, t), 0.0) << "r=" << r << " N=" << n << " T=" << t;
      }
    }
  }
}

TEST(CondensationRate, FreeMolecularLimit)
{
  // For very small particles (Kn → ∞), the Fuchs-Sutugin factor approaches
  // f → 0.75α/Kn, so k_cond → 4π·r·N·D·(0.75α/Kn) = π·r²·N·c̄·α
  // This is the kinetic-theory (free-molecular) limit.
  double r_small = 1.0e-10;  // 0.1 nm — deep free-molecular regime
  auto provider = MakeCondensationRateProvider(D_g, alpha, gas_molecular_weight);

  double c_bar = std::sqrt(8.0 * miam::GAS_CONSTANT * T / (std::numbers::pi * gas_molecular_weight));
  double k_fm_limit = std::numbers::pi * r_small * r_small * N * c_bar * alpha;

  double k_cond = provider.ComputeValue(r_small, N, T);
  EXPECT_NEAR(k_cond / k_fm_limit, 1.0, 1.0e-3);
}

TEST(CondensationRate, DerivativesFreeMolecularLimit)
{
  // In the free-molecular limit k_cond ≈ π·r²·N·c̄·α,
  // so dk/dr ≈ 2π·r·N·c̄·α.
  double r_small = 1.0e-10;  // 0.1 nm — deep free-molecular regime
  auto provider = MakeCondensationRateProvider(D_g, alpha, gas_molecular_weight);

  double c_bar = std::sqrt(8.0 * miam::GAS_CONSTANT * T / (std::numbers::pi * gas_molecular_weight));
  double dk_dr_fm_limit = 2.0 * std::numbers::pi * r_small * N * c_bar * alpha;

  double k_cond, dk_dr, dk_dN;
  provider.ComputeValueAndDerivatives(r_small, N, T, k_cond, dk_dr, dk_dN);
  EXPECT_NEAR(dk_dr / dk_dr_fm_limit, 1.0, 1.0e-3);
}

// ============================================================================
// Limit / Extreme tests (Phase E1)
// ============================================================================

TEST(CondensationRate, ZeroRadiusReturnsZero)
{
  // r_eff = 0: k_cond = 4π·r·N·D·f(Kn,α) → 0 as r → 0 (Kn → ∞, f ~ 0.75α/Kn ~ c̄r/D)
  auto provider = MakeCondensationRateProvider(D_g, alpha, gas_molecular_weight);
  double k = provider.ComputeValue(0.0, N, T);
  EXPECT_NEAR(k, 0.0, 1.0e-30);
  EXPECT_FALSE(std::isnan(k));

  double k_cond, dk_dr, dk_dN;
  provider.ComputeValueAndDerivatives(0.0, N, T, k_cond, dk_dr, dk_dN);
  EXPECT_NEAR(k_cond, 0.0, 1.0e-30);
  EXPECT_FALSE(std::isnan(k_cond));
  EXPECT_FALSE(std::isinf(dk_dr));
}

TEST(CondensationRate, ZeroNumberConcentrationReturnsZero)
{
  // N = 0: k_cond = 4π·r·N·D·f(Kn,α) → 0 linearly in N
  auto provider = MakeCondensationRateProvider(D_g, alpha, gas_molecular_weight);
  double k = provider.ComputeValue(r_eff, 0.0, T);
  EXPECT_NEAR(k, 0.0, 1.0e-30);
  EXPECT_FALSE(std::isnan(k));

  double k_cond, dk_dr, dk_dN;
  provider.ComputeValueAndDerivatives(r_eff, 0.0, T, k_cond, dk_dr, dk_dN);
  EXPECT_NEAR(k_cond, 0.0, 1.0e-30);
  EXPECT_FALSE(std::isnan(k_cond));
  EXPECT_FALSE(std::isnan(dk_dN));
}

TEST(CondensationRate, CAMPZaveri2008FormulaExact)
{
  // Direct formula verification: f(Kn,α) = 0.75α(1+Kn) / (Kn² + (1+0.283α)Kn + 0.75α)
  // at specific Kn = 0.5, α = 0.1.  k_cond = 4π·r·N·D·f(Kn,α).
  constexpr double alpha_test = 0.1;
  auto provider = MakeCondensationRateProvider(D_g, alpha_test, gas_molecular_weight);

  constexpr double r_test = 1.0e-6;  // 1 μm
  constexpr double N_test = 1.0e9;
  constexpr double T_test = 298.15;

  double c_bar = mean_molecular_speed(T_test, gas_molecular_weight);
  double lambda = mean_free_path(D_g, c_bar);
  double Kn = lambda / r_test;
  double f = fuchs_sutugin(Kn, alpha_test);
  double k_expected = 4.0 * std::numbers::pi * r_test * N_test * D_g * f;

  double k_computed = provider.ComputeValue(r_test, N_test, T_test);
  EXPECT_NEAR(k_computed, k_expected, k_expected * 1.0e-12);
}
