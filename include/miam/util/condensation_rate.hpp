// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <micm/Util.hpp>

#include <cmath>
#include <functional>
#include <numbers>
#include <stdexcept>

namespace miam
{
  namespace util
  {
    /// @brief Gas constant [J mol⁻¹ K⁻¹]
    constexpr double R_gas = micm::constants::GAS_CONSTANT;

    /// @brief Provider for condensation rate and its derivatives
    /// @details Encapsulates the Fuchs-Sutugin transition regime calculation:
    ///          k_cond = 4π · r_eff · N · D · f(Kn)
    ///          where f(Kn) = (1 + Kn) / (1 + 2·Kn·(1 + Kn) / α)
    ///          Kn = λ / r_eff, λ = 3·D / c̄, c̄ = √(8·R·T / (π·molecular_weight))
    struct CondensationRateProvider
    {
      /// @brief Compute condensation rate k_cond [s⁻¹]
      /// @param r_eff Effective radius [m]
      /// @param N Number concentration [# m⁻³]
      /// @param T Temperature [K]
      /// @return k_cond [s⁻¹]
      std::function<double(double r_eff, double N, double T)> ComputeValue;

      /// @brief Compute condensation rate and partial derivatives
      /// @param r_eff Effective radius [m]
      /// @param N Number concentration [# m⁻³]
      /// @param T Temperature [K]
      /// @param k_cond Output: condensation rate [s⁻¹]
      /// @param dk_dr Output: ∂k_cond/∂r_eff [s⁻¹ m⁻¹]
      /// @param dk_dN Output: ∂k_cond/∂N [s⁻¹ m³ #⁻¹]
      std::function<void(double r_eff, double N, double T, double& k_cond, double& dk_dr, double& dk_dN)>
          ComputeValueAndDerivatives;
    };

    /// @brief Factory function to create a CondensationRateProvider
    /// @param diffusion_coefficient Gas-phase diffusion coefficient [m² s⁻¹]
    /// @param accommodation_coefficient Mass accommodation coefficient [dimensionless, 0-1]
    /// @param molecular_weight Molecular weight of the gas species [kg mol⁻¹]
    /// @return A CondensationRateProvider with the Fuchs-Sutugin regime correction
    inline CondensationRateProvider MakeCondensationRateProvider(double diffusion_coefficient, double accommodation_coefficient, double molecular_weight)
    {
      if (diffusion_coefficient <= 0)
        throw std::invalid_argument("Diffusion coefficient must be positive.");
      if (accommodation_coefficient <= 0)
        throw std::invalid_argument("Accommodation coefficient must be positive.");
      if (molecular_weight <= 0)
        throw std::invalid_argument("Molecular weight must be positive.");

      CondensationRateProvider provider;

      provider.ComputeValue = [diffusion_coefficient, accommodation_coefficient, molecular_weight](double r_eff, double N, double T) -> double
      {
        if (r_eff <= 0 || N <= 0 || T <= 0)
          return 0.0;
        double c_bar = std::sqrt(8.0 * R_gas * T / (std::numbers::pi * molecular_weight));
        double lambda = 3.0 * diffusion_coefficient / c_bar;
        double Kn = lambda / r_eff;
        double denom = 1.0 + 2.0 * Kn * (1.0 + Kn) / accommodation_coefficient;
        double f = (1.0 + Kn) / denom;
        return 4.0 * std::numbers::pi * r_eff * N * diffusion_coefficient * f;
      };

      provider.ComputeValueAndDerivatives =
          [diffusion_coefficient, accommodation_coefficient, molecular_weight](double r_eff, double N, double T, double& k_cond, double& dk_dr, double& dk_dN)
      {
        if (r_eff <= 0 || N <= 0 || T <= 0)
        {
          k_cond = 0.0;
          dk_dr = 0.0;
          dk_dN = 0.0;
          return;
        }
        double c_bar = std::sqrt(8.0 * R_gas * T / (std::numbers::pi * molecular_weight));
        double lambda = 3.0 * diffusion_coefficient / c_bar;
        double Kn = lambda / r_eff;
        double denom = 1.0 + 2.0 * Kn * (1.0 + Kn) / accommodation_coefficient;
        double f = (1.0 + Kn) / denom;

        k_cond = 4.0 * std::numbers::pi * r_eff * N * diffusion_coefficient * f;

        // df/dKn = (α - 2 - 4·Kn - 2·Kn²) / (α · denom²)
        double df_dKn = (accommodation_coefficient - 2.0 - 4.0 * Kn - 2.0 * Kn * Kn) / (accommodation_coefficient * denom * denom);

        // dk_cond/dr_eff = 4π·N·D · (f + r_eff · df/dKn · dKn/dr)
        // where dKn/dr = -Kn / r_eff
        // simplifies to: 4π·N·D · (f - Kn · df/dKn)
        dk_dr = 4.0 * std::numbers::pi * N * diffusion_coefficient * (f - Kn * df_dKn);

        // dk_cond/dN = k_cond / N  (linear in N)
        dk_dN = k_cond / N;
      };

      return provider;
    }
  }  // namespace util
}  // namespace miam
