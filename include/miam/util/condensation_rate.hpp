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
    ///          k_cond = 4π · r_eff · N · D_g · f(Kn)
    ///          where f(Kn) = (1 + Kn) / (1 + 2·Kn·(1 + Kn) / α)
    ///          Kn = λ / r_eff, λ = 3·D_g / c̄, c̄ = √(8RT / (π·Mw))
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
    /// @param D_g Gas-phase diffusion coefficient [m² s⁻¹]
    /// @param alpha Mass accommodation coefficient [dimensionless, 0-1]
    /// @param Mw_gas Molecular weight of the gas species [kg mol⁻¹]
    /// @return A CondensationRateProvider with the Fuchs-Sutugin regime correction
    inline CondensationRateProvider MakeCondensationRateProvider(double D_g, double alpha, double Mw_gas)
    {
      if (D_g <= 0)
        throw std::invalid_argument("Diffusion coefficient D_g must be positive.");
      if (alpha <= 0)
        throw std::invalid_argument("Accommodation coefficient alpha must be positive.");
      if (Mw_gas <= 0)
        throw std::invalid_argument("Molecular weight Mw_gas must be positive.");

      CondensationRateProvider provider;

      provider.ComputeValue = [D_g, alpha, Mw_gas](double r_eff, double N, double T) -> double
      {
        if (r_eff <= 0 || N <= 0 || T <= 0)
          return 0.0;
        double c_bar = std::sqrt(8.0 * R_gas * T / (std::numbers::pi * Mw_gas));
        double lambda = 3.0 * D_g / c_bar;
        double Kn = lambda / r_eff;
        double denom = 1.0 + 2.0 * Kn * (1.0 + Kn) / alpha;
        double f = (1.0 + Kn) / denom;
        return 4.0 * std::numbers::pi * r_eff * N * D_g * f;
      };

      provider.ComputeValueAndDerivatives =
          [D_g, alpha, Mw_gas](double r_eff, double N, double T, double& k_cond, double& dk_dr, double& dk_dN)
      {
        if (r_eff <= 0 || N <= 0 || T <= 0)
        {
          k_cond = 0.0;
          dk_dr = 0.0;
          dk_dN = 0.0;
          return;
        }
        double c_bar = std::sqrt(8.0 * R_gas * T / (std::numbers::pi * Mw_gas));
        double lambda = 3.0 * D_g / c_bar;
        double Kn = lambda / r_eff;
        double denom = 1.0 + 2.0 * Kn * (1.0 + Kn) / alpha;
        double f = (1.0 + Kn) / denom;

        k_cond = 4.0 * std::numbers::pi * r_eff * N * D_g * f;

        // df/dKn = (α - 2 - 4·Kn - 2·Kn²) / (α · denom²)
        double df_dKn = (alpha - 2.0 - 4.0 * Kn - 2.0 * Kn * Kn) / (alpha * denom * denom);

        // dk_cond/dr_eff = 4π·N·D_g · (f + r_eff · df/dKn · dKn/dr)
        // where dKn/dr = -Kn / r_eff
        // simplifies to: 4π·N·D_g · (f - Kn · df/dKn)
        dk_dr = 4.0 * std::numbers::pi * N * D_g * (f - Kn * df_dKn);

        // dk_cond/dN = k_cond / N  (linear in N)
        dk_dN = k_cond / N;
      };

      return provider;
    }
  }  // namespace util
}  // namespace miam
