// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/util/error.hpp>
#include <miam/util/miam_exception.hpp>

#include <micm/Util.hpp>
#include <micm/util/constants.hpp>

#include <cmath>
#include <numbers>
#include <stdexcept>

namespace miam
{
  /// @brief Provider for condensation rate and its derivatives
  /// @details Encapsulates the Fuchs-Sutugin transition regime calculation \cite Fuchs1971, Zaveri2008:
  ///
  ///          k_cond = 4pi * r_eff * N * D * f(Kn, alpha)                      [s^-1]
  ///
  ///          f(Kn, alpha) = 0.75alpha(1 + Kn) / (Kn^2 + (1 + 0.283alpha)Kn + 0.75alpha)   [dimensionless]
  ///
  ///          Kn = lambda / r_eff                                              [dimensionless]
  ///          lambda  = 3*D / c_bar                                                [m]
  ///          c_bar  = sqrt(8*R*T / (pi*M))                                       [m s^-1]
  ///
  ///          Variable definitions:
  ///            k_cond  First-order condensation rate coefficient         [s^-1]
  ///            r_eff   Effective radius of the aerosol particle          [m]
  ///            N       Aerosol particle number concentration             [# m^-3]
  ///            D       Gas-phase diffusion coefficient                   [m^2 s^-1]
  ///            f       Fuchs-Sutugin transition regime correction factor [dimensionless]
  ///            Kn      Knudsen number                                    [dimensionless]
  ///            lambda       Mean free path of gas molecules                   [m]
  ///            c_bar       Mean molecular speed of gas molecules             [m s^-1]
  ///            R       Ideal gas constant (8.314 J mol^-1 K^-1)            [J mol^-1 K^-1]
  ///            T       Temperature                                       [K]
  ///            M       Molecular weight of the gas species               [kg mol^-1]
  ///            alpha       Mass accommodation coefficient                    [dimensionless, 0-1]
  ///
  ///          The provider is a plain struct with device-callable methods, so device kernels can
  ///          capture it by value. Create it with `MakeCondensationRateProvider`, which checks the inputs.
  struct CondensationRateProvider
  {
    double diffusion_coefficient_ = 0.0;      ///< D [m^2 s^-1]
    double accommodation_coefficient_ = 0.0;  ///< alpha [dimensionless, 0-1]
    double molecular_weight_ = 0.0;           ///< M [kg mol^-1]

    /// @brief Compute condensation rate k_cond [s^-1]
    /// @param r_eff Effective radius [m]
    /// @param N Number concentration [# m^-3]
    /// @param T Temperature [K]
    /// @return k_cond [s^-1]
    MICM_DEVICE_FUNCTION double ComputeValue(double r_eff, double N, double T) const
    {
      if (r_eff <= 0 || N <= 0 || T <= 0)
        return 0.0;
      constexpr double pi = std::numbers::pi;
      const double f = CorrectionFactor(KnudsenNumber(r_eff, T));
      return 4.0 * pi * r_eff * N * diffusion_coefficient_ * f;
    }

    /// @brief Compute condensation rate and partial derivatives
    /// @param r_eff Effective radius [m]
    /// @param N Number concentration [# m^-3]
    /// @param T Temperature [K]
    /// @param k_cond Output: condensation rate [s^-1]
    /// @param dk_dr Output: dk_cond/dr_eff [s^-1 m^-1]
    /// @param dk_dN Output: dk_cond/dN [s^-1 m^3 #^-1]
    MICM_DEVICE_FUNCTION void
    ComputeValueAndDerivatives(double r_eff, double N, double T, double& k_cond, double& dk_dr, double& dk_dN) const
    {
      if (r_eff <= 0 || N <= 0 || T <= 0)
      {
        k_cond = 0.0;
        dk_dr = 0.0;
        dk_dN = 0.0;
        return;
      }
      constexpr double pi = std::numbers::pi;
      const double alpha = accommodation_coefficient_;
      const double Kn = KnudsenNumber(r_eff, T);
      const double denom = Kn * Kn + (1.0 + 0.283 * alpha) * Kn + 0.75 * alpha;
      const double f = 0.75 * alpha * (1.0 + Kn) / denom;

      k_cond = 4.0 * pi * r_eff * N * diffusion_coefficient_ * f;

      // df/dKn = 0.75alpha * (-Kn^2 - 2Kn + (0.467alpha - 1)) / denom^2
      // where 0.467 = 0.75 - 0.283
      const double df_dKn = 0.75 * alpha * (-Kn * Kn - 2.0 * Kn + (0.75 - 0.283) * alpha - 1.0) / (denom * denom);

      // dk_cond/dr_eff = 4pi*N*D * (f + r_eff * df/dKn * dKn/dr)
      // where dKn/dr = -Kn / r_eff
      // simplifies to: 4pi*N*D * (f - Kn * df/dKn)
      dk_dr = 4.0 * pi * N * diffusion_coefficient_ * (f - Kn * df_dKn);

      // dk_cond/dN = k_cond / N  (linear in N)
      dk_dN = k_cond / N;
    }

   private:
    /// @brief Kn = lambda / r_eff, with lambda = 3*D / c_bar and c_bar = sqrt(8*R*T / (pi*M))
    MICM_DEVICE_FUNCTION double KnudsenNumber(double r_eff, double T) const
    {
      constexpr double pi = std::numbers::pi;
      constexpr double gas_constant = micm::constants::GAS_CONSTANT;
      const double c_bar = std::sqrt(8.0 * gas_constant * T / (pi * molecular_weight_));
      const double lambda = 3.0 * diffusion_coefficient_ / c_bar;
      return lambda / r_eff;
    }

    /// @brief f(Kn, alpha) = 0.75alpha(1 + Kn) / (Kn^2 + (1 + 0.283alpha)Kn + 0.75alpha)
    MICM_DEVICE_FUNCTION double CorrectionFactor(double Kn) const
    {
      const double alpha = accommodation_coefficient_;
      const double denom = Kn * Kn + (1.0 + 0.283 * alpha) * Kn + 0.75 * alpha;
      return 0.75 * alpha * (1.0 + Kn) / denom;
    }
  };

  /// @brief Factory function to create a CondensationRateProvider
  /// @param diffusion_coefficient Gas-phase diffusion coefficient [m^2 s^-1]
  /// @param accommodation_coefficient Mass accommodation coefficient [dimensionless, 0-1]
  /// @param molecular_weight Molecular weight of the gas species [kg mol^-1]
  /// @return A CondensationRateProvider with the Fuchs-Sutugin regime correction
  inline CondensationRateProvider
  MakeCondensationRateProvider(double diffusion_coefficient, double accommodation_coefficient, double molecular_weight)
  {
    if (diffusion_coefficient <= 0)
    {
      throw MiamException(
          MIAM_ERROR_CATEGORY_CONFIGURATION,
          MIAM_CONFIGURATION_INVALID_PARAMETER,
          "Diffusion coefficient must be positive.");
    }
    if (accommodation_coefficient <= 0)
    {
      throw MiamException(
          MIAM_ERROR_CATEGORY_CONFIGURATION,
          MIAM_CONFIGURATION_INVALID_PARAMETER,
          "Accommodation coefficient must be positive.");
    }
    if (molecular_weight <= 0)
    {
      throw MiamException(
          MIAM_ERROR_CATEGORY_CONFIGURATION, MIAM_CONFIGURATION_INVALID_PARAMETER, "Molecular weight must be positive.");
    }
    return CondensationRateProvider{ diffusion_coefficient, accommodation_coefficient, molecular_weight };
  }
}  // namespace miam
