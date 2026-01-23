// Copyright (C) 2026 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <miam/model/mode.hpp>

#include <cmath>
#include <format>
#include <stdexcept>
#include <string>

namespace miam
{

  double Mode::GetEffectiveRadius() const
  {
    const double r_g = 0.5 * geometric_mean_diameter_;
    const double ln_sig = std::log(geometric_standard_deviation_);
    return r_g * std::exp(2.5 * ln_sig * ln_sig);
  }

  double Mode::CalculateEffectiveRadius(double mass, double N, double density, double sig_g)
  {
    if (mass < 1e-18 || N < 1e-10)
    {
      throw std::runtime_error(std::format(
          "Cannot calculate effective radius: mass ({}) or number ({}) is below the numerical stability limit for '{}'.",
          mass,
          N,
          name_));
    }

    double V = mass / density;
    double ln_sig = std::log(sig_g);

    // V = N * (4/3) * pi * r_g^3 * exp(9/2 * ln(std_dev) * ln(std_dev))
    double exp_term_vol = std::exp(4.5 * ln_sig * ln_sig);
    double rg = std::pow((3.0 * V) / (4.0 * M_PI * N * exp_term_vol), 1.0 / 3.0);

    // r_eff = (r_g) * exp(5/2 * ln(std_dev) * ln(std_dev))
    return rg * std::exp(2.5 * ln_sig * ln_sig);
  }

}  // namespace miam