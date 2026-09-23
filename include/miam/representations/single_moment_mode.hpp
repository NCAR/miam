// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/representations/aerosol_property.hpp>
#include <miam/representations/phase_volume_fraction_descriptor.hpp>
#include <miam/util/error.hpp>
#include <miam/util/miam_exception.hpp>

#include <micm/system/phase.hpp>

#include <cmath>
#include <map>
#include <numbers>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

namespace miam
{
  /// @brief Effective radius of a single-moment log-normal aerosol mode: `r_eff = GMD * exp(2.5 * ln^2(GSD))`.
  /// @details Depends only on the mode's two shape parameters (GMD, GSD); has no state-variable dependencies.
  template<typename DenseMatrixPolicy>
  class SingleMomentModeEffectiveRadiusDescriptor
  {
   public:
    struct View
    {
      std::size_t gmd_parameter_index_ = 0;
      std::size_t gsd_parameter_index_ = 0;
    };

    SingleMomentModeEffectiveRadiusDescriptor() = default;

    SingleMomentModeEffectiveRadiusDescriptor(std::size_t gmd_parameter_index, std::size_t gsd_parameter_index)
        : view_{ gmd_parameter_index, gsd_parameter_index }
    {
    }

    View GetView() const
    {
      return view_;
    }

    const std::vector<std::size_t>& DependentVariableIndices() const
    {
      return dependent_variable_indices_;
    }

    void Evaluate(
        const DenseMatrixPolicy& state_parameters,
        const DenseMatrixPolicy& state_variables,
        DenseMatrixPolicy& result) const
    {
      const auto view = view_;
      DenseMatrixPolicy::Function(
          MICM_LAMBDA(const typename DenseMatrixPolicy::ConstViewType& params_view, const typename DenseMatrixPolicy::ConstViewType& /*vars_view*/, const typename DenseMatrixPolicy::ViewType& result_view)
          {
            params_view.ForEachRowStrict(
                [](const double& gmd, const double& gsd, double& r_eff)
                {
                  const double ln_gsd = std::log(gsd);
                  r_eff = gmd * std::exp(2.5 * ln_gsd * ln_gsd);
                },
                params_view.GetConstColumnView(view.gmd_parameter_index_),
                params_view.GetConstColumnView(view.gsd_parameter_index_),
                result_view.GetColumnView(0));
          },
          state_parameters,
          state_variables,
          result)(state_parameters, state_variables, result);
    }

    void EvaluateAndDerivatives(
        const DenseMatrixPolicy& state_parameters,
        const DenseMatrixPolicy& state_variables,
        DenseMatrixPolicy& result,
        DenseMatrixPolicy& /*partials*/) const
    {
      Evaluate(state_parameters, state_variables, result);
    }

   private:
    View view_{};
    std::vector<std::size_t> dependent_variable_indices_{};
  };

  /// @brief Number concentration for a single-moment log-normal aerosol mode:
  ///        `N = (Σ_k [species_k] · MW_k / ρ_k) / V_s`, with `V_s = (4/3)·π·GMD³·exp(4.5·ln²(GSD))`.
  /// @details Depends on GMD, GSD, and the aqueous-species concentrations attached to the mode.
  template<typename DenseMatrixPolicy>
  class SingleMomentModeNumberConcentrationDescriptor
  {
   public:
    template<typename U>
    using Vector = typename DenseMatrixPolicy::template VectorType<U>;

    struct View
    {
      std::size_t gmd_parameter_index_ = 0;
      std::size_t gsd_parameter_index_ = 0;
      typename Vector<std::size_t>::ConstViewType species_variable_indices_{};
      typename Vector<double>::ConstViewType species_molar_volumes_{};
    };

    SingleMomentModeNumberConcentrationDescriptor() = default;

    SingleMomentModeNumberConcentrationDescriptor(
        std::size_t gmd_parameter_index,
        std::size_t gsd_parameter_index,
        std::vector<std::size_t> species_variable_indices,
        std::vector<double> species_molar_volumes)
        : gmd_parameter_index_(gmd_parameter_index),
          gsd_parameter_index_(gsd_parameter_index),
          species_variable_indices_(std::move(species_variable_indices)),
          species_molar_volumes_(std::move(species_molar_volumes)),
          dependent_variable_indices_(species_variable_indices_.begin(), species_variable_indices_.end())
    {
      species_variable_indices_.CopyToDevice();
      species_molar_volumes_.CopyToDevice();
    }

    View GetView() const
    {
      return { gmd_parameter_index_, gsd_parameter_index_, species_variable_indices_.GetView(), species_molar_volumes_.GetView() };
    }

    const std::vector<std::size_t>& DependentVariableIndices() const
    {
      return dependent_variable_indices_;
    }

    void Evaluate(
        const DenseMatrixPolicy& state_parameters,
        const DenseMatrixPolicy& state_variables,
        DenseMatrixPolicy& result) const
    {
      const auto view = GetView();
      DenseMatrixPolicy::Function(
          MICM_LAMBDA(const typename DenseMatrixPolicy::ConstViewType& params_view, const typename DenseMatrixPolicy::ConstViewType& vars_view, const typename DenseMatrixPolicy::ViewType& result_view)
          {
            auto N = result_view.GetColumnView(0);
            params_view.ForEachRowStrict([](double& v) { v = 0.0; }, N);
            const std::size_t n = view.species_variable_indices_.size();
            for (std::size_t k = 0; k < n; ++k)
            {
              const double mv = view.species_molar_volumes_[k];
              params_view.ForEachRowStrict(
                  [mv](const double& c, double& V) { V += c * mv; },
                  vars_view.GetConstColumnView(view.species_variable_indices_[k]),
                  N);
            }
            params_view.ForEachRowStrict(
                [](const double& gmd, const double& gsd, double& N_out)
                {
                  const double ln_gsd = std::log(gsd);
                  const double V_s = (4.0 / 3.0) * std::numbers::pi * gmd * gmd * gmd * std::exp(4.5 * ln_gsd * ln_gsd);
                  N_out /= V_s;
                },
                params_view.GetConstColumnView(view.gmd_parameter_index_),
                params_view.GetConstColumnView(view.gsd_parameter_index_),
                N);
          },
          state_parameters,
          state_variables,
          result)(state_parameters, state_variables, result);
    }

    void EvaluateAndDerivatives(
        const DenseMatrixPolicy& state_parameters,
        const DenseMatrixPolicy& state_variables,
        DenseMatrixPolicy& result,
        DenseMatrixPolicy& partials) const
    {
      Evaluate(state_parameters, state_variables, result);
      if (dependent_variable_indices_.empty())
        return;
      const auto view = GetView();
      DenseMatrixPolicy::Function(
          MICM_LAMBDA(const typename DenseMatrixPolicy::ConstViewType& params_view, const typename DenseMatrixPolicy::ConstViewType& /*vars_view*/, const typename DenseMatrixPolicy::ViewType& /*result_view*/, const typename DenseMatrixPolicy::ViewType& partials_view)
          {
            const std::size_t n = view.species_variable_indices_.size();
            for (std::size_t k = 0; k < n; ++k)
            {
              const double mv = view.species_molar_volumes_[k];
              params_view.ForEachRowStrict(
                  [mv](const double& gmd, const double& gsd, double& dN)
                  {
                    const double ln_gsd = std::log(gsd);
                    const double V_s = (4.0 / 3.0) * std::numbers::pi * gmd * gmd * gmd * std::exp(4.5 * ln_gsd * ln_gsd);
                    dN = mv / V_s;
                  },
                  params_view.GetConstColumnView(view.gmd_parameter_index_),
                  params_view.GetConstColumnView(view.gsd_parameter_index_),
                  partials_view.GetColumnView(k));
            }
          },
          state_parameters,
          state_variables,
          result,
          partials)(state_parameters, state_variables, result, partials);
    }

   private:
    std::size_t gmd_parameter_index_ = 0;
    std::size_t gsd_parameter_index_ = 0;
    Vector<std::size_t> species_variable_indices_{};
    Vector<double> species_molar_volumes_{};
    std::vector<std::size_t> dependent_variable_indices_{};
  };

  /// @brief Single moment log-normal particle size distribution representation
  /// @details Represents a single moment log-normal distribution for aerosol or cloud particle size distributions.
  ///          Characterized by a geometric mean radius and geometric standard deviation.
  class SingleMomentMode
  {
   public:
    SingleMomentMode() = delete;

    SingleMomentMode(const std::string& prefix, const std::vector<micm::Phase>& phases)
        : prefix_(prefix),
          phases_(phases),
          default_geometric_mean_radius_(0.0),
          default_geometric_standard_deviation_(1.0)
    {
    }

    SingleMomentMode(
        const std::string& prefix,
        const std::vector<micm::Phase>& phases,
        const double geometric_mean_radius,
        const double geometric_standard_deviation)
        : prefix_(prefix),
          phases_(phases),
          default_geometric_mean_radius_(geometric_mean_radius),
          default_geometric_standard_deviation_(geometric_standard_deviation)
    {
    }

    std::tuple<std::size_t, std::size_t> StateSize() const
    {
      std::size_t size = 0;
      for (const auto& phase : phases_)
      {
        size += phase.StateSize();
      }
      return { size, 2 };  // Two parameters: geometric mean radius and geometric standard deviation
    }

    std::set<std::string> StateVariableNames() const
    {
      std::set<std::string> names;
      for (const auto& phase : phases_)
      {
        for (const auto& species : phase.UniqueNames())
        {
          names.insert(prefix_ + "." + species);
        }
      }
      return names;
    }

    std::set<std::string> StateParameterNames() const
    {
      std::set<std::string> names;
      names.insert(prefix_ + ".GEOMETRIC_MEAN_RADIUS");
      names.insert(prefix_ + ".GEOMETRIC_STANDARD_DEVIATION");
      return names;
    }

    std::string Species(const micm::Phase& phase, const micm::Species& species) const
    {
      return prefix_ + "." + phase.name_ + "." + species.name_;
    }

    std::map<std::string, double> DefaultParameters() const
    {
      return { { prefix_ + ".GEOMETRIC_MEAN_RADIUS", default_geometric_mean_radius_ },
               { prefix_ + ".GEOMETRIC_STANDARD_DEVIATION", default_geometric_standard_deviation_ } };
    }

    std::string GeometricMeanRadius() const
    {
      return prefix_ + ".GEOMETRIC_MEAN_RADIUS";
    }

    std::string GeometricStandardDeviation() const
    {
      return prefix_ + ".GEOMETRIC_STANDARD_DEVIATION";
    }

    void SetDefaultParameters(auto& state) const
    {
      auto gmd_it = state.custom_rate_parameter_map_.find(GeometricMeanRadius());
      if (gmd_it == state.custom_rate_parameter_map_.end())
      {
        throw MiamException(
            MIAM_ERROR_CATEGORY_CONFIGURATION,
            MIAM_CONFIGURATION_MISSING_STATE_PARAMETER,
            "SingleMomentMode::SetDefaultParameters: GEOMETRIC_MEAN_RADIUS parameter not found in state for " + prefix_);
      }
      auto gsd_it = state.custom_rate_parameter_map_.find(GeometricStandardDeviation());
      if (gsd_it == state.custom_rate_parameter_map_.end())
      {
        throw MiamException(
            MIAM_ERROR_CATEGORY_CONFIGURATION,
            MIAM_CONFIGURATION_MISSING_STATE_PARAMETER,
            "SingleMomentMode::SetDefaultParameters: GEOMETRIC_STANDARD_DEVIATION parameter not found in state for " +
                prefix_);
      }
      for (std::size_t i_cell = 0; i_cell < state.variables_.NumRows(); ++i_cell)
      {
        state.custom_rate_parameters_[i_cell][gmd_it->second] = default_geometric_mean_radius_;
        state.custom_rate_parameters_[i_cell][gsd_it->second] = default_geometric_standard_deviation_;
      }
    }

    std::map<std::string, std::size_t> NumPhaseInstances() const
    {
      std::map<std::string, std::size_t> num_instances;
      for (const auto& phase : phases_)
      {
        num_instances[phase.name_] = 1;  // Single moment representation has one instance per phase
      }
      return num_instances;
    }

    /// @brief Returns a map of phase names to sets of state variable prefixes associated with that phase
    ///        The prefix does not include the phase or species names, and each prefix must be unique across all
    ///        representations.
    /// @return Map of phase names to sets of state variable prefixes
    std::map<std::string, std::set<std::string>> PhaseStatePrefixes() const
    {
      std::map<std::string, std::set<std::string>> phase_prefixes;
      for (const auto& phase : phases_)
      {
        phase_prefixes[phase.name_].insert(prefix_);
      }
      return phase_prefixes;
    }

    /// @brief The set of concrete aerosol-property descriptor types this representation can produce.
    template<typename DenseMatrixPolicy>
    using DescriptorVariant = std::variant<
        SingleMomentModeEffectiveRadiusDescriptor<DenseMatrixPolicy>,
        SingleMomentModeNumberConcentrationDescriptor<DenseMatrixPolicy>,
        PhaseVolumeFractionDescriptor<DenseMatrixPolicy>>;

    /// @brief Returns a descriptor for the requested aerosol property.
    /// @tparam DenseMatrixPolicy The dense matrix type used for state data
    /// @param property The aerosol property to provide
    /// @param state_parameter_indices Map of parameter names to column indices
    /// @param state_variable_indices Map of variable names to column indices
    /// @param target_phase_name Phase name for PhaseVolumeFraction (required if multi-phase)
    template<typename DenseMatrixPolicy>
    DescriptorVariant<DenseMatrixPolicy> GetPropertyDescriptor(
        AerosolProperty property,
        const auto& state_parameter_indices,
        const auto& state_variable_indices,
        const std::string& target_phase_name = "") const
    {
      switch (property)
      {
        case AerosolProperty::EffectiveRadius:
        {
          return SingleMomentModeEffectiveRadiusDescriptor<DenseMatrixPolicy>{
            state_parameter_indices.at(GeometricMeanRadius()), state_parameter_indices.at(GeometricStandardDeviation())
          };
        }
        case AerosolProperty::NumberConcentration:
        {
          std::vector<std::size_t> species_indices;
          std::vector<double> molar_volumes;
          for (const auto& phase : phases_)
            for (const auto& ps : phase.phase_species_)
              if (!ps.species_.IsParameterized())
              {
                species_indices.push_back(state_variable_indices.at(prefix_ + "." + phase.name_ + "." + ps.species_.name_));
                molar_volumes.push_back(
                    ps.species_.GetProperty<double>("molecular weight [kg mol-1]") /
                    ps.species_.GetProperty<double>("density [kg m-3]"));
              }
          return SingleMomentModeNumberConcentrationDescriptor<DenseMatrixPolicy>{
            state_parameter_indices.at(GeometricMeanRadius()),
            state_parameter_indices.at(GeometricStandardDeviation()),
            std::move(species_indices),
            std::move(molar_volumes)
          };
        }
        case AerosolProperty::PhaseVolumeFraction:
        {
          if (phases_.size() == 1)
            return PhaseVolumeFractionDescriptor<DenseMatrixPolicy>{ {}, {}, 0, {} };
          if (target_phase_name.empty())
            throw MiamException(
                MIAM_ERROR_CATEGORY_CONFIGURATION,
                MIAM_CONFIGURATION_PHASE_NAME_REQUIRED,
                "SingleMomentMode::GetPropertyDescriptor: target_phase_name required for PhaseVolumeFraction with "
                "multiple phases");
          std::vector<std::size_t> all_species;
          std::vector<double> all_mw_over_rho;
          std::size_t phase_count = 0;
          for (const auto& phase : phases_)
            if (phase.name_ == target_phase_name)
            {
              for (const auto& ps : phase.phase_species_)
                if (!ps.species_.IsParameterized())
                {
                  all_species.push_back(state_variable_indices.at(prefix_ + "." + phase.name_ + "." + ps.species_.name_));
                  all_mw_over_rho.push_back(
                      ps.species_.GetProperty<double>("molecular weight [kg mol-1]") /
                      ps.species_.GetProperty<double>("density [kg m-3]"));
                }
              phase_count = all_species.size();
              break;
            }
          for (const auto& phase : phases_)
          {
            if (phase.name_ == target_phase_name)
              continue;
            for (const auto& ps : phase.phase_species_)
              if (!ps.species_.IsParameterized())
              {
                all_species.push_back(state_variable_indices.at(prefix_ + "." + phase.name_ + "." + ps.species_.name_));
                all_mw_over_rho.push_back(
                    ps.species_.GetProperty<double>("molecular weight [kg mol-1]") /
                    ps.species_.GetProperty<double>("density [kg m-3]"));
              }
          }
          std::vector<std::size_t> deps = all_species;
          return PhaseVolumeFractionDescriptor<DenseMatrixPolicy>{
            std::move(all_species), std::move(all_mw_over_rho), phase_count, std::move(deps)
          };
        }
        default:
          throw MiamException(
              MIAM_ERROR_CATEGORY_CONFIGURATION,
              MIAM_CONFIGURATION_UNSUPPORTED_PROPERTY,
              "SingleMomentMode::GetPropertyDescriptor: unsupported AerosolProperty");
      }
    }

   private:
    std::string prefix_;                           // State name prefix to apply to mode properties
    std::vector<micm::Phase> phases_;              // Phases associated with the mode
    double default_geometric_mean_radius_;         // Default geometric mean radius
    double default_geometric_standard_deviation_;  // Default geometric standard deviation
  };
}  // namespace miam
