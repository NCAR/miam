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
  /// @brief Effective radius for a uniform-section representation: `r_eff = 0.5 * (r_min + r_max)`.
  template<typename DenseMatrixPolicy>
  class UniformSectionEffectiveRadiusDescriptor
  {
   public:
    struct View
    {
      std::size_t rmin_parameter_index_ = 0;
      std::size_t rmax_parameter_index_ = 0;
    };

    UniformSectionEffectiveRadiusDescriptor() = default;

    UniformSectionEffectiveRadiusDescriptor(std::size_t rmin_parameter_index, std::size_t rmax_parameter_index)
        : view_{ rmin_parameter_index, rmax_parameter_index }
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
                [](const double& r_min, const double& r_max, double& r_eff) { r_eff = 0.5 * (r_min + r_max); },
                params_view.GetConstColumnView(view.rmin_parameter_index_),
                params_view.GetConstColumnView(view.rmax_parameter_index_),
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

  /// @brief Number concentration for a uniform-section representation:
  ///        `N = (Σ_k [species_k] · mv_k) / V_s`, with `V_s = (4/3)·π·r_eff³` and `r_eff = 0.5·(r_min + r_max)`.
  template<typename DenseMatrixPolicy>
  class UniformSectionNumberConcentrationDescriptor
  {
   public:
    template<typename U>
    using Vector = typename DenseMatrixPolicy::template VectorType<U>;

    struct View
    {
      std::size_t rmin_parameter_index_ = 0;
      std::size_t rmax_parameter_index_ = 0;
      typename Vector<std::size_t>::ConstViewType species_variable_indices_{};
      typename Vector<double>::ConstViewType species_molar_volumes_{};
    };

    UniformSectionNumberConcentrationDescriptor() = default;

    UniformSectionNumberConcentrationDescriptor(
        std::size_t rmin_parameter_index,
        std::size_t rmax_parameter_index,
        std::vector<std::size_t> species_variable_indices,
        std::vector<double> species_molar_volumes)
        : rmin_parameter_index_(rmin_parameter_index),
          rmax_parameter_index_(rmax_parameter_index),
          species_variable_indices_(std::move(species_variable_indices)),
          species_molar_volumes_(std::move(species_molar_volumes)),
          dependent_variable_indices_(species_variable_indices_.begin(), species_variable_indices_.end())
    {
      species_variable_indices_.CopyToDevice();
      species_molar_volumes_.CopyToDevice();
    }

    View GetView() const
    {
      return { rmin_parameter_index_, rmax_parameter_index_, species_variable_indices_.GetView(), species_molar_volumes_.GetView() };
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
                [](const double& r_min, const double& r_max, double& N_out)
                {
                  const double r_eff = 0.5 * (r_min + r_max);
                  const double V_s = (4.0 / 3.0) * std::numbers::pi * r_eff * r_eff * r_eff;
                  N_out /= V_s;
                },
                params_view.GetConstColumnView(view.rmin_parameter_index_),
                params_view.GetConstColumnView(view.rmax_parameter_index_),
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
                  [mv](const double& r_min, const double& r_max, double& dN)
                  {
                    const double r_eff = 0.5 * (r_min + r_max);
                    const double V_s = (4.0 / 3.0) * std::numbers::pi * r_eff * r_eff * r_eff;
                    dN = mv / V_s;
                  },
                  params_view.GetConstColumnView(view.rmin_parameter_index_),
                  params_view.GetConstColumnView(view.rmax_parameter_index_),
                  partials_view.GetColumnView(k));
            }
          },
          state_parameters,
          state_variables,
          result,
          partials)(state_parameters, state_variables, result, partials);
    }

   private:
    std::size_t rmin_parameter_index_ = 0;
    std::size_t rmax_parameter_index_ = 0;
    Vector<std::size_t> species_variable_indices_{};
    Vector<double> species_molar_volumes_{};
    std::vector<std::size_t> dependent_variable_indices_{};
  };

  /// @brief Sectional particle size distribution representation with uniform sections
  /// @details Represents a sectional distribution with uniform sections for aerosol or cloud particle size distributions.
  ///          Each section is characterized by a fixed size range and variable total volume. Number concentrations are
  ///          derived from the total volume and section size.
  class UniformSection
  {
   public:
    UniformSection() = delete;

    UniformSection(const std::string& prefix, const std::vector<micm::Phase>& phases)
        : prefix_(prefix),
          phases_(phases),
          default_min_radius_(0.0),
          default_max_radius_(0.0)
    {
    }

    UniformSection(
        const std::string& prefix,
        const std::vector<micm::Phase>& phases,
        const double minimum_radius,
        const double maximum_radius)
        : prefix_(prefix),
          phases_(phases),
          default_min_radius_(minimum_radius),
          default_max_radius_(maximum_radius)
    {
    }

    std::tuple<std::size_t, std::size_t> StateSize() const
    {
      std::size_t size = 0;
      for (const auto& phase : phases_)
      {
        size += phase.StateSize();
      }
      return { size, 2 };  // Two parameters: min and max radius
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
      names.insert(prefix_ + ".MIN_RADIUS");
      names.insert(prefix_ + ".MAX_RADIUS");
      return names;
    }

    std::string Species(const micm::Phase& phase, const micm::Species& species) const
    {
      return prefix_ + "." + phase.name_ + "." + species.name_;
    }

    std::map<std::string, double> DefaultParameters() const
    {
      return { { prefix_ + ".MIN_RADIUS", default_min_radius_ }, { prefix_ + ".MAX_RADIUS", default_max_radius_ } };
    }

    std::string MinRadius() const
    {
      return prefix_ + ".MIN_RADIUS";
    }

    std::string MaxRadius() const
    {
      return prefix_ + ".MAX_RADIUS";
    }

    void SetDefaultParameters(auto& state) const
    {
      auto min_radius_it = state.custom_rate_parameter_map_.find(MinRadius());
      if (min_radius_it == state.custom_rate_parameter_map_.end())
      {
        throw MiamException(
            MIAM_ERROR_CATEGORY_CONFIGURATION,
            MIAM_CONFIGURATION_MISSING_STATE_PARAMETER,
            "UniformSection::SetDefaultParameters: MIN_RADIUS parameter not found in state for " + prefix_);
      }
      auto max_radius_it = state.custom_rate_parameter_map_.find(MaxRadius());
      if (max_radius_it == state.custom_rate_parameter_map_.end())
      {
        throw MiamException(
            MIAM_ERROR_CATEGORY_CONFIGURATION,
            MIAM_CONFIGURATION_MISSING_STATE_PARAMETER,
            "UniformSection::SetDefaultParameters: MAX_RADIUS parameter not found in state for " + prefix_);
      }
      for (std::size_t cell = 0; cell < state.variables_.NumRows(); ++cell)
      {
        state.custom_rate_parameters_[cell][min_radius_it->second] = default_min_radius_;
        state.custom_rate_parameters_[cell][max_radius_it->second] = default_max_radius_;
      }
    }

    std::map<std::string, std::size_t> NumPhaseInstances() const
    {
      std::map<std::string, std::size_t> num_instances;
      for (const auto& phase : phases_)
      {
        num_instances[phase.name_] = 1;  // Uniform section representation has one instance per phase
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
        UniformSectionEffectiveRadiusDescriptor<DenseMatrixPolicy>,
        UniformSectionNumberConcentrationDescriptor<DenseMatrixPolicy>,
        PhaseVolumeFractionDescriptor<DenseMatrixPolicy>>;

    /// @brief Returns a descriptor for the requested aerosol property.
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
          return UniformSectionEffectiveRadiusDescriptor<DenseMatrixPolicy>{
            state_parameter_indices.at(MinRadius()), state_parameter_indices.at(MaxRadius())
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
          return UniformSectionNumberConcentrationDescriptor<DenseMatrixPolicy>{
            state_parameter_indices.at(MinRadius()),
            state_parameter_indices.at(MaxRadius()),
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
                "UniformSection::GetPropertyDescriptor: target_phase_name required for PhaseVolumeFraction "
                "with multiple phases");
          std::vector<std::size_t> all_species;
          std::vector<double> all_molar_volumes;
          std::size_t phase_count = 0;
          for (const auto& phase : phases_)
            if (phase.name_ == target_phase_name)
            {
              for (const auto& ps : phase.phase_species_)
                if (!ps.species_.IsParameterized())
                {
                  all_species.push_back(state_variable_indices.at(prefix_ + "." + phase.name_ + "." + ps.species_.name_));
                  all_molar_volumes.push_back(
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
                all_molar_volumes.push_back(
                    ps.species_.GetProperty<double>("molecular weight [kg mol-1]") /
                    ps.species_.GetProperty<double>("density [kg m-3]"));
              }
          }
          std::vector<std::size_t> deps = all_species;
          return PhaseVolumeFractionDescriptor<DenseMatrixPolicy>{
            std::move(all_species), std::move(all_molar_volumes), phase_count, std::move(deps)
          };
        }
        default:
          throw MiamException(
              MIAM_ERROR_CATEGORY_CONFIGURATION,
              MIAM_CONFIGURATION_UNSUPPORTED_PROPERTY,
              "UniformSection::GetPropertyDescriptor: unsupported AerosolProperty");
      }
    }

   private:
    std::string prefix_;               // State name prefix to apply to section properties
    std::vector<micm::Phase> phases_;  // Phases associated with the section
    double default_min_radius_;        // Minimum radius of the section
    double default_max_radius_;        // Maximum radius of the section
  };
}  // namespace miam