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
  /// @brief Effective radius for a two-moment log-normal mode:
  ///        `r_eff = r_mean * exp(2.5·ln²(GSD))`, `r_mean = cbrt(3·V_mean / (4π))`, `V_mean = V_total / N`.
  /// @details Depends on GSD, the mode's number concentration state variable, and the aqueous-species
  ///          concentrations attached to the mode.
  template<typename DenseMatrixPolicy>
  class TwoMomentModeEffectiveRadiusDescriptor
  {
   public:
    template<typename U>
    using Vector = typename DenseMatrixPolicy::template VectorType<U>;

    struct View
    {
      std::size_t gsd_parameter_index_ = 0;
      std::size_t nc_variable_index_ = 0;
      typename Vector<std::size_t>::ConstViewType species_variable_indices_{};
      typename Vector<double>::ConstViewType species_molar_volumes_{};
    };

    TwoMomentModeEffectiveRadiusDescriptor() = default;

    TwoMomentModeEffectiveRadiusDescriptor(
        std::size_t gsd_parameter_index,
        std::size_t nc_variable_index,
        std::vector<std::size_t> species_variable_indices,
        std::vector<double> species_molar_volumes)
        : gsd_parameter_index_(gsd_parameter_index),
          nc_variable_index_(nc_variable_index),
          species_variable_indices_(std::move(species_variable_indices)),
          species_molar_volumes_(std::move(species_molar_volumes))
    {
      dependent_variable_indices_.assign(species_variable_indices_.begin(), species_variable_indices_.end());
      dependent_variable_indices_.push_back(nc_variable_index_);
      species_variable_indices_.CopyToDevice();
      species_molar_volumes_.CopyToDevice();
    }

    View GetView() const
    {
      return { gsd_parameter_index_, nc_variable_index_, species_variable_indices_.GetView(), species_molar_volumes_.GetView() };
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
            auto r = result_view.GetColumnView(0);
            params_view.ForEachRowStrict([](double& v) { v = 0.0; }, r);
            const std::size_t n = view.species_variable_indices_.size();
            for (std::size_t k = 0; k < n; ++k)
            {
              const double mv = view.species_molar_volumes_[k];
              params_view.ForEachRowStrict(
                  [mv](const double& c, double& V) { V += c * mv; },
                  vars_view.GetConstColumnView(view.species_variable_indices_[k]),
                  r);
            }
            params_view.ForEachRowStrict(
                [](const double& gsd, const double& nc, double& r_eff)
                {
                  const double ln_gsd = std::log(gsd);
                  const double V_mean = r_eff / nc;
                  const double r_mean = std::cbrt(3.0 * V_mean / (4.0 * std::numbers::pi));
                  r_eff = r_mean * std::exp(2.5 * ln_gsd * ln_gsd);
                },
                params_view.GetConstColumnView(view.gsd_parameter_index_),
                vars_view.GetConstColumnView(view.nc_variable_index_),
                r);
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
      const auto view = GetView();
      DenseMatrixPolicy::Function(
          MICM_LAMBDA(const typename DenseMatrixPolicy::ConstViewType& params_view, const typename DenseMatrixPolicy::ConstViewType& vars_view, const typename DenseMatrixPolicy::ViewType& result_view, const typename DenseMatrixPolicy::ViewType& partials_view)
          {
            auto V_total = result_view.GetRowVariable();
            params_view.ForEachRowStrict([](double& v) { v = 0.0; }, V_total);
            const std::size_t n = view.species_variable_indices_.size();
            for (std::size_t k = 0; k < n; ++k)
            {
              const double mv = view.species_molar_volumes_[k];
              params_view.ForEachRowStrict(
                  [mv](const double& c, double& V) { V += c * mv; },
                  vars_view.GetConstColumnView(view.species_variable_indices_[k]),
                  V_total);
            }
            for (std::size_t k = 0; k < n; ++k)
            {
              const double mv = view.species_molar_volumes_[k];
              params_view.ForEachRowStrict(
                  [mv](const double& gsd, const double& nc, const double& V_total_row, double& dr)
                  {
                    const double ln_gsd = std::log(gsd);
                    const double V_mean = V_total_row / nc;
                    const double r_mean = std::cbrt(3.0 * V_mean / (4.0 * std::numbers::pi));
                    const double r_eff = r_mean * std::exp(2.5 * ln_gsd * ln_gsd);
                    dr = r_eff * mv / (3.0 * V_total_row);
                  },
                  params_view.GetConstColumnView(view.gsd_parameter_index_),
                  vars_view.GetConstColumnView(view.nc_variable_index_),
                  V_total,
                  partials_view.GetColumnView(k));
            }
            params_view.ForEachRowStrict(
                [](const double& gsd, const double& nc, const double& V_total_row, double& dr_dN)
                {
                  const double ln_gsd = std::log(gsd);
                  const double V_mean = V_total_row / nc;
                  const double r_mean = std::cbrt(3.0 * V_mean / (4.0 * std::numbers::pi));
                  const double r_eff = r_mean * std::exp(2.5 * ln_gsd * ln_gsd);
                  dr_dN = -r_eff / (3.0 * nc);
                },
                params_view.GetConstColumnView(view.gsd_parameter_index_),
                vars_view.GetConstColumnView(view.nc_variable_index_),
                V_total,
                partials_view.GetColumnView(n));
            params_view.ForEachRowStrict(
                [](const double& gsd, const double& nc, const double& V_total_row, double& r_out)
                {
                  const double ln_gsd = std::log(gsd);
                  const double V_mean = V_total_row / nc;
                  const double r_mean = std::cbrt(3.0 * V_mean / (4.0 * std::numbers::pi));
                  r_out = r_mean * std::exp(2.5 * ln_gsd * ln_gsd);
                },
                params_view.GetConstColumnView(view.gsd_parameter_index_),
                vars_view.GetConstColumnView(view.nc_variable_index_),
                V_total,
                result_view.GetColumnView(0));
          },
          state_parameters,
          state_variables,
          result,
          partials)(state_parameters, state_variables, result, partials);
    }

   private:
    std::size_t gsd_parameter_index_ = 0;
    std::size_t nc_variable_index_ = 0;
    Vector<std::size_t> species_variable_indices_{};
    Vector<double> species_molar_volumes_{};
    std::vector<std::size_t> dependent_variable_indices_{};
  };

  /// @brief Number concentration for a two-moment log-normal mode: `N = state_variables[NUMBER_CONCENTRATION]`.
  /// @details Trivially just returns the stored NC state variable and reports `dN/dN = 1`.
  template<typename DenseMatrixPolicy>
  class TwoMomentModeNumberConcentrationDescriptor
  {
   public:
    struct View
    {
      std::size_t nc_variable_index_ = 0;
    };

    TwoMomentModeNumberConcentrationDescriptor() = default;

    explicit TwoMomentModeNumberConcentrationDescriptor(std::size_t nc_variable_index)
        : view_{ nc_variable_index },
          dependent_variable_indices_{ nc_variable_index }
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
          MICM_LAMBDA(const typename DenseMatrixPolicy::ConstViewType& params_view, const typename DenseMatrixPolicy::ConstViewType& vars_view, const typename DenseMatrixPolicy::ViewType& result_view)
          {
            params_view.ForEachRowStrict(
                [](const double& nc, double& N) { N = nc; },
                vars_view.GetConstColumnView(view.nc_variable_index_),
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
        DenseMatrixPolicy& partials) const
    {
      const auto view = view_;
      DenseMatrixPolicy::Function(
          MICM_LAMBDA(const typename DenseMatrixPolicy::ConstViewType& params_view, const typename DenseMatrixPolicy::ConstViewType& vars_view, const typename DenseMatrixPolicy::ViewType& result_view, const typename DenseMatrixPolicy::ViewType& partials_view)
          {
            params_view.ForEachRowStrict(
                [](const double& nc, double& N, double& dN_dN)
                {
                  N = nc;
                  dN_dN = 1.0;
                },
                vars_view.GetConstColumnView(view.nc_variable_index_),
                result_view.GetColumnView(0),
                partials_view.GetColumnView(0));
          },
          state_parameters,
          state_variables,
          result,
          partials)(state_parameters, state_variables, result, partials);
    }

   private:
    View view_{};
    std::vector<std::size_t> dependent_variable_indices_{};
  };

  /// @brief Two moment log-normal particle size distribution representation
  /// @details Represents a two moment log-normal distribution for aerosol or cloud particle size distributions.
  ///          Characterized by number concentration, geometric mean radius, and geometric standard deviation.
  class TwoMomentMode
  {
   public:
    TwoMomentMode() = delete;

    TwoMomentMode(const std::string& prefix, const std::vector<micm::Phase>& phases)
        : prefix_(prefix),
          phases_(phases),
          default_geometric_standard_deviation_(1.0)
    {
    }

    TwoMomentMode(
        const std::string& prefix,
        const std::vector<micm::Phase>& phases,
        const double geometric_standard_deviation)
        : prefix_(prefix),
          phases_(phases),
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
      size++;              // Number concentration
      return { size, 1 };  // One parameter: geometric standard deviation
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
      names.insert(prefix_ + ".NUMBER_CONCENTRATION");
      return names;
    }

    std::set<std::string> StateParameterNames() const
    {
      std::set<std::string> names;
      names.insert(prefix_ + ".GEOMETRIC_STANDARD_DEVIATION");
      return names;
    }

    std::string Species(const micm::Phase& phase, const micm::Species& species) const
    {
      return prefix_ + "." + phase.name_ + "." + species.name_;
    }

    std::map<std::string, double> DefaultParameters() const
    {
      return { { prefix_ + ".GEOMETRIC_STANDARD_DEVIATION", default_geometric_standard_deviation_ } };
    }

    std::string NumberConcentration() const
    {
      return prefix_ + ".NUMBER_CONCENTRATION";
    }

    std::string GeometricStandardDeviation() const
    {
      return prefix_ + ".GEOMETRIC_STANDARD_DEVIATION";
    }

    void SetDefaultParameters(auto& state) const
    {
      auto gsd_it = state.custom_rate_parameter_map_.find(GeometricStandardDeviation());
      if (gsd_it == state.custom_rate_parameter_map_.end())
      {
        throw MiamException(
            MIAM_ERROR_CATEGORY_CONFIGURATION,
            MIAM_CONFIGURATION_MISSING_STATE_PARAMETER,
            "TwoMomentMode::SetDefaultParameters: Geometric standard deviation parameter not found in state.");
      }
      for (std::size_t i_cell = 0; i_cell < state.variables_.NumRows(); ++i_cell)
      {
        state.custom_rate_parameters_[i_cell][gsd_it->second] = default_geometric_standard_deviation_;
      }
    }

    std::map<std::string, std::size_t> NumPhaseInstances() const
    {
      std::map<std::string, std::size_t> num_instances;
      for (const auto& phase : phases_)
      {
        num_instances[phase.name_] = 1;  // Two moment representation has one instance per phase
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
        TwoMomentModeEffectiveRadiusDescriptor<DenseMatrixPolicy>,
        TwoMomentModeNumberConcentrationDescriptor<DenseMatrixPolicy>,
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
          return TwoMomentModeEffectiveRadiusDescriptor<DenseMatrixPolicy>{
            state_parameter_indices.at(GeometricStandardDeviation()),
            state_variable_indices.at(NumberConcentration()),
            std::move(species_indices),
            std::move(molar_volumes)
          };
        }
        case AerosolProperty::NumberConcentration:
        {
          return TwoMomentModeNumberConcentrationDescriptor<DenseMatrixPolicy>{
            state_variable_indices.at(NumberConcentration())
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
                "TwoMomentMode::GetPropertyDescriptor: target_phase_name required for PhaseVolumeFraction "
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
              "TwoMomentMode::GetPropertyDescriptor: unsupported AerosolProperty");
      }
    }

   private:
    std::string prefix_;                           // State name prefix to apply to shape properties
    std::vector<micm::Phase> phases_;              // Phases associated with the mode
    double default_geometric_standard_deviation_;  // Default geometric standard deviation
  };
}  // namespace miam