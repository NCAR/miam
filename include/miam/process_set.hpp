// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/aerosol_property.hpp>

#include <micm/system/conditions.hpp>

#include <functional>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace miam
{
  /// @brief Type-erased wrapper around any MIAM process
  /// @details Follows the same pattern MICM uses for ExternalModelProcessSet. Allows
  ///          the Model to store all process types in a single collection and iterate
  ///          them generically through a common interface.
  /// @tparam DenseMatrixPolicy The dense matrix type used for state data
  /// @tparam SparseMatrixPolicy The sparse matrix type used for Jacobian data
  template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
  struct MiamProcessSet
  {
    using PhaseMap = std::map<std::string, std::set<std::string>>;
    using IndexMap = std::unordered_map<std::string, std::size_t>;
    using ProviderMap = std::map<std::string, std::map<AerosolProperty, AerosolPropertyProvider<DenseMatrixPolicy>>>;

    std::function<std::set<std::string>(const PhaseMap&)> process_parameter_names_;
    std::function<std::set<std::string>(const PhaseMap&)> species_used_;
    std::function<std::map<std::string, std::vector<AerosolProperty>>()> required_aerosol_properties_;
    std::function<std::set<std::pair<std::size_t, std::size_t>>(const PhaseMap&, const IndexMap&, const ProviderMap&)>
        non_zero_jacobian_elements_;
    std::function<
        std::function<void(const std::vector<micm::Conditions>&, DenseMatrixPolicy&)>(const PhaseMap&, const IndexMap&)>
        update_state_parameters_function_;
    std::function<std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)>(
        const PhaseMap&,
        const IndexMap&,
        const IndexMap&,
        ProviderMap)>
        get_forcing_function_;
    std::function<std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)>(
        const PhaseMap&,
        const IndexMap&,
        const IndexMap&,
        const SparseMatrixPolicy&,
        ProviderMap)>
        get_jacobian_function_;

    /// @brief Construct a MiamProcessSet from any process type that satisfies the common interface
    /// @tparam ProcessType The concrete process type
    /// @param process The process instance to wrap
    template<typename ProcessType>
    MiamProcessSet(ProcessType&& process)
    {
      auto shared = std::make_shared<std::decay_t<ProcessType>>(std::forward<ProcessType>(process));

      process_parameter_names_ = [shared](const PhaseMap& pp) { return shared->ProcessParameterNames(pp); };

      species_used_ = [shared](const PhaseMap& pp) { return shared->SpeciesUsed(pp); };

      required_aerosol_properties_ = [shared]() { return shared->RequiredAerosolProperties(); };

      non_zero_jacobian_elements_ = [shared](const PhaseMap& pp, const IndexMap& vi, const ProviderMap& prov)
      { return shared->NonZeroJacobianElements(pp, vi, prov); };

      update_state_parameters_function_ = [shared](const PhaseMap& pp, const IndexMap& pi)
      { return shared->template UpdateStateParametersFunction<DenseMatrixPolicy>(pp, pi); };

      get_forcing_function_ = [shared](const PhaseMap& pp, const IndexMap& pi, const IndexMap& vi, ProviderMap prov)
      { return shared->template ForcingFunction<DenseMatrixPolicy>(pp, pi, vi, std::move(prov)); };

      get_jacobian_function_ =
          [shared](
              const PhaseMap& pp, const IndexMap& pi, const IndexMap& vi, const SparseMatrixPolicy& jac, ProviderMap prov)
      { return shared->template JacobianFunction<DenseMatrixPolicy, SparseMatrixPolicy>(pp, pi, vi, jac, std::move(prov)); };
    }
  };
}  // namespace miam
