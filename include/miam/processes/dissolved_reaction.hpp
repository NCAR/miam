// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/processes/constants/rate_constant.hpp>
#include <miam/representations/aerosol_property.hpp>
#include <miam/util/error.hpp>
#include <miam/util/miam_exception.hpp>
#include <miam/util/uuid.hpp>

#include <micm/system/conditions.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/types.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <variant>
#include <vector>

namespace miam
{
  /// @brief A dissolved (irreversible) reaction
  /// @details Dissolved reactions involve reactants and products in solution, and
  ///          are characterized by a single rate constant.
  ///
  ///          The reaction takes the form:
  ///          \f$ \mathrm{Reactants} \rightarrow \mathrm{Products} \f$
  ///
  ///          with rate constant \f$ k \f$.
  ///
  ///          The rate expression (regularized to prevent singularity as [S] → 0) is:
  ///          \f[
  ///            r = k \cdot \frac{[S]}{([S] + \delta)^{n_r}} \prod_i [R_i]
  ///          \f]
  ///
  ///          where [S] is the solvent concentration (mol m⁻³ air),
  ///          \f$ n_r \f$ is the number of reactants, and \f$ \delta \f$ =
  ///          `solvent_floor_` is a small regularization constant (default 10⁻²⁰).
  ///          In the limit \f$ \delta \to 0 \f$ this reduces to the clean form:
  ///          \f$ r = k / [S]^{n_r - 1} \cdot \prod_i [R_i] \f$.
  ///
  ///          The solvent-normalization factor absorbs concentration dimensions,
  ///          so \f$ k \f$ always has units of s⁻¹. To convert from a literature
  ///          rate constant \f$ k_{lit} \f$ in molar units:
  ///          \f$ k = k_{lit} \times c_{H_2O}^{n_r - 1} \f$
  ///          where \f$ c_{H_2O} = 55.51 \f$ mol/L (1000 g/L ÷ 18.015 g/mol).
  ///
  ///          Partial derivatives used in the Jacobian:
  ///          \f[
  ///            \frac{\partial r}{\partial [R_j]} = k \cdot \frac{[S]}{([S]+\delta)^{n_r}}
  ///              \prod_{i \neq j} [R_i]
  ///          \f]
  ///          \f[
  ///            \frac{\partial r}{\partial [S]} = k \cdot
  ///              \frac{\delta + (1 - n_r)[S]}{([S]+\delta)^{n_r+1}} \prod_i [R_i]
  ///          \f]
  class DissolvedReaction
  {
   public:
    RateConstant rate_constant_;  ///< Rate constant
    std::vector<micm::Species> reactants_;                                     ///< Reactant species
    std::vector<micm::Species> products_;                                      ///< Product species
    micm::Species solvent_;                                                    ///< Solvent species
    micm::Phase phase_;                                                        ///< Phase in which the reaction occurs
    std::string uuid_;                                                         ///< Unique identifier for the reaction
    double solvent_floor_{
      1.0e-20
    };  ///< Floor [mol m⁻³] added to [S] in ([S]+δ)^n denominator to prevent singularity as [S] → 0
    double min_halflife_{
      0.0
    };  ///< When > 0, caps the reaction rate so no reactant is depleted faster than this half-life [s]

    DissolvedReaction() = delete;

    /// @brief Constructor
    DissolvedReaction(
        RateConstant rate_constant,
        const std::vector<micm::Species>& reactants,
        const std::vector<micm::Species>& products,
        micm::Species solvent,
        micm::Phase phase,
        double solvent_floor = 1.0e-20,
        double min_halflife = 0.0)
        : rate_constant_(rate_constant),
          reactants_(reactants),
          products_(products),
          solvent_(solvent),
          phase_(phase),
          uuid_(GenerateUuid()),
          solvent_floor_(solvent_floor),
          min_halflife_(min_halflife)
    {
    }

    /// @brief Create a copy of this reaction with a new UUID
    /// @return A new DissolvedReaction with the same properties but a unique UUID
    DissolvedReaction CopyWithNewUuid() const
    {
      return DissolvedReaction(rate_constant_, reactants_, products_, solvent_, phase_, solvent_floor_, min_halflife_);
    }

    /// @brief Returns a set of unique parameter names for this process
    /// @param phase_prefixes Map of phase names to sets of state variable prefixes (prefix does not include phase or
    /// species names)
    /// @return Set of unique parameter names for this process
    std::set<std::string> ProcessParameterNames(const std::map<std::string, std::set<std::string>>& phase_prefixes) const
    {
      std::set<std::string> parameter_names;
      // The conditions are shared by the whole system, so we just need one value for the rate
      // constant. We can use the phase name and uuid to create a unique parameter name.
      parameter_names.insert(phase_.name_ + "." + uuid_ + ".k");
      return parameter_names;
    }

    /// @brief Returns participating species' unique state names
    /// @param phase_prefixes Map of phase names to sets of state variable prefixes (prefix does not include phase or
    /// species names)
    /// @return Set of unique state variable names for all species involved in the reaction
    std::set<std::string> SpeciesUsed(const std::map<std::string, std::set<std::string>>& phase_prefixes) const
    {
      std::set<std::string> species_names;
      auto phase_it = phase_prefixes.find(phase_.name_);
      if (phase_it != phase_prefixes.end())
      {
        const auto& prefixes = phase_it->second;
        for (const auto& prefix : prefixes)
        {
          for (const auto& reactant : reactants_)
          {
            species_names.insert(prefix + "." + phase_.name_ + "." + reactant.name_);
          }
          for (const auto& product : products_)
          {
            species_names.insert(prefix + "." + phase_.name_ + "." + product.name_);
          }
          species_names.insert(prefix + "." + phase_.name_ + "." + solvent_.name_);
        }
      }
      else
      {
        throw MiamException(
            MIAM_ERROR_CATEGORY_INTERNAL,
            MIAM_INTERNAL_MISSING_PHASE_PREFIX,
            "Internal Error: Phase " + phase_.name_ + " not found in phase_prefixes map for process " + uuid_);
      }
      return species_names;
    }

    /// @brief Returns the aerosol properties required by this process
    /// @details DissolvedReaction does not depend on aerosol properties.
    /// @return Empty map
    std::map<std::string, std::vector<AerosolProperty>> RequiredAerosolProperties() const
    {
      return {};
    }

    /// @brief Returns a set of Jacobian index pairs for this process
    /// @param phase_prefixes Map of phase names to sets of state variable prefixes (prefix does not include phase or
    /// species names)
    /// @param state_variable_indices Map of state variable names to their corresponding indices in the Jacobian
    /// @return Set of index pairs representing the Jacobian entries affected by this process (dependent, independent)
    std::set<std::pair<std::size_t, std::size_t>> NonZeroJacobianElements(
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const auto& state_variable_indices  // acts like std::unordered_map<std::string, std::size_t>
    ) const
    {
      std::set<std::pair<std::size_t, std::size_t>> jacobian_indices;

      StateVariableIndices variable_indices = GetStateVariableIndices(phase_prefixes, state_variable_indices);
      // Get pairs for each phase instance
      for (std::size_t i_phase = 0; i_phase < variable_indices.number_of_phase_instances_; ++i_phase)
      {
        // For a forward-only reaction, the rate depends only on reactants and solvent.
        // Products do NOT appear as independent variables.

        // Reactant contributions (each reactant is an independent variable)
        for (std::size_t r = 0; r < variable_indices.reactant_indices_.NumColumns(); ++r)
        {
          std::size_t independent_index = variable_indices.reactant_indices_[i_phase][r];
          // Each reactant affects all reactants
          for (std::size_t r2 = 0; r2 < variable_indices.reactant_indices_.NumColumns(); ++r2)
          {
            std::size_t dependent_index = variable_indices.reactant_indices_[i_phase][r2];
            jacobian_indices.insert({ dependent_index, independent_index });
          }
          // Each reactant affects all products
          for (std::size_t p = 0; p < variable_indices.product_indices_.NumColumns(); ++p)
          {
            std::size_t dependent_index = variable_indices.product_indices_[i_phase][p];
            jacobian_indices.insert({ dependent_index, independent_index });
          }
        }
        // Solvent contributions (affects all reactants and products)
        std::size_t independent_index = variable_indices.solvent_indices_[i_phase];
        for (std::size_t r = 0; r < variable_indices.reactant_indices_.NumColumns(); ++r)
        {
          std::size_t dependent_index = variable_indices.reactant_indices_[i_phase][r];
          jacobian_indices.insert({ dependent_index, independent_index });
        }
        for (std::size_t p = 0; p < variable_indices.product_indices_.NumColumns(); ++p)
        {
          std::size_t dependent_index = variable_indices.product_indices_[i_phase][p];
          jacobian_indices.insert({ dependent_index, independent_index });
        }
      }

      return jacobian_indices;
    }

    /// @brief Returns non-zero Jacobian elements (common interface overload accepting providers)
    /// @details Delegates to the existing two-argument version; providers are unused.
    template<typename DenseMatrixPolicy>
    std::set<std::pair<std::size_t, std::size_t>> NonZeroJacobianElements(
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const std::unordered_map<std::string, std::size_t>& state_variable_indices,
        const std::map<std::string, std::map<AerosolProperty, AerosolPropertyProvider<DenseMatrixPolicy>>>& /* providers */)
        const
    {
      return NonZeroJacobianElements(phase_prefixes, state_variable_indices);
    }

    /// @brief Returns a function that updates state parameters for this process
    /// @param phase_prefixes Map of phase names to sets of state variable prefixes (prefix does not include phase or
    /// species names)
    /// @param state_parameter_indices Map of state parameter names to their corresponding indices in the state parameter
    /// vector
    /// @return Function that updates state parameters for this process
    template<typename DenseMatrixPolicy>
    std::function<void(const typename DenseMatrixPolicy::template VectorType<micm::Conditions>&, DenseMatrixPolicy&)>
    UpdateStateParametersFunction(
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const auto& state_parameter_indices  // acts like std::unordered_map<std::string, std::size_t>
    ) const
    {
      std::string k_param = phase_.name_ + "." + uuid_ + ".k";
      if (state_parameter_indices.find(k_param) == state_parameter_indices.end())
      {
        throw MiamException(
            MIAM_ERROR_CATEGORY_INTERNAL,
            MIAM_INTERNAL_MISSING_STATE_PARAMETER,
            "Internal Error: UpdateStateParametersFunction: Rate constant parameter " + k_param +
                " not found in state_parameter_indices");
      }
      std::size_t k_index = state_parameter_indices.at(k_param);
      std::size_t num_params = state_parameter_indices.size();

      // NVCC forbids an extended __host__ __device__ lambda inside the generic std::visit lambda.
      return std::visit(
          [k_index, num_params](const auto& expr) { return MakeRateUpdateFn<DenseMatrixPolicy>(expr, k_index, num_params); },
          rate_constant_);
    }

    /// @brief Builds the rate constant update function for one expression type
    /// @details NVCC requires this function to be public.
    template<typename DenseMatrixPolicy, typename ExprT>
    static std::function<void(const typename DenseMatrixPolicy::template VectorType<micm::Conditions>&, DenseMatrixPolicy&)>
    MakeRateUpdateFn(const ExprT& expr, std::size_t k_index, std::size_t num_params)
    {
      const ExprT expr_copy = expr;
      // Set up dummy arguments to build the function
      DenseMatrixPolicy state_parameters{ 1, num_params, 0.0 };
      typename DenseMatrixPolicy::template VectorType<micm::Conditions> conditions_vector;

      // return a function that updates the rate constant parameter based on the current conditions
      return DenseMatrixPolicy::Function(
          MICM_LAMBDA(
              const typename DenseMatrixPolicy::template VectorType<micm::Conditions>::ConstViewType& conditions_view,
              const typename DenseMatrixPolicy::ViewType& params_view)
          {
            params_view.ForEachRowStrict(
                [expr_copy](const micm::Conditions& condition, micm::Real& parameter)
                { parameter = Calculate(expr_copy, condition); },
                conditions_view,
                params_view.GetColumnView(k_index));
          },
          conditions_vector,
          state_parameters);
    }

    /// @brief Returns a function that calculates the forcing terms for this process
    /// @param phase_prefixes Map of phase names to sets of state variable prefixes (prefix does not include phase or
    /// species names)
    /// @param state_parameter_indices Map of state parameter names to their corresponding indices in the state parameter
    /// vector
    /// @param state_variable_indices Map of state variable names to their corresponding indices in the state variable
    /// vector
    /// @return Function that calculates the forcing terms for this process
    template<typename DenseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ForcingFunction(
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const auto& state_parameter_indices,  // acts like std::unordered_map<std::string, std::size_t>
        const auto& state_variable_indices,   // acts like std::unordered_map<std::string, std::size_t>
        std::map<std::string, std::map<AerosolProperty, AerosolPropertyProvider<DenseMatrixPolicy>>> /* providers */
    ) const
    {
      return ForcingFunction<DenseMatrixPolicy>(phase_prefixes, state_parameter_indices, state_variable_indices);
    }

    /// @brief Returns a function that calculates the forcing terms for this process (original)
    template<typename DenseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ForcingFunction(
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const auto& state_parameter_indices,  // acts like std::unordered_map<std::string, std::size_t>
        const auto& state_variable_indices    // acts like std::unordered_map<std::string, std::size_t>
    ) const
    {
      StateVariableIndices variable_indices = GetStateVariableIndices(phase_prefixes, state_variable_indices);
      std::size_t k_index = GetParameterIndex(state_parameter_indices);

      DenseMatrixPolicy dummy_state_parameters{ 1, state_parameter_indices.size(), 0.0 };
      DenseMatrixPolicy dummy_state_variables{ 1, state_variable_indices.size(), 0.0 };

      if (min_halflife_ > 0.0)
      {
        return ForcingFunctionCapped<DenseMatrixPolicy>(
            variable_indices, k_index, dummy_state_parameters, dummy_state_variables);
      }

      auto storage = std::make_shared<IndexStorage<DenseMatrixPolicy>>(variable_indices);
      auto reactant_view = storage->reactant_indices_.GetView();
      auto product_view = storage->product_indices_.GetView();
      auto solvent_view = storage->solvent_indices_.GetView();
      const std::size_t num_phases = variable_indices.number_of_phase_instances_;
      const std::size_t num_reactants = reactants_.size();
      const std::size_t num_products = products_.size();
      const double eps = solvent_floor_;

      auto function = DenseMatrixPolicy::Function(
            MICM_LAMBDA(
                const typename DenseMatrixPolicy::ConstViewType& params_view,
                const typename DenseMatrixPolicy::ConstViewType& state_view,
                const typename DenseMatrixPolicy::ViewType& forcing_view)
          {
            // For each phase instance, calculate the reaction rate and update the forcing terms
              for (std::size_t phase = 0; phase < num_phases; ++phase)
            {
                const std::size_t solvent_idx = solvent_view[phase];
                auto rate = forcing_view.GetRowVariable();
              // Calculate the damped rate: k * [S] / ([S] + eps)^n_r * prod([reactants])
                forcing_view.ForEachRowStrict(
                    [num_reactants, eps](const micm::Real& k, const micm::Real& solvent, micm::Real& out)
                    { out = k * solvent / std::pow(solvent + eps, num_reactants); },
                    params_view.GetConstColumnView(k_index),
                    state_view.GetConstColumnView(solvent_idx),
                  rate);
                for (std::size_t r = 0; r < num_reactants; ++r)
              {
                  forcing_view.ForEachRowStrict(
                      [](const micm::Real& reactant, micm::Real& out) { out *= reactant; },
                      state_view.GetConstColumnView(reactant_view[phase * num_reactants + r]),
                    rate);
              }

              // Apply the reaction rate to the forcing terms for reactants and products
                for (std::size_t r = 0; r < num_reactants; ++r)
              {
                  forcing_view.ForEachRowStrict(
                      [](const micm::Real& rate, micm::Real& forcing) { forcing -= rate; },
                    rate,
                      forcing_view.GetColumnView(reactant_view[phase * num_reactants + r]));
              }
                for (std::size_t p = 0; p < num_products; ++p)
              {
                  forcing_view.ForEachRowStrict(
                      [](const micm::Real& rate, micm::Real& forcing) { forcing += rate; },
                    rate,
                      forcing_view.GetColumnView(product_view[phase * num_products + p]));
              }
            }
          },
          dummy_state_parameters,
          dummy_state_variables,
          dummy_state_variables);

      return [storage, function](
                 const DenseMatrixPolicy& state_parameters,
                 const DenseMatrixPolicy& state_variables,
                 DenseMatrixPolicy& forcing_terms) mutable { function(state_parameters, state_variables, forcing_terms); };
    }

    /// @brief Returns a function that calculates the Jacobian contributions for this process (common interface overload)
    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> JacobianFunction(
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const auto& state_parameter_indices,  // acts like std::unordered_map<std::string, std::size_t>
        const auto& state_variable_indices,   // acts like std::unordered_map<std::string, std::size_t>
        const SparseMatrixPolicy& jacobian,
        std::map<std::string, std::map<AerosolProperty, AerosolPropertyProvider<DenseMatrixPolicy>>> /* providers */) const
    {
      return JacobianFunction<DenseMatrixPolicy, SparseMatrixPolicy>(
          phase_prefixes, state_parameter_indices, state_variable_indices, jacobian);
    }

    /// @brief Returns a function that calculates the Jacobian contributions for this process (original)
    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> JacobianFunction(
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const auto& state_parameter_indices,  // acts like std::unordered_map<std::string, std::size_t>
        const auto& state_variable_indices,   // acts like std::unordered_map<std::string, std::size_t>
        const SparseMatrixPolicy& jacobian) const
    {
      StateVariableIndices variable_indices = GetStateVariableIndices(phase_prefixes, state_variable_indices);
      JacobianIndices jacobian_indices = GetJacobianIndices(variable_indices, jacobian);
      std::size_t k_index = GetParameterIndex(state_parameter_indices);

      DenseMatrixPolicy dummy_state_parameters{ 1, state_parameter_indices.size(), 0.0 };
      DenseMatrixPolicy dummy_state_variables{ 1, state_variable_indices.size(), 0.0 };

      if (min_halflife_ > 0.0)
      {
        return JacobianFunctionCapped<DenseMatrixPolicy, SparseMatrixPolicy>(
            variable_indices, jacobian_indices, k_index, dummy_state_parameters, dummy_state_variables, jacobian);
      }

      auto storage = std::make_shared<IndexStorage<SparseMatrixPolicy>>(variable_indices, jacobian_indices);
      auto reactant_view = storage->reactant_indices_.GetView();
      auto solvent_view = storage->solvent_indices_.GetView();
      auto jac_id_view = storage->jacobian_flat_ids_.GetView();
      const std::size_t num_phases = variable_indices.number_of_phase_instances_;
      const std::size_t num_reactants = reactants_.size();
      const std::size_t num_products = products_.size();
      const double eps = solvent_floor_;

      auto function = SparseMatrixPolicy::Function(
            MICM_LAMBDA(
                const typename DenseMatrixPolicy::ConstViewType& params_view,
                const typename DenseMatrixPolicy::ConstViewType& state_view,
                const typename SparseMatrixPolicy::ViewType& jac_view)
          {
              const std::size_t pairs_per_phase = (num_reactants + 1) * (num_reactants + num_products);

            // For each phase instance, calculate the partial derivatives for the Jacobian entries
              for (std::size_t phase = 0; phase < num_phases; ++phase)
            {
                const std::size_t solvent_idx = solvent_view[phase];
                auto d_rate_d_ind = jac_view.GetBlockVariable();
                std::size_t pair = phase * pairs_per_phase;

              // Calculate partials for independent reactants
                for (std::size_t i_ind = 0; i_ind < num_reactants; ++i_ind)
              {
                // dr/d[R_i] = k * [S] / ([S]+eps)^n_r * prod(R_j, j!=i)
                  jac_view.ForEachBlockStrict(
                      [num_reactants, eps](const micm::Real& k, const micm::Real& solvent, micm::Real& partial)
                      { partial = k * solvent / std::pow(solvent + eps, num_reactants); },
                      params_view.GetConstColumnView(k_index),
                      state_view.GetConstColumnView(solvent_idx),
                    d_rate_d_ind);
                // add contributions to the partial from the other reactants
                  for (std::size_t r = 0; r < num_reactants; ++r)
                {
                  if (r == i_ind)
                    continue;  // Skip the variable we're taking the derivative with respect to
                    jac_view.ForEachBlockStrict(
                        [](const micm::Real& reactant, micm::Real& partial) { partial *= reactant; },
                        state_view.GetConstColumnView(reactant_view[phase * num_reactants + r]),
                      d_rate_d_ind);
                }
                // apply partial to dependent reactants (subtract: -J convention)
                  for (std::size_t i_dep = 0; i_dep < num_reactants; ++i_dep)
                {
                    jac_view.ForEachBlockStrict(
                        [](const micm::Real& partial, micm::Real& jac) { jac += partial; },
                      d_rate_d_ind,
                        jac_view.GetBlockView(jac_id_view[pair++]));
                }
                // apply partial to dependent products (subtract: -J convention)
                  for (std::size_t i_dep = 0; i_dep < num_products; ++i_dep)
                {
                    jac_view.ForEachBlockStrict(
                        [](const micm::Real& partial, micm::Real& jac) { jac -= partial; },
                      d_rate_d_ind,
                        jac_view.GetBlockView(jac_id_view[pair++]));
                }
              }
              // Calculate partials for independent solvent
              // dr/d[S] = k * (eps + (1-n_r)*[S]) / ([S]+eps)^(n_r+1) * prod([R_i])
                jac_view.ForEachBlockStrict(
                    [num_reactants, eps](const micm::Real& k, const micm::Real& solvent, micm::Real& partial)
                    {
                      partial = k * (eps + (1.0 - static_cast<micm::Real>(num_reactants)) * solvent) /
                                std::pow(solvent + eps, num_reactants + 1);
                  },
                    params_view.GetConstColumnView(k_index),
                    state_view.GetConstColumnView(solvent_idx),
                  d_rate_d_ind);
              // add contributions to the partial from the reactants
                for (std::size_t r = 0; r < num_reactants; ++r)
              {
                  jac_view.ForEachBlockStrict(
                      [](const micm::Real& reactant, micm::Real& partial) { partial *= reactant; },
                      state_view.GetConstColumnView(reactant_view[phase * num_reactants + r]),
                    d_rate_d_ind);
              }
              // apply partials to dependent reactants (subtract: -J convention)
                for (std::size_t i_dep = 0; i_dep < num_reactants; ++i_dep)
              {
                  jac_view.ForEachBlockStrict(
                      [](const micm::Real& partial, micm::Real& jac) { jac += partial; },
                    d_rate_d_ind,
                      jac_view.GetBlockView(jac_id_view[pair++]));
              }
              // apply partials to dependent products (subtract: -J convention)
                for (std::size_t i_dep = 0; i_dep < num_products; ++i_dep)
              {
                  jac_view.ForEachBlockStrict(
                      [](const micm::Real& partial, micm::Real& jac) { jac -= partial; },
                    d_rate_d_ind,
                      jac_view.GetBlockView(jac_id_view[pair++]));
              }
            }
          },
          dummy_state_parameters,
          dummy_state_variables,
          jacobian);

      return [storage, function](
                 const DenseMatrixPolicy& state_parameters,
                 const DenseMatrixPolicy& state_variables,
                 SparseMatrixPolicy& jacobian_values) mutable { function(state_parameters, state_variables, jacobian_values); };
    }

   private:
    /// @brief Helper struct for keeping track of state variable indices for reactants, products, and solvent across
    /// multiple phase instances (e.g. grid cells)
    struct StateVariableIndices
    {
      std::size_t number_of_phase_instances_;  ///< Number of instances of the phase in the system
      micm::Matrix<std::size_t>
          reactant_indices_;  ///< Matrix of state variable indices for reactants (num_prefixes x num_reactants)
      micm::Matrix<std::size_t>
          product_indices_;  ///< Matrix of state variable indices for products (num_prefixes x num_products)
      std::vector<std::size_t> solvent_indices_;  ///< Vector of state variable indices for solvent (num_prefixes)
    };

    /// @brief Helper struct for keeping track of Jacobian sparse matrix elements
    struct JacobianIndices
    {
      micm::Matrix<std::size_t>
          indices_;  // Index in sparse matrix for each dependent/independent pair (num_pairs x num_prefixes)
    };

    /// @brief Device-ready copies of the indices that the solve-time functions read
    /// @details The solve-time functions hold a shared pointer to this storage. Their kernels capture views into it.
    template<typename MatrixPolicy>
    struct IndexStorage
    {
      using Vector = typename MatrixPolicy::template VectorType<std::size_t>;
      Vector reactant_indices_;
      Vector product_indices_;
      Vector solvent_indices_;
      Vector jacobian_flat_ids_;

      IndexStorage(const StateVariableIndices& variable_indices, const JacobianIndices& jacobian_indices = {})
          : reactant_indices_(variable_indices.reactant_indices_.AsVector()),
            product_indices_(variable_indices.product_indices_.AsVector()),
            solvent_indices_(variable_indices.solvent_indices_),
            jacobian_flat_ids_(jacobian_indices.indices_.AsVector())
      {
        reactant_indices_.CopyToDevice();
        product_indices_.CopyToDevice();
        solvent_indices_.CopyToDevice();
        jacobian_flat_ids_.CopyToDevice();
      }
    };

   public:
    /// @brief Returns the capped forcing function (called only when min_halflife_ > 0)
    /// @details The soft-min uses the exponent p = 10 and the floor 1e-300. NVCC requires this function to be public.
    template<typename DenseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)> ForcingFunctionCapped(
        const StateVariableIndices& variable_indices,
        const std::size_t k_index,
        DenseMatrixPolicy& dummy_state_parameters,
        DenseMatrixPolicy& dummy_state_variables) const
    {
      auto storage = std::make_shared<IndexStorage<DenseMatrixPolicy>>(variable_indices);
      auto reactant_view = storage->reactant_indices_.GetView();
      auto product_view = storage->product_indices_.GetView();
      auto solvent_view = storage->solvent_indices_.GetView();
      const std::size_t num_phases = variable_indices.number_of_phase_instances_;
      const std::size_t num_reactants = reactants_.size();
      const std::size_t num_products = products_.size();
            const double eps = solvent_floor_;
            const double t_half = min_halflife_;

      auto function = DenseMatrixPolicy::Function(
            MICM_LAMBDA(
                const typename DenseMatrixPolicy::ConstViewType& params_view,
                const typename DenseMatrixPolicy::ConstViewType& state_view,
                const typename DenseMatrixPolicy::ViewType& forcing_view)
            {
              for (std::size_t phase = 0; phase < num_phases; ++phase)
              {
                const std::size_t solvent_idx = solvent_view[phase];
                auto rate = forcing_view.GetRowVariable();
                auto accum = forcing_view.GetRowVariable();

              // 1. Compute raw rate: k * [S] / ([S] + eps)^n_r * prod([R_i])
                forcing_view.ForEachRowStrict(
                    [num_reactants, eps](const micm::Real& k, const micm::Real& solvent, micm::Real& out)
                    { out = k * solvent / std::pow(solvent + eps, num_reactants); },
                    params_view.GetConstColumnView(k_index),
                    state_view.GetConstColumnView(solvent_idx),
                  rate);
                for (std::size_t r = 0; r < num_reactants; ++r)
              {
                  forcing_view.ForEachRowStrict(
                      [](const micm::Real& reactant, micm::Real& out) { out *= reactant; },
                      state_view.GetConstColumnView(reactant_view[phase * num_reactants + r]),
                    rate);
              }

              // 2. Compute soft-min of reactant concentrations: C_min = (sum R_i^{-p})^{-1/p}
                forcing_view.ForEachRowStrict(
                    [](const micm::Real& R, micm::Real& acc)
                    { acc = std::pow(R > micm::Real(1.0e-300) ? R : micm::Real(1.0e-300), -micm::Real(10.0)); },
                    state_view.GetConstColumnView(reactant_view[phase * num_reactants + 0]),
                  accum);
                for (std::size_t r = 1; r < num_reactants; ++r)
              {
                  forcing_view.ForEachRowStrict(
                      [](const micm::Real& R, micm::Real& acc)
                      { acc += std::pow(R > micm::Real(1.0e-300) ? R : micm::Real(1.0e-300), -micm::Real(10.0)); },
                      state_view.GetConstColumnView(reactant_view[phase * num_reactants + r]),
                    accum);
              }

              // 3. Apply tanh cap: rate = r_max * tanh(rate / r_max), where r_max = C_min / t_half
                forcing_view.ForEachRowStrict(
                    [t_half](micm::Real& out, micm::Real& acc)
                  {
                      const micm::Real c_min = std::pow(acc, -1.0 / micm::Real(10.0));
                      const micm::Real r_max = c_min / t_half;
                      if (r_max > micm::Real(1.0e-300))
                        out = r_max * std::tanh(out / r_max);
                  },
                  rate,
                  accum);

              // 4. Apply capped rate to forcing
                for (std::size_t r = 0; r < num_reactants; ++r)
              {
                  forcing_view.ForEachRowStrict(
                      [](const micm::Real& rate, micm::Real& forcing) { forcing -= rate; },
                    rate,
                      forcing_view.GetColumnView(reactant_view[phase * num_reactants + r]));
              }
                for (std::size_t p = 0; p < num_products; ++p)
              {
                  forcing_view.ForEachRowStrict(
                      [](const micm::Real& rate, micm::Real& forcing) { forcing += rate; },
                    rate,
                      forcing_view.GetColumnView(product_view[phase * num_products + p]));
              }
            }
          },
          dummy_state_parameters,
          dummy_state_variables,
          dummy_state_variables);

      return [storage, function](
                 const DenseMatrixPolicy& state_parameters,
                 const DenseMatrixPolicy& state_variables,
                 DenseMatrixPolicy& forcing_terms) mutable { function(state_parameters, state_variables, forcing_terms); };
    }

    /// @brief Returns the capped Jacobian function (called only when min_halflife_ > 0)
    /// @details The capped rate is r_c = r_max * tanh(r / r_max), where r_max = C_min / t_half
    ///          and C_min = (sum R_i^{-p})^{-1/p} is a smooth approximation to min(R_i).
    ///
    ///          For independent reactant R_j:
    ///            dr_c/dR_j = sech^2(u) * dr/dR_j + [tanh(u) - u*sech^2(u)] * (1/t_half) * (C_min/R_j)^{p+1}
    ///
    ///          For independent solvent S (r_max doesn't depend on S):
    ///            dr_c/dS = sech^2(u) * dr/dS
    ///
    ///          NVCC requires this function to be public.
    template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
    std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)> JacobianFunctionCapped(
        const StateVariableIndices& variable_indices,
        const JacobianIndices& jacobian_indices,
        const std::size_t k_index,
        DenseMatrixPolicy& dummy_state_parameters,
        DenseMatrixPolicy& dummy_state_variables,
        const SparseMatrixPolicy& jacobian) const
    {
      auto storage = std::make_shared<IndexStorage<SparseMatrixPolicy>>(variable_indices, jacobian_indices);
      auto reactant_view = storage->reactant_indices_.GetView();
      auto solvent_view = storage->solvent_indices_.GetView();
      auto jac_id_view = storage->jacobian_flat_ids_.GetView();
      const std::size_t num_phases = variable_indices.number_of_phase_instances_;
      const std::size_t num_reactants = reactants_.size();
      const std::size_t num_products = products_.size();
            const double eps = solvent_floor_;
            const double t_half = min_halflife_;

      auto function = SparseMatrixPolicy::Function(
            MICM_LAMBDA(
                const typename DenseMatrixPolicy::ConstViewType& params_view,
                const typename DenseMatrixPolicy::ConstViewType& state_view,
                const typename SparseMatrixPolicy::ViewType& jac_view)
            {
              const std::size_t pairs_per_phase = (num_reactants + 1) * (num_reactants + num_products);

              for (std::size_t phase = 0; phase < num_phases; ++phase)
              {
                const std::size_t solvent_idx = solvent_view[phase];
                auto d_rate_d_ind = jac_view.GetBlockVariable();
                auto raw_rate = jac_view.GetBlockVariable();
                auto sech2_var = jac_view.GetBlockVariable();
                auto corr_var = jac_view.GetBlockVariable();
                auto c_min_var = jac_view.GetBlockVariable();
                std::size_t pair = phase * pairs_per_phase;

              // Step A: Compute raw rate into raw_rate
                jac_view.ForEachBlockStrict(
                    [num_reactants, eps](const micm::Real& k, const micm::Real& solvent, micm::Real& rr)
                    { rr = k * solvent / std::pow(solvent + eps, num_reactants); },
                    params_view.GetConstColumnView(k_index),
                    state_view.GetConstColumnView(solvent_idx),
                  raw_rate);
                for (std::size_t r = 0; r < num_reactants; ++r)
              {
                  jac_view.ForEachBlockStrict(
                      [](const micm::Real& reactant, micm::Real& rr) { rr *= reactant; },
                      state_view.GetConstColumnView(reactant_view[phase * num_reactants + r]),
                    raw_rate);
              }

              // Step B: Compute soft-min sum into c_min_var
                jac_view.ForEachBlockStrict(
                    [](const micm::Real& R, micm::Real& cm)
                    { cm = std::pow(R > micm::Real(1.0e-300) ? R : micm::Real(1.0e-300), -micm::Real(10.0)); },
                    state_view.GetConstColumnView(reactant_view[phase * num_reactants + 0]),
                  c_min_var);
                for (std::size_t r = 1; r < num_reactants; ++r)
              {
                  jac_view.ForEachBlockStrict(
                      [](const micm::Real& R, micm::Real& cm)
                      { cm += std::pow(R > micm::Real(1.0e-300) ? R : micm::Real(1.0e-300), -micm::Real(10.0)); },
                      state_view.GetConstColumnView(reactant_view[phase * num_reactants + r]),
                    c_min_var);
              }

              // Step C: Convert to C_min, compute sech^2(u) and correction factor
                jac_view.ForEachBlockStrict(
                    [t_half](micm::Real& rr, micm::Real& cm, micm::Real& s2, micm::Real& cr)
                  {
                      cm = std::pow(cm, -1.0 / micm::Real(10.0));
                      const micm::Real r_max = cm / t_half;
                      if (r_max > micm::Real(1.0e-300))
                    {
                        const micm::Real u = rr / r_max;
                        const micm::Real th = std::tanh(u);
                      s2 = 1.0 - th * th;  // sech^2(u)
                      cr = (th - u * s2) / t_half;
                    }
                    else
                    {
                      s2 = 1.0;  // sech^2(0) = 1: uncapped Jacobian
                      cr = 0.0;  // no cap correction when r_max ≈ 0
                    }
                  },
                  raw_rate,
                  c_min_var,
                  sech2_var,
                  corr_var);

              // Step D: Partials for independent reactants
                for (std::size_t i_ind = 0; i_ind < num_reactants; ++i_ind)
              {
                // D1: Compute raw partial dr/dR_{i_ind}
                  jac_view.ForEachBlockStrict(
                      [num_reactants, eps](const micm::Real& k, const micm::Real& solvent, micm::Real& partial)
                      { partial = k * solvent / std::pow(solvent + eps, num_reactants); },
                      params_view.GetConstColumnView(k_index),
                      state_view.GetConstColumnView(solvent_idx),
                    d_rate_d_ind);
                  for (std::size_t r = 0; r < num_reactants; ++r)
                {
                  if (r == i_ind)
                    continue;
                    jac_view.ForEachBlockStrict(
                        [](const micm::Real& reactant, micm::Real& partial) { partial *= reactant; },
                        state_view.GetConstColumnView(reactant_view[phase * num_reactants + r]),
                      d_rate_d_ind);
                }

                // D2: Apply capping: partial = sech2 * partial + corr * (c_min/R_j)^{p+1}
                  jac_view.ForEachBlockStrict(
                      [](const micm::Real& s2,
                         const micm::Real& cr,
                         const micm::Real& cm,
                         const micm::Real& R,
                         micm::Real& partial)
                    {
                        const micm::Real ratio = cm / (R > micm::Real(1.0e-300) ? R : micm::Real(1.0e-300));
                        partial = s2 * partial + cr * std::pow(ratio, micm::Real(10.0) + 1.0);
                    },
                    sech2_var,
                    corr_var,
                    c_min_var,
                      state_view.GetConstColumnView(reactant_view[phase * num_reactants + i_ind]),
                    d_rate_d_ind);

                // D3: Apply to dependent reactants (-J convention)
                  for (std::size_t i_dep = 0; i_dep < num_reactants; ++i_dep)
                {
                    jac_view.ForEachBlockStrict(
                        [](const micm::Real& partial, micm::Real& jac) { jac += partial; },
                      d_rate_d_ind,
                        jac_view.GetBlockView(jac_id_view[pair++]));
                }
                // D4: Apply to dependent products (-J convention)
                  for (std::size_t i_dep = 0; i_dep < num_products; ++i_dep)
                {
                    jac_view.ForEachBlockStrict(
                        [](const micm::Real& partial, micm::Real& jac) { jac -= partial; },
                      d_rate_d_ind,
                        jac_view.GetBlockView(jac_id_view[pair++]));
                }
              }

              // Step E: Partial for independent solvent
              // dr/dS = k * (eps + (1-n_r)*S) / (S+eps)^{n_r+1} * prod(R_i)
                jac_view.ForEachBlockStrict(
                    [num_reactants, eps](const micm::Real& k, const micm::Real& solvent, micm::Real& partial)
                    {
                      partial = k * (eps + (1.0 - static_cast<micm::Real>(num_reactants)) * solvent) /
                                std::pow(solvent + eps, num_reactants + 1);
                  },
                    params_view.GetConstColumnView(k_index),
                    state_view.GetConstColumnView(solvent_idx),
                  d_rate_d_ind);
                for (std::size_t r = 0; r < num_reactants; ++r)
              {
                  jac_view.ForEachBlockStrict(
                      [](const micm::Real& reactant, micm::Real& partial) { partial *= reactant; },
                      state_view.GetConstColumnView(reactant_view[phase * num_reactants + r]),
                    d_rate_d_ind);
              }

              // Solvent capping: partial *= sech^2(u) (r_max doesn't depend on S)
                jac_view.ForEachBlockStrict(
                    [](const micm::Real& s2, micm::Real& partial) { partial *= s2; }, sech2_var, d_rate_d_ind);

              // Apply to dependent reactants
                for (std::size_t i_dep = 0; i_dep < num_reactants; ++i_dep)
              {
                  jac_view.ForEachBlockStrict(
                      [](const micm::Real& partial, micm::Real& jac) { jac += partial; },
                    d_rate_d_ind,
                      jac_view.GetBlockView(jac_id_view[pair++]));
              }
              // Apply to dependent products
                for (std::size_t i_dep = 0; i_dep < num_products; ++i_dep)
              {
                  jac_view.ForEachBlockStrict(
                      [](const micm::Real& partial, micm::Real& jac) { jac -= partial; },
                    d_rate_d_ind,
                      jac_view.GetBlockView(jac_id_view[pair++]));
              }
            }
          },
          dummy_state_parameters,
          dummy_state_variables,
          jacobian);

      return [storage, function](
                 const DenseMatrixPolicy& state_parameters,
                 const DenseMatrixPolicy& state_variables,
                 SparseMatrixPolicy& jacobian_values) mutable { function(state_parameters, state_variables, jacobian_values); };
    }

   private:
    /// @brief Helper function to return parameter index for the rate constant
    std::size_t GetParameterIndex(
        const auto& state_parameter_indices  // acts like std::unordered_map<std::string, std::size_t>
    ) const
    {
      std::string k_param = phase_.name_ + "." + uuid_ + ".k";
      if (state_parameter_indices.find(k_param) == state_parameter_indices.end())
      {
        throw MiamException(
            MIAM_ERROR_CATEGORY_INTERNAL,
            MIAM_INTERNAL_MISSING_STATE_PARAMETER,
            "Internal Error: GetParameterIndex: Rate constant parameter " + k_param +
                " not found in state_parameter_indices");
      }
      return state_parameter_indices.at(k_param);
    }

    /// @brief Helper function to return variable indices for all species involved in the reaction
    StateVariableIndices GetStateVariableIndices(
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const auto& state_variable_indices  // acts like std::unordered_map<std::string, std::size_t>
    ) const
    {
      StateVariableIndices indices;
      auto phase_it = phase_prefixes.find(phase_.name_);
      const auto& prefixes = phase_it->second;
      indices.number_of_phase_instances_ = prefixes.size();
      indices.reactant_indices_ = micm::Matrix<std::size_t>(prefixes.size(), reactants_.size());
      indices.product_indices_ = micm::Matrix<std::size_t>(prefixes.size(), products_.size());
      indices.solvent_indices_ = std::vector<std::size_t>(prefixes.size());
      std::size_t i_phase = 0;
      for (const auto& prefix : prefixes)
      {
        for (std::size_t i_reactant = 0; i_reactant < reactants_.size(); ++i_reactant)
        {
          std::string reactant_var = prefix + "." + phase_.name_ + "." + reactants_[i_reactant].name_;
          if (state_variable_indices.find(reactant_var) == state_variable_indices.end())
          {
            throw MiamException(
                MIAM_ERROR_CATEGORY_INTERNAL,
                MIAM_INTERNAL_MISSING_STATE_VARIABLE,
                "Internal Error: GetStateVariableIndices: Reactant variable " + reactant_var +
                    " not found in state_variable_indices");
          }
          indices.reactant_indices_[i_phase][i_reactant] = state_variable_indices.at(reactant_var);
        }
        for (std::size_t i_product = 0; i_product < products_.size(); ++i_product)
        {
          std::string product_var = prefix + "." + phase_.name_ + "." + products_[i_product].name_;
          if (state_variable_indices.find(product_var) == state_variable_indices.end())
          {
            throw MiamException(
                MIAM_ERROR_CATEGORY_INTERNAL,
                MIAM_INTERNAL_MISSING_STATE_VARIABLE,
                "Internal Error: GetStateVariableIndices: Product variable " + product_var +
                    " not found in state_variable_indices");
          }
          indices.product_indices_[i_phase][i_product] = state_variable_indices.at(product_var);
        }
        std::string solvent_var = prefix + "." + phase_.name_ + "." + solvent_.name_;
        if (state_variable_indices.find(solvent_var) == state_variable_indices.end())
        {
          throw MiamException(
              MIAM_ERROR_CATEGORY_INTERNAL,
              MIAM_INTERNAL_MISSING_STATE_VARIABLE,
              "Internal Error: GetStateVariableIndices: Solvent variable " + solvent_var +
                  " not found in state_variable_indices");
        }
        indices.solvent_indices_[i_phase] = state_variable_indices.at(solvent_var);
        ++i_phase;
      }
      return indices;
    }

    /// @brief Helper function to return Jacobian sparse matrix indices for all pairs of species involved in the reaction
    JacobianIndices GetJacobianIndices(
        const StateVariableIndices& variable_indices,
        const auto& jacobian  // sparse matrix policy object for the Jacobian structure
    ) const
    {
      // For a forward-only reaction, independent variables are only reactants and solvent (not products).
      // Each reactant and each product depends on all reactants and the solvent.
      std::size_t num_pairs = (reactants_.size() + products_.size()) * (reactants_.size() + 1);  // +1 for solvent
      JacobianIndices jacobian_indices;
      jacobian_indices.indices_ = micm::Matrix<std::size_t>(variable_indices.number_of_phase_instances_, num_pairs);
      for (std::size_t i_phase = 0; i_phase < variable_indices.number_of_phase_instances_; ++i_phase)
      {
        std::size_t pair_index = 0;
        // add terms for independent reactants
        for (std::size_t i_ind = 0; i_ind < reactants_.size(); ++i_ind)
        {
          // ... and dependent reactants
          for (std::size_t i_dep = 0; i_dep < reactants_.size(); ++i_dep)
          {
            jacobian_indices.indices_[i_phase][pair_index++] = jacobian.VectorIndex(
                0, variable_indices.reactant_indices_[i_phase][i_dep], variable_indices.reactant_indices_[i_phase][i_ind]);
          }
          // ... and dependent products
          for (std::size_t i_dep = 0; i_dep < products_.size(); ++i_dep)
          {
            jacobian_indices.indices_[i_phase][pair_index++] = jacobian.VectorIndex(
                0, variable_indices.product_indices_[i_phase][i_dep], variable_indices.reactant_indices_[i_phase][i_ind]);
          }
        }
        // add terms for independent solvent
        // ... and dependent reactants
        for (std::size_t i_dep = 0; i_dep < reactants_.size(); ++i_dep)
        {
          jacobian_indices.indices_[i_phase][pair_index++] = jacobian.VectorIndex(
              0, variable_indices.reactant_indices_[i_phase][i_dep], variable_indices.solvent_indices_[i_phase]);
        }
        // ... and dependent products
        for (std::size_t i_dep = 0; i_dep < products_.size(); ++i_dep)
        {
          jacobian_indices.indices_[i_phase][pair_index++] = jacobian.VectorIndex(
              0, variable_indices.product_indices_[i_phase][i_dep], variable_indices.solvent_indices_[i_phase]);
        }
      }
      return jacobian_indices;
    }
  };
}  // namespace miam
