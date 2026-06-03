// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <miam/constraints/henry_law_equilibrium_constraint.hpp>
#include <miam/constraints/henry_law_equilibrium_constraint_builder.hpp>
#include <miam/math/condensation_rate.hpp>

#include <micm/system/conditions.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
#include <micm/util/jacobian_verification.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_standard_ordering.hpp>
#include <micm/util/sparse_matrix_vector_ordering.hpp>
#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

using namespace miam;

namespace
{
  using DMP = micm::Matrix<double>;

  constexpr double water_molecular_weight = 0.018;
  constexpr double water_density = 1000.0;
  auto h2o =
      micm::Species{ "H2O",
                     { { "molecular weight [kg mol-1]", water_molecular_weight }, { "density [kg m-3]", water_density } } };
  auto A_g = micm::Species{ "A_g" };
  auto A_aq = micm::Species{ "A_aq" };
  auto aqueous_phase = micm::Phase{ "AQUEOUS", { { h2o }, { A_aq } } };

  /// @brief Build state parameter index map from a constraint's ConstraintStateParameterNames
  template<typename ConstraintT>
  std::unordered_map<std::string, std::size_t> BuildParamIndices(
      const ConstraintT& constraint,
      const std::map<std::string, std::set<std::string>>& phase_prefixes)
  {
    std::unordered_map<std::string, std::size_t> indices;
    std::size_t idx = 0;
    for (const auto& name : constraint.ConstraintStateParameterNames(phase_prefixes))
      indices[name] = idx++;
    return indices;
  }
}  // namespace

// ── ConstraintAlgebraicVariableNames ──

TEST(HenryLawEquilibriumConstraint, AlgebraicVariableNamesSinglePrefix)
{
  auto hlc = [](const micm::Conditions&) { return 5.0e3; };
  auto constraint = HenryLawEquilibriumConstraintBuilder()
                        .SetGasSpecies(A_g)
                        .SetCondensedSpecies(A_aq)
                        .SetSolvent(h2o)
                        .SetCondensedPhase(aqueous_phase)
                        .SetHenryLawConstant(hlc)
                        .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  auto names = constraint.ConstraintAlgebraicVariableNames(phase_prefixes);
  EXPECT_EQ(names.size(), 1);
  EXPECT_TRUE(names.count("DROP.AQUEOUS.A_aq"));
}

TEST(HenryLawEquilibriumConstraint, AlgebraicVariableNamesMultiplePrefixes)
{
  auto hlc = [](const micm::Conditions&) { return 5.0e3; };
  auto constraint = HenryLawEquilibriumConstraintBuilder()
                        .SetGasSpecies(A_g)
                        .SetCondensedSpecies(A_aq)
                        .SetSolvent(h2o)
                        .SetCondensedPhase(aqueous_phase)
                        .SetHenryLawConstant(hlc)
                        .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("LARGE");
  phase_prefixes["AQUEOUS"].insert("SMALL");

  auto names = constraint.ConstraintAlgebraicVariableNames(phase_prefixes);
  EXPECT_EQ(names.size(), 2);
  EXPECT_TRUE(names.count("SMALL.AQUEOUS.A_aq"));
  EXPECT_TRUE(names.count("LARGE.AQUEOUS.A_aq"));
}

// ── ConstraintSpeciesDependencies ──

TEST(HenryLawEquilibriumConstraint, SpeciesDependencies)
{
  auto hlc = [](const micm::Conditions&) { return 5.0e3; };
  auto constraint = HenryLawEquilibriumConstraintBuilder()
                        .SetGasSpecies(A_g)
                        .SetCondensedSpecies(A_aq)
                        .SetSolvent(h2o)
                        .SetCondensedPhase(aqueous_phase)
                        .SetHenryLawConstant(hlc)
                        .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  auto deps = constraint.ConstraintSpeciesDependencies(phase_prefixes);
  // Gas species (standalone) + condensed species + solvent
  EXPECT_EQ(deps.size(), 3);
  EXPECT_TRUE(deps.count("A_g"));                // gas, no prefix
  EXPECT_TRUE(deps.count("DROP.AQUEOUS.A_aq"));  // condensed
  EXPECT_TRUE(deps.count("DROP.AQUEOUS.H2O"));   // solvent
}

TEST(HenryLawEquilibriumConstraint, SpeciesDependenciesMultipleInstances)
{
  auto hlc = [](const micm::Conditions&) { return 5.0e3; };
  auto constraint = HenryLawEquilibriumConstraintBuilder()
                        .SetGasSpecies(A_g)
                        .SetCondensedSpecies(A_aq)
                        .SetSolvent(h2o)
                        .SetCondensedPhase(aqueous_phase)
                        .SetHenryLawConstant(hlc)
                        .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("SMALL");
  phase_prefixes["AQUEOUS"].insert("LARGE");

  auto deps = constraint.ConstraintSpeciesDependencies(phase_prefixes);
  // 1 gas + 2*2 condensed = 5
  EXPECT_EQ(deps.size(), 5);
  EXPECT_TRUE(deps.count("A_g"));
  EXPECT_TRUE(deps.count("SMALL.AQUEOUS.A_aq"));
  EXPECT_TRUE(deps.count("SMALL.AQUEOUS.H2O"));
  EXPECT_TRUE(deps.count("LARGE.AQUEOUS.A_aq"));
  EXPECT_TRUE(deps.count("LARGE.AQUEOUS.H2O"));
}

// ── NonZeroConstraintJacobianElements ──

TEST(HenryLawEquilibriumConstraint, NonZeroJacobianElements)
{
  auto hlc = [](const micm::Conditions&) { return 5.0e3; };
  auto constraint = HenryLawEquilibriumConstraintBuilder()
                        .SetGasSpecies(A_g)
                        .SetCondensedSpecies(A_aq)
                        .SetSolvent(h2o)
                        .SetCondensedPhase(aqueous_phase)
                        .SetHenryLawConstant(hlc)
                        .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  std::unordered_map<std::string, std::size_t> state_indices;
  state_indices["A_g"] = 0;
  state_indices["DROP.AQUEOUS.A_aq"] = 1;
  state_indices["DROP.AQUEOUS.H2O"] = 2;

  auto elements = constraint.NonZeroConstraintJacobianElements(phase_prefixes, state_indices);
  // Algebraic row = A_aq (1). Depends on: A_g (0), A_aq (1), H2O (2)
  EXPECT_EQ(elements.size(), 3);
  EXPECT_TRUE(elements.count({ 1, 0 }));  // dG/d[A_g]
  EXPECT_TRUE(elements.count({ 1, 1 }));  // dG/d[A_aq]
  EXPECT_TRUE(elements.count({ 1, 2 }));  // dG/d[H2O]
}

// ── ConstraintResidualFunction ──

TEST(HenryLawEquilibriumConstraint, ResidualSingleInstance)
{
  // G = HLC * R * T * f_v * [A_g] - [A_aq]
  // f_v = [H2O] * Mw / rho
  double HLC = 5.0e3;
  double T = 298.15;
  auto hlc = [HLC](const micm::Conditions&) { return HLC; };
  auto constraint = HenryLawEquilibriumConstraintBuilder()
                        .SetGasSpecies(A_g)
                        .SetCondensedSpecies(A_aq)
                        .SetSolvent(h2o)
                        .SetCondensedPhase(aqueous_phase)
                        .SetHenryLawConstant(hlc)
                        .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  std::unordered_map<std::string, std::size_t> state_indices;
  state_indices["A_g"] = 0;
  state_indices["DROP.AQUEOUS.A_aq"] = 1;
  state_indices["DROP.AQUEOUS.H2O"] = 2;

  // Initialize HLC*R*T into state parameters
  using DMP = micm::Matrix<double>;
  auto pi = BuildParamIndices(constraint, phase_prefixes);
  DMP state_params{ 1, std::max(pi.size(), std::size_t{ 1 }), 0.0 };
  std::vector<micm::Conditions> conditions(1);
  conditions[0].temperature_ = T;
  auto update_fn = constraint.UpdateConstraintParametersFunction<DMP>(phase_prefixes, pi);
  update_fn(conditions, state_params);

  auto residual_fn = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, pi, state_indices);

  double gas_conc = 1.0e-6;
  double aq_conc = 0.5;
  double solvent_conc = 300.0;

  DMP state_variables{ 1, 3, 0.0 };
  state_variables[0][0] = gas_conc;
  state_variables[0][1] = aq_conc;
  state_variables[0][2] = solvent_conc;

  DMP residual{ 1, 3, 0.0 };
  residual_fn(state_variables, state_params, residual);

  double f_v = solvent_conc * water_molecular_weight / water_density;
  double hlc_rt = HLC * micm::constants::GAS_CONSTANT * T;
  double expected = hlc_rt * f_v * gas_conc - aq_conc;
  EXPECT_NEAR(residual[0][1], expected, std::abs(expected) * 1.0e-12);
}

// ── ConstraintResidualFunction — multiple instances ──

TEST(HenryLawEquilibriumConstraint, ResidualMultipleInstances)
{
  double HLC = 3.0e3;
  double T = 300.0;
  auto hlc = [HLC](const micm::Conditions&) { return HLC; };
  auto constraint = HenryLawEquilibriumConstraintBuilder()
                        .SetGasSpecies(A_g)
                        .SetCondensedSpecies(A_aq)
                        .SetSolvent(h2o)
                        .SetCondensedPhase(aqueous_phase)
                        .SetHenryLawConstant(hlc)
                        .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("LARGE");
  phase_prefixes["AQUEOUS"].insert("SMALL");

  std::unordered_map<std::string, std::size_t> state_indices;
  state_indices["A_g"] = 0;
  state_indices["LARGE.AQUEOUS.A_aq"] = 1;
  state_indices["LARGE.AQUEOUS.H2O"] = 2;
  state_indices["SMALL.AQUEOUS.A_aq"] = 3;
  state_indices["SMALL.AQUEOUS.H2O"] = 4;

  using DMP = micm::Matrix<double>;
  auto pi = BuildParamIndices(constraint, phase_prefixes);
  DMP state_params{ 1, std::max(pi.size(), std::size_t{ 1 }), 0.0 };
  std::vector<micm::Conditions> conditions(1);
  conditions[0].temperature_ = T;
  auto update_fn = constraint.UpdateConstraintParametersFunction<DMP>(phase_prefixes, pi);
  update_fn(conditions, state_params);

  auto residual_fn = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, pi, state_indices);

  DMP state_variables{ 1, 5, 0.0 };
  state_variables[0][0] = 1.0e-5;  // [A_g]
  state_variables[0][1] = 0.3;     // LARGE [A_aq]
  state_variables[0][2] = 400.0;   // LARGE [H2O]
  state_variables[0][3] = 0.1;     // SMALL [A_aq]
  state_variables[0][4] = 200.0;   // SMALL [H2O]

  DMP residual{ 1, 5, 0.0 };
  residual_fn(state_variables, state_params, residual);

  double hlc_rt = HLC * micm::constants::GAS_CONSTANT * T;
  double gas_conc = 1.0e-5;

  // LARGE: f_v = 400.0 * 0.018 / 1000.0 = 0.0072
  double fv_large = 400.0 * water_molecular_weight / water_density;
  EXPECT_NEAR(residual[0][1], hlc_rt * fv_large * gas_conc - 0.3, std::abs(hlc_rt * fv_large * gas_conc) * 1.0e-10);

  // SMALL: f_v = 200.0 * 0.018 / 1000.0 = 0.0036
  double fv_small = 200.0 * water_molecular_weight / water_density;
  EXPECT_NEAR(residual[0][3], hlc_rt * fv_small * gas_conc - 0.1, std::abs(hlc_rt * fv_small * gas_conc) * 1.0e-10);
}

// ── ConstraintJacobianFunction ──

TEST(HenryLawEquilibriumConstraint, JacobianSingleInstance)
{
  double HLC = 5.0e3;
  double T = 298.15;
  auto hlc = [HLC](const micm::Conditions&) { return HLC; };
  auto constraint = HenryLawEquilibriumConstraintBuilder()
                        .SetGasSpecies(A_g)
                        .SetCondensedSpecies(A_aq)
                        .SetSolvent(h2o)
                        .SetCondensedPhase(aqueous_phase)
                        .SetHenryLawConstant(hlc)
                        .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  std::unordered_map<std::string, std::size_t> state_indices;
  state_indices["A_g"] = 0;
  state_indices["DROP.AQUEOUS.A_aq"] = 1;
  state_indices["DROP.AQUEOUS.H2O"] = 2;

  using DMP = micm::Matrix<double>;
  using SMP = micm::SparseMatrix<double, micm::SparseMatrixStandardOrderingCompressedSparseRow>;
  auto pi = BuildParamIndices(constraint, phase_prefixes);
  DMP state_params{ 1, std::max(pi.size(), std::size_t{ 1 }), 0.0 };
  std::vector<micm::Conditions> conditions(1);
  conditions[0].temperature_ = T;
  auto update_fn = constraint.UpdateConstraintParametersFunction<DMP>(phase_prefixes, pi);
  update_fn(conditions, state_params);

  auto nz_elements = constraint.NonZeroConstraintJacobianElements(phase_prefixes, state_indices);
  auto builder = SMP::Create(3).SetNumberOfBlocks(1);
  for (const auto& [row, col] : nz_elements)
    builder.WithElement(row, col);
  SMP jacobian(builder);

  auto jac_fn = constraint.ConstraintJacobianFunction<DMP, SMP>(phase_prefixes, pi, state_indices, jacobian);

  double gas_conc = 1.0e-6;
  double solvent_conc = 300.0;

  DMP state_variables{ 1, 3, 0.0 };
  state_variables[0][0] = gas_conc;
  state_variables[0][1] = 0.5;           // [A_aq]
  state_variables[0][2] = solvent_conc;  // [H2O]

  for (auto& v : jacobian.AsVector())
    v = 0.0;
  jac_fn(state_variables, state_params, jacobian);

  double hlc_rt = HLC * micm::constants::GAS_CONSTANT * T;
  double f_v = solvent_conc * water_molecular_weight / water_density;
  double Mw_rho = water_molecular_weight / water_density;

  // jac -= dG/d[A_g] = HLC*R*T*f_v → jac[1,0] = -HLC*R*T*f_v
  EXPECT_NEAR(jacobian.AsVector()[jacobian.VectorIndex(0, 1, 0)], -hlc_rt * f_v, std::abs(hlc_rt * f_v) * 1.0e-10);

  // jac -= dG/d[A_aq] = -1 → jac[1,1] = +1
  EXPECT_NEAR(jacobian.AsVector()[jacobian.VectorIndex(0, 1, 1)], 1.0, 1.0e-12);

  // jac -= dG/d[H2O] = HLC*R*T*Mw/rho*[A_g] → jac[1,2] = -HLC*R*T*Mw/rho*[A_g]
  EXPECT_NEAR(
      jacobian.AsVector()[jacobian.VectorIndex(0, 1, 2)],
      -hlc_rt * Mw_rho * gas_conc,
      std::abs(hlc_rt * Mw_rho * gas_conc) * 1.0e-10);
}

// ── UpdateConstraintParameters — temperature-dependent HLC ──

TEST(HenryLawEquilibriumConstraint, UpdateConstraintParametersTemperatureDep)
{
  // HLC(T) = 1000.0 / T  =>  HLC*R*T = 1000 * R  (temperature-independent)
  auto hlc = [](const micm::Conditions& c) { return 1000.0 / c.temperature_; };
  auto constraint = HenryLawEquilibriumConstraintBuilder()
                        .SetGasSpecies(A_g)
                        .SetCondensedSpecies(A_aq)
                        .SetSolvent(h2o)
                        .SetCondensedPhase(aqueous_phase)
                        .SetHenryLawConstant(hlc)
                        .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  using DMP = micm::Matrix<double>;
  auto pi = BuildParamIndices(constraint, phase_prefixes);
  ASSERT_EQ(pi.size(), 1u);

  std::vector<micm::Conditions> conditions(2);
  conditions[0].temperature_ = 250.0;
  conditions[1].temperature_ = 350.0;

  DMP state_params{ 2, 1, 0.0 };
  auto update_fn = constraint.UpdateConstraintParametersFunction<DMP>(phase_prefixes, pi);
  update_fn(conditions, state_params);

  // HLC*R*T = (1000/T) * R * T = 1000 * R (temperature-independent)
  double expected = 1000.0 * micm::constants::GAS_CONSTANT;
  std::size_t hlc_rt_col = pi.begin()->second;
  EXPECT_NEAR(state_params[0][hlc_rt_col], expected, expected * 1.0e-12);
  EXPECT_NEAR(state_params[1][hlc_rt_col], expected, expected * 1.0e-12);
}

// ── Builder ──

TEST(HenryLawEquilibriumConstraint, BuilderValidation)
{
  struct FakeHLC
  {
    double Calculate(const micm::Conditions&) const
    {
      return 5.0e3;
    }
  };

  // Missing gas species
  EXPECT_THROW(
      HenryLawEquilibriumConstraintBuilder()
          .SetCondensedSpecies(A_aq)
          .SetSolvent(h2o)
          .SetCondensedPhase(aqueous_phase)
          .SetHenryLawConstant(FakeHLC{})
          .Build(),
      std::runtime_error);

  // Valid build
  auto constraint = HenryLawEquilibriumConstraintBuilder()
                        .SetGasSpecies(A_g)
                        .SetCondensedSpecies(A_aq)
                        .SetSolvent(h2o)
                        .SetCondensedPhase(aqueous_phase)
                        .SetHenryLawConstant(FakeHLC{})
                        .Build();
  EXPECT_EQ(constraint.gas_species_.name_, "A_g");
  EXPECT_EQ(constraint.condensed_species_.name_, "A_aq");
}

// ═══════════════════════════════════════════════════════════════════
// Expanded Scenarios
// ═══════════════════════════════════════════════════════════════════

namespace
{
  using DMP = micm::Matrix<double>;
  using SMP = micm::SparseMatrix<double, micm::SparseMatrixStandardOrderingCompressedSparseRow>;

  SMP BuildConstraintJacobian(
      const HenryLawEquilibriumConstraint& constraint,
      const std::map<std::string, std::set<std::string>>& phase_prefixes,
      const std::unordered_map<std::string, std::size_t>& state_indices,
      std::size_t num_blocks)
  {
    auto elements = constraint.NonZeroConstraintJacobianElements(phase_prefixes, state_indices);
    auto builder = SMP::Create(state_indices.size()).SetNumberOfBlocks(num_blocks).InitialValue(0.0);
    for (const auto& [row, col] : elements)
      builder = builder.WithElement(row, col);
    return SMP(builder);
  }

  /// @brief Initialize HLC*R*T values into a state parameter matrix; returns the matrix
  DMP InitHlcRt(
      const HenryLawEquilibriumConstraint& constraint,
      const std::map<std::string, std::set<std::string>>& phase_prefixes,
      const std::unordered_map<std::string, std::size_t>& param_idx,
      std::size_t num_cells,
      double T = 298.15)
  {
    std::vector<micm::Conditions> conditions(num_cells);
    for (auto& c : conditions)
      c.temperature_ = T;
    DMP state_params{ num_cells, std::max(param_idx.size(), std::size_t{ 1 }), 0.0 };
    auto fn = constraint.UpdateConstraintParametersFunction<DMP>(phase_prefixes, param_idx);
    fn(conditions, state_params);
    return state_params;
  }

  DMP InitHlcRt(
      const HenryLawEquilibriumConstraint& constraint,
      const std::map<std::string, std::set<std::string>>& phase_prefixes,
      const std::unordered_map<std::string, std::size_t>& param_idx,
      const std::vector<micm::Conditions>& conditions)
  {
    DMP state_params{ conditions.size(), std::max(param_idx.size(), std::size_t{ 1 }), 0.0 };
    auto fn = constraint.UpdateConstraintParametersFunction<DMP>(phase_prefixes, param_idx);
    fn(conditions, state_params);
    return state_params;
  }

  /// @brief Finite-difference check for the constraint Jacobian
  /// @brief Compare analytical constraint Jacobian against central finite-difference approximation
  ///        using MICM's FiniteDifferenceJacobian / CompareJacobianToFiniteDifference utilities.
  void CheckConstraintFDJacobian(
      const HenryLawEquilibriumConstraint& constraint,
      const std::map<std::string, std::set<std::string>>& phase_prefixes,
      const std::unordered_map<std::string, std::size_t>& param_idx,
      const std::unordered_map<std::string, std::size_t>& state_indices,
      const DMP& state_variables,
      const DMP& state_params)
  {
    const std::size_t num_blocks = state_variables.NumRows();
    const std::size_t num_vars = state_indices.size();

    // Build sparse Jacobian structure and compute analytical Jacobian
    auto jacobian = BuildConstraintJacobian(constraint, phase_prefixes, state_indices, num_blocks);
    auto jac_fn = constraint.ConstraintJacobianFunction<DMP, SMP>(phase_prefixes, param_idx, state_indices, jacobian);
    jac_fn(state_variables, state_params, jacobian);

    // Build FD Jacobian — bind state_params into the residual callable
    // Use perturbation=1e-7 to match old central-difference scheme: h = max(|x|, 1) * 1e-7.
    // Use atol=rtol=1e-5 to match old rel_tol tolerance.
    auto rf = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, param_idx, state_indices);
    auto fd_jac = micm::FiniteDifferenceJacobian<DMP>(
        [&](const DMP& vars, DMP& out)
        {
          out.Fill(0.0);
          rf(vars, state_params, out);
        },
        state_variables,
        num_vars,
        1e-7);

    // Compare analytical vs FD
    auto cmp = micm::CompareJacobianToFiniteDifference(jacobian, fd_jac, num_vars, 1e-5, 1e-5);
    EXPECT_TRUE(cmp.passed_) << "FD mismatch: block=" << cmp.worst_block_ << " row=" << cmp.worst_row_
                             << " col=" << cmp.worst_col_ << " +J(analytical)=" << cmp.worst_analytical_
                             << " +J(fd)=" << cmp.worst_fd_;

    // Verify no significant FD signal outside the declared sparsity pattern
    auto spc = micm::CheckJacobianSparsityCompleteness(jacobian, fd_jac, num_vars);
    EXPECT_TRUE(spc.passed_) << "Missing sparsity entry: block=" << spc.worst_block_ << " row=" << spc.worst_row_
                             << " col=" << spc.worst_col_ << " fd=" << spc.worst_fd_;
  }
}  // namespace

// ── FD Jacobian: single instance ──

TEST(HenryLawEquilibriumConstraint, JacobianFDSingleInstance)
{
  double HLC = 5.0e3;
  double T = 298.15;
  auto hlc = [HLC](const micm::Conditions&) { return HLC; };
  auto constraint = HenryLawEquilibriumConstraintBuilder()
                        .SetGasSpecies(A_g)
                        .SetCondensedSpecies(A_aq)
                        .SetSolvent(h2o)
                        .SetCondensedPhase(aqueous_phase)
                        .SetHenryLawConstant(hlc)
                        .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  std::unordered_map<std::string, std::size_t> si;
  si["A_g"] = 0;
  si["DROP.AQUEOUS.A_aq"] = 1;
  si["DROP.AQUEOUS.H2O"] = 2;

  auto pi = BuildParamIndices(constraint, phase_prefixes);
  auto sp = InitHlcRt(constraint, phase_prefixes, pi, 1, T);

  DMP sv{ 1, 3, 0.0 };
  sv[0][0] = 1.0e-6;
  sv[0][1] = 0.5;
  sv[0][2] = 300.0;

  CheckConstraintFDJacobian(constraint, phase_prefixes, pi, si, sv, sp);
}

// ── FD Jacobian: multiple instances ──

TEST(HenryLawEquilibriumConstraint, JacobianFDMultipleInstances)
{
  double HLC = 3.0e3;
  double T = 300.0;
  auto hlc = [HLC](const micm::Conditions&) { return HLC; };
  auto constraint = HenryLawEquilibriumConstraintBuilder()
                        .SetGasSpecies(A_g)
                        .SetCondensedSpecies(A_aq)
                        .SetSolvent(h2o)
                        .SetCondensedPhase(aqueous_phase)
                        .SetHenryLawConstant(hlc)
                        .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("LARGE");
  phase_prefixes["AQUEOUS"].insert("SMALL");

  std::unordered_map<std::string, std::size_t> si;
  si["A_g"] = 0;
  si["LARGE.AQUEOUS.A_aq"] = 1;
  si["LARGE.AQUEOUS.H2O"] = 2;
  si["SMALL.AQUEOUS.A_aq"] = 3;
  si["SMALL.AQUEOUS.H2O"] = 4;

  auto pi = BuildParamIndices(constraint, phase_prefixes);
  auto sp = InitHlcRt(constraint, phase_prefixes, pi, 1, T);

  DMP sv{ 1, 5, 0.0 };
  sv[0][0] = 1.0e-5;
  sv[0][1] = 0.3;
  sv[0][2] = 400.0;
  sv[0][3] = 0.1;
  sv[0][4] = 200.0;

  CheckConstraintFDJacobian(constraint, phase_prefixes, pi, si, sv, sp);
}

// ── Multiple grid cells: residual ──

TEST(HenryLawEquilibriumConstraint, ResidualMultipleCells)
{
  double HLC = 5.0e3;
  double T = 298.15;
  auto hlc = [HLC](const micm::Conditions&) { return HLC; };
  auto constraint = HenryLawEquilibriumConstraintBuilder()
                        .SetGasSpecies(A_g)
                        .SetCondensedSpecies(A_aq)
                        .SetSolvent(h2o)
                        .SetCondensedPhase(aqueous_phase)
                        .SetHenryLawConstant(hlc)
                        .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  std::unordered_map<std::string, std::size_t> si;
  si["A_g"] = 0;
  si["DROP.AQUEOUS.A_aq"] = 1;
  si["DROP.AQUEOUS.H2O"] = 2;

  std::size_t nc = 3;
  auto pi = BuildParamIndices(constraint, phase_prefixes);
  auto sp = InitHlcRt(constraint, phase_prefixes, pi, nc, T);

  double hlc_rt = HLC * micm::constants::GAS_CONSTANT * T;

  DMP sv{ nc, 3, 0.0 };
  sv[0][0] = 1.0e-6;
  sv[0][1] = 0.5;
  sv[0][2] = 300.0;
  sv[1][0] = 2.0e-5;
  sv[1][1] = 1.0;
  sv[1][2] = 200.0;
  sv[2][0] = 5.0e-7;
  sv[2][1] = 0.01;
  sv[2][2] = 0.017;

  auto rf = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, pi, si);
  DMP residual{ nc, 3, 0.0 };
  rf(sv, sp, residual);

  for (std::size_t c = 0; c < nc; ++c)
  {
    double fv = sv[c][2] * water_molecular_weight / water_density;
    double expected = hlc_rt * fv * sv[c][0] - sv[c][1];
    EXPECT_NEAR(residual[c][1], expected, std::max(std::abs(expected) * 1e-12, 1e-20));
    // Non-algebraic rows untouched
    EXPECT_NEAR(residual[c][0], 0.0, 1e-30);
    EXPECT_NEAR(residual[c][2], 0.0, 1e-30);
  }
}

// ── Multiple cells + multiple instances + FD ──

TEST(HenryLawEquilibriumConstraint, MultiInstanceMultiCellFD)
{
  double HLC = 4.0e3;
  double T = 290.0;
  auto hlc = [HLC](const micm::Conditions&) { return HLC; };
  auto constraint = HenryLawEquilibriumConstraintBuilder()
                        .SetGasSpecies(A_g)
                        .SetCondensedSpecies(A_aq)
                        .SetSolvent(h2o)
                        .SetCondensedPhase(aqueous_phase)
                        .SetHenryLawConstant(hlc)
                        .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("LARGE");
  phase_prefixes["AQUEOUS"].insert("SMALL");

  std::unordered_map<std::string, std::size_t> si;
  si["A_g"] = 0;
  si["LARGE.AQUEOUS.A_aq"] = 1;
  si["LARGE.AQUEOUS.H2O"] = 2;
  si["SMALL.AQUEOUS.A_aq"] = 3;
  si["SMALL.AQUEOUS.H2O"] = 4;

  std::size_t nc = 4;
  auto pi = BuildParamIndices(constraint, phase_prefixes);
  auto sp = InitHlcRt(constraint, phase_prefixes, pi, nc, T);

  DMP sv{ nc, 5, 0.0 };
  sv[0][0] = 1e-6;
  sv[0][1] = 0.3;
  sv[0][2] = 300.0;
  sv[0][3] = 0.05;
  sv[0][4] = 100.0;
  sv[1][0] = 5e-5;
  sv[1][1] = 2.0;
  sv[1][2] = 400.0;
  sv[1][3] = 0.8;
  sv[1][4] = 250.0;
  sv[2][0] = 1e-7;
  sv[2][1] = 0.01;
  sv[2][2] = 0.017;
  sv[2][3] = 0.005;
  sv[2][4] = 0.017;
  sv[3][0] = 1e-4;
  sv[3][1] = 10.0;
  sv[3][2] = 500.0;
  sv[3][3] = 5.0;
  sv[3][4] = 350.0;

  CheckConstraintFDJacobian(constraint, phase_prefixes, pi, si, sv, sp);
}

// ── Three instances ──

TEST(HenryLawEquilibriumConstraint, ThreeInstancesFD)
{
  double HLC = 2.0e3;
  double T = 310.0;
  auto hlc = [HLC](const micm::Conditions&) { return HLC; };
  auto constraint = HenryLawEquilibriumConstraintBuilder()
                        .SetGasSpecies(A_g)
                        .SetCondensedSpecies(A_aq)
                        .SetSolvent(h2o)
                        .SetCondensedPhase(aqueous_phase)
                        .SetHenryLawConstant(hlc)
                        .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("A");
  phase_prefixes["AQUEOUS"].insert("B");
  phase_prefixes["AQUEOUS"].insert("C");

  std::unordered_map<std::string, std::size_t> si;
  si["A_g"] = 0;
  si["A.AQUEOUS.A_aq"] = 1;
  si["A.AQUEOUS.H2O"] = 2;
  si["B.AQUEOUS.A_aq"] = 3;
  si["B.AQUEOUS.H2O"] = 4;
  si["C.AQUEOUS.A_aq"] = 5;
  si["C.AQUEOUS.H2O"] = 6;

  std::size_t nc = 2;
  auto pi = BuildParamIndices(constraint, phase_prefixes);
  auto sp = InitHlcRt(constraint, phase_prefixes, pi, nc, T);

  DMP sv{ nc, 7, 0.0 };
  sv[0][0] = 1e-5;
  sv[0][1] = 0.1;
  sv[0][2] = 300.0;
  sv[0][3] = 0.2;
  sv[0][4] = 200.0;
  sv[0][5] = 0.3;
  sv[0][6] = 100.0;
  sv[1][0] = 3e-5;
  sv[1][1] = 0.5;
  sv[1][2] = 400.0;
  sv[1][3] = 1.0;
  sv[1][4] = 350.0;
  sv[1][5] = 1.5;
  sv[1][6] = 250.0;

  CheckConstraintFDJacobian(constraint, phase_prefixes, pi, si, sv, sp);
}

// ── Jacobian accumulates ──

TEST(HenryLawEquilibriumConstraint, JacobianAccumulates)
{
  double HLC = 5.0e3;
  double T = 298.15;
  auto hlc = [HLC](const micm::Conditions&) { return HLC; };
  auto constraint = HenryLawEquilibriumConstraintBuilder()
                        .SetGasSpecies(A_g)
                        .SetCondensedSpecies(A_aq)
                        .SetSolvent(h2o)
                        .SetCondensedPhase(aqueous_phase)
                        .SetHenryLawConstant(hlc)
                        .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  std::unordered_map<std::string, std::size_t> si;
  si["A_g"] = 0;
  si["DROP.AQUEOUS.A_aq"] = 1;
  si["DROP.AQUEOUS.H2O"] = 2;

  auto pi = BuildParamIndices(constraint, phase_prefixes);
  auto sp = InitHlcRt(constraint, phase_prefixes, pi, 1, T);

  auto jacobian = BuildConstraintJacobian(constraint, phase_prefixes, si, 1);
  auto jac_fn = constraint.ConstraintJacobianFunction<DMP, SMP>(phase_prefixes, pi, si, jacobian);

  DMP sv{ 1, 3, 0.0 };
  sv[0][0] = 1.0e-6;
  sv[0][1] = 0.5;
  sv[0][2] = 300.0;

  // First call
  jac_fn(sv, sp, jacobian);
  double j10_once = jacobian.AsVector()[jacobian.VectorIndex(0, 1, 0)];
  double j11_once = jacobian.AsVector()[jacobian.VectorIndex(0, 1, 1)];
  double j12_once = jacobian.AsVector()[jacobian.VectorIndex(0, 1, 2)];

  // Second call accumulates
  jac_fn(sv, sp, jacobian);
  EXPECT_NEAR(jacobian.AsVector()[jacobian.VectorIndex(0, 1, 0)], 2.0 * j10_once, std::abs(j10_once) * 1e-12);
  EXPECT_NEAR(jacobian.AsVector()[jacobian.VectorIndex(0, 1, 1)], 2.0 * j11_once, 1e-12);
  EXPECT_NEAR(jacobian.AsVector()[jacobian.VectorIndex(0, 1, 2)], 2.0 * j12_once, std::abs(j12_once) * 1e-12);
}

// ── Residual sets (does not accumulate) ──

TEST(HenryLawEquilibriumConstraint, ResidualSetsNotAccumulates)
{
  double HLC = 5.0e3;
  double T = 298.15;
  auto hlc = [HLC](const micm::Conditions&) { return HLC; };
  auto constraint = HenryLawEquilibriumConstraintBuilder()
                        .SetGasSpecies(A_g)
                        .SetCondensedSpecies(A_aq)
                        .SetSolvent(h2o)
                        .SetCondensedPhase(aqueous_phase)
                        .SetHenryLawConstant(hlc)
                        .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  std::unordered_map<std::string, std::size_t> si;
  si["A_g"] = 0;
  si["DROP.AQUEOUS.A_aq"] = 1;
  si["DROP.AQUEOUS.H2O"] = 2;

  auto pi = BuildParamIndices(constraint, phase_prefixes);
  auto sp = InitHlcRt(constraint, phase_prefixes, pi, 1, T);

  auto rf = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, pi, si);

  DMP sv{ 1, 3, 0.0 };
  sv[0][0] = 1.0e-6;
  sv[0][1] = 0.5;
  sv[0][2] = 300.0;

  DMP residual{ 1, 3, 999.0 };
  rf(sv, sp, residual);
  double val1 = residual[0][1];

  rf(sv, sp, residual);
  EXPECT_NEAR(residual[0][1], val1, 1e-15);
}

// ── Temperature-dependent HLC with multi-cell ──

TEST(HenryLawEquilibriumConstraint, TemperatureDependentHlcMultiCell)
{
  // HLC(T) = 5000 * exp(2400 * (1/T - 1/298.15))
  auto hlc = [](const micm::Conditions& c) { return 5000.0 * std::exp(2400.0 * (1.0 / c.temperature_ - 1.0 / 298.15)); };
  auto constraint = HenryLawEquilibriumConstraintBuilder()
                        .SetGasSpecies(A_g)
                        .SetCondensedSpecies(A_aq)
                        .SetSolvent(h2o)
                        .SetCondensedPhase(aqueous_phase)
                        .SetHenryLawConstant(hlc)
                        .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  std::unordered_map<std::string, std::size_t> si;
  si["A_g"] = 0;
  si["DROP.AQUEOUS.A_aq"] = 1;
  si["DROP.AQUEOUS.H2O"] = 2;

  std::size_t nc = 3;
  std::vector<micm::Conditions> conditions(nc);
  conditions[0].temperature_ = 270.0;
  conditions[1].temperature_ = 298.15;
  conditions[2].temperature_ = 320.0;
  auto pi = BuildParamIndices(constraint, phase_prefixes);
  auto sp = InitHlcRt(constraint, phase_prefixes, pi, conditions);

  DMP sv{ nc, 3, 0.0 };
  sv[0][0] = 1e-6;
  sv[0][1] = 0.5;
  sv[0][2] = 300.0;
  sv[1][0] = 2e-6;
  sv[1][1] = 1.0;
  sv[1][2] = 250.0;
  sv[2][0] = 5e-7;
  sv[2][1] = 0.1;
  sv[2][2] = 0.017;

  // Check residuals with per-cell HLC*R*T
  auto rf = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, pi, si);
  DMP residual{ nc, 3, 0.0 };
  rf(sv, sp, residual);

  for (std::size_t c = 0; c < nc; ++c)
  {
    double T = conditions[c].temperature_;
    double hlc_val = 5000.0 * std::exp(2400.0 * (1.0 / T - 1.0 / 298.15));
    double hlc_rt = hlc_val * micm::constants::GAS_CONSTANT * T;
    double fv = sv[c][2] * water_molecular_weight / water_density;
    double expected = hlc_rt * fv * sv[c][0] - sv[c][1];
    EXPECT_NEAR(residual[c][1], expected, std::max(std::abs(expected) * 1e-10, 1e-20)) << "Cell " << c;
  }

  // FD Jacobian
  CheckConstraintFDJacobian(constraint, phase_prefixes, pi, si, sv, sp);
}

// ── Cross-instance isolation ──

TEST(HenryLawEquilibriumConstraint, CrossInstanceJacobianIsolation)
{
  double HLC = 5.0e3;
  auto hlc = [HLC](const micm::Conditions&) { return HLC; };
  auto constraint = HenryLawEquilibriumConstraintBuilder()
                        .SetGasSpecies(A_g)
                        .SetCondensedSpecies(A_aq)
                        .SetSolvent(h2o)
                        .SetCondensedPhase(aqueous_phase)
                        .SetHenryLawConstant(hlc)
                        .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("MODE1");
  phase_prefixes["AQUEOUS"].insert("MODE2");

  std::unordered_map<std::string, std::size_t> si;
  si["A_g"] = 0;
  si["MODE1.AQUEOUS.A_aq"] = 1;
  si["MODE1.AQUEOUS.H2O"] = 2;
  si["MODE2.AQUEOUS.A_aq"] = 3;
  si["MODE2.AQUEOUS.H2O"] = 4;

  auto elements = constraint.NonZeroConstraintJacobianElements(phase_prefixes, si);

  // MODE1 algebraic row (1) should only depend on cols 0 (gas), 1 (aq), 2 (solvent)
  for (const auto& [row, col] : elements)
  {
    if (row == 1)
    {
      EXPECT_TRUE(col == 0 || col == 1 || col == 2) << "MODE1 algebraic depends on unexpected col " << col;
    }
    // MODE2 algebraic row (3) should only depend on cols 0 (gas), 3 (aq), 4 (solvent)
    if (row == 3)
    {
      EXPECT_TRUE(col == 0 || col == 3 || col == 4) << "MODE2 algebraic depends on unexpected col " << col;
    }
  }

  // Both instances share dependency on gas column
  bool mode1_has_gas = elements.count({ 1, 0 }) > 0;
  bool mode2_has_gas = elements.count({ 3, 0 }) > 0;
  EXPECT_TRUE(mode1_has_gas);
  EXPECT_TRUE(mode2_has_gas);
}

// ── Analytical Jacobian: multiple instances, verify all entries ──

TEST(HenryLawEquilibriumConstraint, JacobianMultipleInstancesAnalytical)
{
  double HLC = 3.0e3;
  double T = 300.0;
  auto hlc = [HLC](const micm::Conditions&) { return HLC; };
  auto constraint = HenryLawEquilibriumConstraintBuilder()
                        .SetGasSpecies(A_g)
                        .SetCondensedSpecies(A_aq)
                        .SetSolvent(h2o)
                        .SetCondensedPhase(aqueous_phase)
                        .SetHenryLawConstant(hlc)
                        .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("LARGE");
  phase_prefixes["AQUEOUS"].insert("SMALL");

  std::unordered_map<std::string, std::size_t> si;
  si["A_g"] = 0;
  si["LARGE.AQUEOUS.A_aq"] = 1;
  si["LARGE.AQUEOUS.H2O"] = 2;
  si["SMALL.AQUEOUS.A_aq"] = 3;
  si["SMALL.AQUEOUS.H2O"] = 4;

  auto pi = BuildParamIndices(constraint, phase_prefixes);
  auto sp = InitHlcRt(constraint, phase_prefixes, pi, 1, T);

  auto jacobian = BuildConstraintJacobian(constraint, phase_prefixes, si, 1);
  auto jac_fn = constraint.ConstraintJacobianFunction<DMP, SMP>(phase_prefixes, pi, si, jacobian);

  DMP sv{ 1, 5, 0.0 };
  sv[0][0] = 1.0e-5;  // gas
  sv[0][1] = 0.3;     // LARGE aq
  sv[0][2] = 400.0;   // LARGE solvent
  sv[0][3] = 0.1;     // SMALL aq
  sv[0][4] = 200.0;   // SMALL solvent

  jac_fn(sv, sp, jacobian);

  double hlc_rt = HLC * micm::constants::GAS_CONSTANT * T;
  double Mw_rho = water_molecular_weight / water_density;

  // LARGE instance (row 1):
  double fv_large = 400.0 * Mw_rho;
  EXPECT_NEAR(jacobian[0][1][0], -hlc_rt * fv_large, std::abs(hlc_rt * fv_large) * 1e-12);  // dG/d[A_g]
  EXPECT_NEAR(jacobian[0][1][1], 1.0, 1e-12);                                               // dG/d[A_aq] = -1
  EXPECT_NEAR(jacobian[0][1][2], -hlc_rt * Mw_rho * 1.0e-5, std::abs(hlc_rt * Mw_rho * 1e-5) * 1e-12);

  // SMALL instance (row 3):
  double fv_small = 200.0 * Mw_rho;
  EXPECT_NEAR(jacobian[0][3][0], -hlc_rt * fv_small, std::abs(hlc_rt * fv_small) * 1e-12);
  EXPECT_NEAR(jacobian[0][3][3], 1.0, 1e-12);
  EXPECT_NEAR(jacobian[0][3][4], -hlc_rt * Mw_rho * 1.0e-5, std::abs(hlc_rt * Mw_rho * 1e-5) * 1e-12);

  // Cross-instance: LARGE row shouldn't have entries for SMALL columns and vice versa
  EXPECT_THROW(jacobian[0][1][3], std::exception);
  EXPECT_THROW(jacobian[0][1][4], std::exception);
  EXPECT_THROW(jacobian[0][3][1], std::exception);
  EXPECT_THROW(jacobian[0][3][2], std::exception);
}

// ── Large HLC ──

TEST(HenryLawEquilibriumConstraint, LargeHLC)
{
  double HLC = 1.0e8;
  double T = 298.15;
  auto hlc = [HLC](const micm::Conditions&) { return HLC; };
  auto constraint = HenryLawEquilibriumConstraintBuilder()
                        .SetGasSpecies(A_g)
                        .SetCondensedSpecies(A_aq)
                        .SetSolvent(h2o)
                        .SetCondensedPhase(aqueous_phase)
                        .SetHenryLawConstant(hlc)
                        .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  std::unordered_map<std::string, std::size_t> si;
  si["A_g"] = 0;
  si["DROP.AQUEOUS.A_aq"] = 1;
  si["DROP.AQUEOUS.H2O"] = 2;

  auto pi = BuildParamIndices(constraint, phase_prefixes);
  auto sp = InitHlcRt(constraint, phase_prefixes, pi, 1, T);

  DMP sv{ 1, 3, 0.0 };
  sv[0][0] = 1.0e-9;
  sv[0][1] = 0.1;
  sv[0][2] = 0.017;

  auto rf = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, pi, si);
  DMP residual{ 1, 3, 0.0 };
  rf(sv, sp, residual);

  double hlc_rt = HLC * micm::constants::GAS_CONSTANT * T;
  double fv = 0.017 * water_molecular_weight / water_density;
  double expected = hlc_rt * fv * 1.0e-9 - 0.1;
  EXPECT_NEAR(residual[0][1], expected, std::abs(expected) * 1e-10);

  CheckConstraintFDJacobian(constraint, phase_prefixes, pi, si, sv, sp);
}

// ── CopyWithNewUuid preserves behavior ──

TEST(HenryLawEquilibriumConstraint, CopiedConstraintProducesSameResults)
{
  double HLC = 5.0e3;
  double T = 298.15;
  auto hlc = [HLC](const micm::Conditions&) { return HLC; };
  auto original = HenryLawEquilibriumConstraintBuilder()
                      .SetGasSpecies(A_g)
                      .SetCondensedSpecies(A_aq)
                      .SetSolvent(h2o)
                      .SetCondensedPhase(aqueous_phase)
                      .SetHenryLawConstant(hlc)
                      .Build();

  auto copy = original.CopyWithNewUuid();
  EXPECT_NE(copy.uuid_, original.uuid_);

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  std::unordered_map<std::string, std::size_t> si;
  si["A_g"] = 0;
  si["DROP.AQUEOUS.A_aq"] = 1;
  si["DROP.AQUEOUS.H2O"] = 2;

  auto pi_orig = BuildParamIndices(original, phase_prefixes);
  auto pi_copy = BuildParamIndices(copy, phase_prefixes);
  auto sp_orig = InitHlcRt(original, phase_prefixes, pi_orig, 1, T);
  auto sp_copy = InitHlcRt(copy, phase_prefixes, pi_copy, 1, T);

  DMP sv{ 1, 3, 0.0 };
  sv[0][0] = 1.0e-6;
  sv[0][1] = 0.5;
  sv[0][2] = 300.0;

  auto rf_orig = original.ConstraintResidualFunction<DMP>(phase_prefixes, pi_orig, si);
  auto rf_copy = copy.ConstraintResidualFunction<DMP>(phase_prefixes, pi_copy, si);

  DMP res_orig{ 1, 3, 0.0 };
  DMP res_copy{ 1, 3, 0.0 };
  rf_orig(sv, sp_orig, res_orig);
  rf_copy(sv, sp_copy, res_copy);

  EXPECT_NEAR(res_orig[0][1], res_copy[0][1], 1e-15);
}

// ── Equilibrium check: G=0 when at equilibrium ──

TEST(HenryLawEquilibriumConstraint, ResidualZeroAtEquilibrium)
{
  double HLC = 5.0e3;
  double T = 298.15;
  auto hlc = [HLC](const micm::Conditions&) { return HLC; };
  auto constraint = HenryLawEquilibriumConstraintBuilder()
                        .SetGasSpecies(A_g)
                        .SetCondensedSpecies(A_aq)
                        .SetSolvent(h2o)
                        .SetCondensedPhase(aqueous_phase)
                        .SetHenryLawConstant(hlc)
                        .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  std::unordered_map<std::string, std::size_t> si;
  si["A_g"] = 0;
  si["DROP.AQUEOUS.A_aq"] = 1;
  si["DROP.AQUEOUS.H2O"] = 2;

  auto pi = BuildParamIndices(constraint, phase_prefixes);
  auto sp = InitHlcRt(constraint, phase_prefixes, pi, 1, T);

  double hlc_rt = HLC * micm::constants::GAS_CONSTANT * T;
  double solvent_conc = 0.017;
  double fv = solvent_conc * water_molecular_weight / water_density;
  double gas_conc = 1.0e-6;
  double aq_conc_eq = hlc_rt * fv * gas_conc;  // equilibrium: G = 0

  DMP sv{ 1, 3, 0.0 };
  sv[0][0] = gas_conc;
  sv[0][1] = aq_conc_eq;
  sv[0][2] = solvent_conc;

  auto rf = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, pi, si);
  DMP residual{ 1, 3, 0.0 };
  rf(sv, sp, residual);

  EXPECT_NEAR(residual[0][1], 0.0, 1e-15);
}

// ── Kitchen-sink: multi-instance, multi-cell, T-dependent, FD ──

TEST(HenryLawEquilibriumConstraint, KitchenSinkFD)
{
  auto hlc = [](const micm::Conditions& c) { return 3000.0 * std::exp(1500.0 * (1.0 / c.temperature_ - 1.0 / 298.15)); };
  auto constraint = HenryLawEquilibriumConstraintBuilder()
                        .SetGasSpecies(A_g)
                        .SetCondensedSpecies(A_aq)
                        .SetSolvent(h2o)
                        .SetCondensedPhase(aqueous_phase)
                        .SetHenryLawConstant(hlc)
                        .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("M1");
  phase_prefixes["AQUEOUS"].insert("M2");
  phase_prefixes["AQUEOUS"].insert("M3");

  std::unordered_map<std::string, std::size_t> si;
  si["A_g"] = 0;
  si["M1.AQUEOUS.A_aq"] = 1;
  si["M1.AQUEOUS.H2O"] = 2;
  si["M2.AQUEOUS.A_aq"] = 3;
  si["M2.AQUEOUS.H2O"] = 4;
  si["M3.AQUEOUS.A_aq"] = 5;
  si["M3.AQUEOUS.H2O"] = 6;

  std::size_t nc = 3;
  std::vector<micm::Conditions> conditions(nc);
  conditions[0].temperature_ = 275.0;
  conditions[1].temperature_ = 298.15;
  conditions[2].temperature_ = 315.0;
  auto pi = BuildParamIndices(constraint, phase_prefixes);
  auto sp = InitHlcRt(constraint, phase_prefixes, pi, conditions);

  DMP sv{ nc, 7, 0.0 };
  // Cell 0
  sv[0][0] = 1e-6;
  sv[0][1] = 0.1;
  sv[0][2] = 300.0;
  sv[0][3] = 0.2;
  sv[0][4] = 200.0;
  sv[0][5] = 0.05;
  sv[0][6] = 100.0;
  // Cell 1
  sv[1][0] = 5e-5;
  sv[1][1] = 5.0;
  sv[1][2] = 400.0;
  sv[1][3] = 2.0;
  sv[1][4] = 350.0;
  sv[1][5] = 0.5;
  sv[1][6] = 0.017;
  // Cell 2
  sv[2][0] = 1e-7;
  sv[2][1] = 0.01;
  sv[2][2] = 0.017;
  sv[2][3] = 0.005;
  sv[2][4] = 0.017;
  sv[2][5] = 0.001;
  sv[2][6] = 0.017;

  CheckConstraintFDJacobian(constraint, phase_prefixes, pi, si, sv, sp);
}

// ── VectorMatrix typed tests ──
// These tests exercise the fix for the per-cell state-parameter bug:
// with L>1, each cell's HLC*R*T must be stored correctly in the block-ordered
// DenseMatrixPolicy, not in a flat shared_ptr<vector>.

namespace
{
  /// @brief Test HL constraint with a templated DenseMatrixPolicy, num_cells cells,
  ///        with each cell having a distinct temperature (varying_temperature=true
  ///        is the critical regression test for the VectorMatrix bug).
  template<typename VDM, typename VSM>
  void TestHLConstraintVectorMatrix(std::size_t num_cells, double T_base, bool varying_temperature)
  {
    double HLC = 5.0e3;
    auto hlc_fn = [HLC](const micm::Conditions&) { return HLC; };
    auto constraint = HenryLawEquilibriumConstraintBuilder()
                          .SetGasSpecies(A_g)
                          .SetCondensedSpecies(A_aq)
                          .SetSolvent(h2o)
                          .SetCondensedPhase(aqueous_phase)
                          .SetHenryLawConstant(hlc_fn)
                          .Build();

    std::map<std::string, std::set<std::string>> phase_prefixes;
    phase_prefixes["AQUEOUS"].insert("DROP");

    std::unordered_map<std::string, std::size_t> si;
    si["A_g"] = 0;
    si["DROP.AQUEOUS.A_aq"] = 1;
    si["DROP.AQUEOUS.H2O"] = 2;

    // Build param indices using the templated constraint
    std::unordered_map<std::string, std::size_t> pi;
    std::size_t idx = 0;
    for (const auto& name : constraint.ConstraintStateParameterNames(phase_prefixes))
      pi[name] = idx++;
    ASSERT_EQ(pi.size(), 1u);
    std::size_t hlc_rt_col = pi.begin()->second;

    // Build conditions with per-cell temperatures
    std::vector<micm::Conditions> conditions(num_cells);
    for (std::size_t c = 0; c < num_cells; ++c)
      conditions[c].temperature_ = T_base + (varying_temperature ? static_cast<double>(c) * 10.0 : 0.0);

    // Update state parameters
    VDM state_params{ num_cells, std::max(pi.size(), std::size_t{ 1 }), 0.0 };
    auto update_fn = constraint.template UpdateConstraintParametersFunction<VDM>(phase_prefixes, pi);
    update_fn(conditions, state_params);

    // Verify each cell's HLC*R*T value
    for (std::size_t c = 0; c < num_cells; ++c)
    {
      double T = conditions[c].temperature_;
      double expected = HLC * micm::constants::GAS_CONSTANT * T;
      EXPECT_NEAR(state_params[c][hlc_rt_col], expected, expected * 1e-12) << "Cell " << c << " T=" << T;
    }

    // Build state variables with per-cell values
    VDM state_variables{ num_cells, 3, 0.0 };
    for (std::size_t c = 0; c < num_cells; ++c)
    {
      state_variables[c][0] = 1.0e-6 * static_cast<double>(c + 1);    // gas
      state_variables[c][1] = 0.5 * static_cast<double>(c + 1);       // aq
      state_variables[c][2] = 300.0 + static_cast<double>(c) * 10.0;  // solvent
    }

    // Residual: each cell independently
    VDM residual{ num_cells, 3, 0.0 };
    auto rf = constraint.template ConstraintResidualFunction<VDM>(phase_prefixes, pi, si);
    rf(state_variables, state_params, residual);

    for (std::size_t c = 0; c < num_cells; ++c)
    {
      double hlc_rt = state_params[c][hlc_rt_col];
      double solvent = state_variables[c][2];
      double fv = solvent * water_molecular_weight / water_density;
      double gas = state_variables[c][0];
      double aq = state_variables[c][1];
      double expected_res = hlc_rt * fv * gas - aq;
      EXPECT_NEAR(residual[c][1], expected_res, std::max(std::abs(expected_res) * 1e-10, 1e-20))
          << "Residual mismatch cell " << c;
    }

    // Jacobian: build sparse and check via FD
    auto nz = constraint.NonZeroConstraintJacobianElements(phase_prefixes, si);
    auto jac_builder = VSM::Create(si.size()).SetNumberOfBlocks(num_cells).InitialValue(0.0);
    for (const auto& [row, col] : nz)
      jac_builder = jac_builder.WithElement(row, col);
    VSM jacobian(jac_builder);

    auto jac_fn = constraint.template ConstraintJacobianFunction<VDM, VSM>(phase_prefixes, pi, si, jacobian);
    jac_fn(state_variables, state_params, jacobian);

    // FD check per cell
    double eps = 1e-7;
    for (std::size_t c = 0; c < num_cells; ++c)
    {
      for (std::size_t j = 0; j < si.size(); ++j)
      {
        VDM vp(state_variables);
        VDM vm(state_variables);
        double h = std::max(std::abs(state_variables[c][j]) * eps, eps);
        vp[c][j] += h;
        vm[c][j] -= h;

        VDM rp(num_cells, 3, 0.0);
        VDM rm(num_cells, 3, 0.0);
        rf(vp, state_params, rp);
        rf(vm, state_params, rm);

        for (std::size_t i = 0; i < si.size(); ++i)
        {
          double fd = (rp[c][i] - rm[c][i]) / (2.0 * h);
          double analytical;
          try
          {
            analytical = jacobian[c][i][j];
          }
          catch (...)
          {
            if (std::abs(fd) > 1e-10)
              ADD_FAILURE() << "Missing Jacobian cell=" << c << " row=" << i << " col=" << j << " fd=" << fd;
            continue;
          }
          double scale = std::max(std::abs(analytical), std::abs(fd));
          if (scale > 1e-15)
          {
            EXPECT_NEAR(analytical + fd, 0.0, scale * 1e-5) << "FD Jac mismatch cell=" << c << " row=" << i << " col=" << j;
          }
        }
      }
    }
  }
}  // namespace

TEST(HenryLawEquilibriumConstraint, VectorMatrix_L1_4cells)
{
  using VDM = micm::VectorMatrix<double, 1>;
  using VSM = micm::SparseMatrix<double, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>;
  TestHLConstraintVectorMatrix<VDM, VSM>(4, 298.15, true);
}

TEST(HenryLawEquilibriumConstraint, VectorMatrix_L2_4cells)
{
  using VDM = micm::VectorMatrix<double, 2>;
  using VSM = micm::SparseMatrix<double, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>;
  TestHLConstraintVectorMatrix<VDM, VSM>(4, 298.15, true);
}

TEST(HenryLawEquilibriumConstraint, VectorMatrix_L4_4cells)
{
  using VDM = micm::VectorMatrix<double, 4>;
  using VSM = micm::SparseMatrix<double, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>;
  TestHLConstraintVectorMatrix<VDM, VSM>(4, 298.15, true);
}

// ============================================================================
// Limit / Extreme tests (Phase D2)
// ============================================================================

TEST(HenryLawEquilibriumConstraint, ResidualZeroGasConcentration)
{
  // G = HLC*R*T*f_v*[A_g] - [A_aq]; with [A_g]=0: G = -[A_aq]
  double HLC = 5.0e3;
  double T = 298.15;
  auto hlc = [HLC](const micm::Conditions&) { return HLC; };
  auto constraint = HenryLawEquilibriumConstraintBuilder()
                        .SetGasSpecies(A_g)
                        .SetCondensedSpecies(A_aq)
                        .SetSolvent(h2o)
                        .SetCondensedPhase(aqueous_phase)
                        .SetHenryLawConstant(hlc)
                        .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  std::unordered_map<std::string, std::size_t> si;
  si["A_g"] = 0;
  si["DROP.AQUEOUS.A_aq"] = 1;
  si["DROP.AQUEOUS.H2O"] = 2;

  auto pi = BuildParamIndices(constraint, phase_prefixes);
  auto sp = InitHlcRt(constraint, phase_prefixes, pi, 1, T);

  double A_aq_conc = 0.5;
  DMP sv{ 1, 3, 0.0 };
  sv[0][0] = 0.0;  // [A_g] = 0
  sv[0][1] = A_aq_conc;
  sv[0][2] = 300.0;

  auto rf = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, pi, si);
  DMP residual{ 1, 3, 0.0 };
  rf(sv, sp, residual);

  EXPECT_NEAR(residual[0][1], -A_aq_conc, 1.0e-12);
}

TEST(HenryLawEquilibriumConstraint, ResidualZeroAqueousConcentration)
{
  // G = HLC*R*T*f_v*[A_g] - [A_aq]; with [A_aq]=0: G = HLC*R*T*f_v*[A_g]
  double HLC = 5.0e3;
  double T = 298.15;
  auto hlc = [HLC](const micm::Conditions&) { return HLC; };
  auto constraint = HenryLawEquilibriumConstraintBuilder()
                        .SetGasSpecies(A_g)
                        .SetCondensedSpecies(A_aq)
                        .SetSolvent(h2o)
                        .SetCondensedPhase(aqueous_phase)
                        .SetHenryLawConstant(hlc)
                        .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  std::unordered_map<std::string, std::size_t> si;
  si["A_g"] = 0;
  si["DROP.AQUEOUS.A_aq"] = 1;
  si["DROP.AQUEOUS.H2O"] = 2;

  auto pi = BuildParamIndices(constraint, phase_prefixes);
  auto sp = InitHlcRt(constraint, phase_prefixes, pi, 1, T);

  double A_g_conc = 1.0e-6;
  double H2O_conc = 300.0;
  DMP sv{ 1, 3, 0.0 };
  sv[0][0] = A_g_conc;
  sv[0][1] = 0.0;  // [A_aq] = 0
  sv[0][2] = H2O_conc;

  auto rf = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, pi, si);
  DMP residual{ 1, 3, 0.0 };
  rf(sv, sp, residual);

  double f_v = H2O_conc * water_molecular_weight / water_density;
  double hlc_rt = HLC * micm::constants::GAS_CONSTANT * T;
  EXPECT_NEAR(residual[0][1], hlc_rt * f_v * A_g_conc, 1.0e-10);
  EXPECT_GT(residual[0][1], 0.0);  // positive: equilibrium favors aqueous phase
}

TEST(HenryLawEquilibriumConstraint, JacobianFDZeroConcentrations)
{
  double HLC = 5.0e3;
  double T = 298.15;
  auto hlc = [HLC](const micm::Conditions&) { return HLC; };
  auto constraint = HenryLawEquilibriumConstraintBuilder()
                        .SetGasSpecies(A_g)
                        .SetCondensedSpecies(A_aq)
                        .SetSolvent(h2o)
                        .SetCondensedPhase(aqueous_phase)
                        .SetHenryLawConstant(hlc)
                        .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  std::unordered_map<std::string, std::size_t> si;
  si["A_g"] = 0;
  si["DROP.AQUEOUS.A_aq"] = 1;
  si["DROP.AQUEOUS.H2O"] = 2;

  auto pi = BuildParamIndices(constraint, phase_prefixes);
  auto sp = InitHlcRt(constraint, phase_prefixes, pi, 1, T);

  // [A_g]=0, [A_aq] non-zero
  {
    DMP sv{ 1, 3, 0.0 };
    sv[0][0] = 0.0;
    sv[0][1] = 0.5;
    sv[0][2] = 300.0;
    CheckConstraintFDJacobian(constraint, phase_prefixes, pi, si, sv, sp);
  }

  // [A_aq]=0, [A_g] non-zero
  {
    DMP sv{ 1, 3, 0.0 };
    sv[0][0] = 1.0e-6;
    sv[0][1] = 0.0;
    sv[0][2] = 300.0;
    CheckConstraintFDJacobian(constraint, phase_prefixes, pi, si, sv, sp);
  }
}

TEST(HenryLawEquilibriumConstraint, JacobianFDTemperatureExtremes)
{
  double HLC = 5.0e3;
  auto hlc = [HLC](const micm::Conditions&) { return HLC; };
  auto constraint = HenryLawEquilibriumConstraintBuilder()
                        .SetGasSpecies(A_g)
                        .SetCondensedSpecies(A_aq)
                        .SetSolvent(h2o)
                        .SetCondensedPhase(aqueous_phase)
                        .SetHenryLawConstant(hlc)
                        .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  std::unordered_map<std::string, std::size_t> si;
  si["A_g"] = 0;
  si["DROP.AQUEOUS.A_aq"] = 1;
  si["DROP.AQUEOUS.H2O"] = 2;

  DMP sv{ 1, 3, 0.0 };
  sv[0][0] = 1.0e-6;
  sv[0][1] = 0.5;
  sv[0][2] = 300.0;

  for (double T : { 200.0, 298.15, 350.0 })
  {
    auto pi = BuildParamIndices(constraint, phase_prefixes);
    auto sp = InitHlcRt(constraint, phase_prefixes, pi, 1, T);
    CheckConstraintFDJacobian(constraint, phase_prefixes, pi, si, sv, sp);
  }
}
