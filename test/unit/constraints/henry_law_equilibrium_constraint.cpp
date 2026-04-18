// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <miam/constraints/henry_law_equilibrium_constraint.hpp>
#include <miam/constraints/henry_law_equilibrium_constraint_builder.hpp>
#include <miam/util/condensation_rate.hpp>
#include <micm/system/conditions.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_standard_ordering.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

using namespace miam::constraint;

namespace
{
  using DMP = micm::Matrix<double>;
  std::unordered_map<std::string, std::size_t> param_indices;
  const DMP no_params{};

  auto h2o = micm::Species{ "H2O" };
  auto A_g = micm::Species{ "A_g" };
  auto A_aq = micm::Species{ "A_aq" };
  auto aqueous_phase = micm::Phase{ "AQUEOUS", { { h2o }, { A_aq } } };
  constexpr double Mw_water = 0.018;
  constexpr double rho_water = 1000.0;
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
      .SetMwSolvent(Mw_water)
      .SetRhoSolvent(rho_water)
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
      .SetMwSolvent(Mw_water)
      .SetRhoSolvent(rho_water)
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
      .SetMwSolvent(Mw_water)
      .SetRhoSolvent(rho_water)
      .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  auto deps = constraint.ConstraintSpeciesDependencies(phase_prefixes);
  // Gas species (standalone) + condensed species + solvent
  EXPECT_EQ(deps.size(), 3);
  EXPECT_TRUE(deps.count("A_g"));                     // gas, no prefix
  EXPECT_TRUE(deps.count("DROP.AQUEOUS.A_aq"));       // condensed
  EXPECT_TRUE(deps.count("DROP.AQUEOUS.H2O"));        // solvent
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
      .SetMwSolvent(Mw_water)
      .SetRhoSolvent(rho_water)
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
      .SetMwSolvent(Mw_water)
      .SetRhoSolvent(rho_water)
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
      .SetMwSolvent(Mw_water)
      .SetRhoSolvent(rho_water)
      .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  std::unordered_map<std::string, std::size_t> state_indices;
  state_indices["A_g"] = 0;
  state_indices["DROP.AQUEOUS.A_aq"] = 1;
  state_indices["DROP.AQUEOUS.H2O"] = 2;

  // Initialize HLC*R*T
  std::vector<micm::Conditions> conditions(1);
  conditions[0].temperature_ = T;
  auto update_fn = constraint.UpdateConstraintParametersFunction<micm::Matrix<double>>();
  update_fn(conditions);

  using DMP = micm::Matrix<double>;
  auto residual_fn = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, param_indices, state_indices);

  double gas_conc = 1.0e-6;
  double aq_conc = 0.5;
  double solvent_conc = 300.0;

  DMP state_variables{ 1, 3, 0.0 };
  state_variables[0][0] = gas_conc;
  state_variables[0][1] = aq_conc;
  state_variables[0][2] = solvent_conc;

  DMP residual{ 1, 3, 0.0 };
  residual_fn(state_variables, no_params, residual);

  double f_v = solvent_conc * Mw_water / rho_water;
  double hlc_rt = HLC * miam::util::R_gas * T;
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
      .SetMwSolvent(Mw_water)
      .SetRhoSolvent(rho_water)
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

  std::vector<micm::Conditions> conditions(1);
  conditions[0].temperature_ = T;
  auto update_fn = constraint.UpdateConstraintParametersFunction<micm::Matrix<double>>();
  update_fn(conditions);

  using DMP = micm::Matrix<double>;
  auto residual_fn = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, param_indices, state_indices);

  DMP state_variables{ 1, 5, 0.0 };
  state_variables[0][0] = 1.0e-5;  // [A_g]
  state_variables[0][1] = 0.3;     // LARGE [A_aq]
  state_variables[0][2] = 400.0;   // LARGE [H2O]
  state_variables[0][3] = 0.1;     // SMALL [A_aq]
  state_variables[0][4] = 200.0;   // SMALL [H2O]

  DMP residual{ 1, 5, 0.0 };
  residual_fn(state_variables, no_params, residual);

  double hlc_rt = HLC * miam::util::R_gas * T;
  double gas_conc = 1.0e-5;

  // LARGE: f_v = 400.0 * 0.018 / 1000.0 = 0.0072
  double fv_large = 400.0 * Mw_water / rho_water;
  EXPECT_NEAR(residual[0][1], hlc_rt * fv_large * gas_conc - 0.3, std::abs(hlc_rt * fv_large * gas_conc) * 1.0e-10);

  // SMALL: f_v = 200.0 * 0.018 / 1000.0 = 0.0036
  double fv_small = 200.0 * Mw_water / rho_water;
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
      .SetMwSolvent(Mw_water)
      .SetRhoSolvent(rho_water)
      .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  std::unordered_map<std::string, std::size_t> state_indices;
  state_indices["A_g"] = 0;
  state_indices["DROP.AQUEOUS.A_aq"] = 1;
  state_indices["DROP.AQUEOUS.H2O"] = 2;

  std::vector<micm::Conditions> conditions(1);
  conditions[0].temperature_ = T;
  auto update_fn = constraint.UpdateConstraintParametersFunction<micm::Matrix<double>>();
  update_fn(conditions);

  using DMP = micm::Matrix<double>;
  using SMP = micm::SparseMatrix<double, micm::SparseMatrixStandardOrderingCompressedSparseRow>;
  auto nz_elements = constraint.NonZeroConstraintJacobianElements(phase_prefixes, state_indices);
  auto builder = SMP::Create(3).SetNumberOfBlocks(1);
  for (const auto& [row, col] : nz_elements)
    builder.WithElement(row, col);
  SMP jacobian(builder);

  auto jac_fn = constraint.ConstraintJacobianFunction<DMP, SMP>(phase_prefixes, state_indices, jacobian);

  double gas_conc = 1.0e-6;
  double solvent_conc = 300.0;

  DMP state_variables{ 1, 3, 0.0 };
  state_variables[0][0] = gas_conc;
  state_variables[0][1] = 0.5;          // [A_aq]
  state_variables[0][2] = solvent_conc;  // [H2O]

  for (auto& v : jacobian.AsVector())
    v = 0.0;
  jac_fn(state_variables, jacobian);

  double hlc_rt = HLC * miam::util::R_gas * T;
  double f_v = solvent_conc * Mw_water / rho_water;
  double Mw_rho = Mw_water / rho_water;

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
  // HLC(T) = 1000.0 / T
  auto hlc = [](const micm::Conditions& c) { return 1000.0 / c.temperature_; };
  auto constraint = HenryLawEquilibriumConstraintBuilder()
      .SetGasSpecies(A_g)
      .SetCondensedSpecies(A_aq)
      .SetSolvent(h2o)
      .SetCondensedPhase(aqueous_phase)
      .SetHenryLawConstant(hlc)
      .SetMwSolvent(Mw_water)
      .SetRhoSolvent(rho_water)
      .Build();

  auto update_fn = constraint.UpdateConstraintParametersFunction<micm::Matrix<double>>();

  std::vector<micm::Conditions> conditions(2);
  conditions[0].temperature_ = 250.0;
  conditions[1].temperature_ = 350.0;
  update_fn(conditions);

  EXPECT_EQ(constraint.hlc_rt_values_->size(), 2);
  // HLC*R*T = (1000/T) * R * T = 1000 * R
  double expected = 1000.0 * miam::util::R_gas;
  EXPECT_NEAR((*constraint.hlc_rt_values_)[0], expected, expected * 1.0e-12);
  EXPECT_NEAR((*constraint.hlc_rt_values_)[1], expected, expected * 1.0e-12);
}

// ── Builder ──

TEST(HenryLawEquilibriumConstraint, BuilderValidation)
{
  struct FakeHLC
  {
    double Calculate(const micm::Conditions&) const { return 5.0e3; }
  };

  // Missing gas species
  EXPECT_THROW(
      HenryLawEquilibriumConstraintBuilder()
          .SetCondensedSpecies(A_aq)
          .SetSolvent(h2o)
          .SetCondensedPhase(aqueous_phase)
          .SetHenryLawConstant(FakeHLC{})
          .SetMwSolvent(Mw_water)
          .SetRhoSolvent(rho_water)
          .Build(),
      std::runtime_error);

  // Valid build
  auto constraint = HenryLawEquilibriumConstraintBuilder()
                        .SetGasSpecies(A_g)
                        .SetCondensedSpecies(A_aq)
                        .SetSolvent(h2o)
                        .SetCondensedPhase(aqueous_phase)
                        .SetHenryLawConstant(FakeHLC{})
                        .SetMwSolvent(Mw_water)
                        .SetRhoSolvent(rho_water)
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

  void InitHlcRt(HenryLawEquilibriumConstraint& constraint, std::size_t num_cells, double T = 298.15)
  {
    std::vector<micm::Conditions> conditions(num_cells);
    for (auto& c : conditions)
      c.temperature_ = T;
    auto update_fn = constraint.UpdateConstraintParametersFunction<DMP>();
    update_fn(conditions);
  }

  void InitHlcRt(
      HenryLawEquilibriumConstraint& constraint,
      const std::vector<micm::Conditions>& conditions)
  {
    auto update_fn = constraint.UpdateConstraintParametersFunction<DMP>();
    update_fn(conditions);
  }

  /// @brief Finite-difference check for the constraint Jacobian
  void CheckConstraintFDJacobian(
      HenryLawEquilibriumConstraint& constraint,
      const std::map<std::string, std::set<std::string>>& phase_prefixes,
      const std::unordered_map<std::string, std::size_t>& state_indices,
      const DMP& state_variables,
      double rel_tol = 1e-5,
      double abs_tol = 1e-8)
  {
    std::size_t num_blocks = state_variables.NumRows();
    std::size_t num_vars = state_indices.size();

    // Analytical Jacobian
    auto jacobian = BuildConstraintJacobian(constraint, phase_prefixes, state_indices, num_blocks);
    auto jac_fn = constraint.ConstraintJacobianFunction<DMP, SMP>(phase_prefixes, state_indices, jacobian);
    jac_fn(state_variables, jacobian);

    // Central-difference FD
    double eps = 1e-7;
    for (std::size_t j = 0; j < num_vars; ++j)
    {
      DMP vars_plus(state_variables);
      DMP vars_minus(state_variables);
      for (std::size_t b = 0; b < num_blocks; ++b)
      {
        double h = std::max(std::abs(state_variables[b][j]) * eps, eps);
        vars_plus[b][j] += h;
        vars_minus[b][j] -= h;
      }

      DMP res_plus(num_blocks, num_vars, 0.0);
      DMP res_minus(num_blocks, num_vars, 0.0);
      auto rf_plus = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, param_indices, state_indices);
      auto rf_minus = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, param_indices, state_indices);
      rf_plus(vars_plus, no_params, res_plus);
      rf_minus(vars_minus, no_params, res_minus);

      for (std::size_t b = 0; b < num_blocks; ++b)
      {
        double h = std::max(std::abs(state_variables[b][j]) * eps, eps);
        for (std::size_t i = 0; i < num_vars; ++i)
        {
          double fd = (res_plus[b][i] - res_minus[b][i]) / (2.0 * h);
          double analytical;
          try
          {
            analytical = jacobian[b][i][j];
          }
          catch (...)
          {
            if (std::abs(fd) > 1e-10)
              ADD_FAILURE() << "Missing Jacobian at block=" << b << " row=" << i << " col=" << j << " fd=" << fd;
            continue;
          }

          double scale = std::max(std::abs(analytical), std::abs(fd));
          if (scale > 1e-15)
          {
            double tol = std::max(scale * rel_tol, abs_tol);
            EXPECT_NEAR(analytical + fd, 0.0, tol)
                << "FD mismatch: block=" << b << " row=" << i << " col=" << j
                << " analytical(-dG/dy)=" << analytical << " fd(dG/dy)=" << fd;
          }
        }
      }
    }
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
      .SetMwSolvent(Mw_water)
      .SetRhoSolvent(rho_water)
      .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  std::unordered_map<std::string, std::size_t> si;
  si["A_g"] = 0;
  si["DROP.AQUEOUS.A_aq"] = 1;
  si["DROP.AQUEOUS.H2O"] = 2;

  InitHlcRt(constraint, 1, T);

  DMP sv{ 1, 3, 0.0 };
  sv[0][0] = 1.0e-6;
  sv[0][1] = 0.5;
  sv[0][2] = 300.0;

  CheckConstraintFDJacobian(constraint, phase_prefixes, si, sv);
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
      .SetMwSolvent(Mw_water)
      .SetRhoSolvent(rho_water)
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

  InitHlcRt(constraint, 1, T);

  DMP sv{ 1, 5, 0.0 };
  sv[0][0] = 1.0e-5;
  sv[0][1] = 0.3;
  sv[0][2] = 400.0;
  sv[0][3] = 0.1;
  sv[0][4] = 200.0;

  CheckConstraintFDJacobian(constraint, phase_prefixes, si, sv);
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
      .SetMwSolvent(Mw_water)
      .SetRhoSolvent(rho_water)
      .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  std::unordered_map<std::string, std::size_t> si;
  si["A_g"] = 0;
  si["DROP.AQUEOUS.A_aq"] = 1;
  si["DROP.AQUEOUS.H2O"] = 2;

  std::size_t nc = 3;
  InitHlcRt(constraint, nc, T);

  double hlc_rt = HLC * miam::util::R_gas * T;

  DMP sv{ nc, 3, 0.0 };
  sv[0][0] = 1.0e-6;  sv[0][1] = 0.5;  sv[0][2] = 300.0;
  sv[1][0] = 2.0e-5;  sv[1][1] = 1.0;  sv[1][2] = 200.0;
  sv[2][0] = 5.0e-7;  sv[2][1] = 0.01; sv[2][2] = 0.017;

  auto rf = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, param_indices, si);
  DMP residual{ nc, 3, 0.0 };
  rf(sv, no_params, residual);

  for (std::size_t c = 0; c < nc; ++c)
  {
    double fv = sv[c][2] * Mw_water / rho_water;
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
      .SetMwSolvent(Mw_water)
      .SetRhoSolvent(rho_water)
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
  InitHlcRt(constraint, nc, T);

  DMP sv{ nc, 5, 0.0 };
  sv[0][0] = 1e-6;   sv[0][1] = 0.3;   sv[0][2] = 300.0;  sv[0][3] = 0.05;  sv[0][4] = 100.0;
  sv[1][0] = 5e-5;   sv[1][1] = 2.0;   sv[1][2] = 400.0;  sv[1][3] = 0.8;   sv[1][4] = 250.0;
  sv[2][0] = 1e-7;   sv[2][1] = 0.01;  sv[2][2] = 0.017;   sv[2][3] = 0.005; sv[2][4] = 0.017;
  sv[3][0] = 1e-4;   sv[3][1] = 10.0;  sv[3][2] = 500.0;  sv[3][3] = 5.0;   sv[3][4] = 350.0;

  CheckConstraintFDJacobian(constraint, phase_prefixes, si, sv);
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
      .SetMwSolvent(Mw_water)
      .SetRhoSolvent(rho_water)
      .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("A");
  phase_prefixes["AQUEOUS"].insert("B");
  phase_prefixes["AQUEOUS"].insert("C");

  std::unordered_map<std::string, std::size_t> si;
  si["A_g"] = 0;
  si["A.AQUEOUS.A_aq"] = 1;  si["A.AQUEOUS.H2O"] = 2;
  si["B.AQUEOUS.A_aq"] = 3;  si["B.AQUEOUS.H2O"] = 4;
  si["C.AQUEOUS.A_aq"] = 5;  si["C.AQUEOUS.H2O"] = 6;

  std::size_t nc = 2;
  InitHlcRt(constraint, nc, T);

  DMP sv{ nc, 7, 0.0 };
  sv[0][0] = 1e-5;
  sv[0][1] = 0.1;  sv[0][2] = 300.0;
  sv[0][3] = 0.2;  sv[0][4] = 200.0;
  sv[0][5] = 0.3;  sv[0][6] = 100.0;
  sv[1][0] = 3e-5;
  sv[1][1] = 0.5;  sv[1][2] = 400.0;
  sv[1][3] = 1.0;  sv[1][4] = 350.0;
  sv[1][5] = 1.5;  sv[1][6] = 250.0;

  CheckConstraintFDJacobian(constraint, phase_prefixes, si, sv);
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
      .SetMwSolvent(Mw_water)
      .SetRhoSolvent(rho_water)
      .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  std::unordered_map<std::string, std::size_t> si;
  si["A_g"] = 0;
  si["DROP.AQUEOUS.A_aq"] = 1;
  si["DROP.AQUEOUS.H2O"] = 2;

  InitHlcRt(constraint, 1, T);

  auto jacobian = BuildConstraintJacobian(constraint, phase_prefixes, si, 1);
  auto jac_fn = constraint.ConstraintJacobianFunction<DMP, SMP>(phase_prefixes, si, jacobian);

  DMP sv{ 1, 3, 0.0 };
  sv[0][0] = 1.0e-6;  sv[0][1] = 0.5;  sv[0][2] = 300.0;

  // First call
  jac_fn(sv, jacobian);
  double j10_once = jacobian.AsVector()[jacobian.VectorIndex(0, 1, 0)];
  double j11_once = jacobian.AsVector()[jacobian.VectorIndex(0, 1, 1)];
  double j12_once = jacobian.AsVector()[jacobian.VectorIndex(0, 1, 2)];

  // Second call accumulates
  jac_fn(sv, jacobian);
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
      .SetMwSolvent(Mw_water)
      .SetRhoSolvent(rho_water)
      .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  std::unordered_map<std::string, std::size_t> si;
  si["A_g"] = 0;
  si["DROP.AQUEOUS.A_aq"] = 1;
  si["DROP.AQUEOUS.H2O"] = 2;

  InitHlcRt(constraint, 1, T);

  auto rf = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, param_indices, si);

  DMP sv{ 1, 3, 0.0 };
  sv[0][0] = 1.0e-6;  sv[0][1] = 0.5;  sv[0][2] = 300.0;

  DMP residual{ 1, 3, 999.0 };
  rf(sv, no_params, residual);
  double val1 = residual[0][1];

  rf(sv, no_params, residual);
  EXPECT_NEAR(residual[0][1], val1, 1e-15);
}

// ── Temperature-dependent HLC with multi-cell ──

TEST(HenryLawEquilibriumConstraint, TemperatureDependentHlcMultiCell)
{
  // HLC(T) = 5000 * exp(2400 * (1/T - 1/298.15))
  auto hlc = [](const micm::Conditions& c)
  { return 5000.0 * std::exp(2400.0 * (1.0 / c.temperature_ - 1.0 / 298.15)); };
  auto constraint = HenryLawEquilibriumConstraintBuilder()
      .SetGasSpecies(A_g)
      .SetCondensedSpecies(A_aq)
      .SetSolvent(h2o)
      .SetCondensedPhase(aqueous_phase)
      .SetHenryLawConstant(hlc)
      .SetMwSolvent(Mw_water)
      .SetRhoSolvent(rho_water)
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
  InitHlcRt(constraint, conditions);

  DMP sv{ nc, 3, 0.0 };
  sv[0][0] = 1e-6;  sv[0][1] = 0.5;  sv[0][2] = 300.0;
  sv[1][0] = 2e-6;  sv[1][1] = 1.0;  sv[1][2] = 250.0;
  sv[2][0] = 5e-7;  sv[2][1] = 0.1;  sv[2][2] = 0.017;

  // Check residuals with per-cell HLC*R*T
  auto rf = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, param_indices, si);
  DMP residual{ nc, 3, 0.0 };
  rf(sv, no_params, residual);

  for (std::size_t c = 0; c < nc; ++c)
  {
    double T = conditions[c].temperature_;
    double hlc_val = 5000.0 * std::exp(2400.0 * (1.0 / T - 1.0 / 298.15));
    double hlc_rt = hlc_val * miam::util::R_gas * T;
    double fv = sv[c][2] * Mw_water / rho_water;
    double expected = hlc_rt * fv * sv[c][0] - sv[c][1];
    EXPECT_NEAR(residual[c][1], expected, std::max(std::abs(expected) * 1e-10, 1e-20))
        << "Cell " << c;
  }

  // FD Jacobian
  CheckConstraintFDJacobian(constraint, phase_prefixes, si, sv);
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
      .SetMwSolvent(Mw_water)
      .SetRhoSolvent(rho_water)
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
      EXPECT_TRUE(col == 0 || col == 1 || col == 2)
          << "MODE1 algebraic depends on unexpected col " << col;
    }
    // MODE2 algebraic row (3) should only depend on cols 0 (gas), 3 (aq), 4 (solvent)
    if (row == 3)
    {
      EXPECT_TRUE(col == 0 || col == 3 || col == 4)
          << "MODE2 algebraic depends on unexpected col " << col;
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
      .SetMwSolvent(Mw_water)
      .SetRhoSolvent(rho_water)
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

  InitHlcRt(constraint, 1, T);

  auto jacobian = BuildConstraintJacobian(constraint, phase_prefixes, si, 1);
  auto jac_fn = constraint.ConstraintJacobianFunction<DMP, SMP>(phase_prefixes, si, jacobian);

  DMP sv{ 1, 5, 0.0 };
  sv[0][0] = 1.0e-5;   // gas
  sv[0][1] = 0.3;      // LARGE aq
  sv[0][2] = 400.0;    // LARGE solvent
  sv[0][3] = 0.1;      // SMALL aq
  sv[0][4] = 200.0;    // SMALL solvent

  jac_fn(sv, jacobian);

  double hlc_rt = HLC * miam::util::R_gas * T;
  double Mw_rho = Mw_water / rho_water;

  // LARGE instance (row 1):
  double fv_large = 400.0 * Mw_rho;
  EXPECT_NEAR(jacobian[0][1][0], -hlc_rt * fv_large, std::abs(hlc_rt * fv_large) * 1e-12);  // dG/d[A_g]
  EXPECT_NEAR(jacobian[0][1][1], 1.0, 1e-12);                                                // dG/d[A_aq] = -1
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
      .SetMwSolvent(Mw_water)
      .SetRhoSolvent(rho_water)
      .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  std::unordered_map<std::string, std::size_t> si;
  si["A_g"] = 0;
  si["DROP.AQUEOUS.A_aq"] = 1;
  si["DROP.AQUEOUS.H2O"] = 2;

  InitHlcRt(constraint, 1, T);

  DMP sv{ 1, 3, 0.0 };
  sv[0][0] = 1.0e-9;
  sv[0][1] = 0.1;
  sv[0][2] = 0.017;

  auto rf = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, param_indices, si);
  DMP residual{ 1, 3, 0.0 };
  rf(sv, no_params, residual);

  double hlc_rt = HLC * miam::util::R_gas * T;
  double fv = 0.017 * Mw_water / rho_water;
  double expected = hlc_rt * fv * 1.0e-9 - 0.1;
  EXPECT_NEAR(residual[0][1], expected, std::abs(expected) * 1e-10);

  CheckConstraintFDJacobian(constraint, phase_prefixes, si, sv);
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
      .SetMwSolvent(Mw_water)
      .SetRhoSolvent(rho_water)
      .Build();

  auto copy = original.CopyWithNewUuid();
  EXPECT_NE(copy.uuid_, original.uuid_);

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  std::unordered_map<std::string, std::size_t> si;
  si["A_g"] = 0;
  si["DROP.AQUEOUS.A_aq"] = 1;
  si["DROP.AQUEOUS.H2O"] = 2;

  InitHlcRt(original, 1, T);
  InitHlcRt(copy, 1, T);

  DMP sv{ 1, 3, 0.0 };
  sv[0][0] = 1.0e-6;  sv[0][1] = 0.5;  sv[0][2] = 300.0;

  auto rf_orig = original.ConstraintResidualFunction<DMP>(phase_prefixes, param_indices, si);
  auto rf_copy = copy.ConstraintResidualFunction<DMP>(phase_prefixes, param_indices, si);

  DMP res_orig{ 1, 3, 0.0 };
  DMP res_copy{ 1, 3, 0.0 };
  rf_orig(sv, no_params, res_orig);
  rf_copy(sv, no_params, res_copy);

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
      .SetMwSolvent(Mw_water)
      .SetRhoSolvent(rho_water)
      .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  std::unordered_map<std::string, std::size_t> si;
  si["A_g"] = 0;
  si["DROP.AQUEOUS.A_aq"] = 1;
  si["DROP.AQUEOUS.H2O"] = 2;

  InitHlcRt(constraint, 1, T);

  double hlc_rt = HLC * miam::util::R_gas * T;
  double solvent_conc = 0.017;
  double fv = solvent_conc * Mw_water / rho_water;
  double gas_conc = 1.0e-6;
  double aq_conc_eq = hlc_rt * fv * gas_conc;  // equilibrium: G = 0

  DMP sv{ 1, 3, 0.0 };
  sv[0][0] = gas_conc;
  sv[0][1] = aq_conc_eq;
  sv[0][2] = solvent_conc;

  auto rf = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, param_indices, si);
  DMP residual{ 1, 3, 0.0 };
  rf(sv, no_params, residual);

  EXPECT_NEAR(residual[0][1], 0.0, 1e-15);
}

// ── Kitchen-sink: multi-instance, multi-cell, T-dependent, FD ──

TEST(HenryLawEquilibriumConstraint, KitchenSinkFD)
{
  auto hlc = [](const micm::Conditions& c)
  { return 3000.0 * std::exp(1500.0 * (1.0 / c.temperature_ - 1.0 / 298.15)); };
  auto constraint = HenryLawEquilibriumConstraintBuilder()
      .SetGasSpecies(A_g)
      .SetCondensedSpecies(A_aq)
      .SetSolvent(h2o)
      .SetCondensedPhase(aqueous_phase)
      .SetHenryLawConstant(hlc)
      .SetMwSolvent(Mw_water)
      .SetRhoSolvent(rho_water)
      .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("M1");
  phase_prefixes["AQUEOUS"].insert("M2");
  phase_prefixes["AQUEOUS"].insert("M3");

  std::unordered_map<std::string, std::size_t> si;
  si["A_g"] = 0;
  si["M1.AQUEOUS.A_aq"] = 1;  si["M1.AQUEOUS.H2O"] = 2;
  si["M2.AQUEOUS.A_aq"] = 3;  si["M2.AQUEOUS.H2O"] = 4;
  si["M3.AQUEOUS.A_aq"] = 5;  si["M3.AQUEOUS.H2O"] = 6;

  std::size_t nc = 3;
  std::vector<micm::Conditions> conditions(nc);
  conditions[0].temperature_ = 275.0;
  conditions[1].temperature_ = 298.15;
  conditions[2].temperature_ = 315.0;
  InitHlcRt(constraint, conditions);

  DMP sv{ nc, 7, 0.0 };
  // Cell 0
  sv[0][0] = 1e-6;
  sv[0][1] = 0.1;  sv[0][2] = 300.0;
  sv[0][3] = 0.2;  sv[0][4] = 200.0;
  sv[0][5] = 0.05; sv[0][6] = 100.0;
  // Cell 1
  sv[1][0] = 5e-5;
  sv[1][1] = 5.0;  sv[1][2] = 400.0;
  sv[1][3] = 2.0;  sv[1][4] = 350.0;
  sv[1][5] = 0.5;  sv[1][6] = 0.017;
  // Cell 2
  sv[2][0] = 1e-7;
  sv[2][1] = 0.01; sv[2][2] = 0.017;
  sv[2][3] = 0.005; sv[2][4] = 0.017;
  sv[2][5] = 0.001; sv[2][6] = 0.017;

  CheckConstraintFDJacobian(constraint, phase_prefixes, si, sv);
}
