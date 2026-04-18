// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <miam/constraints/dissolved_equilibrium_constraint.hpp>
#include <miam/constraints/dissolved_equilibrium_constraint_builder.hpp>
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
  auto A = micm::Species{ "A" };
  auto B = micm::Species{ "B" };
  auto C = micm::Species{ "C" };
  auto hp = micm::Species{ "H+" };
  auto aqueous_phase = micm::Phase{ "AQUEOUS", { { h2o }, { A }, { B }, { C }, { hp } } };
}  // namespace

// ── ConstraintAlgebraicVariableNames ──

TEST(DissolvedEquilibriumConstraint, AlgebraicVariableNamesSinglePrefix)
{
  auto keq = [](const micm::Conditions&) { return 10.0; };
  auto constraint = DissolvedEquilibriumConstraintBuilder()
      .SetPhase(aqueous_phase)
      .SetReactants({ A })
      .SetProducts({ B })
      .SetAlgebraicSpecies(B)
      .SetSolvent(h2o)
      .SetEquilibriumConstant(keq)
      .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("SMALL");

  auto names = constraint.ConstraintAlgebraicVariableNames(phase_prefixes);
  EXPECT_EQ(names.size(), 1);
  EXPECT_TRUE(names.count("SMALL.AQUEOUS.B"));
}

TEST(DissolvedEquilibriumConstraint, AlgebraicVariableNamesMultiplePrefixes)
{
  auto keq = [](const micm::Conditions&) { return 10.0; };
  auto constraint = DissolvedEquilibriumConstraintBuilder()
      .SetPhase(aqueous_phase)
      .SetReactants({ A })
      .SetProducts({ B })
      .SetAlgebraicSpecies(B)
      .SetSolvent(h2o)
      .SetEquilibriumConstant(keq)
      .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("SMALL");
  phase_prefixes["AQUEOUS"].insert("LARGE");

  auto names = constraint.ConstraintAlgebraicVariableNames(phase_prefixes);
  EXPECT_EQ(names.size(), 2);
  EXPECT_TRUE(names.count("SMALL.AQUEOUS.B"));
  EXPECT_TRUE(names.count("LARGE.AQUEOUS.B"));
}

TEST(DissolvedEquilibriumConstraint, AlgebraicVariableNamesNoMatchingPhase)
{
  auto keq = [](const micm::Conditions&) { return 10.0; };
  auto constraint = DissolvedEquilibriumConstraintBuilder()
      .SetPhase(aqueous_phase)
      .SetReactants({ A })
      .SetProducts({ B })
      .SetAlgebraicSpecies(B)
      .SetSolvent(h2o)
      .SetEquilibriumConstant(keq)
      .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["OTHER"].insert("PREFIX");

  auto names = constraint.ConstraintAlgebraicVariableNames(phase_prefixes);
  EXPECT_EQ(names.size(), 0);
}

// ── ConstraintSpeciesDependencies ──

TEST(DissolvedEquilibriumConstraint, SpeciesDependencies)
{
  auto keq = [](const micm::Conditions&) { return 10.0; };
  auto constraint = DissolvedEquilibriumConstraintBuilder()
      .SetPhase(aqueous_phase)
      .SetReactants({ A })
      .SetProducts({ B, hp })
      .SetAlgebraicSpecies(B)
      .SetSolvent(h2o)
      .SetEquilibriumConstant(keq)
      .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("SMALL");

  auto deps = constraint.ConstraintSpeciesDependencies(phase_prefixes);
  // A (reactant), B (product), H+ (product), H2O (solvent)
  EXPECT_EQ(deps.size(), 4);
  EXPECT_TRUE(deps.count("SMALL.AQUEOUS.A"));
  EXPECT_TRUE(deps.count("SMALL.AQUEOUS.B"));
  EXPECT_TRUE(deps.count("SMALL.AQUEOUS.H+"));
  EXPECT_TRUE(deps.count("SMALL.AQUEOUS.H2O"));
}

// ── NonZeroConstraintJacobianElements ──

TEST(DissolvedEquilibriumConstraint, NonZeroJacobianElements)
{
  auto keq = [](const micm::Conditions&) { return 10.0; };
  auto constraint = DissolvedEquilibriumConstraintBuilder()
      .SetPhase(aqueous_phase)
      .SetReactants({ A })
      .SetProducts({ B, hp })
      .SetAlgebraicSpecies(B)
      .SetSolvent(h2o)
      .SetEquilibriumConstant(keq)
      .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("SMALL");

  std::unordered_map<std::string, std::size_t> state_indices;
  state_indices["SMALL.AQUEOUS.A"] = 0;
  state_indices["SMALL.AQUEOUS.B"] = 1;
  state_indices["SMALL.AQUEOUS.H+"] = 2;
  state_indices["SMALL.AQUEOUS.H2O"] = 3;

  auto elements = constraint.NonZeroConstraintJacobianElements(phase_prefixes, state_indices);
  // Algebraic row = B (idx 1). Depends on: A (0), B (1), H+ (2), H2O (3)
  EXPECT_EQ(elements.size(), 4);
  EXPECT_TRUE(elements.count({ 1, 0 }));  // dG/dA
  EXPECT_TRUE(elements.count({ 1, 1 }));  // dG/dB
  EXPECT_TRUE(elements.count({ 1, 2 }));  // dG/dH+
  EXPECT_TRUE(elements.count({ 1, 3 }));  // dG/dH2O (solvent)
}

// ── ConstraintResidualFunction — simple A <-> B ──

TEST(DissolvedEquilibriumConstraint, ResidualSimpleAB)
{
  // G = K_eq * [A] - [B] = 0
  double K_eq = 10.0;
  auto keq = [K_eq](const micm::Conditions&) { return K_eq; };
  auto constraint = DissolvedEquilibriumConstraintBuilder()
      .SetPhase(aqueous_phase)
      .SetReactants({ A })
      .SetProducts({ B })
      .SetAlgebraicSpecies(B)
      .SetSolvent(h2o)
      .SetEquilibriumConstant(keq)
      .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("SMALL");

  std::unordered_map<std::string, std::size_t> state_indices;
  state_indices["SMALL.AQUEOUS.A"] = 0;
  state_indices["SMALL.AQUEOUS.B"] = 1;
  state_indices["SMALL.AQUEOUS.H2O"] = 2;

  // Initialize K_eq values (simulating UpdateConstraintParameters)
  std::vector<micm::Conditions> conditions(1);
  auto update_fn = constraint.UpdateConstraintParametersFunction<micm::Matrix<double>>();
  update_fn(conditions);

  using DMP = micm::Matrix<double>;
  auto residual_fn = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, param_indices, state_indices);

  DMP state_variables{ 1, 3, 0.0 };
  state_variables[0][0] = 1.0;    // [A]
  state_variables[0][1] = 5.0;    // [B]
  state_variables[0][2] = 300.0;  // [H2O]

  DMP residual{ 1, 3, 0.0 };
  residual_fn(state_variables, no_params, residual);

  // G = K_eq * [A] / [S]^(1-1) - [B] / [S]^(1-1) = K_eq * [A] - [B]
  // = 10.0 * 1.0 - 5.0 = 5.0
  EXPECT_NEAR(residual[0][1], 5.0, 1.0e-12);  // residual written to algebraic species (B) row
}

// ── ConstraintResidualFunction — A + B <-> C with solvent ──

TEST(DissolvedEquilibriumConstraint, ResidualMultiReactant)
{
  // G = K_eq * [A]*[B] / [S] - [C] / [S]^0 = K_eq * [A]*[B]/[S] - [C]
  double K_eq = 2.0;
  auto keq = [K_eq](const micm::Conditions&) { return K_eq; };
  auto constraint = DissolvedEquilibriumConstraintBuilder()
      .SetPhase(aqueous_phase)
      .SetReactants({ A, B })
      .SetProducts({ C })
      .SetAlgebraicSpecies(C)
      .SetSolvent(h2o)
      .SetEquilibriumConstant(keq)
      .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("SMALL");

  std::unordered_map<std::string, std::size_t> state_indices;
  state_indices["SMALL.AQUEOUS.A"] = 0;
  state_indices["SMALL.AQUEOUS.B"] = 1;
  state_indices["SMALL.AQUEOUS.C"] = 2;
  state_indices["SMALL.AQUEOUS.H2O"] = 3;

  std::vector<micm::Conditions> conditions(1);
  auto update_fn = constraint.UpdateConstraintParametersFunction<micm::Matrix<double>>();
  update_fn(conditions);

  using DMP = micm::Matrix<double>;
  auto residual_fn = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, param_indices, state_indices);

  DMP state_variables{ 1, 4, 0.0 };
  state_variables[0][0] = 3.0;    // [A]
  state_variables[0][1] = 4.0;    // [B]
  state_variables[0][2] = 20.0;   // [C]
  state_variables[0][3] = 300.0;  // [H2O]

  DMP residual{ 1, 4, 0.0 };
  residual_fn(state_variables, no_params, residual);

  // G = K_eq * [A]*[B]/[S]^(2-1) - [C]/[S]^(1-1)
  // = 2.0 * 3.0*4.0/300.0 - 20.0 = 0.08 - 20.0 = -19.92
  double expected = K_eq * 3.0 * 4.0 / 300.0 - 20.0;
  EXPECT_NEAR(residual[0][2], expected, 1.0e-12);
}

// ── ConstraintResidualFunction — multiple phase instances ──

TEST(DissolvedEquilibriumConstraint, ResidualMultipleInstances)
{
  double K_eq = 5.0;
  auto keq = [K_eq](const micm::Conditions&) { return K_eq; };
  auto constraint = DissolvedEquilibriumConstraintBuilder()
      .SetPhase(aqueous_phase)
      .SetReactants({ A })
      .SetProducts({ B })
      .SetAlgebraicSpecies(B)
      .SetSolvent(h2o)
      .SetEquilibriumConstant(keq)
      .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("LARGE");
  phase_prefixes["AQUEOUS"].insert("SMALL");

  std::unordered_map<std::string, std::size_t> state_indices;
  state_indices["LARGE.AQUEOUS.A"] = 0;
  state_indices["LARGE.AQUEOUS.B"] = 1;
  state_indices["LARGE.AQUEOUS.H2O"] = 2;
  state_indices["SMALL.AQUEOUS.A"] = 3;
  state_indices["SMALL.AQUEOUS.B"] = 4;
  state_indices["SMALL.AQUEOUS.H2O"] = 5;

  std::vector<micm::Conditions> conditions(1);
  auto update_fn = constraint.UpdateConstraintParametersFunction<micm::Matrix<double>>();
  update_fn(conditions);

  using DMP = micm::Matrix<double>;
  auto residual_fn = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, param_indices, state_indices);

  DMP state_variables{ 1, 6, 0.0 };
  state_variables[0][0] = 2.0;    // LARGE [A]
  state_variables[0][1] = 8.0;    // LARGE [B]
  state_variables[0][2] = 300.0;  // LARGE [H2O]
  state_variables[0][3] = 1.0;    // SMALL [A]
  state_variables[0][4] = 3.0;    // SMALL [B]
  state_variables[0][5] = 300.0;  // SMALL [H2O]

  DMP residual{ 1, 6, 0.0 };
  residual_fn(state_variables, no_params, residual);

  // LARGE: G = 5.0 * 2.0 - 8.0 = 2.0
  EXPECT_NEAR(residual[0][1], 2.0, 1.0e-12);
  // SMALL: G = 5.0 * 1.0 - 3.0 = 2.0
  EXPECT_NEAR(residual[0][4], 2.0, 1.0e-12);
}

// ── ConstraintResidualFunction — multiple grid cells ──

TEST(DissolvedEquilibriumConstraint, ResidualMultipleCells)
{
  double K_eq = 4.0;
  auto keq = [K_eq](const micm::Conditions&) { return K_eq; };
  auto constraint = DissolvedEquilibriumConstraintBuilder()
      .SetPhase(aqueous_phase)
      .SetReactants({ A })
      .SetProducts({ B })
      .SetAlgebraicSpecies(B)
      .SetSolvent(h2o)
      .SetEquilibriumConstant(keq)
      .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  std::unordered_map<std::string, std::size_t> state_indices;
  state_indices["DROP.AQUEOUS.A"] = 0;
  state_indices["DROP.AQUEOUS.B"] = 1;
  state_indices["DROP.AQUEOUS.H2O"] = 2;

  std::vector<micm::Conditions> conditions(2);
  auto update_fn = constraint.UpdateConstraintParametersFunction<micm::Matrix<double>>();
  update_fn(conditions);

  using DMP = micm::Matrix<double>;
  auto residual_fn = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, param_indices, state_indices);

  DMP state_variables{ 2, 3, 0.0 };
  state_variables[0][0] = 1.0;    // cell 0 [A]
  state_variables[0][1] = 2.0;    // cell 0 [B]
  state_variables[0][2] = 300.0;  // cell 0 [H2O]
  state_variables[1][0] = 3.0;    // cell 1 [A]
  state_variables[1][1] = 10.0;   // cell 1 [B]
  state_variables[1][2] = 300.0;  // cell 1 [H2O]

  DMP residual{ 2, 3, 0.0 };
  residual_fn(state_variables, no_params, residual);

  // cell 0: G = 4*1 - 2 = 2
  EXPECT_NEAR(residual[0][1], 2.0, 1.0e-12);
  // cell 1: G = 4*3 - 10 = 2
  EXPECT_NEAR(residual[1][1], 2.0, 1.0e-12);
}

// ── ConstraintJacobianFunction — simple A <-> B ──

TEST(DissolvedEquilibriumConstraint, JacobianSimpleAB)
{
  // G = K_eq * [A] - [B]
  // dG/dA = K_eq, dG/dB = -1
  double K_eq = 10.0;
  auto keq = [K_eq](const micm::Conditions&) { return K_eq; };
  auto constraint = DissolvedEquilibriumConstraintBuilder()
      .SetPhase(aqueous_phase)
      .SetReactants({ A })
      .SetProducts({ B })
      .SetAlgebraicSpecies(B)
      .SetSolvent(h2o)
      .SetEquilibriumConstant(keq)
      .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("SMALL");

  std::unordered_map<std::string, std::size_t> state_indices;
  state_indices["SMALL.AQUEOUS.A"] = 0;
  state_indices["SMALL.AQUEOUS.B"] = 1;
  state_indices["SMALL.AQUEOUS.H2O"] = 2;

  std::vector<micm::Conditions> conditions(1);
  auto update_fn = constraint.UpdateConstraintParametersFunction<micm::Matrix<double>>();
  update_fn(conditions);

  // Build sparse Jacobian with the right sparsity
  using DMP = micm::Matrix<double>;
  using SMP = micm::SparseMatrix<double, micm::SparseMatrixStandardOrderingCompressedSparseRow>;
  auto nz_elements = constraint.NonZeroConstraintJacobianElements(phase_prefixes, state_indices);
  auto builder = SMP::Create(3).SetNumberOfBlocks(1);
  for (const auto& [row, col] : nz_elements)
    builder.WithElement(row, col);
  SMP jacobian(builder);

  auto jac_fn = constraint.ConstraintJacobianFunction<DMP, SMP>(phase_prefixes, state_indices, jacobian);

  DMP state_variables{ 1, 3, 0.0 };
  state_variables[0][0] = 1.0;    // [A]
  state_variables[0][1] = 5.0;    // [B]
  state_variables[0][2] = 300.0;  // [H2O]

  // Zero the Jacobian
  for (auto& v : jacobian.AsVector())
    v = 0.0;

  jac_fn(state_variables, jacobian);

  // MICM convention: jac -= dG/dy
  // jac[B,A] -= K_eq → jac[B,A] = -K_eq = -10.0
  EXPECT_NEAR(jacobian.AsVector()[jacobian.VectorIndex(0, 1, 0)], -K_eq, 1.0e-12);
  // jac[B,B] -= (-1) → jac[B,B] = +1.0
  EXPECT_NEAR(jacobian.AsVector()[jacobian.VectorIndex(0, 1, 1)], 1.0, 1.0e-12);
}

// ── ConstraintJacobianFunction — A + B <-> C with solvent ──

TEST(DissolvedEquilibriumConstraint, JacobianMultiReactant)
{
  // G = K_eq * [A]*[B]/[S] - [C]
  double K_eq = 2.0;
  auto keq = [K_eq](const micm::Conditions&) { return K_eq; };
  auto constraint = DissolvedEquilibriumConstraintBuilder()
      .SetPhase(aqueous_phase)
      .SetReactants({ A, B })
      .SetProducts({ C })
      .SetAlgebraicSpecies(C)
      .SetSolvent(h2o)
      .SetEquilibriumConstant(keq)
      .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("SMALL");

  std::unordered_map<std::string, std::size_t> state_indices;
  state_indices["SMALL.AQUEOUS.A"] = 0;
  state_indices["SMALL.AQUEOUS.B"] = 1;
  state_indices["SMALL.AQUEOUS.C"] = 2;
  state_indices["SMALL.AQUEOUS.H2O"] = 3;

  std::vector<micm::Conditions> conditions(1);
  auto update_fn = constraint.UpdateConstraintParametersFunction<micm::Matrix<double>>();
  update_fn(conditions);

  using DMP = micm::Matrix<double>;
  using SMP = micm::SparseMatrix<double, micm::SparseMatrixStandardOrderingCompressedSparseRow>;
  auto nz_elements = constraint.NonZeroConstraintJacobianElements(phase_prefixes, state_indices);
  auto builder = SMP::Create(4).SetNumberOfBlocks(1);
  for (const auto& [row, col] : nz_elements)
    builder.WithElement(row, col);
  SMP jacobian(builder);

  auto jac_fn = constraint.ConstraintJacobianFunction<DMP, SMP>(phase_prefixes, state_indices, jacobian);

  DMP state_variables{ 1, 4, 0.0 };
  state_variables[0][0] = 3.0;    // [A]
  state_variables[0][1] = 4.0;    // [B]
  state_variables[0][2] = 20.0;   // [C]
  state_variables[0][3] = 300.0;  // [S]

  for (auto& v : jacobian.AsVector())
    v = 0.0;
  jac_fn(state_variables, jacobian);

  // dG/dA = K_eq * [B] / [S] = 2.0 * 4.0 / 300.0
  double dG_dA = K_eq * 4.0 / 300.0;
  EXPECT_NEAR(jacobian.AsVector()[jacobian.VectorIndex(0, 2, 0)], -dG_dA, 1.0e-12);

  // dG/dB = K_eq * [A] / [S] = 2.0 * 3.0 / 300.0
  double dG_dB = K_eq * 3.0 / 300.0;
  EXPECT_NEAR(jacobian.AsVector()[jacobian.VectorIndex(0, 2, 1)], -dG_dB, 1.0e-12);

  // dG/dC = -1
  EXPECT_NEAR(jacobian.AsVector()[jacobian.VectorIndex(0, 2, 2)], 1.0, 1.0e-12);

  // dG/dS = K_eq * [A] * [B] * (1 - 2) / [S]^2 - [C] * (1 - 1) / [S]^1
  //       = -K_eq * [A] * [B] / [S]^2
  double dG_dS = -K_eq * 3.0 * 4.0 / (300.0 * 300.0);
  EXPECT_NEAR(jacobian.AsVector()[jacobian.VectorIndex(0, 2, 3)], -dG_dS, 1.0e-12);
}

// ── UpdateConstraintParameters — temperature-dependent K_eq ──

TEST(DissolvedEquilibriumConstraint, UpdateConstraintParametersTemperatureDep)
{
  // K_eq(T) = 1000.0 / T
  auto keq = [](const micm::Conditions& c) { return 1000.0 / c.temperature_; };
  auto constraint = DissolvedEquilibriumConstraintBuilder()
      .SetPhase(aqueous_phase)
      .SetReactants({ A })
      .SetProducts({ B })
      .SetAlgebraicSpecies(B)
      .SetSolvent(h2o)
      .SetEquilibriumConstant(keq)
      .Build();

  auto update_fn = constraint.UpdateConstraintParametersFunction<micm::Matrix<double>>();

  std::vector<micm::Conditions> conditions(3);
  conditions[0].temperature_ = 250.0;
  conditions[1].temperature_ = 300.0;
  conditions[2].temperature_ = 350.0;
  update_fn(conditions);

  EXPECT_EQ(constraint.k_eq_values_->size(), 3);
  EXPECT_NEAR((*constraint.k_eq_values_)[0], 4.0, 1.0e-12);
  EXPECT_NEAR((*constraint.k_eq_values_)[1], 1000.0 / 300.0, 1.0e-12);
  EXPECT_NEAR((*constraint.k_eq_values_)[2], 1000.0 / 350.0, 1.0e-12);
}

// ── Builder ──

TEST(DissolvedEquilibriumConstraint, BuilderValidation)
{
  // Requires that algebraic species is one of the products
  struct FakeConstant
  {
    double Calculate(const micm::Conditions&) const { return 1.0; }
  };

  EXPECT_THROW(
      DissolvedEquilibriumConstraintBuilder()
          .SetPhase(aqueous_phase)
          .SetReactants({ A })
          .SetProducts({ B })
          .SetAlgebraicSpecies(A)  // A is a reactant, not a product
          .SetSolvent(h2o)
          .SetEquilibriumConstant(FakeConstant{})
          .Build(),
      std::invalid_argument);

  // Missing phase
  EXPECT_THROW(
      DissolvedEquilibriumConstraintBuilder().SetReactants({ A }).SetProducts({ B }).SetAlgebraicSpecies(B).SetSolvent(h2o).SetEquilibriumConstant(FakeConstant{}).Build(),
      std::runtime_error);

  // Valid build
  auto constraint = DissolvedEquilibriumConstraintBuilder()
                        .SetPhase(aqueous_phase)
                        .SetReactants({ A })
                        .SetProducts({ B })
                        .SetAlgebraicSpecies(B)
                        .SetSolvent(h2o)
                        .SetEquilibriumConstant(FakeConstant{})
                        .Build();
  EXPECT_EQ(constraint.algebraic_species_.name_, "B");
}

// ═══════════════════════════════════════════════════════════════════
// Expanded Scenarios
// ═══════════════════════════════════════════════════════════════════

namespace
{
  using DMP = micm::Matrix<double>;
  using SMP = micm::SparseMatrix<double, micm::SparseMatrixStandardOrderingCompressedSparseRow>;

  SMP BuildConstraintJacobian(
      const DissolvedEquilibriumConstraint& constraint,
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

  /// @brief Finite-difference check for the constraint Jacobian
  void CheckConstraintFDJacobian(
      DissolvedEquilibriumConstraint& constraint,
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

  /// @brief Helper to initialize K_eq values for a constraint
  void InitKeq(DissolvedEquilibriumConstraint& constraint, std::size_t num_cells)
  {
    std::vector<micm::Conditions> conditions(num_cells);
    auto update_fn = constraint.UpdateConstraintParametersFunction<DMP>();
    update_fn(conditions);
  }
}  // namespace

// ── FD Jacobian: simple A <-> B ──

TEST(DissolvedEquilibriumConstraint, JacobianFDSimpleAB)
{
  double K_eq = 10.0;
  auto keq = [K_eq](const micm::Conditions&) { return K_eq; };
  auto constraint = DissolvedEquilibriumConstraintBuilder()
      .SetPhase(aqueous_phase)
      .SetReactants({ A })
      .SetProducts({ B })
      .SetAlgebraicSpecies(B)
      .SetSolvent(h2o)
      .SetEquilibriumConstant(keq)
      .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("SMALL");

  std::unordered_map<std::string, std::size_t> si;
  si["SMALL.AQUEOUS.A"] = 0;
  si["SMALL.AQUEOUS.B"] = 1;
  si["SMALL.AQUEOUS.H2O"] = 2;

  InitKeq(constraint, 1);

  DMP sv{ 1, 3, 0.0 };
  sv[0][0] = 1.0;
  sv[0][1] = 5.0;
  sv[0][2] = 300.0;

  CheckConstraintFDJacobian(constraint, phase_prefixes, si, sv);
}

// ── FD Jacobian: multi-reactant A + B <-> C ──

TEST(DissolvedEquilibriumConstraint, JacobianFDMultiReactant)
{
  double K_eq = 2.0;
  auto keq = [K_eq](const micm::Conditions&) { return K_eq; };
  auto constraint = DissolvedEquilibriumConstraintBuilder()
      .SetPhase(aqueous_phase)
      .SetReactants({ A, B })
      .SetProducts({ C })
      .SetAlgebraicSpecies(C)
      .SetSolvent(h2o)
      .SetEquilibriumConstant(keq)
      .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("SMALL");

  std::unordered_map<std::string, std::size_t> si;
  si["SMALL.AQUEOUS.A"] = 0;
  si["SMALL.AQUEOUS.B"] = 1;
  si["SMALL.AQUEOUS.C"] = 2;
  si["SMALL.AQUEOUS.H2O"] = 3;

  InitKeq(constraint, 1);

  DMP sv{ 1, 4, 0.0 };
  sv[0][0] = 3.0;
  sv[0][1] = 4.0;
  sv[0][2] = 20.0;
  sv[0][3] = 300.0;

  CheckConstraintFDJacobian(constraint, phase_prefixes, si, sv);
}

// ── FD Jacobian: multi-product A <-> B + H+ ──

TEST(DissolvedEquilibriumConstraint, JacobianFDMultiProduct)
{
  double K_eq = 0.5;
  auto keq = [K_eq](const micm::Conditions&) { return K_eq; };
  auto constraint = DissolvedEquilibriumConstraintBuilder()
      .SetPhase(aqueous_phase)
      .SetReactants({ A })
      .SetProducts({ B, hp })
      .SetAlgebraicSpecies(B)
      .SetSolvent(h2o)
      .SetEquilibriumConstant(keq)
      .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  std::unordered_map<std::string, std::size_t> si;
  si["DROP.AQUEOUS.A"] = 0;
  si["DROP.AQUEOUS.B"] = 1;
  si["DROP.AQUEOUS.H+"] = 2;
  si["DROP.AQUEOUS.H2O"] = 3;

  InitKeq(constraint, 1);

  DMP sv{ 1, 4, 0.0 };
  sv[0][0] = 0.1;
  sv[0][1] = 0.01;
  sv[0][2] = 1.0e-4;
  sv[0][3] = 0.017;

  CheckConstraintFDJacobian(constraint, phase_prefixes, si, sv);
}

// ── FD Jacobian: multi-reactant + multi-product  A + B <-> C + H+ ──

TEST(DissolvedEquilibriumConstraint, JacobianFDMultiReactantMultiProduct)
{
  double K_eq = 3.0;
  auto keq = [K_eq](const micm::Conditions&) { return K_eq; };
  auto constraint = DissolvedEquilibriumConstraintBuilder()
      .SetPhase(aqueous_phase)
      .SetReactants({ A, B })
      .SetProducts({ C, hp })
      .SetAlgebraicSpecies(C)
      .SetSolvent(h2o)
      .SetEquilibriumConstant(keq)
      .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("M1");

  std::unordered_map<std::string, std::size_t> si;
  si["M1.AQUEOUS.A"] = 0;
  si["M1.AQUEOUS.B"] = 1;
  si["M1.AQUEOUS.C"] = 2;
  si["M1.AQUEOUS.H+"] = 3;
  si["M1.AQUEOUS.H2O"] = 4;

  InitKeq(constraint, 1);

  DMP sv{ 1, 5, 0.0 };
  sv[0][0] = 0.5;
  sv[0][1] = 0.3;
  sv[0][2] = 0.1;
  sv[0][3] = 1.0e-3;
  sv[0][4] = 0.017;

  CheckConstraintFDJacobian(constraint, phase_prefixes, si, sv);
}

// ── Multiple cells: residual + FD ──

TEST(DissolvedEquilibriumConstraint, ResidualAndJacobianFDMultipleCells)
{
  double K_eq = 5.0;
  auto keq = [K_eq](const micm::Conditions&) { return K_eq; };
  auto constraint = DissolvedEquilibriumConstraintBuilder()
      .SetPhase(aqueous_phase)
      .SetReactants({ A })
      .SetProducts({ B })
      .SetAlgebraicSpecies(B)
      .SetSolvent(h2o)
      .SetEquilibriumConstant(keq)
      .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  std::unordered_map<std::string, std::size_t> si;
  si["DROP.AQUEOUS.A"] = 0;
  si["DROP.AQUEOUS.B"] = 1;
  si["DROP.AQUEOUS.H2O"] = 2;

  std::size_t nc = 4;
  InitKeq(constraint, nc);

  DMP sv{ nc, 3, 0.0 };
  sv[0][0] = 1.0;   sv[0][1] = 5.0;   sv[0][2] = 300.0;
  sv[1][0] = 0.5;   sv[1][1] = 2.5;   sv[1][2] = 0.017;
  sv[2][0] = 10.0;  sv[2][1] = 40.0;  sv[2][2] = 100.0;
  sv[3][0] = 0.01;  sv[3][1] = 0.1;   sv[3][2] = 0.017;

  // Check residuals
  auto rf = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, param_indices, si);
  DMP residual{ nc, 3, 0.0 };
  rf(sv, no_params, residual);

  // G = K_eq * [A] - [B]
  EXPECT_NEAR(residual[0][1], 5.0 * 1.0 - 5.0, 1e-12);
  EXPECT_NEAR(residual[1][1], 5.0 * 0.5 - 2.5, 1e-12);
  EXPECT_NEAR(residual[2][1], 5.0 * 10.0 - 40.0, 1e-12);
  EXPECT_NEAR(residual[3][1], 5.0 * 0.01 - 0.1, 1e-12);

  // FD Jacobian check
  CheckConstraintFDJacobian(constraint, phase_prefixes, si, sv);
}

// ── Multiple instances + multiple cells + FD ──

TEST(DissolvedEquilibriumConstraint, MultiInstanceMultiCellFD)
{
  double K_eq = 3.0;
  auto keq = [K_eq](const micm::Conditions&) { return K_eq; };
  auto constraint = DissolvedEquilibriumConstraintBuilder()
      .SetPhase(aqueous_phase)
      .SetReactants({ A, B })
      .SetProducts({ C })
      .SetAlgebraicSpecies(C)
      .SetSolvent(h2o)
      .SetEquilibriumConstant(keq)
      .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("LARGE");
  phase_prefixes["AQUEOUS"].insert("SMALL");

  std::unordered_map<std::string, std::size_t> si;
  si["LARGE.AQUEOUS.A"] = 0;
  si["LARGE.AQUEOUS.B"] = 1;
  si["LARGE.AQUEOUS.C"] = 2;
  si["LARGE.AQUEOUS.H2O"] = 3;
  si["SMALL.AQUEOUS.A"] = 4;
  si["SMALL.AQUEOUS.B"] = 5;
  si["SMALL.AQUEOUS.C"] = 6;
  si["SMALL.AQUEOUS.H2O"] = 7;

  std::size_t nc = 3;
  InitKeq(constraint, nc);

  DMP sv{ nc, 8, 0.0 };
  // Cell 0
  sv[0][0] = 0.5;  sv[0][1] = 0.3;  sv[0][2] = 0.1;  sv[0][3] = 0.017;
  sv[0][4] = 0.2;  sv[0][5] = 0.1;  sv[0][6] = 0.05; sv[0][7] = 0.017;
  // Cell 1
  sv[1][0] = 1.0;  sv[1][1] = 2.0;  sv[1][2] = 3.0;  sv[1][3] = 100.0;
  sv[1][4] = 0.5;  sv[1][5] = 1.0;  sv[1][6] = 1.5;  sv[1][7] = 80.0;
  // Cell 2
  sv[2][0] = 5.0;  sv[2][1] = 3.0;  sv[2][2] = 10.0; sv[2][3] = 0.017;
  sv[2][4] = 2.0;  sv[2][5] = 1.5;  sv[2][6] = 4.0;  sv[2][7] = 0.017;

  // Residuals: G = K_eq * [A]*[B]/[S] - [C]
  auto rf = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, param_indices, si);
  DMP residual{ nc, 8, 0.0 };
  rf(sv, no_params, residual);

  // LARGE instance, cell 0: G = 3.0 * 0.5*0.3/0.017 - 0.1
  double expected_L0 = 3.0 * 0.5 * 0.3 / 0.017 - 0.1;
  EXPECT_NEAR(residual[0][2], expected_L0, 1e-12);

  // SMALL instance, cell 0: G = 3.0 * 0.2*0.1/0.017 - 0.05
  double expected_S0 = 3.0 * 0.2 * 0.1 / 0.017 - 0.05;
  EXPECT_NEAR(residual[0][6], expected_S0, 1e-12);

  // Non-algebraic rows untouched
  EXPECT_NEAR(residual[0][0], 0.0, 1e-30);
  EXPECT_NEAR(residual[0][1], 0.0, 1e-30);
  EXPECT_NEAR(residual[0][3], 0.0, 1e-30);

  // FD check
  CheckConstraintFDJacobian(constraint, phase_prefixes, si, sv);
}

// ── Three instances + FD ──

TEST(DissolvedEquilibriumConstraint, ThreeInstancesFD)
{
  double K_eq = 7.0;
  auto keq = [K_eq](const micm::Conditions&) { return K_eq; };
  auto constraint = DissolvedEquilibriumConstraintBuilder()
      .SetPhase(aqueous_phase)
      .SetReactants({ A })
      .SetProducts({ B, hp })
      .SetAlgebraicSpecies(B)
      .SetSolvent(h2o)
      .SetEquilibriumConstant(keq)
      .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("A");
  phase_prefixes["AQUEOUS"].insert("B");
  phase_prefixes["AQUEOUS"].insert("C");

  std::unordered_map<std::string, std::size_t> si;
  si["A.AQUEOUS.A"] = 0;   si["A.AQUEOUS.B"] = 1;   si["A.AQUEOUS.H+"] = 2;   si["A.AQUEOUS.H2O"] = 3;
  si["B.AQUEOUS.A"] = 4;   si["B.AQUEOUS.B"] = 5;   si["B.AQUEOUS.H+"] = 6;   si["B.AQUEOUS.H2O"] = 7;
  si["C.AQUEOUS.A"] = 8;   si["C.AQUEOUS.B"] = 9;   si["C.AQUEOUS.H+"] = 10;  si["C.AQUEOUS.H2O"] = 11;

  std::size_t nc = 2;
  InitKeq(constraint, nc);

  DMP sv{ nc, 12, 0.0 };
  // Cell 0
  sv[0][0] = 0.5;  sv[0][1] = 0.1;  sv[0][2] = 1e-4;  sv[0][3] = 0.017;
  sv[0][4] = 1.0;  sv[0][5] = 0.2;  sv[0][6] = 2e-4;  sv[0][7] = 0.017;
  sv[0][8] = 2.0;  sv[0][9] = 0.5;  sv[0][10] = 5e-4; sv[0][11] = 0.017;
  // Cell 1
  sv[1][0] = 3.0;  sv[1][1] = 1.0;  sv[1][2] = 1e-3;  sv[1][3] = 40.0;
  sv[1][4] = 4.0;  sv[1][5] = 2.0;  sv[1][6] = 2e-3;  sv[1][7] = 50.0;
  sv[1][8] = 5.0;  sv[1][9] = 3.0;  sv[1][10] = 3e-3; sv[1][11] = 60.0;

  CheckConstraintFDJacobian(constraint, phase_prefixes, si, sv);
}

// ── Jacobian accumulates (subtracts) ──

TEST(DissolvedEquilibriumConstraint, JacobianAccumulates)
{
  double K_eq = 10.0;
  auto keq = [K_eq](const micm::Conditions&) { return K_eq; };
  auto constraint = DissolvedEquilibriumConstraintBuilder()
      .SetPhase(aqueous_phase)
      .SetReactants({ A })
      .SetProducts({ B })
      .SetAlgebraicSpecies(B)
      .SetSolvent(h2o)
      .SetEquilibriumConstant(keq)
      .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  std::unordered_map<std::string, std::size_t> si;
  si["DROP.AQUEOUS.A"] = 0;
  si["DROP.AQUEOUS.B"] = 1;
  si["DROP.AQUEOUS.H2O"] = 2;

  InitKeq(constraint, 1);

  auto jacobian = BuildConstraintJacobian(constraint, phase_prefixes, si, 1);
  auto jac_fn = constraint.ConstraintJacobianFunction<DMP, SMP>(phase_prefixes, si, jacobian);

  DMP sv{ 1, 3, 0.0 };
  sv[0][0] = 1.0;  sv[0][1] = 5.0;  sv[0][2] = 300.0;

  // First call
  jac_fn(sv, jacobian);
  double j_BA_once = jacobian.AsVector()[jacobian.VectorIndex(0, 1, 0)];
  double j_BB_once = jacobian.AsVector()[jacobian.VectorIndex(0, 1, 1)];

  // Second call should accumulate
  jac_fn(sv, jacobian);
  EXPECT_NEAR(jacobian.AsVector()[jacobian.VectorIndex(0, 1, 0)], 2.0 * j_BA_once, 1e-12);
  EXPECT_NEAR(jacobian.AsVector()[jacobian.VectorIndex(0, 1, 1)], 2.0 * j_BB_once, 1e-12);
}

// ── Residual sets (does not accumulate) ──

TEST(DissolvedEquilibriumConstraint, ResidualSetsNotAccumulates)
{
  double K_eq = 10.0;
  auto keq = [K_eq](const micm::Conditions&) { return K_eq; };
  auto constraint = DissolvedEquilibriumConstraintBuilder()
      .SetPhase(aqueous_phase)
      .SetReactants({ A })
      .SetProducts({ B })
      .SetAlgebraicSpecies(B)
      .SetSolvent(h2o)
      .SetEquilibriumConstant(keq)
      .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  std::unordered_map<std::string, std::size_t> si;
  si["DROP.AQUEOUS.A"] = 0;
  si["DROP.AQUEOUS.B"] = 1;
  si["DROP.AQUEOUS.H2O"] = 2;

  InitKeq(constraint, 1);

  auto rf = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, param_indices, si);

  DMP sv{ 1, 3, 0.0 };
  sv[0][0] = 1.0;  sv[0][1] = 5.0;  sv[0][2] = 300.0;

  DMP residual{ 1, 3, 999.0 };
  rf(sv, no_params, residual);
  double val1 = residual[0][1];

  rf(sv, no_params, residual);
  EXPECT_NEAR(residual[0][1], val1, 1e-15);
}

// ── Temperature-dependent K_eq with multi-cell FD ──

TEST(DissolvedEquilibriumConstraint, TemperatureDependentKeqFD)
{
  // K_eq(T) = 1000.0 / T
  auto keq = [](const micm::Conditions& c) { return 1000.0 / c.temperature_; };
  auto constraint = DissolvedEquilibriumConstraintBuilder()
      .SetPhase(aqueous_phase)
      .SetReactants({ A })
      .SetProducts({ B })
      .SetAlgebraicSpecies(B)
      .SetSolvent(h2o)
      .SetEquilibriumConstant(keq)
      .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  std::unordered_map<std::string, std::size_t> si;
  si["DROP.AQUEOUS.A"] = 0;
  si["DROP.AQUEOUS.B"] = 1;
  si["DROP.AQUEOUS.H2O"] = 2;

  std::size_t nc = 3;
  std::vector<micm::Conditions> conditions(nc);
  conditions[0].temperature_ = 250.0;
  conditions[1].temperature_ = 300.0;
  conditions[2].temperature_ = 350.0;
  auto update_fn = constraint.UpdateConstraintParametersFunction<DMP>();
  update_fn(conditions);

  DMP sv{ nc, 3, 0.0 };
  sv[0][0] = 1.0;  sv[0][1] = 2.0;  sv[0][2] = 0.017;
  sv[1][0] = 0.5;  sv[1][1] = 1.0;  sv[1][2] = 0.017;
  sv[2][0] = 2.0;  sv[2][1] = 5.0;  sv[2][2] = 0.017;

  // Check residuals with different K_eq per cell
  auto rf = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, param_indices, si);
  DMP residual{ nc, 3, 0.0 };
  rf(sv, no_params, residual);

  EXPECT_NEAR(residual[0][1], (1000.0 / 250.0) * 1.0 - 2.0, 1e-12);
  EXPECT_NEAR(residual[1][1], (1000.0 / 300.0) * 0.5 - 1.0, 1e-12);
  EXPECT_NEAR(residual[2][1], (1000.0 / 350.0) * 2.0 - 5.0, 1e-12);

  // FD check
  CheckConstraintFDJacobian(constraint, phase_prefixes, si, sv);
}

// ── Large K_eq range ──

TEST(DissolvedEquilibriumConstraint, LargeKeqRange)
{
  // Very large K_eq
  double K_eq = 1.0e8;
  auto keq = [K_eq](const micm::Conditions&) { return K_eq; };
  auto constraint = DissolvedEquilibriumConstraintBuilder()
      .SetPhase(aqueous_phase)
      .SetReactants({ A })
      .SetProducts({ B })
      .SetAlgebraicSpecies(B)
      .SetSolvent(h2o)
      .SetEquilibriumConstant(keq)
      .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  std::unordered_map<std::string, std::size_t> si;
  si["DROP.AQUEOUS.A"] = 0;
  si["DROP.AQUEOUS.B"] = 1;
  si["DROP.AQUEOUS.H2O"] = 2;

  InitKeq(constraint, 1);

  DMP sv{ 1, 3, 0.0 };
  sv[0][0] = 1.0e-6;   // small [A]
  sv[0][1] = 100.0;    // [B] = K_eq * A
  sv[0][2] = 0.017;

  auto rf = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, param_indices, si);
  DMP residual{ 1, 3, 0.0 };
  rf(sv, no_params, residual);

  // G = 1e8 * 1e-6 - 100 = 100 - 100 = 0
  EXPECT_NEAR(residual[0][1], 0.0, 1e-6);

  CheckConstraintFDJacobian(constraint, phase_prefixes, si, sv);
}

// ── CopyWithNewUuid preserves behavior ──

TEST(DissolvedEquilibriumConstraint, CopiedConstraintProducesSameResults)
{
  double K_eq = 10.0;
  auto keq = [K_eq](const micm::Conditions&) { return K_eq; };
  auto original = DissolvedEquilibriumConstraintBuilder()
      .SetPhase(aqueous_phase)
      .SetReactants({ A, B })
      .SetProducts({ C })
      .SetAlgebraicSpecies(C)
      .SetSolvent(h2o)
      .SetEquilibriumConstant(keq)
      .Build();

  auto copy = original.CopyWithNewUuid();
  EXPECT_NE(copy.uuid_, original.uuid_);

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("MODE1");

  std::unordered_map<std::string, std::size_t> si;
  si["MODE1.AQUEOUS.A"] = 0;
  si["MODE1.AQUEOUS.B"] = 1;
  si["MODE1.AQUEOUS.C"] = 2;
  si["MODE1.AQUEOUS.H2O"] = 3;

  // Initialize both k_eq_values
  InitKeq(original, 1);
  InitKeq(copy, 1);

  DMP sv{ 1, 4, 0.0 };
  sv[0][0] = 3.0;  sv[0][1] = 4.0;  sv[0][2] = 20.0;  sv[0][3] = 300.0;

  auto rf_orig = original.ConstraintResidualFunction<DMP>(phase_prefixes, param_indices, si);
  auto rf_copy = copy.ConstraintResidualFunction<DMP>(phase_prefixes, param_indices, si);

  DMP res_orig{ 1, 4, 0.0 };
  DMP res_copy{ 1, 4, 0.0 };
  rf_orig(sv, no_params, res_orig);
  rf_copy(sv, no_params, res_copy);

  EXPECT_NEAR(res_orig[0][2], res_copy[0][2], 1e-15);
}

// ── Solvent dependence: verify dG/dS for multi-reactant case ──

TEST(DissolvedEquilibriumConstraint, SolventJacobianMultiReactant)
{
  // G = K_eq * [A]*[B] / [S] - [C]
  // dG/dS = K_eq * [A]*[B] * (-1) / [S]^2 = -K_eq*[A]*[B]/[S]^2
  double K_eq = 2.0;
  auto keq = [K_eq](const micm::Conditions&) { return K_eq; };
  auto constraint = DissolvedEquilibriumConstraintBuilder()
      .SetPhase(aqueous_phase)
      .SetReactants({ A, B })
      .SetProducts({ C })
      .SetAlgebraicSpecies(C)
      .SetSolvent(h2o)
      .SetEquilibriumConstant(keq)
      .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  std::unordered_map<std::string, std::size_t> si;
  si["DROP.AQUEOUS.A"] = 0;
  si["DROP.AQUEOUS.B"] = 1;
  si["DROP.AQUEOUS.C"] = 2;
  si["DROP.AQUEOUS.H2O"] = 3;

  InitKeq(constraint, 1);

  auto jacobian = BuildConstraintJacobian(constraint, phase_prefixes, si, 1);
  auto jac_fn = constraint.ConstraintJacobianFunction<DMP, SMP>(phase_prefixes, si, jacobian);

  DMP sv{ 1, 4, 0.0 };
  sv[0][0] = 3.0;
  sv[0][1] = 4.0;
  sv[0][2] = 20.0;
  sv[0][3] = 300.0;

  jac_fn(sv, jacobian);

  // dG/dS = K_eq * [A]*[B] * (1-2)/[S]^2 - [C] * (1-1)/[S]^1
  //       = -K_eq * [A]*[B] / [S]^2 - 0
  double dG_dS = -K_eq * 3.0 * 4.0 / (300.0 * 300.0);
  // jac -= dG/dS
  EXPECT_NEAR(jacobian.AsVector()[jacobian.VectorIndex(0, 2, 3)], -dG_dS, 1e-12);

  // Also FD check
  CheckConstraintFDJacobian(constraint, phase_prefixes, si, sv);
}

// ── Kitchen-sink: multi-reactant, multi-product, multi-instance, multi-cell, FD ──

TEST(DissolvedEquilibriumConstraint, KitchenSinkFD)
{
  // A + B <-> C + H+,  K_eq(T) = 500.0 / T
  auto keq = [](const micm::Conditions& c) { return 500.0 / c.temperature_; };
  auto constraint = DissolvedEquilibriumConstraintBuilder()
      .SetPhase(aqueous_phase)
      .SetReactants({ A, B })
      .SetProducts({ C, hp })
      .SetAlgebraicSpecies(C)
      .SetSolvent(h2o)
      .SetEquilibriumConstant(keq)
      .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("M1");
  phase_prefixes["AQUEOUS"].insert("M2");

  std::unordered_map<std::string, std::size_t> si;
  si["M1.AQUEOUS.A"] = 0;   si["M1.AQUEOUS.B"] = 1;   si["M1.AQUEOUS.C"] = 2;
  si["M1.AQUEOUS.H+"] = 3;  si["M1.AQUEOUS.H2O"] = 4;
  si["M2.AQUEOUS.A"] = 5;   si["M2.AQUEOUS.B"] = 6;   si["M2.AQUEOUS.C"] = 7;
  si["M2.AQUEOUS.H+"] = 8;  si["M2.AQUEOUS.H2O"] = 9;

  std::size_t nc = 3;
  std::vector<micm::Conditions> conditions(nc);
  conditions[0].temperature_ = 280.0;
  conditions[1].temperature_ = 300.0;
  conditions[2].temperature_ = 320.0;
  auto update_fn = constraint.UpdateConstraintParametersFunction<DMP>();
  update_fn(conditions);

  DMP sv{ nc, 10, 0.0 };
  // Cell 0
  sv[0][0] = 0.5;  sv[0][1] = 0.3;  sv[0][2] = 0.1;  sv[0][3] = 1e-4;  sv[0][4] = 0.017;
  sv[0][5] = 0.2;  sv[0][6] = 0.4;  sv[0][7] = 0.05; sv[0][8] = 2e-4;  sv[0][9] = 40.0;
  // Cell 1
  sv[1][0] = 1.0;  sv[1][1] = 2.0;  sv[1][2] = 0.5;  sv[1][3] = 5e-3;  sv[1][4] = 0.017;
  sv[1][5] = 3.0;  sv[1][6] = 1.0;  sv[1][7] = 1.0;  sv[1][8] = 1e-3;  sv[1][9] = 0.017;
  // Cell 2
  sv[2][0] = 0.1;  sv[2][1] = 0.1;  sv[2][2] = 0.01; sv[2][3] = 1e-5;  sv[2][4] = 0.017;
  sv[2][5] = 5.0;  sv[2][6] = 2.0;  sv[2][7] = 3.0;  sv[2][8] = 5e-3;  sv[2][9] = 0.017;

  CheckConstraintFDJacobian(constraint, phase_prefixes, si, sv);
}

// ── Cross-instance isolation: Jacobian entries only in own instance block ──

TEST(DissolvedEquilibriumConstraint, CrossInstanceIsolation)
{
  double K_eq = 5.0;
  auto keq = [K_eq](const micm::Conditions&) { return K_eq; };
  auto constraint = DissolvedEquilibriumConstraintBuilder()
      .SetPhase(aqueous_phase)
      .SetReactants({ A })
      .SetProducts({ B })
      .SetAlgebraicSpecies(B)
      .SetSolvent(h2o)
      .SetEquilibriumConstant(keq)
      .Build();

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("MODE1");
  phase_prefixes["AQUEOUS"].insert("MODE2");

  std::unordered_map<std::string, std::size_t> si;
  si["MODE1.AQUEOUS.A"] = 0;
  si["MODE1.AQUEOUS.B"] = 1;
  si["MODE1.AQUEOUS.H2O"] = 2;
  si["MODE2.AQUEOUS.A"] = 3;
  si["MODE2.AQUEOUS.B"] = 4;
  si["MODE2.AQUEOUS.H2O"] = 5;

  // Sparsity: MODE1's algebraic row (1) shouldn't have entries in MODE2 columns (3,4,5) and vice versa
  auto elements = constraint.NonZeroConstraintJacobianElements(phase_prefixes, si);
  for (const auto& [row, col] : elements)
  {
    if (row == 1)
    {
      EXPECT_LE(col, 2u) << "MODE1 algebraic row depends on MODE2 variable " << col;
    }
    else if (row == 4)
    {
      EXPECT_GE(col, 3u) << "MODE2 algebraic row depends on MODE1 variable " << col;
    }
  }
}
