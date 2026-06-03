// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <miam/constraints/dissolved_equilibrium_constraint.hpp>
#include <miam/constraints/dissolved_equilibrium_constraint_builder.hpp>
#include <micm/util/jacobian_verification.hpp>
#include <miam/util/miam_exception.hpp>

#include <micm/system/conditions.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
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

  auto h2o = micm::Species{ "H2O" };
  auto A = micm::Species{ "A" };
  auto B = micm::Species{ "B" };
  auto C = micm::Species{ "C" };
  auto hp = micm::Species{ "H+" };
  auto aqueous_phase = micm::Phase{ "AQUEOUS", { { h2o }, { A }, { B }, { C }, { hp } } };

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

  // Initialize K_eq values into state parameters
  using DMP = micm::Matrix<double>;
  auto pi = BuildParamIndices(constraint, phase_prefixes);
  DMP state_params{ 1, std::max(pi.size(), std::size_t{ 1 }), 0.0 };
  std::vector<micm::Conditions> conditions(1);
  auto update_fn = constraint.UpdateConstraintParametersFunction<DMP>(phase_prefixes, pi);
  update_fn(conditions, state_params);

  auto residual_fn = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, pi, state_indices);

  DMP state_variables{ 1, 3, 0.0 };
  state_variables[0][0] = 1.0;    // [A]
  state_variables[0][1] = 5.0;    // [B]
  state_variables[0][2] = 300.0;  // [H2O]

  DMP residual{ 1, 3, 0.0 };
  residual_fn(state_variables, state_params, residual);

  // G = K_eq * [A] / [S]^(1-1) - [B] / [S]^(1-1) = K_eq * [A] - [B]
  // = 10.0 * 1.0 - 5.0 = 5.0
  EXPECT_NEAR(residual[0][1], 5.0, 1e-6);  // residual written to algebraic species (B) row
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

  using DMP = micm::Matrix<double>;
  auto pi = BuildParamIndices(constraint, phase_prefixes);
  DMP state_params{ 1, std::max(pi.size(), std::size_t{ 1 }), 0.0 };
  std::vector<micm::Conditions> conditions(1);
  auto update_fn = constraint.UpdateConstraintParametersFunction<DMP>(phase_prefixes, pi);
  update_fn(conditions, state_params);

  auto residual_fn = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, pi, state_indices);

  DMP state_variables{ 1, 4, 0.0 };
  state_variables[0][0] = 3.0;    // [A]
  state_variables[0][1] = 4.0;    // [B]
  state_variables[0][2] = 20.0;   // [C]
  state_variables[0][3] = 300.0;  // [H2O]

  DMP residual{ 1, 4, 0.0 };
  residual_fn(state_variables, state_params, residual);

  // G = K_eq * [A]*[B]/[S]^(2-1) - [C]/[S]^(1-1)
  // = 2.0 * 3.0*4.0/300.0 - 20.0 = 0.08 - 20.0 = -19.92
  double expected = K_eq * 3.0 * 4.0 / 300.0 - 20.0;
  EXPECT_NEAR(residual[0][2], expected, 1e-6);
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

  using DMP = micm::Matrix<double>;
  auto pi = BuildParamIndices(constraint, phase_prefixes);
  DMP state_params{ 1, std::max(pi.size(), std::size_t{ 1 }), 0.0 };
  std::vector<micm::Conditions> conditions(1);
  auto update_fn = constraint.UpdateConstraintParametersFunction<DMP>(phase_prefixes, pi);
  update_fn(conditions, state_params);

  auto residual_fn = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, pi, state_indices);

  DMP state_variables{ 1, 6, 0.0 };
  state_variables[0][0] = 2.0;    // LARGE [A]
  state_variables[0][1] = 8.0;    // LARGE [B]
  state_variables[0][2] = 300.0;  // LARGE [H2O]
  state_variables[0][3] = 1.0;    // SMALL [A]
  state_variables[0][4] = 3.0;    // SMALL [B]
  state_variables[0][5] = 300.0;  // SMALL [H2O]

  DMP residual{ 1, 6, 0.0 };
  residual_fn(state_variables, state_params, residual);

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

  using DMP = micm::Matrix<double>;
  auto pi = BuildParamIndices(constraint, phase_prefixes);
  DMP state_params{ 2, std::max(pi.size(), std::size_t{ 1 }), 0.0 };
  std::vector<micm::Conditions> conditions(2);
  auto update_fn = constraint.UpdateConstraintParametersFunction<DMP>(phase_prefixes, pi);
  update_fn(conditions, state_params);

  auto residual_fn = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, pi, state_indices);

  DMP state_variables{ 2, 3, 0.0 };
  state_variables[0][0] = 1.0;    // cell 0 [A]
  state_variables[0][1] = 2.0;    // cell 0 [B]
  state_variables[0][2] = 300.0;  // cell 0 [H2O]
  state_variables[1][0] = 3.0;    // cell 1 [A]
  state_variables[1][1] = 10.0;   // cell 1 [B]
  state_variables[1][2] = 300.0;  // cell 1 [H2O]

  DMP residual{ 2, 3, 0.0 };
  residual_fn(state_variables, state_params, residual);

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

  using DMP = micm::Matrix<double>;
  using SMP = micm::SparseMatrix<double, micm::SparseMatrixStandardOrderingCompressedSparseRow>;
  auto pi = BuildParamIndices(constraint, phase_prefixes);
  DMP state_params{ 1, std::max(pi.size(), std::size_t{ 1 }), 0.0 };
  std::vector<micm::Conditions> conditions(1);
  auto update_fn = constraint.UpdateConstraintParametersFunction<DMP>(phase_prefixes, pi);
  update_fn(conditions, state_params);

  // Build sparse Jacobian with the right sparsity
  auto nz_elements = constraint.NonZeroConstraintJacobianElements(phase_prefixes, state_indices);
  auto builder = SMP::Create(3).SetNumberOfBlocks(1);
  for (const auto& [row, col] : nz_elements)
    builder.WithElement(row, col);
  SMP jacobian(builder);

  auto jac_fn = constraint.ConstraintJacobianFunction<DMP, SMP>(phase_prefixes, pi, state_indices, jacobian);

  DMP state_variables{ 1, 3, 0.0 };
  state_variables[0][0] = 1.0;    // [A]
  state_variables[0][1] = 5.0;    // [B]
  state_variables[0][2] = 300.0;  // [H2O]

  // Zero the Jacobian
  for (auto& v : jacobian.AsVector())
    v = 0.0;

  jac_fn(state_variables, state_params, jacobian);

  // MICM convention: jac -= dG/dy
  // jac[B,A] -= K_eq → jac[B,A] = -K_eq = -10.0
  EXPECT_NEAR(jacobian.AsVector()[jacobian.VectorIndex(0, 1, 0)], -K_eq, 1e-6);
  // jac[B,B] -= (-1) → jac[B,B] = +1.0
  EXPECT_NEAR(jacobian.AsVector()[jacobian.VectorIndex(0, 1, 1)], 1.0, 1e-6);
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

  using DMP = micm::Matrix<double>;
  using SMP = micm::SparseMatrix<double, micm::SparseMatrixStandardOrderingCompressedSparseRow>;
  auto pi = BuildParamIndices(constraint, phase_prefixes);
  DMP state_params{ 1, std::max(pi.size(), std::size_t{ 1 }), 0.0 };
  std::vector<micm::Conditions> conditions(1);
  auto update_fn = constraint.UpdateConstraintParametersFunction<DMP>(phase_prefixes, pi);
  update_fn(conditions, state_params);

  auto nz_elements = constraint.NonZeroConstraintJacobianElements(phase_prefixes, state_indices);
  auto builder = SMP::Create(4).SetNumberOfBlocks(1);
  for (const auto& [row, col] : nz_elements)
    builder.WithElement(row, col);
  SMP jacobian(builder);

  auto jac_fn = constraint.ConstraintJacobianFunction<DMP, SMP>(phase_prefixes, pi, state_indices, jacobian);

  DMP state_variables{ 1, 4, 0.0 };
  state_variables[0][0] = 3.0;    // [A]
  state_variables[0][1] = 4.0;    // [B]
  state_variables[0][2] = 20.0;   // [C]
  state_variables[0][3] = 300.0;  // [S]

  for (auto& v : jacobian.AsVector())
    v = 0.0;
  jac_fn(state_variables, state_params, jacobian);

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

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  using DMP = micm::Matrix<double>;
  auto pi = BuildParamIndices(constraint, phase_prefixes);
  ASSERT_EQ(pi.size(), 1u);

  std::vector<micm::Conditions> conditions(3);
  conditions[0].temperature_ = 250.0;
  conditions[1].temperature_ = 300.0;
  conditions[2].temperature_ = 350.0;

  DMP state_params{ 3, 1, 0.0 };
  auto update_fn = constraint.UpdateConstraintParametersFunction<DMP>(phase_prefixes, pi);
  update_fn(conditions, state_params);

  std::size_t keq_col = pi.begin()->second;
  EXPECT_NEAR(state_params[0][keq_col], 1000.0 / 250.0, 1.0e-12);
  EXPECT_NEAR(state_params[1][keq_col], 1000.0 / 300.0, 1.0e-12);
  EXPECT_NEAR(state_params[2][keq_col], 1000.0 / 350.0, 1.0e-12);
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
      miam::MiamException);

  // Missing phase
  EXPECT_THROW(
      DissolvedEquilibriumConstraintBuilder().SetReactants({ A }).SetProducts({ B }).SetAlgebraicSpecies(B).SetSolvent(h2o).SetEquilibriumConstant(FakeConstant{}).Build(),
      miam::MiamException);

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

  /// @brief Initialize K_eq values into a state parameter matrix; returns the matrix
  DMP InitKeq(
      const DissolvedEquilibriumConstraint& constraint,
      const std::map<std::string, std::set<std::string>>& phase_prefixes,
      const std::unordered_map<std::string, std::size_t>& param_idx,
      std::size_t num_cells)
  {
    std::vector<micm::Conditions> conditions(num_cells);
    DMP state_params{ num_cells, std::max(param_idx.size(), std::size_t{ 1 }), 0.0 };
    auto fn = constraint.UpdateConstraintParametersFunction<DMP>(phase_prefixes, param_idx);
    fn(conditions, state_params);
    return state_params;
  }

  DMP InitKeq(
      const DissolvedEquilibriumConstraint& constraint,
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
      const DissolvedEquilibriumConstraint& constraint,
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
    EXPECT_TRUE(spc.passed_) << "Missing sparsity entry: block=" << spc.worst_block_
                            << " row=" << spc.worst_row_ << " col=" << spc.worst_col_
                            << " fd=" << spc.worst_fd_;
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

  auto pi = BuildParamIndices(constraint, phase_prefixes);
  auto sp = InitKeq(constraint, phase_prefixes, pi, 1);

  DMP sv{ 1, 3, 0.0 };
  sv[0][0] = 1.0;
  sv[0][1] = 5.0;
  sv[0][2] = 300.0;

  CheckConstraintFDJacobian(constraint, phase_prefixes, pi, si, sv, sp);
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

  auto pi = BuildParamIndices(constraint, phase_prefixes);
  auto sp = InitKeq(constraint, phase_prefixes, pi, 1);

  DMP sv{ 1, 4, 0.0 };
  sv[0][0] = 3.0;
  sv[0][1] = 4.0;
  sv[0][2] = 20.0;
  sv[0][3] = 300.0;

  CheckConstraintFDJacobian(constraint, phase_prefixes, pi, si, sv, sp);
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

  auto pi = BuildParamIndices(constraint, phase_prefixes);
  auto sp = InitKeq(constraint, phase_prefixes, pi, 1);

  DMP sv{ 1, 4, 0.0 };
  sv[0][0] = 0.1;
  sv[0][1] = 0.01;
  sv[0][2] = 1.0e-4;
  sv[0][3] = 0.017;

  CheckConstraintFDJacobian(constraint, phase_prefixes, pi, si, sv, sp);
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

  auto pi = BuildParamIndices(constraint, phase_prefixes);
  auto sp = InitKeq(constraint, phase_prefixes, pi, 1);

  DMP sv{ 1, 5, 0.0 };
  sv[0][0] = 0.5;
  sv[0][1] = 0.3;
  sv[0][2] = 0.1;
  sv[0][3] = 1.0e-3;
  sv[0][4] = 0.017;

  CheckConstraintFDJacobian(constraint, phase_prefixes, pi, si, sv, sp);
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
  auto pi = BuildParamIndices(constraint, phase_prefixes);
  auto sp = InitKeq(constraint, phase_prefixes, pi, nc);

  DMP sv{ nc, 3, 0.0 };
  sv[0][0] = 1.0;   sv[0][1] = 5.0;   sv[0][2] = 300.0;
  sv[1][0] = 0.5;   sv[1][1] = 2.5;   sv[1][2] = 0.017;
  sv[2][0] = 10.0;  sv[2][1] = 40.0;  sv[2][2] = 100.0;
  sv[3][0] = 0.01;  sv[3][1] = 0.1;   sv[3][2] = 0.017;

  // Check residuals
  auto rf = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, pi, si);
  DMP residual{ nc, 3, 0.0 };
  rf(sv, sp, residual);

  // G = K_eq * [A] - [B]
  EXPECT_NEAR(residual[0][1], 5.0 * 1.0 - 5.0, 1e-6);
  EXPECT_NEAR(residual[1][1], 5.0 * 0.5 - 2.5, 1e-6);
  EXPECT_NEAR(residual[2][1], 5.0 * 10.0 - 40.0, 1e-6);
  EXPECT_NEAR(residual[3][1], 5.0 * 0.01 - 0.1, 1e-6);

  // FD Jacobian check
  CheckConstraintFDJacobian(constraint, phase_prefixes, pi, si, sv, sp);
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
  auto pi = BuildParamIndices(constraint, phase_prefixes);
  auto sp = InitKeq(constraint, phase_prefixes, pi, nc);

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
  auto rf = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, pi, si);
  DMP residual{ nc, 8, 0.0 };
  rf(sv, sp, residual);

  // LARGE instance, cell 0: G = 3.0 * 0.5*0.3/0.017 - 0.1
  double expected_L0 = 3.0 * 0.5 * 0.3 / 0.017 - 0.1;
  EXPECT_NEAR(residual[0][2], expected_L0, 1e-6);

  // SMALL instance, cell 0: G = 3.0 * 0.2*0.1/0.017 - 0.05
  double expected_S0 = 3.0 * 0.2 * 0.1 / 0.017 - 0.05;
  EXPECT_NEAR(residual[0][6], expected_S0, 1e-6);

  // Non-algebraic rows untouched
  EXPECT_NEAR(residual[0][0], 0.0, 1e-30);
  EXPECT_NEAR(residual[0][1], 0.0, 1e-30);
  EXPECT_NEAR(residual[0][3], 0.0, 1e-30);

  // FD check
  CheckConstraintFDJacobian(constraint, phase_prefixes, pi, si, sv, sp);
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
  auto pi = BuildParamIndices(constraint, phase_prefixes);
  auto sp = InitKeq(constraint, phase_prefixes, pi, nc);

  DMP sv{ nc, 12, 0.0 };
  // Cell 0
  sv[0][0] = 0.5;  sv[0][1] = 0.1;  sv[0][2] = 1e-4;  sv[0][3] = 0.017;
  sv[0][4] = 1.0;  sv[0][5] = 0.2;  sv[0][6] = 2e-4;  sv[0][7] = 0.017;
  sv[0][8] = 2.0;  sv[0][9] = 0.5;  sv[0][10] = 5e-4; sv[0][11] = 0.017;
  // Cell 1
  sv[1][0] = 3.0;  sv[1][1] = 1.0;  sv[1][2] = 1e-3;  sv[1][3] = 40.0;
  sv[1][4] = 4.0;  sv[1][5] = 2.0;  sv[1][6] = 2e-3;  sv[1][7] = 50.0;
  sv[1][8] = 5.0;  sv[1][9] = 3.0;  sv[1][10] = 3e-3; sv[1][11] = 60.0;

  CheckConstraintFDJacobian(constraint, phase_prefixes, pi, si, sv, sp);
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

  auto pi = BuildParamIndices(constraint, phase_prefixes);
  auto sp = InitKeq(constraint, phase_prefixes, pi, 1);

  auto jacobian = BuildConstraintJacobian(constraint, phase_prefixes, si, 1);
  auto jac_fn = constraint.ConstraintJacobianFunction<DMP, SMP>(phase_prefixes, pi, si, jacobian);

  DMP sv{ 1, 3, 0.0 };
  sv[0][0] = 1.0;  sv[0][1] = 5.0;  sv[0][2] = 300.0;

  // First call
  jac_fn(sv, sp, jacobian);
  double j_BA_once = jacobian.AsVector()[jacobian.VectorIndex(0, 1, 0)];
  double j_BB_once = jacobian.AsVector()[jacobian.VectorIndex(0, 1, 1)];

  // Second call should accumulate
  jac_fn(sv, sp, jacobian);
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

  auto pi = BuildParamIndices(constraint, phase_prefixes);
  auto sp = InitKeq(constraint, phase_prefixes, pi, 1);

  auto rf = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, pi, si);

  DMP sv{ 1, 3, 0.0 };
  sv[0][0] = 1.0;  sv[0][1] = 5.0;  sv[0][2] = 300.0;

  DMP residual{ 1, 3, 999.0 };
  rf(sv, sp, residual);
  double val1 = residual[0][1];

  rf(sv, sp, residual);
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
  auto pi = BuildParamIndices(constraint, phase_prefixes);
  auto sp = InitKeq(constraint, phase_prefixes, pi, conditions);

  DMP sv{ nc, 3, 0.0 };
  sv[0][0] = 1.0;  sv[0][1] = 2.0;  sv[0][2] = 0.017;
  sv[1][0] = 0.5;  sv[1][1] = 1.0;  sv[1][2] = 0.017;
  sv[2][0] = 2.0;  sv[2][1] = 5.0;  sv[2][2] = 0.017;

  // Check residuals with different K_eq per cell
  auto rf = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, pi, si);
  DMP residual{ nc, 3, 0.0 };
  rf(sv, sp, residual);

  EXPECT_NEAR(residual[0][1], (1000.0 / 250.0) * 1.0 - 2.0, 1e-6);
  EXPECT_NEAR(residual[1][1], (1000.0 / 300.0) * 0.5 - 1.0, 1e-6);
  EXPECT_NEAR(residual[2][1], (1000.0 / 350.0) * 2.0 - 5.0, 1e-6);

  // FD check
  CheckConstraintFDJacobian(constraint, phase_prefixes, pi, si, sv, sp);
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

  auto pi = BuildParamIndices(constraint, phase_prefixes);
  auto sp = InitKeq(constraint, phase_prefixes, pi, 1);

  DMP sv{ 1, 3, 0.0 };
  sv[0][0] = 1.0e-6;   // small [A]
  sv[0][1] = 100.0;    // [B] = K_eq * A
  sv[0][2] = 0.017;

  auto rf = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, pi, si);
  DMP residual{ 1, 3, 0.0 };
  rf(sv, sp, residual);

  // G = 1e8 * 1e-6 - 100 = 100 - 100 = 0
  EXPECT_NEAR(residual[0][1], 0.0, 1e-6);

  CheckConstraintFDJacobian(constraint, phase_prefixes, pi, si, sv, sp);
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

  // Initialize both k_eq state params (different UUIDs, so different param names)
  auto pi_orig = BuildParamIndices(original, phase_prefixes);
  auto pi_copy = BuildParamIndices(copy, phase_prefixes);
  auto sp_orig = InitKeq(original, phase_prefixes, pi_orig, 1);
  auto sp_copy = InitKeq(copy, phase_prefixes, pi_copy, 1);

  DMP sv{ 1, 4, 0.0 };
  sv[0][0] = 3.0;  sv[0][1] = 4.0;  sv[0][2] = 20.0;  sv[0][3] = 300.0;

  auto rf_orig = original.ConstraintResidualFunction<DMP>(phase_prefixes, pi_orig, si);
  auto rf_copy = copy.ConstraintResidualFunction<DMP>(phase_prefixes, pi_copy, si);

  DMP res_orig{ 1, 4, 0.0 };
  DMP res_copy{ 1, 4, 0.0 };
  rf_orig(sv, sp_orig, res_orig);
  rf_copy(sv, sp_copy, res_copy);

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

  auto pi = BuildParamIndices(constraint, phase_prefixes);
  auto sp = InitKeq(constraint, phase_prefixes, pi, 1);

  auto jacobian = BuildConstraintJacobian(constraint, phase_prefixes, si, 1);
  auto jac_fn = constraint.ConstraintJacobianFunction<DMP, SMP>(phase_prefixes, pi, si, jacobian);

  DMP sv{ 1, 4, 0.0 };
  sv[0][0] = 3.0;
  sv[0][1] = 4.0;
  sv[0][2] = 20.0;
  sv[0][3] = 300.0;

  jac_fn(sv, sp, jacobian);

  // dG/dS = K_eq * [A]*[B] * (1-2)/[S]^2 - [C] * (1-1)/[S]^1
  //       = -K_eq * [A]*[B] / [S]^2 - 0
  double dG_dS = -K_eq * 3.0 * 4.0 / (300.0 * 300.0);
  // jac -= dG/dS
  EXPECT_NEAR(jacobian.AsVector()[jacobian.VectorIndex(0, 2, 3)], -dG_dS, 1e-12);

  // Also FD check
  CheckConstraintFDJacobian(constraint, phase_prefixes, pi, si, sv, sp);
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
  auto pi = BuildParamIndices(constraint, phase_prefixes);
  auto sp = InitKeq(constraint, phase_prefixes, pi, conditions);

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

  CheckConstraintFDJacobian(constraint, phase_prefixes, pi, si, sv, sp);
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

// ── VectorMatrix typed tests ──
// These tests verify the fix for the per-cell state-parameter bug with L>1.
// DissolvedEquilibriumConstraint K_eq is temperature-independent in this test,
// but the important check is that state_params indexing is correct across blocks.

namespace
{
  template<typename VDM, typename VSM>
  void TestDissolvedConstraintVectorMatrix(std::size_t num_cells, double K_eq_val)
  {
    auto keq_fn = [K_eq_val](const micm::Conditions&) { return K_eq_val; };
    auto constraint = DissolvedEquilibriumConstraintBuilder()
        .SetPhase(aqueous_phase)
        .SetReactants({ A })
        .SetProducts({ B })
        .SetAlgebraicSpecies(B)
        .SetSolvent(h2o)
        .SetEquilibriumConstant(keq_fn)
        .Build();

    std::map<std::string, std::set<std::string>> phase_prefixes;
    phase_prefixes["AQUEOUS"].insert("DROP");

    std::unordered_map<std::string, std::size_t> si;
    si["DROP.AQUEOUS.A"] = 0;
    si["DROP.AQUEOUS.B"] = 1;
    si["DROP.AQUEOUS.H2O"] = 2;

    // Build param indices
    std::unordered_map<std::string, std::size_t> pi;
    std::size_t idx = 0;
    for (const auto& name : constraint.ConstraintStateParameterNames(phase_prefixes))
      pi[name] = idx++;
    ASSERT_EQ(pi.size(), 1u);
    std::size_t keq_col = pi.begin()->second;

    // Update state parameters with distinct per-cell conditions
    std::vector<micm::Conditions> conditions(num_cells);
    for (std::size_t c = 0; c < num_cells; ++c)
      conditions[c].temperature_ = 298.15 + static_cast<double>(c) * 10.0;

    VDM state_params{ num_cells, std::max(pi.size(), std::size_t{ 1 }), 0.0 };
    auto update_fn = constraint.template UpdateConstraintParametersFunction<VDM>(phase_prefixes, pi);
    update_fn(conditions, state_params);

    // K_eq is temperature-independent in this test, so all cells have same value
    for (std::size_t c = 0; c < num_cells; ++c)
    {
      EXPECT_NEAR(state_params[c][keq_col], K_eq_val, K_eq_val * 1e-12)
          << "Cell " << c;
    }

    // Build state variables with per-cell values
    VDM state_variables{ num_cells, 3, 0.0 };
    for (std::size_t c = 0; c < num_cells; ++c)
    {
      state_variables[c][0] = 1.0 * static_cast<double>(c + 1);   // [A]
      state_variables[c][1] = 3.0 * static_cast<double>(c + 1);   // [B]
      state_variables[c][2] = 300.0;                               // [H2O]
    }

    // Residual check: G = K_eq * [A] - [B]
    VDM residual{ num_cells, 3, 0.0 };
    auto rf = constraint.template ConstraintResidualFunction<VDM>(phase_prefixes, pi, si);
    rf(state_variables, state_params, residual);

    for (std::size_t c = 0; c < num_cells; ++c)
    {
      double A_val = state_variables[c][0];
      double B_val = state_variables[c][1];
      double expected = K_eq_val * A_val - B_val;
      EXPECT_NEAR(residual[c][1], expected, std::max(std::abs(expected) * 1e-10, 1e-20))
          << "Residual mismatch cell " << c;
    }

    // Jacobian FD check
    auto nz = constraint.NonZeroConstraintJacobianElements(phase_prefixes, si);
    auto jac_builder = VSM::Create(si.size()).SetNumberOfBlocks(num_cells).InitialValue(0.0);
    for (const auto& [row, col] : nz)
      jac_builder = jac_builder.WithElement(row, col);
    VSM jacobian(jac_builder);

    auto jac_fn = constraint.template ConstraintJacobianFunction<VDM, VSM>(phase_prefixes, pi, si, jacobian);
    jac_fn(state_variables, state_params, jacobian);

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
            EXPECT_NEAR(analytical + fd, 0.0, scale * 1e-5)
                << "FD Jac mismatch cell=" << c << " row=" << i << " col=" << j;
          }
        }
      }
    }
  }
}  // namespace

TEST(DissolvedEquilibriumConstraint, VectorMatrix_L1_4cells)
{
  using VDM = micm::VectorMatrix<double, 1>;
  using VSM = micm::SparseMatrix<double, micm::SparseMatrixVectorOrderingCompressedSparseRow<1>>;
  TestDissolvedConstraintVectorMatrix<VDM, VSM>(4, 10.0);
}

TEST(DissolvedEquilibriumConstraint, VectorMatrix_L2_4cells)
{
  using VDM = micm::VectorMatrix<double, 2>;
  using VSM = micm::SparseMatrix<double, micm::SparseMatrixVectorOrderingCompressedSparseRow<2>>;
  TestDissolvedConstraintVectorMatrix<VDM, VSM>(4, 10.0);
}

TEST(DissolvedEquilibriumConstraint, VectorMatrix_L4_4cells)
{
  using VDM = micm::VectorMatrix<double, 4>;
  using VSM = micm::SparseMatrix<double, micm::SparseMatrixVectorOrderingCompressedSparseRow<4>>;
  TestDissolvedConstraintVectorMatrix<VDM, VSM>(4, 10.0);
}

// ============================================================================
// Limit / Extreme tests (Phase D1)
// ============================================================================

TEST(DissolvedEquilibriumConstraint, ResidualZeroReactant)
{
  // G = K_eq*[A]/[H2O]^0 - [B]/[H2O]^0 (n_r=1, n_p=1, same powers)
  // With [A]=0: G = -[B] (forward term vanishes, algebraic variable forced to 0)
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
  si["SMALL.AQUEOUS.A"] = 0; si["SMALL.AQUEOUS.B"] = 1; si["SMALL.AQUEOUS.H2O"] = 2;

  auto pi = BuildParamIndices(constraint, phase_prefixes);
  auto sp = InitKeq(constraint, phase_prefixes, pi, 1);

  DMP sv{ 1, 3, 0.0 };
  sv[0][0] = 0.0;   // [A] = 0
  sv[0][1] = 0.5;   // [B]
  sv[0][2] = 55.0;  // [H2O]

  auto rf = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, pi, si);
  DMP residual{ 1, 3, 0.0 };
  rf(sv, sp, residual);

  // G = K_eq * [A] / [H2O]^0 - [B] / [H2O]^0 = K_eq*0 - [B] = -0.5
  EXPECT_NEAR(residual[0][1], -0.5, 1.0e-12);
}

TEST(DissolvedEquilibriumConstraint, ResidualExtremeKeq)
{
  // Very large and very small K_eq should produce numerically stable residuals
  auto make_and_eval = [&](double K_eq_val, double A_conc, double B_conc) -> double
  {
    auto keq = [K_eq_val](const micm::Conditions&) { return K_eq_val; };
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
    si["SMALL.AQUEOUS.A"] = 0; si["SMALL.AQUEOUS.B"] = 1; si["SMALL.AQUEOUS.H2O"] = 2;

    auto pi = BuildParamIndices(constraint, phase_prefixes);
    auto sp = InitKeq(constraint, phase_prefixes, pi, 1);

    DMP sv{ 1, 3, 0.0 };
    sv[0][0] = A_conc; sv[0][1] = B_conc; sv[0][2] = 55.0;

    auto rf = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, pi, si);
    DMP residual{ 1, 3, 0.0 };
    rf(sv, sp, residual);
    return residual[0][1];
  };

  // K_eq = 1e15: at equilibrium [B] = K_eq * [A] → G = 0
  double large_K = 1.0e15;
  double A_large = 1.0e-10;
  double B_large_eq = large_K * A_large;  // 1e5
  EXPECT_NEAR(make_and_eval(large_K, A_large, B_large_eq), 0.0, 1.0e-3);  // relative precision limited by large values

  // K_eq = 1e-15: at equilibrium [B] = K_eq * [A] → G = 0
  double small_K = 1.0e-15;
  double A_small = 1.0;
  double B_small_eq = small_K * A_small;  // 1e-15
  EXPECT_NEAR(make_and_eval(small_K, A_small, B_small_eq), 0.0, 1.0e-27);

  // Results should be finite (no NaN/Inf)
  EXPECT_FALSE(std::isnan(make_and_eval(large_K, A_large, B_large_eq)));
  EXPECT_FALSE(std::isnan(make_and_eval(small_K, A_small, B_small_eq)));
}

TEST(DissolvedEquilibriumConstraint, JacobianFDZeroReactant)
{
  // FD Jacobian check with [A]=0 (forward term vanishes — Jacobian still well-defined)
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
  si["SMALL.AQUEOUS.A"] = 0; si["SMALL.AQUEOUS.B"] = 1; si["SMALL.AQUEOUS.H2O"] = 2;

  auto pi = BuildParamIndices(constraint, phase_prefixes);
  auto sp = InitKeq(constraint, phase_prefixes, pi, 1);

  DMP sv{ 1, 3, 0.0 };
  sv[0][0] = 0.0;   // [A] = 0
  sv[0][1] = 0.5;
  sv[0][2] = 55.0;

  CheckConstraintFDJacobian(constraint, phase_prefixes, pi, si, sv, sp);
}

TEST(DissolvedEquilibriumConstraint, JacobianFDExtremeKeq)
{
  // FD Jacobian check at extreme K_eq values
  for (double K_eq_val : { 1.0e-10, 1.0, 1.0e10 })
  {
    auto keq = [K_eq_val](const micm::Conditions&) { return K_eq_val; };
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
    si["SMALL.AQUEOUS.A"] = 0; si["SMALL.AQUEOUS.B"] = 1; si["SMALL.AQUEOUS.H2O"] = 2;

    auto pi = BuildParamIndices(constraint, phase_prefixes);
    auto sp = InitKeq(constraint, phase_prefixes, pi, 1);

    DMP sv{ 1, 3, 0.0 };
    sv[0][0] = 1.0; sv[0][1] = K_eq_val;  // at equilibrium
    sv[0][2] = 0.017;  // small [H2O] to stress-test

    CheckConstraintFDJacobian(constraint, phase_prefixes, pi, si, sv, sp);
  }
}
