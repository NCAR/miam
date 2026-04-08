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
  DissolvedEquilibriumConstraint constraint(keq, { A }, { B }, B, h2o, aqueous_phase);

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("SMALL");

  auto names = constraint.ConstraintAlgebraicVariableNames(phase_prefixes);
  EXPECT_EQ(names.size(), 1);
  EXPECT_TRUE(names.count("SMALL.AQUEOUS.B"));
}

TEST(DissolvedEquilibriumConstraint, AlgebraicVariableNamesMultiplePrefixes)
{
  auto keq = [](const micm::Conditions&) { return 10.0; };
  DissolvedEquilibriumConstraint constraint(keq, { A }, { B }, B, h2o, aqueous_phase);

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
  DissolvedEquilibriumConstraint constraint(keq, { A }, { B }, B, h2o, aqueous_phase);

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["OTHER"].insert("PREFIX");

  auto names = constraint.ConstraintAlgebraicVariableNames(phase_prefixes);
  EXPECT_EQ(names.size(), 0);
}

// ── ConstraintSpeciesDependencies ──

TEST(DissolvedEquilibriumConstraint, SpeciesDependencies)
{
  auto keq = [](const micm::Conditions&) { return 10.0; };
  DissolvedEquilibriumConstraint constraint(keq, { A }, { B, hp }, B, h2o, aqueous_phase);

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
  DissolvedEquilibriumConstraint constraint(keq, { A }, { B, hp }, B, h2o, aqueous_phase);

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
  DissolvedEquilibriumConstraint constraint(keq, { A }, { B }, B, h2o, aqueous_phase);

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
  auto residual_fn = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, state_indices);

  DMP state_variables{ 1, 3, 0.0 };
  state_variables[0][0] = 1.0;    // [A]
  state_variables[0][1] = 5.0;    // [B]
  state_variables[0][2] = 300.0;  // [H2O]

  DMP residual{ 1, 3, 0.0 };
  residual_fn(state_variables, residual);

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
  DissolvedEquilibriumConstraint constraint(keq, { A, B }, { C }, C, h2o, aqueous_phase);

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
  auto residual_fn = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, state_indices);

  DMP state_variables{ 1, 4, 0.0 };
  state_variables[0][0] = 3.0;    // [A]
  state_variables[0][1] = 4.0;    // [B]
  state_variables[0][2] = 20.0;   // [C]
  state_variables[0][3] = 300.0;  // [H2O]

  DMP residual{ 1, 4, 0.0 };
  residual_fn(state_variables, residual);

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
  DissolvedEquilibriumConstraint constraint(keq, { A }, { B }, B, h2o, aqueous_phase);

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
  auto residual_fn = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, state_indices);

  DMP state_variables{ 1, 6, 0.0 };
  state_variables[0][0] = 2.0;    // LARGE [A]
  state_variables[0][1] = 8.0;    // LARGE [B]
  state_variables[0][2] = 300.0;  // LARGE [H2O]
  state_variables[0][3] = 1.0;    // SMALL [A]
  state_variables[0][4] = 3.0;    // SMALL [B]
  state_variables[0][5] = 300.0;  // SMALL [H2O]

  DMP residual{ 1, 6, 0.0 };
  residual_fn(state_variables, residual);

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
  DissolvedEquilibriumConstraint constraint(keq, { A }, { B }, B, h2o, aqueous_phase);

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
  auto residual_fn = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, state_indices);

  DMP state_variables{ 2, 3, 0.0 };
  state_variables[0][0] = 1.0;    // cell 0 [A]
  state_variables[0][1] = 2.0;    // cell 0 [B]
  state_variables[0][2] = 300.0;  // cell 0 [H2O]
  state_variables[1][0] = 3.0;    // cell 1 [A]
  state_variables[1][1] = 10.0;   // cell 1 [B]
  state_variables[1][2] = 300.0;  // cell 1 [H2O]

  DMP residual{ 2, 3, 0.0 };
  residual_fn(state_variables, residual);

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
  DissolvedEquilibriumConstraint constraint(keq, { A }, { B }, B, h2o, aqueous_phase);

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
  DissolvedEquilibriumConstraint constraint(keq, { A, B }, { C }, C, h2o, aqueous_phase);

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
  DissolvedEquilibriumConstraint constraint(keq, { A }, { B }, B, h2o, aqueous_phase);

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
