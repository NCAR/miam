// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <miam/constraints/linear_constraint.hpp>
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
  auto A_g = micm::Species{ "A_g" };
  auto A_aq = micm::Species{ "A_aq" };
  auto B_aq = micm::Species{ "B_aq" };
  auto hp = micm::Species{ "H+" };
  auto am = micm::Species{ "A-" };
  auto bm = micm::Species{ "B-" };
  auto gas_phase = micm::Phase{ "GAS", { { A_g } } };
  auto aqueous_phase = micm::Phase{ "AQUEOUS", { { h2o }, { A_aq }, { B_aq }, { hp }, { am }, { bm } } };
}  // namespace

// ── Global constraint: algebraic in non-instanced (gas) phase ──

TEST(LinearConstraint, AlgebraicVariableNamesGlobal)
{
  // Mass conservation: [A_g] + Σ [A_aq_i] = total
  LinearConstraint constraint(
      gas_phase,
      A_g,
      { { gas_phase, A_g, 1.0 }, { aqueous_phase, A_aq, 1.0 } },
      100.0);

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("SMALL");
  phase_prefixes["AQUEOUS"].insert("LARGE");
  // GAS not in phase_prefixes → non-instanced

  auto names = constraint.ConstraintAlgebraicVariableNames(phase_prefixes);
  EXPECT_EQ(names.size(), 1);
  EXPECT_TRUE(names.count("A_g"));  // gas species, no prefix
}

TEST(LinearConstraint, SpeciesDependenciesGlobal)
{
  LinearConstraint constraint(
      gas_phase,
      A_g,
      { { gas_phase, A_g, 1.0 }, { aqueous_phase, A_aq, 1.0 } },
      100.0);

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("SMALL");
  phase_prefixes["AQUEOUS"].insert("LARGE");

  auto deps = constraint.ConstraintSpeciesDependencies(phase_prefixes);
  // A_g (gas) + SMALL.AQUEOUS.A_aq + LARGE.AQUEOUS.A_aq = 3
  EXPECT_EQ(deps.size(), 3);
  EXPECT_TRUE(deps.count("A_g"));
  EXPECT_TRUE(deps.count("SMALL.AQUEOUS.A_aq"));
  EXPECT_TRUE(deps.count("LARGE.AQUEOUS.A_aq"));
}

TEST(LinearConstraint, NonZeroJacobianGlobal)
{
  LinearConstraint constraint(
      gas_phase,
      A_g,
      { { gas_phase, A_g, 1.0 }, { aqueous_phase, A_aq, 1.0 } },
      100.0);

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("LARGE");
  phase_prefixes["AQUEOUS"].insert("SMALL");

  std::unordered_map<std::string, std::size_t> state_indices;
  state_indices["A_g"] = 0;
  state_indices["LARGE.AQUEOUS.A_aq"] = 1;
  state_indices["SMALL.AQUEOUS.A_aq"] = 2;

  auto elements = constraint.NonZeroConstraintJacobianElements(phase_prefixes, state_indices);
  // Single algebraic row (A_g = 0). Columns: A_g (0), LARGE.A_aq (1), SMALL.A_aq (2)
  EXPECT_EQ(elements.size(), 3);
  EXPECT_TRUE(elements.count({ 0, 0 }));
  EXPECT_TRUE(elements.count({ 0, 1 }));
  EXPECT_TRUE(elements.count({ 0, 2 }));
}

TEST(LinearConstraint, ResidualGlobal)
{
  // G = 1.0*[A_g] + 1.0*[A_aq_LARGE] + 1.0*[A_aq_SMALL] - 100.0
  LinearConstraint constraint(
      gas_phase,
      A_g,
      { { gas_phase, A_g, 1.0 }, { aqueous_phase, A_aq, 1.0 } },
      100.0);

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("LARGE");
  phase_prefixes["AQUEOUS"].insert("SMALL");

  std::unordered_map<std::string, std::size_t> state_indices;
  state_indices["A_g"] = 0;
  state_indices["LARGE.AQUEOUS.A_aq"] = 1;
  state_indices["SMALL.AQUEOUS.A_aq"] = 2;

  using DMP = micm::Matrix<double>;
  auto residual_fn = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, state_indices);

  DMP state_variables{ 1, 3, 0.0 };
  state_variables[0][0] = 50.0;  // [A_g]
  state_variables[0][1] = 30.0;  // LARGE [A_aq]
  state_variables[0][2] = 20.0;  // SMALL [A_aq]

  DMP residual{ 1, 3, 0.0 };
  residual_fn(state_variables, residual);

  // G = 50 + 30 + 20 - 100 = 0
  EXPECT_NEAR(residual[0][0], 0.0, 1.0e-12);

  // Non-equilibrium case
  state_variables[0][0] = 60.0;
  for (auto& v : residual.AsVector())
    v = 0.0;
  residual_fn(state_variables, residual);
  // G = 60 + 30 + 20 - 100 = 10
  EXPECT_NEAR(residual[0][0], 10.0, 1.0e-12);
}

TEST(LinearConstraint, JacobianGlobal)
{
  LinearConstraint constraint(
      gas_phase,
      A_g,
      { { gas_phase, A_g, 1.0 }, { aqueous_phase, A_aq, 1.0 } },
      100.0);

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("LARGE");
  phase_prefixes["AQUEOUS"].insert("SMALL");

  std::unordered_map<std::string, std::size_t> state_indices;
  state_indices["A_g"] = 0;
  state_indices["LARGE.AQUEOUS.A_aq"] = 1;
  state_indices["SMALL.AQUEOUS.A_aq"] = 2;

  using DMP = micm::Matrix<double>;
  using SMP = micm::SparseMatrix<double, micm::SparseMatrixStandardOrderingCompressedSparseRow>;
  auto nz_elements = constraint.NonZeroConstraintJacobianElements(phase_prefixes, state_indices);
  auto builder = SMP::Create(3).SetNumberOfBlocks(1);
  for (const auto& [row, col] : nz_elements)
    builder.WithElement(row, col);
  SMP jacobian(builder);

  auto jac_fn = constraint.ConstraintJacobianFunction<DMP, SMP>(phase_prefixes, state_indices, jacobian);

  DMP state_variables{ 1, 3, 0.0 };  // values don't matter for linear Jacobian

  for (auto& v : jacobian.AsVector())
    v = 0.0;
  jac_fn(state_variables, jacobian);

  // All coefficients are 1.0; jac -= 1.0 → jac = -1.0
  EXPECT_NEAR(jacobian.AsVector()[jacobian.VectorIndex(0, 0, 0)], -1.0, 1.0e-12);
  EXPECT_NEAR(jacobian.AsVector()[jacobian.VectorIndex(0, 0, 1)], -1.0, 1.0e-12);
  EXPECT_NEAR(jacobian.AsVector()[jacobian.VectorIndex(0, 0, 2)], -1.0, 1.0e-12);
}

// ── Per-instance constraint: algebraic in instanced phase ──

TEST(LinearConstraint, AlgebraicVariableNamesPerInstance)
{
  // Charge balance: [H+] - [A-] - [B-] = 0, H+ algebraic
  LinearConstraint constraint(
      aqueous_phase,
      hp,
      { { aqueous_phase, hp, 1.0 }, { aqueous_phase, am, -1.0 }, { aqueous_phase, bm, -1.0 } },
      0.0);

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("SMALL");
  phase_prefixes["AQUEOUS"].insert("LARGE");

  auto names = constraint.ConstraintAlgebraicVariableNames(phase_prefixes);
  EXPECT_EQ(names.size(), 2);
  EXPECT_TRUE(names.count("SMALL.AQUEOUS.H+"));
  EXPECT_TRUE(names.count("LARGE.AQUEOUS.H+"));
}

TEST(LinearConstraint, SpeciesDependenciesPerInstance)
{
  LinearConstraint constraint(
      aqueous_phase,
      hp,
      { { aqueous_phase, hp, 1.0 }, { aqueous_phase, am, -1.0 }, { aqueous_phase, bm, -1.0 } },
      0.0);

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("SMALL");
  phase_prefixes["AQUEOUS"].insert("LARGE");

  auto deps = constraint.ConstraintSpeciesDependencies(phase_prefixes);
  // 3 species * 2 instances = 6
  EXPECT_EQ(deps.size(), 6);
  EXPECT_TRUE(deps.count("SMALL.AQUEOUS.H+"));
  EXPECT_TRUE(deps.count("SMALL.AQUEOUS.A-"));
  EXPECT_TRUE(deps.count("SMALL.AQUEOUS.B-"));
  EXPECT_TRUE(deps.count("LARGE.AQUEOUS.H+"));
  EXPECT_TRUE(deps.count("LARGE.AQUEOUS.A-"));
  EXPECT_TRUE(deps.count("LARGE.AQUEOUS.B-"));
}

TEST(LinearConstraint, ResidualPerInstance)
{
  // Charge balance: G_i = [H+_i] - [A-_i] - [B-_i] = 0
  LinearConstraint constraint(
      aqueous_phase,
      hp,
      { { aqueous_phase, hp, 1.0 }, { aqueous_phase, am, -1.0 }, { aqueous_phase, bm, -1.0 } },
      0.0);

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("LARGE");
  phase_prefixes["AQUEOUS"].insert("SMALL");

  std::unordered_map<std::string, std::size_t> state_indices;
  state_indices["LARGE.AQUEOUS.H+"] = 0;
  state_indices["LARGE.AQUEOUS.A-"] = 1;
  state_indices["LARGE.AQUEOUS.B-"] = 2;
  state_indices["SMALL.AQUEOUS.H+"] = 3;
  state_indices["SMALL.AQUEOUS.A-"] = 4;
  state_indices["SMALL.AQUEOUS.B-"] = 5;

  using DMP = micm::Matrix<double>;
  auto residual_fn = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, state_indices);

  DMP state_variables{ 1, 6, 0.0 };
  state_variables[0][0] = 1.0e-4;  // LARGE H+
  state_variables[0][1] = 5.0e-5;  // LARGE A-
  state_variables[0][2] = 3.0e-5;  // LARGE B-
  state_variables[0][3] = 2.0e-4;  // SMALL H+
  state_variables[0][4] = 1.2e-4;  // SMALL A-
  state_variables[0][5] = 8.0e-5;  // SMALL B-

  DMP residual{ 1, 6, 0.0 };
  residual_fn(state_variables, residual);

  // LARGE: G = 1e-4 - 5e-5 - 3e-5 = 2e-5
  EXPECT_NEAR(residual[0][0], 2.0e-5, 1.0e-18);
  // SMALL: G = 2e-4 - 1.2e-4 - 8e-5 = 0
  EXPECT_NEAR(residual[0][3], 0.0, 1.0e-18);
}

TEST(LinearConstraint, JacobianPerInstance)
{
  LinearConstraint constraint(
      aqueous_phase,
      hp,
      { { aqueous_phase, hp, 1.0 }, { aqueous_phase, am, -1.0 }, { aqueous_phase, bm, -1.0 } },
      0.0);

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("LARGE");
  phase_prefixes["AQUEOUS"].insert("SMALL");

  std::unordered_map<std::string, std::size_t> state_indices;
  state_indices["LARGE.AQUEOUS.H+"] = 0;
  state_indices["LARGE.AQUEOUS.A-"] = 1;
  state_indices["LARGE.AQUEOUS.B-"] = 2;
  state_indices["SMALL.AQUEOUS.H+"] = 3;
  state_indices["SMALL.AQUEOUS.A-"] = 4;
  state_indices["SMALL.AQUEOUS.B-"] = 5;

  auto nz_elements = constraint.NonZeroConstraintJacobianElements(phase_prefixes, state_indices);
  // LARGE: (0,0), (0,1), (0,2); SMALL: (3,3), (3,4), (3,5)
  EXPECT_EQ(nz_elements.size(), 6);

  using DMP = micm::Matrix<double>;
  using SMP = micm::SparseMatrix<double, micm::SparseMatrixStandardOrderingCompressedSparseRow>;
  auto builder = SMP::Create(6).SetNumberOfBlocks(1);
  for (const auto& [row, col] : nz_elements)
    builder.WithElement(row, col);
  SMP jacobian(builder);

  auto jac_fn = constraint.ConstraintJacobianFunction<DMP, SMP>(phase_prefixes, state_indices, jacobian);

  DMP state_variables{ 1, 6, 0.0 };
  for (auto& v : jacobian.AsVector())
    v = 0.0;
  jac_fn(state_variables, jacobian);

  // LARGE algebraic row (0): jac -= coeff => jac[0,0] = -1, jac[0,1] = +1, jac[0,2] = +1
  EXPECT_NEAR(jacobian.AsVector()[jacobian.VectorIndex(0, 0, 0)], -1.0, 1.0e-12);
  EXPECT_NEAR(jacobian.AsVector()[jacobian.VectorIndex(0, 0, 1)], 1.0, 1.0e-12);
  EXPECT_NEAR(jacobian.AsVector()[jacobian.VectorIndex(0, 0, 2)], 1.0, 1.0e-12);

  // SMALL algebraic row (3): same pattern
  EXPECT_NEAR(jacobian.AsVector()[jacobian.VectorIndex(0, 3, 3)], -1.0, 1.0e-12);
  EXPECT_NEAR(jacobian.AsVector()[jacobian.VectorIndex(0, 3, 4)], 1.0, 1.0e-12);
  EXPECT_NEAR(jacobian.AsVector()[jacobian.VectorIndex(0, 3, 5)], 1.0, 1.0e-12);
}

// ── Builder ──

TEST(LinearConstraint, BuilderValidation)
{
  // Missing algebraic species
  EXPECT_THROW(
      LinearConstraintBuilder().AddTerm(aqueous_phase, hp, 1.0).Build(),
      std::runtime_error);

  // Missing terms
  EXPECT_THROW(
      LinearConstraintBuilder().SetAlgebraicSpecies(aqueous_phase, hp).Build(),
      std::runtime_error);

  // Valid build with constant
  auto constraint = LinearConstraintBuilder()
                        .SetAlgebraicSpecies(gas_phase, A_g)
                        .AddTerm(gas_phase, A_g, 1.0)
                        .AddTerm(aqueous_phase, A_aq, 1.0)
                        .SetConstant(42.0)
                        .Build();
  EXPECT_EQ(constraint.algebraic_species_.name_, "A_g");
  EXPECT_DOUBLE_EQ(constraint.constant_, 42.0);
}

// ── Update constraint parameters (no-op) ──

TEST(LinearConstraint, UpdateConstraintParametersNoOp)
{
  LinearConstraint constraint(
      gas_phase,
      A_g,
      { { gas_phase, A_g, 1.0 } },
      0.0);

  auto update_fn = constraint.UpdateConstraintParametersFunction<micm::Matrix<double>>();
  std::vector<micm::Conditions> conditions(3);
  // Should not throw
  update_fn(conditions);
}

// ── Mixed coefficients with non-unit values ──

TEST(LinearConstraint, ResidualWithNonUnitCoefficients)
{
  // G = 2.0*[A_g] + 0.5*[A_aq] - 10.0
  LinearConstraint constraint(
      gas_phase,
      A_g,
      { { gas_phase, A_g, 2.0 }, { aqueous_phase, A_aq, 0.5 } },
      10.0);

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  std::unordered_map<std::string, std::size_t> state_indices;
  state_indices["A_g"] = 0;
  state_indices["DROP.AQUEOUS.A_aq"] = 1;

  using DMP = micm::Matrix<double>;
  auto residual_fn = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, state_indices);

  DMP state_variables{ 1, 2, 0.0 };
  state_variables[0][0] = 3.0;  // [A_g]
  state_variables[0][1] = 8.0;  // [A_aq]

  DMP residual{ 1, 2, 0.0 };
  residual_fn(state_variables, residual);

  // G = 2*3 + 0.5*8 - 10 = 6 + 4 - 10 = 0
  EXPECT_NEAR(residual[0][0], 0.0, 1.0e-12);
}
