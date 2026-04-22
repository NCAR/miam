// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <miam/constraints/linear_constraint.hpp>
#include <miam/constraints/linear_constraint_builder.hpp>
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
  auto residual_fn = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, param_indices, state_indices);

  DMP state_variables{ 1, 3, 0.0 };
  state_variables[0][0] = 50.0;  // [A_g]
  state_variables[0][1] = 30.0;  // LARGE [A_aq]
  state_variables[0][2] = 20.0;  // SMALL [A_aq]

  DMP residual{ 1, 3, 0.0 };
  residual_fn(state_variables, no_params, residual);

  // G = 50 + 30 + 20 - 100 = 0
  EXPECT_NEAR(residual[0][0], 0.0, 1.0e-12);

  // Non-equilibrium case
  state_variables[0][0] = 60.0;
  for (auto& v : residual.AsVector())
    v = 0.0;
  residual_fn(state_variables, no_params, residual);
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
  auto residual_fn = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, param_indices, state_indices);

  DMP state_variables{ 1, 6, 0.0 };
  state_variables[0][0] = 1.0e-4;  // LARGE H+
  state_variables[0][1] = 5.0e-5;  // LARGE A-
  state_variables[0][2] = 3.0e-5;  // LARGE B-
  state_variables[0][3] = 2.0e-4;  // SMALL H+
  state_variables[0][4] = 1.2e-4;  // SMALL A-
  state_variables[0][5] = 8.0e-5;  // SMALL B-

  DMP residual{ 1, 6, 0.0 };
  residual_fn(state_variables, no_params, residual);

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
  auto residual_fn = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, param_indices, state_indices);

  DMP state_variables{ 1, 2, 0.0 };
  state_variables[0][0] = 3.0;  // [A_g]
  state_variables[0][1] = 8.0;  // [A_aq]

  DMP residual{ 1, 2, 0.0 };
  residual_fn(state_variables, no_params, residual);

  // G = 2*3 + 0.5*8 - 10 = 6 + 4 - 10 = 0
  EXPECT_NEAR(residual[0][0], 0.0, 1.0e-12);
}

// ═══════════════════════════════════════════════════════════════════
// Expanded Scenarios
// ═══════════════════════════════════════════════════════════════════

namespace
{
  using DMP = micm::Matrix<double>;
  using SMP = micm::SparseMatrix<double, micm::SparseMatrixStandardOrderingCompressedSparseRow>;

  SMP BuildConstraintJacobian(
      const LinearConstraint& constraint,
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

  SMP BuildConstraintJacobian(
      const std::vector<std::reference_wrapper<const LinearConstraint>>& constraints,
      const std::map<std::string, std::set<std::string>>& phase_prefixes,
      const std::unordered_map<std::string, std::size_t>& state_indices,
      std::size_t num_blocks)
  {
    std::set<std::pair<std::size_t, std::size_t>> elements;
    for (const auto& c : constraints)
    {
      auto elts = c.get().NonZeroConstraintJacobianElements(phase_prefixes, state_indices);
      elements.insert(elts.begin(), elts.end());
    }
    auto builder = SMP::Create(state_indices.size()).SetNumberOfBlocks(num_blocks).InitialValue(0.0);
    for (const auto& [row, col] : elements)
      builder = builder.WithElement(row, col);
    return SMP(builder);
  }

  /// @brief FD check for constraint Jacobian vs residual
  /// @details For linear constraints, dG/dy is constant, so the FD check should be exact
  ///          up to floating point precision.
  void CheckConstraintFDJacobian(
      const LinearConstraint& constraint,
      const std::map<std::string, std::set<std::string>>& phase_prefixes,
      const std::unordered_map<std::string, std::size_t>& state_indices,
      const DMP& state_variables,
      double rel_tol = 1e-8)
  {
    std::size_t num_blocks = state_variables.NumRows();
    std::size_t num_vars = state_indices.size();

    // Build analytical Jacobian
    auto jacobian = BuildConstraintJacobian(constraint, phase_prefixes, state_indices, num_blocks);
    auto jac_fn = constraint.ConstraintJacobianFunction<DMP, SMP>(phase_prefixes, state_indices, jacobian);
    jac_fn(state_variables, jacobian);

    // Build residual function for FD
    // Use a relatively large step: for linear constraints the FD is exact
    // regardless of step size, and a larger step reduces floating-point noise.
    double eps = 1e-4;
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
          double fd = (res_plus[b][i] - res_minus[b][i]) / (2.0 * h);  // dG_i/dy_j
          // Jacobian stores -dG/dy, so analytical + fd ≈ 0
          double analytical;
          try
          {
            analytical = jacobian[b][i][j];
          }
          catch (...)
          {
            // Zero element in sparse matrix
            if (std::abs(fd) > 1e-12)
              ADD_FAILURE() << "Missing Jacobian at block=" << b << " row=" << i << " col=" << j << " fd=" << fd;
            continue;
          }

          double scale = std::max(std::abs(analytical), std::abs(fd));
          if (scale > 1e-20)
          {
            EXPECT_NEAR(analytical + fd, 0.0, scale * rel_tol)
                << "FD mismatch: block=" << b << " row=" << i << " col=" << j
                << " analytical(-dG/dy)=" << analytical << " fd(dG/dy)=" << fd;
          }
        }
      }
    }
  }
}  // namespace

// ── Global constraint: multiple grid cells ──

TEST(LinearConstraint, ResidualGlobalMultipleCells)
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

  std::size_t num_cells = 4;
  DMP state_variables{ num_cells, 3, 0.0 };
  // Cell 0: balanced (sum = 100)
  state_variables[0][0] = 50.0;  state_variables[0][1] = 30.0;  state_variables[0][2] = 20.0;
  // Cell 1: excess (sum = 110)
  state_variables[1][0] = 60.0;  state_variables[1][1] = 30.0;  state_variables[1][2] = 20.0;
  // Cell 2: deficit (sum = 80)
  state_variables[2][0] = 40.0;  state_variables[2][1] = 25.0;  state_variables[2][2] = 15.0;
  // Cell 3: zero concentrations (sum = 0)
  state_variables[3][0] = 0.0;   state_variables[3][1] = 0.0;   state_variables[3][2] = 0.0;

  auto residual_fn = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, param_indices, state_indices);
  DMP residual{ num_cells, 3, 0.0 };
  residual_fn(state_variables, no_params, residual);

  EXPECT_NEAR(residual[0][0], 0.0, 1e-12);
  EXPECT_NEAR(residual[1][0], 10.0, 1e-12);
  EXPECT_NEAR(residual[2][0], -20.0, 1e-12);
  EXPECT_NEAR(residual[3][0], -100.0, 1e-12);

  // Non-algebraic rows should be untouched
  for (std::size_t c = 0; c < num_cells; ++c)
  {
    EXPECT_NEAR(residual[c][1], 0.0, 1e-30);
    EXPECT_NEAR(residual[c][2], 0.0, 1e-30);
  }
}

// ── Global constraint: Jacobian multi-cell + FD ──

TEST(LinearConstraint, JacobianFDGlobalMultipleCells)
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

  std::size_t num_cells = 3;
  DMP state_variables{ num_cells, 3, 0.0 };
  state_variables[0][0] = 50.0;  state_variables[0][1] = 30.0;  state_variables[0][2] = 20.0;
  state_variables[1][0] = 10.0;  state_variables[1][1] = 80.0;  state_variables[1][2] = 5.0;
  state_variables[2][0] = 1.0;   state_variables[2][1] = 2.0;   state_variables[2][2] = 3.0;

  CheckConstraintFDJacobian(constraint, phase_prefixes, state_indices, state_variables);
}

// ── Global constraint: non-unit coefficients, Jacobian analytical + FD ──

TEST(LinearConstraint, JacobianGlobalNonUnitCoefficients)
{
  // G = 2.0*[A_g] + 0.5*[A_aq_LARGE] + 3.0*[A_aq_SMALL] - 42.0
  LinearConstraint constraint(
      gas_phase,
      A_g,
      { { gas_phase, A_g, 2.0 }, { aqueous_phase, A_aq, 0.5 }, { aqueous_phase, B_aq, 3.0 } },
      42.0);

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("MODE1");

  std::unordered_map<std::string, std::size_t> state_indices;
  state_indices["A_g"] = 0;
  state_indices["MODE1.AQUEOUS.A_aq"] = 1;
  state_indices["MODE1.AQUEOUS.B_aq"] = 2;

  auto jacobian = BuildConstraintJacobian(constraint, phase_prefixes, state_indices, 1);
  auto jac_fn = constraint.ConstraintJacobianFunction<DMP, SMP>(phase_prefixes, state_indices, jacobian);

  DMP state_variables{ 1, 3, 0.0 };
  state_variables[0][0] = 5.0;  state_variables[0][1] = 10.0;  state_variables[0][2] = 7.0;
  jac_fn(state_variables, jacobian);

  // jac -= dG/dy, so jac[0,0] = -2.0, jac[0,1] = -0.5, jac[0,2] = -3.0
  EXPECT_NEAR(jacobian[0][0][0], -2.0, 1e-12);
  EXPECT_NEAR(jacobian[0][0][1], -0.5, 1e-12);
  EXPECT_NEAR(jacobian[0][0][2], -3.0, 1e-12);

  // FD check
  CheckConstraintFDJacobian(constraint, phase_prefixes, state_indices, state_variables);
}

// ── Per-instance: multiple grid cells ──

TEST(LinearConstraint, ResidualPerInstanceMultipleCells)
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

  std::size_t num_cells = 3;
  DMP state_variables{ num_cells, 6, 0.0 };

  // Cell 0: LARGE balanced, SMALL imbalanced
  state_variables[0][0] = 1.0e-4;  state_variables[0][1] = 5.0e-5;  state_variables[0][2] = 5.0e-5;
  state_variables[0][3] = 2.0e-4;  state_variables[0][4] = 1.0e-4;  state_variables[0][5] = 5.0e-5;
  // Cell 1: both balanced
  state_variables[1][0] = 3.0e-3;  state_variables[1][1] = 1.0e-3;  state_variables[1][2] = 2.0e-3;
  state_variables[1][3] = 5.0e-4;  state_variables[1][4] = 3.0e-4;  state_variables[1][5] = 2.0e-4;
  // Cell 2: both imbalanced
  state_variables[2][0] = 1.0e-5;  state_variables[2][1] = 2.0e-5;  state_variables[2][2] = 3.0e-5;
  state_variables[2][3] = 4.0e-5;  state_variables[2][4] = 1.0e-5;  state_variables[2][5] = 1.0e-5;

  auto residual_fn = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, param_indices, state_indices);
  DMP residual{ num_cells, 6, 0.0 };
  residual_fn(state_variables, no_params, residual);

  // Cell 0: LARGE: 1e-4 - 5e-5 - 5e-5 = 0; SMALL: 2e-4 - 1e-4 - 5e-5 = 5e-5
  EXPECT_NEAR(residual[0][0], 0.0, 1e-18);
  EXPECT_NEAR(residual[0][3], 5.0e-5, 1e-18);
  // Cell 1: LARGE: 3e-3 - 1e-3 - 2e-3 = 0; SMALL: 5e-4 - 3e-4 - 2e-4 = 0
  EXPECT_NEAR(residual[1][0], 0.0, 1e-18);
  EXPECT_NEAR(residual[1][3], 0.0, 1e-18);
  // Cell 2: LARGE: 1e-5 - 2e-5 - 3e-5 = -4e-5; SMALL: 4e-5 - 1e-5 - 1e-5 = 2e-5
  EXPECT_NEAR(residual[2][0], -4.0e-5, 1e-18);
  EXPECT_NEAR(residual[2][3], 2.0e-5, 1e-18);
}

// ── Per-instance: Jacobian multi-cell + FD ──

TEST(LinearConstraint, JacobianFDPerInstanceMultipleCells)
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

  std::size_t num_cells = 3;
  DMP state_variables{ num_cells, 6, 0.0 };
  state_variables[0][0] = 1.0e-4;  state_variables[0][1] = 5.0e-5;  state_variables[0][2] = 5.0e-5;
  state_variables[0][3] = 2.0e-4;  state_variables[0][4] = 1.0e-4;  state_variables[0][5] = 5.0e-5;
  state_variables[1][0] = 3.0e-3;  state_variables[1][1] = 1.0e-3;  state_variables[1][2] = 2.0e-3;
  state_variables[1][3] = 5.0e-4;  state_variables[1][4] = 3.0e-4;  state_variables[1][5] = 2.0e-4;
  state_variables[2][0] = 1.0e-5;  state_variables[2][1] = 2.0e-5;  state_variables[2][2] = 3.0e-5;
  state_variables[2][3] = 4.0e-5;  state_variables[2][4] = 1.0e-5;  state_variables[2][5] = 1.0e-5;

  CheckConstraintFDJacobian(constraint, phase_prefixes, state_indices, state_variables);
}

// ── Per-instance: non-unit mixed coefficients + analytical Jacobian ──

TEST(LinearConstraint, JacobianPerInstanceNonUnitCoefficients)
{
  // G_i = 2.0*[H+_i] - 0.5*[A-_i] + 3.0*[B-_i] - 1.0
  LinearConstraint constraint(
      aqueous_phase,
      hp,
      { { aqueous_phase, hp, 2.0 }, { aqueous_phase, am, -0.5 }, { aqueous_phase, bm, 3.0 } },
      1.0);

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("MODE1");
  phase_prefixes["AQUEOUS"].insert("MODE2");

  std::unordered_map<std::string, std::size_t> si;
  si["MODE1.AQUEOUS.H+"] = 0;
  si["MODE1.AQUEOUS.A-"] = 1;
  si["MODE1.AQUEOUS.B-"] = 2;
  si["MODE2.AQUEOUS.H+"] = 3;
  si["MODE2.AQUEOUS.A-"] = 4;
  si["MODE2.AQUEOUS.B-"] = 5;

  // Check analytical Jacobian entries
  auto jacobian = BuildConstraintJacobian(constraint, phase_prefixes, si, 1);
  auto jac_fn = constraint.ConstraintJacobianFunction<DMP, SMP>(phase_prefixes, si, jacobian);

  DMP state_variables{ 1, 6, 0.0 };
  state_variables[0][0] = 1.0;  state_variables[0][1] = 2.0;  state_variables[0][2] = 3.0;
  state_variables[0][3] = 4.0;  state_variables[0][4] = 5.0;  state_variables[0][5] = 6.0;

  jac_fn(state_variables, jacobian);

  // MODE1 algebraic row (0): jac -= {2.0, -0.5, 3.0} => {-2.0, 0.5, -3.0}
  EXPECT_NEAR(jacobian[0][0][0], -2.0, 1e-12);
  EXPECT_NEAR(jacobian[0][0][1], 0.5, 1e-12);
  EXPECT_NEAR(jacobian[0][0][2], -3.0, 1e-12);

  // MODE2 algebraic row (3): same coefficients
  EXPECT_NEAR(jacobian[0][3][3], -2.0, 1e-12);
  EXPECT_NEAR(jacobian[0][3][4], 0.5, 1e-12);
  EXPECT_NEAR(jacobian[0][3][5], -3.0, 1e-12);

  // Cross-instance entries should not exist
  EXPECT_THROW(jacobian[0][0][3], std::exception);
  EXPECT_THROW(jacobian[0][3][0], std::exception);

  // FD check
  CheckConstraintFDJacobian(constraint, phase_prefixes, si, state_variables);
}

// ── Three instances ──

TEST(LinearConstraint, ThreeInstances)
{
  LinearConstraint constraint(
      aqueous_phase,
      hp,
      { { aqueous_phase, hp, 1.0 }, { aqueous_phase, am, -1.0 } },
      0.0);

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("A");
  phase_prefixes["AQUEOUS"].insert("B");
  phase_prefixes["AQUEOUS"].insert("C");

  std::unordered_map<std::string, std::size_t> si;
  si["A.AQUEOUS.H+"] = 0;  si["A.AQUEOUS.A-"] = 1;
  si["B.AQUEOUS.H+"] = 2;  si["B.AQUEOUS.A-"] = 3;
  si["C.AQUEOUS.H+"] = 4;  si["C.AQUEOUS.A-"] = 5;

  // Check residual
  DMP state_variables{ 2, 6, 0.0 };
  // Cell 0
  state_variables[0][0] = 1.0;  state_variables[0][1] = 0.5;
  state_variables[0][2] = 2.0;  state_variables[0][3] = 2.0;
  state_variables[0][4] = 3.0;  state_variables[0][5] = 4.0;
  // Cell 1
  state_variables[1][0] = 10.0;  state_variables[1][1] = 10.0;
  state_variables[1][2] = 20.0;  state_variables[1][3] = 15.0;
  state_variables[1][4] = 30.0;  state_variables[1][5] = 25.0;

  auto residual_fn = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, param_indices, si);
  DMP residual{ 2, 6, 0.0 };
  residual_fn(state_variables, no_params, residual);

  // Cell 0: A: 1.0 - 0.5 = 0.5; B: 2.0 - 2.0 = 0.0; C: 3.0 - 4.0 = -1.0
  EXPECT_NEAR(residual[0][0], 0.5, 1e-12);
  EXPECT_NEAR(residual[0][2], 0.0, 1e-12);
  EXPECT_NEAR(residual[0][4], -1.0, 1e-12);

  // Cell 1: A: 10 - 10 = 0; B: 20 - 15 = 5; C: 30 - 25 = 5
  EXPECT_NEAR(residual[1][0], 0.0, 1e-12);
  EXPECT_NEAR(residual[1][2], 5.0, 1e-12);
  EXPECT_NEAR(residual[1][4], 5.0, 1e-12);

  // FD Jacobian check
  CheckConstraintFDJacobian(constraint, phase_prefixes, si, state_variables);
}

// ── Global constraint with many instances ──

TEST(LinearConstraint, GlobalWithManyInstances)
{
  // Gas-phase algebraic variable = total of gas + all aq instances
  LinearConstraint constraint(
      gas_phase,
      A_g,
      { { gas_phase, A_g, 1.0 }, { aqueous_phase, A_aq, 1.0 } },
      50.0);

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("M1");
  phase_prefixes["AQUEOUS"].insert("M2");
  phase_prefixes["AQUEOUS"].insert("M3");
  phase_prefixes["AQUEOUS"].insert("M4");

  std::unordered_map<std::string, std::size_t> si;
  si["A_g"] = 0;
  si["M1.AQUEOUS.A_aq"] = 1;
  si["M2.AQUEOUS.A_aq"] = 2;
  si["M3.AQUEOUS.A_aq"] = 3;
  si["M4.AQUEOUS.A_aq"] = 4;

  DMP sv{ 1, 5, 0.0 };
  sv[0][0] = 10.0;  sv[0][1] = 10.0;  sv[0][2] = 10.0;  sv[0][3] = 10.0;  sv[0][4] = 10.0;

  auto residual_fn = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, param_indices, si);
  DMP residual{ 1, 5, 0.0 };
  residual_fn(sv, no_params, residual);
  // G = 10+10+10+10+10 - 50 = 0
  EXPECT_NEAR(residual[0][0], 0.0, 1e-12);

  // Jacobian: all columns in row 0 should have -1
  auto jacobian = BuildConstraintJacobian(constraint, phase_prefixes, si, 1);
  auto jac_fn = constraint.ConstraintJacobianFunction<DMP, SMP>(phase_prefixes, si, jacobian);
  jac_fn(sv, jacobian);
  for (std::size_t j = 0; j < 5; ++j)
    EXPECT_NEAR(jacobian[0][0][j], -1.0, 1e-12);

  // FD check
  CheckConstraintFDJacobian(constraint, phase_prefixes, si, sv);
}

// ── Jacobian accumulation ──

TEST(LinearConstraint, JacobianAccumulates)
{
  LinearConstraint constraint(
      gas_phase,
      A_g,
      { { gas_phase, A_g, 2.0 }, { aqueous_phase, A_aq, 0.5 } },
      10.0);

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("MODE1");

  std::unordered_map<std::string, std::size_t> si;
  si["A_g"] = 0;
  si["MODE1.AQUEOUS.A_aq"] = 1;

  auto jacobian = BuildConstraintJacobian(constraint, phase_prefixes, si, 1);
  auto jac_fn = constraint.ConstraintJacobianFunction<DMP, SMP>(phase_prefixes, si, jacobian);

  DMP sv{ 1, 2, 0.0 };
  sv[0][0] = 1.0;  sv[0][1] = 2.0;

  // First call
  jac_fn(sv, jacobian);
  double j00_once = jacobian[0][0][0];
  double j01_once = jacobian[0][0][1];

  EXPECT_NEAR(j00_once, -2.0, 1e-12);
  EXPECT_NEAR(j01_once, -0.5, 1e-12);

  // Second call accumulates
  jac_fn(sv, jacobian);
  EXPECT_NEAR(jacobian[0][0][0], 2.0 * j00_once, 1e-12);
  EXPECT_NEAR(jacobian[0][0][1], 2.0 * j01_once, 1e-12);
}

// ── Residual sets (does not accumulate) ──

TEST(LinearConstraint, ResidualSetsNotAccumulates)
{
  LinearConstraint constraint(
      gas_phase,
      A_g,
      { { gas_phase, A_g, 1.0 } },
      10.0);

  std::map<std::string, std::set<std::string>> phase_prefixes;

  std::unordered_map<std::string, std::size_t> si;
  si["A_g"] = 0;

  auto residual_fn = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, param_indices, si);

  DMP sv{ 1, 1, 0.0 };
  sv[0][0] = 7.0;

  DMP residual{ 1, 1, 999.0 };  // pre-filled with junk
  residual_fn(sv, no_params, residual);
  EXPECT_NEAR(residual[0][0], -3.0, 1e-12);  // 7 - 10 = -3

  // Second call still gives same result (assignment, not accumulation)
  residual_fn(sv, no_params, residual);
  EXPECT_NEAR(residual[0][0], -3.0, 1e-12);
}

// ── Multiple constraints on same variable set ──

TEST(LinearConstraint, MultipleConstraintsCombined)
{
  // Constraint 1: mass conservation [A_g] + [A_aq] = 100
  LinearConstraint mass_cons(
      gas_phase,
      A_g,
      { { gas_phase, A_g, 1.0 }, { aqueous_phase, A_aq, 1.0 } },
      100.0);

  // Constraint 2: charge balance [H+] - [A-] = 0, per instance
  LinearConstraint charge_bal(
      aqueous_phase,
      hp,
      { { aqueous_phase, hp, 1.0 }, { aqueous_phase, am, -1.0 } },
      0.0);

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("MODE1");

  std::unordered_map<std::string, std::size_t> si;
  si["A_g"] = 0;
  si["MODE1.AQUEOUS.A_aq"] = 1;
  si["MODE1.AQUEOUS.H+"] = 2;
  si["MODE1.AQUEOUS.A-"] = 3;

  // Residuals from both constraints
  auto rf1 = mass_cons.ConstraintResidualFunction<DMP>(phase_prefixes, param_indices, si);
  auto rf2 = charge_bal.ConstraintResidualFunction<DMP>(phase_prefixes, param_indices, si);

  DMP sv{ 1, 4, 0.0 };
  sv[0][0] = 60.0;   // A_g
  sv[0][1] = 40.0;   // A_aq
  sv[0][2] = 1.0e-3; // H+
  sv[0][3] = 1.0e-3; // A-

  DMP residual{ 1, 4, 0.0 };
  rf1(sv, no_params, residual);
  rf2(sv, no_params, residual);

  // Mass: 60 + 40 - 100 = 0 → residual at A_g (idx 0)
  EXPECT_NEAR(residual[0][0], 0.0, 1e-12);
  // Charge: 1e-3 - 1e-3 = 0 → residual at H+ (idx 2)
  EXPECT_NEAR(residual[0][2], 0.0, 1e-15);

  // Combined Jacobian
  auto jacobian = BuildConstraintJacobian(
      { std::cref(mass_cons), std::cref(charge_bal) }, phase_prefixes, si, 1);

  auto jf1 = mass_cons.ConstraintJacobianFunction<DMP, SMP>(phase_prefixes, si, jacobian);
  auto jf2 = charge_bal.ConstraintJacobianFunction<DMP, SMP>(phase_prefixes, si, jacobian);
  jf1(sv, jacobian);
  jf2(sv, jacobian);

  // Mass row (0): dG/d[A_g] = 1, dG/d[A_aq] = 1 → jac = -1, -1
  EXPECT_NEAR(jacobian[0][0][0], -1.0, 1e-12);
  EXPECT_NEAR(jacobian[0][0][1], -1.0, 1e-12);
  // Charge row (2): dG/d[H+] = 1, dG/d[A-] = -1 → jac = -1, +1
  EXPECT_NEAR(jacobian[0][2][2], -1.0, 1e-12);
  EXPECT_NEAR(jacobian[0][2][3], 1.0, 1e-12);
}

// ── Mixed: global algebraic + instanced terms, multi-cell, FD ──

TEST(LinearConstraint, JacobianFDGlobalAlgWithInstancedTerms)
{
  // G = [A_g] + [A_aq_M1] + [A_aq_M2] + [A_aq_M3] - 200.0
  LinearConstraint constraint(
      gas_phase,
      A_g,
      { { gas_phase, A_g, 1.0 }, { aqueous_phase, A_aq, 1.0 } },
      200.0);

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("M1");
  phase_prefixes["AQUEOUS"].insert("M2");
  phase_prefixes["AQUEOUS"].insert("M3");

  std::unordered_map<std::string, std::size_t> si;
  si["A_g"] = 0;
  si["M1.AQUEOUS.A_aq"] = 1;
  si["M2.AQUEOUS.A_aq"] = 2;
  si["M3.AQUEOUS.A_aq"] = 3;

  std::size_t nc = 3;
  DMP sv{ nc, 4, 0.0 };
  for (std::size_t c = 0; c < nc; ++c)
  {
    sv[c][0] = 50.0 + 10.0 * c;
    sv[c][1] = 40.0 - 5.0 * c;
    sv[c][2] = 30.0 + 3.0 * c;
    sv[c][3] = 20.0 + 2.0 * c;
  }

  CheckConstraintFDJacobian(constraint, phase_prefixes, si, sv);
}

// ── Per-instance with non-instanced (gas) terms ──

TEST(LinearConstraint, PerInstanceWithGasTerms)
{
  // G_i = [A_aq_i] + [A_g] - 10.0   (algebraic in aqueous, depends on gas)
  LinearConstraint constraint(
      aqueous_phase,
      A_aq,
      { { aqueous_phase, A_aq, 1.0 }, { gas_phase, A_g, 1.0 } },
      10.0);

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("SMALL");
  phase_prefixes["AQUEOUS"].insert("LARGE");

  std::unordered_map<std::string, std::size_t> si;
  si["A_g"] = 0;
  si["LARGE.AQUEOUS.A_aq"] = 1;
  si["SMALL.AQUEOUS.A_aq"] = 2;

  // Residual
  DMP sv{ 1, 3, 0.0 };
  sv[0][0] = 3.0;  // gas
  sv[0][1] = 4.0;  // LARGE aq
  sv[0][2] = 8.0;  // SMALL aq

  auto rf = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, param_indices, si);
  DMP residual{ 1, 3, 0.0 };
  rf(sv, no_params, residual);

  // LARGE: 4 + 3 - 10 = -3
  EXPECT_NEAR(residual[0][1], -3.0, 1e-12);
  // SMALL: 8 + 3 - 10 = 1
  EXPECT_NEAR(residual[0][2], 1.0, 1e-12);

  // Jacobian: each instance row has entries for (A_g, own A_aq)
  auto jacobian = BuildConstraintJacobian(constraint, phase_prefixes, si, 1);
  auto jf = constraint.ConstraintJacobianFunction<DMP, SMP>(phase_prefixes, si, jacobian);
  jf(sv, jacobian);

  // LARGE row (1): dG/d[A_g] = 1, dG/d[LARGE.A_aq] = 1
  EXPECT_NEAR(jacobian[0][1][0], -1.0, 1e-12);  // jac -= 1
  EXPECT_NEAR(jacobian[0][1][1], -1.0, 1e-12);
  // SMALL row (2): dG/d[A_g] = 1, dG/d[SMALL.A_aq] = 1
  EXPECT_NEAR(jacobian[0][2][0], -1.0, 1e-12);
  EXPECT_NEAR(jacobian[0][2][2], -1.0, 1e-12);

  // Cross-instance isolation: LARGE doesn't depend on SMALL
  EXPECT_THROW(jacobian[0][1][2], std::exception);
  EXPECT_THROW(jacobian[0][2][1], std::exception);

  // FD check multi-cell
  std::size_t nc = 3;
  DMP sv_mc{ nc, 3, 0.0 };
  for (std::size_t c = 0; c < nc; ++c)
  {
    sv_mc[c][0] = 3.0 + c;
    sv_mc[c][1] = 4.0 + 2.0 * c;
    sv_mc[c][2] = 8.0 - c;
  }
  CheckConstraintFDJacobian(constraint, phase_prefixes, si, sv_mc);
}

// ── Kitchen-sink: multi-constraint, multi-instance, multi-cell, FD ──

TEST(LinearConstraint, KitchenSinkFDCheck)
{
  // Mass: [A_g] + [A_aq_M1] + [A_aq_M2] = 100
  LinearConstraint mass_cons(
      gas_phase,
      A_g,
      { { gas_phase, A_g, 1.0 }, { aqueous_phase, A_aq, 1.0 } },
      100.0);

  // Charge: [H+_i] - 0.5*[A-_i] - 2.0*[B-_i] = 0
  LinearConstraint charge_bal(
      aqueous_phase,
      hp,
      { { aqueous_phase, hp, 1.0 }, { aqueous_phase, am, -0.5 }, { aqueous_phase, bm, -2.0 } },
      0.0);

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("M1");
  phase_prefixes["AQUEOUS"].insert("M2");

  // Variables: A_g, M1.{A_aq, H+, A-, B-}, M2.{A_aq, H+, A-, B-}
  std::unordered_map<std::string, std::size_t> si;
  si["A_g"] = 0;
  si["M1.AQUEOUS.A_aq"] = 1;
  si["M1.AQUEOUS.H+"] = 2;
  si["M1.AQUEOUS.A-"] = 3;
  si["M1.AQUEOUS.B-"] = 4;
  si["M2.AQUEOUS.A_aq"] = 5;
  si["M2.AQUEOUS.H+"] = 6;
  si["M2.AQUEOUS.A-"] = 7;
  si["M2.AQUEOUS.B-"] = 8;

  std::size_t nc = 2;
  DMP sv{ nc, 9, 0.0 };
  sv[0][0] = 40.0;
  sv[0][1] = 30.0;  sv[0][2] = 1.0e-3;  sv[0][3] = 1.5e-3;  sv[0][4] = 0.5e-4;
  sv[0][5] = 30.0;  sv[0][6] = 2.0e-3;  sv[0][7] = 3.0e-3;  sv[0][8] = 2.5e-4;
  sv[1][0] = 50.0;
  sv[1][1] = 25.0;  sv[1][2] = 5.0e-4;  sv[1][3] = 8.0e-4;  sv[1][4] = 1.0e-4;
  sv[1][5] = 25.0;  sv[1][6] = 3.0e-4;  sv[1][7] = 5.0e-4;  sv[1][8] = 4.0e-5;

  // FD check each constraint individually
  CheckConstraintFDJacobian(mass_cons, phase_prefixes, si, sv);
  CheckConstraintFDJacobian(charge_bal, phase_prefixes, si, sv);
}

// ── Single term (trivial) ──

TEST(LinearConstraint, SingleTerm)
{
  // G = 3.0 * [A_g] - 15.0, algebraic [A_g]
  LinearConstraint constraint(gas_phase, A_g, { { gas_phase, A_g, 3.0 } }, 15.0);

  std::map<std::string, std::set<std::string>> phase_prefixes;
  std::unordered_map<std::string, std::size_t> si;
  si["A_g"] = 0;

  DMP sv{ 1, 1, 0.0 };
  sv[0][0] = 5.0;

  auto rf = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, param_indices, si);
  DMP residual{ 1, 1, 0.0 };
  rf(sv, no_params, residual);
  // 3*5 - 15 = 0
  EXPECT_NEAR(residual[0][0], 0.0, 1e-12);

  sv[0][0] = 4.0;
  rf(sv, no_params, residual);
  // 3*4 - 15 = -3
  EXPECT_NEAR(residual[0][0], -3.0, 1e-12);

  // Jacobian: -3.0
  auto jacobian = BuildConstraintJacobian(constraint, phase_prefixes, si, 1);
  auto jf = constraint.ConstraintJacobianFunction<DMP, SMP>(phase_prefixes, si, jacobian);
  jf(sv, jacobian);
  EXPECT_NEAR(jacobian[0][0][0], -3.0, 1e-12);

  CheckConstraintFDJacobian(constraint, phase_prefixes, si, sv);
}

// ── Large coefficient range ──

TEST(LinearConstraint, LargeCoefficientRange)
{
  // G = 1e6 * [A_g] + 1e-6 * [A_aq] - 500
  LinearConstraint constraint(
      gas_phase,
      A_g,
      { { gas_phase, A_g, 1.0e6 }, { aqueous_phase, A_aq, 1.0e-6 } },
      500.0);

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  std::unordered_map<std::string, std::size_t> si;
  si["A_g"] = 0;
  si["DROP.AQUEOUS.A_aq"] = 1;

  DMP sv{ 1, 2, 0.0 };
  sv[0][0] = 5.0e-4;   // 1e6 * 5e-4 = 500
  sv[0][1] = 1.0e6;    // 1e-6 * 1e6 = 1

  auto rf = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, param_indices, si);
  DMP residual{ 1, 2, 0.0 };
  rf(sv, no_params, residual);
  // 500 + 1 - 500 = 1
  EXPECT_NEAR(residual[0][0], 1.0, 1e-9);

  CheckConstraintFDJacobian(constraint, phase_prefixes, si, sv);
}

// ── CopyWithNewUuid preserves behavior ──

TEST(LinearConstraint, CopiedConstraintProducesSameResults)
{
  LinearConstraint original(
      gas_phase,
      A_g,
      { { gas_phase, A_g, 2.0 }, { aqueous_phase, A_aq, 0.5 } },
      10.0);

  auto copy = original.CopyWithNewUuid();
  EXPECT_NE(copy.uuid_, original.uuid_);

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("MODE1");

  std::unordered_map<std::string, std::size_t> si;
  si["A_g"] = 0;
  si["MODE1.AQUEOUS.A_aq"] = 1;

  DMP sv{ 1, 2, 0.0 };
  sv[0][0] = 3.0;  sv[0][1] = 8.0;

  auto rf_orig = original.ConstraintResidualFunction<DMP>(phase_prefixes, param_indices, si);
  auto rf_copy = copy.ConstraintResidualFunction<DMP>(phase_prefixes, param_indices, si);

  DMP res_orig{ 1, 2, 0.0 };
  DMP res_copy{ 1, 2, 0.0 };
  rf_orig(sv, no_params, res_orig);
  rf_copy(sv, no_params, res_copy);

  EXPECT_NEAR(res_orig[0][0], res_copy[0][0], 1e-15);
}
