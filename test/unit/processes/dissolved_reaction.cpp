// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <miam/processes/dissolved_reaction.hpp>

#include <micm/system/conditions.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
#include <micm/util/jacobian_verification.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <unordered_map>

using namespace miam;

namespace
{
  using MatrixPolicy = micm::Matrix<double>;
  using SparseMatrixPolicy = micm::SparseMatrix<double, micm::SparseMatrixStandardOrderingCompressedSparseRow>;

  /// @brief Compare analytical Jacobian against central finite-difference approximation
  ///        using MICM's FiniteDifferenceJacobian / CompareJacobianToFiniteDifference utilities.
  void CheckFiniteDifferenceJacobian(
      const DissolvedReaction& reaction,
      const std::map<std::string, std::set<std::string>>& phase_prefixes,
      const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
      const std::unordered_map<std::string, std::size_t>& state_variable_indices,
      const MatrixPolicy& state_parameters,
      const MatrixPolicy& state_variables)
  {
    const std::size_t num_blocks = state_parameters.NumRows();
    const std::size_t num_vars = state_variable_indices.size();

    // Build sparse Jacobian structure and compute analytical Jacobian
    auto jac_elements = reaction.NonZeroJacobianElements(phase_prefixes, state_variable_indices);
    auto jac_builder = SparseMatrixPolicy::Create(num_vars).SetNumberOfBlocks(num_blocks);
    for (const auto& elem : jac_elements)
      jac_builder.WithElement(elem.first, elem.second);
    SparseMatrixPolicy jacobian(jac_builder);
    jacobian.Fill(0.0);

    auto jac_func = reaction.JacobianFunction<MatrixPolicy, SparseMatrixPolicy>(
        phase_prefixes, state_parameter_indices, state_variable_indices, jacobian);
    jac_func(state_parameters, state_variables, jacobian);

    // Build FD Jacobian — bind state_parameters into the forcing callable
    auto ff = reaction.ForcingFunction<MatrixPolicy>(phase_prefixes, state_parameter_indices, state_variable_indices);
    auto fd_jac = micm::FiniteDifferenceJacobian<MatrixPolicy>(
        [&](const MatrixPolicy& vars, MatrixPolicy& out)
        {
          out.Fill(0.0);
          ff(state_parameters, vars, out);
        },
        state_variables,
        num_vars);

    // Compare analytical vs FD (MICM defaults: atol=1e-7, rtol=1e-7)
    auto cmp = micm::CompareJacobianToFiniteDifference(jacobian, fd_jac, num_vars);
    EXPECT_TRUE(cmp.passed_) << "FD mismatch: block=" << cmp.worst_block_ << " row=" << cmp.worst_row_
                             << " col=" << cmp.worst_col_ << " +J(analytical)=" << cmp.worst_analytical_
                             << " +J(fd)=" << cmp.worst_fd_;

    // Verify no significant FD signal outside the declared sparsity pattern
    auto spc = micm::CheckJacobianSparsityCompleteness(jacobian, fd_jac, num_vars);
    EXPECT_TRUE(spc.passed_) << "Missing sparsity entry: block=" << spc.worst_block_ << " row=" << spc.worst_row_
                             << " col=" << spc.worst_col_ << " fd=" << spc.worst_fd_;
  }
}  // namespace

// ============================================================================
// SpeciesUsed Tests
// ============================================================================

TEST(DissolvedReaction, SpeciesUsedWithSinglePrefix)
{
  auto a = micm::Species{ "A" };
  auto b = micm::Species{ "B" };
  auto solvent = micm::Species{ "S" };

  auto phase = micm::Phase{ "AQUEOUS", { { a }, { b }, { solvent } } };

  auto rate = [](const micm::Conditions& conditions) { return 0.1; };

  DissolvedReaction reaction{ { { "SMALL_DROP", rate } },
                              { a },    // reactants
                              { b },    // products
                              solvent,  // solvent
                              phase };

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("SMALL_DROP");

  auto species_used = reaction.SpeciesUsed(phase_prefixes);

  // Should have 3 species: A (reactant), B (product), S (solvent)
  EXPECT_EQ(species_used.size(), 3);
  EXPECT_TRUE(species_used.find("SMALL_DROP.AQUEOUS.A") != species_used.end());
  EXPECT_TRUE(species_used.find("SMALL_DROP.AQUEOUS.B") != species_used.end());
  EXPECT_TRUE(species_used.find("SMALL_DROP.AQUEOUS.S") != species_used.end());
}

TEST(DissolvedReaction, SpeciesUsedWithMultiplePrefixes)
{
  auto a = micm::Species{ "A" };
  auto b = micm::Species{ "B" };
  auto c = micm::Species{ "C" };
  auto solvent = micm::Species{ "S" };

  auto phase = micm::Phase{ "AQUEOUS", { { a }, { b }, { c }, { solvent } } };

  auto rate = [](const micm::Conditions& conditions) { return 0.1; };

  DissolvedReaction reaction{ { { "SMALL_DROP", rate }, { "LARGE_DROP", rate } },
                              { a, b },  // reactants
                              { c },     // products
                              solvent,   // solvent
                              phase };

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("SMALL_DROP");
  phase_prefixes["AQUEOUS"].insert("LARGE_DROP");

  auto species_used = reaction.SpeciesUsed(phase_prefixes);

  // 4 species (A, B, C, S) × 2 representations = 8 total
  EXPECT_EQ(species_used.size(), 8);

  for (const auto& prefix : { "SMALL_DROP", "LARGE_DROP" })
  {
    EXPECT_TRUE(species_used.find(std::string(prefix) + ".AQUEOUS.A") != species_used.end());
    EXPECT_TRUE(species_used.find(std::string(prefix) + ".AQUEOUS.B") != species_used.end());
    EXPECT_TRUE(species_used.find(std::string(prefix) + ".AQUEOUS.C") != species_used.end());
    EXPECT_TRUE(species_used.find(std::string(prefix) + ".AQUEOUS.S") != species_used.end());
  }
}

TEST(DissolvedReaction, SpeciesUsedWithNoMatchingPhase)
{
  auto a = micm::Species{ "A" };
  auto b = micm::Species{ "B" };
  auto solvent = micm::Species{ "S" };

  auto phase = micm::Phase{ "AQUEOUS", { { a }, { b }, { solvent } } };

  auto rate = [](const micm::Conditions& conditions) { return 0.1; };

  DissolvedReaction reaction{ { { "DROP", rate } }, { a }, { b }, solvent, phase };

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["ORGANIC"].insert("AEROSOL_MODE");

  EXPECT_THROW(reaction.SpeciesUsed(phase_prefixes), std::runtime_error);
}

TEST(DissolvedReaction, SpeciesUsedDuplicateHandling)
{
  // Test case where solvent is the same as a reactant
  auto a = micm::Species{ "A" };
  auto b = micm::Species{ "B" };

  auto phase = micm::Phase{ "AQUEOUS", { { a }, { b } } };

  auto rate = [](const micm::Conditions& conditions) { return 0.1; };

  DissolvedReaction reaction{ { { "MODE1", rate } },
                              { a },  // reactants
                              { b },  // products
                              a,      // solvent (same as reactant)
                              phase };

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("MODE1");

  auto species_used = reaction.SpeciesUsed(phase_prefixes);

  // Should have 2 unique species (A is both reactant and solvent, no duplicate)
  EXPECT_EQ(species_used.size(), 2);
  EXPECT_TRUE(species_used.find("MODE1.AQUEOUS.A") != species_used.end());
  EXPECT_TRUE(species_used.find("MODE1.AQUEOUS.B") != species_used.end());
}

// ============================================================================
// NonZeroJacobianElements Tests
// ============================================================================

TEST(DissolvedReaction, NonZeroJacobianElementsBasic)
{
  // A -> B + C with solvent S
  // Independent: A (reactant), S (solvent) = 2
  // Dependent: A, B, C = 3
  // Total: 3 × 2 = 6 elements
  auto a = micm::Species{ "A" };
  auto b = micm::Species{ "B" };
  auto c = micm::Species{ "C" };
  auto solvent = micm::Species{ "S" };

  auto phase = micm::Phase{ "AQUEOUS", { { a }, { b }, { c }, { solvent } } };

  auto rate = [](const micm::Conditions& conditions) { return 0.1; };

  DissolvedReaction reaction{ { { "MODE1", rate } },
                              { a },     // 1 reactant
                              { b, c },  // 2 products
                              solvent,
                              phase };

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("MODE1");

  std::unordered_map<std::string, std::size_t> state_indices;
  state_indices["MODE1.AQUEOUS.A"] = 0;
  state_indices["MODE1.AQUEOUS.B"] = 1;
  state_indices["MODE1.AQUEOUS.C"] = 2;
  state_indices["MODE1.AQUEOUS.S"] = 3;

  auto jacobian_elements = reaction.NonZeroJacobianElements(phase_prefixes, state_indices);

  // 3 dependent (A, B, C) × 2 independent (A, S) = 6 elements
  EXPECT_EQ(jacobian_elements.size(), 6);

  // Verify specific elements exist (dependent, independent)
  EXPECT_TRUE(jacobian_elements.find({ 0, 0 }) != jacobian_elements.end());  // d[A]/d[A]
  EXPECT_TRUE(jacobian_elements.find({ 1, 0 }) != jacobian_elements.end());  // d[B]/d[A]
  EXPECT_TRUE(jacobian_elements.find({ 2, 0 }) != jacobian_elements.end());  // d[C]/d[A]
  EXPECT_TRUE(jacobian_elements.find({ 0, 3 }) != jacobian_elements.end());  // d[A]/d[S]
  EXPECT_TRUE(jacobian_elements.find({ 1, 3 }) != jacobian_elements.end());  // d[B]/d[S]
  EXPECT_TRUE(jacobian_elements.find({ 2, 3 }) != jacobian_elements.end());  // d[C]/d[S]

  // Verify NO entries where products are independent variables
  EXPECT_TRUE(jacobian_elements.find({ 0, 1 }) == jacobian_elements.end());  // d[A]/d[B] should NOT exist
  EXPECT_TRUE(jacobian_elements.find({ 0, 2 }) == jacobian_elements.end());  // d[A]/d[C] should NOT exist
  EXPECT_TRUE(jacobian_elements.find({ 1, 1 }) == jacobian_elements.end());  // d[B]/d[B] should NOT exist
  EXPECT_TRUE(jacobian_elements.find({ 1, 2 }) == jacobian_elements.end());  // d[B]/d[C] should NOT exist
  EXPECT_TRUE(jacobian_elements.find({ 2, 1 }) == jacobian_elements.end());  // d[C]/d[B] should NOT exist
  EXPECT_TRUE(jacobian_elements.find({ 2, 2 }) == jacobian_elements.end());  // d[C]/d[C] should NOT exist
}

TEST(DissolvedReaction, NonZeroJacobianElementsMultipleReactants)
{
  // A + B -> C with solvent S
  // Independent: A, B (reactants), S (solvent) = 3
  // Dependent: A, B, C = 3
  // Total: 3 × 3 = 9 elements
  auto a = micm::Species{ "A" };
  auto b = micm::Species{ "B" };
  auto c = micm::Species{ "C" };
  auto solvent = micm::Species{ "S" };

  auto phase = micm::Phase{ "AQUEOUS", { { a }, { b }, { c }, { solvent } } };

  auto rate = [](const micm::Conditions& conditions) { return 0.1; };

  DissolvedReaction reaction{ { { "DROP", rate } },
                              { a, b },  // 2 reactants
                              { c },     // 1 product
                              solvent,
                              phase };

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  std::unordered_map<std::string, std::size_t> state_indices;
  state_indices["DROP.AQUEOUS.A"] = 0;
  state_indices["DROP.AQUEOUS.B"] = 1;
  state_indices["DROP.AQUEOUS.C"] = 2;
  state_indices["DROP.AQUEOUS.S"] = 3;

  auto jacobian_elements = reaction.NonZeroJacobianElements(phase_prefixes, state_indices);

  // 3 dependent (A, B, C) × 3 independent (A, B, S) = 9 elements
  EXPECT_EQ(jacobian_elements.size(), 9);

  // Verify NO entries where product C is independent
  EXPECT_TRUE(jacobian_elements.find({ 0, 2 }) == jacobian_elements.end());  // d[A]/d[C]
  EXPECT_TRUE(jacobian_elements.find({ 1, 2 }) == jacobian_elements.end());  // d[B]/d[C]
  EXPECT_TRUE(jacobian_elements.find({ 2, 2 }) == jacobian_elements.end());  // d[C]/d[C]
}

TEST(DissolvedReaction, NonZeroJacobianElementsMultiplePrefixes)
{
  auto a = micm::Species{ "A" };
  auto b = micm::Species{ "B" };
  auto solvent = micm::Species{ "S" };

  auto phase = micm::Phase{ "AQUEOUS", { { a }, { b }, { solvent } } };

  auto rate = [](const micm::Conditions& conditions) { return 0.1; };

  // A -> B with solvent S
  DissolvedReaction reaction{ { { "SMALL_DROP", rate }, { "LARGE_DROP", rate } }, { a }, { b }, solvent, phase };

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("SMALL_DROP");
  phase_prefixes["AQUEOUS"].insert("LARGE_DROP");

  std::unordered_map<std::string, std::size_t> state_indices;
  state_indices["SMALL_DROP.AQUEOUS.A"] = 0;
  state_indices["SMALL_DROP.AQUEOUS.B"] = 1;
  state_indices["SMALL_DROP.AQUEOUS.S"] = 2;
  state_indices["LARGE_DROP.AQUEOUS.A"] = 3;
  state_indices["LARGE_DROP.AQUEOUS.B"] = 4;
  state_indices["LARGE_DROP.AQUEOUS.S"] = 5;

  auto jacobian_elements = reaction.NonZeroJacobianElements(phase_prefixes, state_indices);

  // A -> B: 2 dependent (A, B) × 2 independent (A, S) = 4 per prefix × 2 prefixes = 8
  EXPECT_EQ(jacobian_elements.size(), 8);
}

// ============================================================================
// UpdateStateParametersFunction Tests
// ============================================================================

TEST(DissolvedReaction, UpdateStateParametersFunctionBasic)
{
  using MatrixPolicy = micm::VectorMatrix<double>;

  auto a = micm::Species{ "A" };
  auto b = micm::Species{ "B" };
  auto solvent = micm::Species{ "S" };

  auto phase = micm::Phase{ "AQUEOUS", { { a }, { b }, { solvent } } };

  double k = 0.1;
  auto rate = [k](const micm::Conditions& conditions) { return k; };

  DissolvedReaction reaction{ { { "MODE1", rate } }, { a }, { b }, solvent, phase };

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("MODE1");

  std::string k_param = "MODE1." + phase.name_ + "." + reaction.uuid_ + ".k";

  std::unordered_map<std::string, std::size_t> state_parameter_indices;
  state_parameter_indices[k_param] = 0;

  auto update_func = reaction.UpdateStateParametersFunction<MatrixPolicy>(phase_prefixes, state_parameter_indices);

  MatrixPolicy state_parameters(3, 1, 0.0);

  std::vector<micm::Conditions> conditions(3);
  conditions[0].temperature_ = 298.15;
  conditions[1].temperature_ = 310.0;
  conditions[2].temperature_ = 285.0;

  update_func(conditions, state_parameters);

  for (std::size_t i_cell = 0; i_cell < 3; ++i_cell)
  {
    EXPECT_EQ(state_parameters[i_cell][0], k);
  }
}

TEST(DissolvedReaction, UpdateStateParametersFunctionTemperatureDependent)
{
  using MatrixPolicy = micm::VectorMatrix<double>;

  auto a = micm::Species{ "A" };
  auto b = micm::Species{ "B" };
  auto solvent = micm::Species{ "S" };

  auto phase = micm::Phase{ "AQUEOUS", { { a }, { b }, { solvent } } };

  auto rate = [](const micm::Conditions& conditions) { return 0.1 * std::exp(-3000.0 / conditions.temperature_); };

  DissolvedReaction reaction{ { { "DROPLET", rate } }, { a }, { b }, solvent, phase };

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROPLET");

  std::string k_param = "DROPLET." + phase.name_ + "." + reaction.uuid_ + ".k";

  std::unordered_map<std::string, std::size_t> state_parameter_indices;
  state_parameter_indices[k_param] = 0;

  auto update_func = reaction.UpdateStateParametersFunction<MatrixPolicy>(phase_prefixes, state_parameter_indices);

  MatrixPolicy state_parameters(3, 1, 0.0);

  std::vector<micm::Conditions> conditions(3);
  conditions[0].temperature_ = 298.15;
  conditions[1].temperature_ = 310.0;
  conditions[2].temperature_ = 273.15;

  update_func(conditions, state_parameters);

  double k_298 = 0.1 * std::exp(-3000.0 / 298.15);
  double k_310 = 0.1 * std::exp(-3000.0 / 310.0);
  double k_273 = 0.1 * std::exp(-3000.0 / 273.15);

  EXPECT_NEAR(state_parameters[0][0], k_298, 1e-15);
  EXPECT_NEAR(state_parameters[1][0], k_310, 1e-15);
  EXPECT_NEAR(state_parameters[2][0], k_273, 1e-15);
}

TEST(DissolvedReaction, UpdateStateParametersFunctionMissingParameter)
{
  using MatrixPolicy = micm::VectorMatrix<double>;

  auto a = micm::Species{ "A" };
  auto b = micm::Species{ "B" };
  auto solvent = micm::Species{ "S" };

  auto phase = micm::Phase{ "AQUEOUS", { { a }, { b }, { solvent } } };

  auto rate = [](const micm::Conditions& conditions) { return 0.1; };

  DissolvedReaction reaction{ { { "MODE1", rate } }, { a }, { b }, solvent, phase };

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("MODE1");

  // Empty parameter indices - should throw
  std::unordered_map<std::string, std::size_t> state_parameter_indices;

  EXPECT_THROW(
      reaction.UpdateStateParametersFunction<MatrixPolicy>(phase_prefixes, state_parameter_indices), std::runtime_error);
}

// ============================================================================
// ForcingFunction Tests
// ============================================================================

TEST(DissolvedReaction, ForcingFunctionBasicRates)
{
  using MatrixPolicy = micm::VectorMatrix<double>;

  auto a = micm::Species{ "A" };
  auto b = micm::Species{ "B" };
  auto solvent = micm::Species{ "S" };

  auto phase = micm::Phase{ "AQUEOUS", { { a }, { b }, { solvent } } };

  double k = 0.1;
  auto rate = [k](const micm::Conditions& conditions) { return k; };

  // A -> B with solvent S
  DissolvedReaction reaction{ { { "MODE1", rate } }, { a }, { b }, solvent, phase };

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("MODE1");

  std::string k_param = "MODE1." + phase.name_ + "." + reaction.uuid_ + ".k";

  std::unordered_map<std::string, std::size_t> state_parameter_indices;
  state_parameter_indices[k_param] = 0;

  std::unordered_map<std::string, std::size_t> state_variable_indices;
  state_variable_indices["MODE1.AQUEOUS.A"] = 0;
  state_variable_indices["MODE1.AQUEOUS.B"] = 1;
  state_variable_indices["MODE1.AQUEOUS.S"] = 2;

  auto forcing_func =
      reaction.ForcingFunction<MatrixPolicy>(phase_prefixes, state_parameter_indices, state_variable_indices);

  MatrixPolicy state_parameters(1, 1);
  state_parameters[0][0] = k;

  MatrixPolicy state_variables(1, 3);
  state_variables[0][0] = 2.0;   // [A]
  state_variables[0][1] = 0.5;   // [B]
  state_variables[0][2] = 55.0;  // [S] (solvent)

  MatrixPolicy forcing_terms(1, 3, 0.0);

  forcing_func(state_parameters, state_variables, forcing_terms);

  // rate = k / [S]^(1-1) * [A] = k * [A] = 0.1 * 2.0 = 0.2
  double expected_rate = k * 2.0;

  // A (reactant): -rate
  EXPECT_NEAR(forcing_terms[0][0], -expected_rate, 1e-10);
  // B (product): +rate
  EXPECT_NEAR(forcing_terms[0][1], expected_rate, 1e-10);
  // S (solvent, not a reactant/product): unchanged
  EXPECT_EQ(forcing_terms[0][2], 0.0);
}

TEST(DissolvedReaction, ForcingFunctionSolventNormalization)
{
  using MatrixPolicy = micm::VectorMatrix<double>;

  auto a = micm::Species{ "A" };
  auto b = micm::Species{ "B" };
  auto c = micm::Species{ "C" };
  auto solvent = micm::Species{ "S" };

  auto phase = micm::Phase{ "AQUEOUS", { { a }, { b }, { c }, { solvent } } };

  double k = 1.0;
  auto rate = [k](const micm::Conditions& conditions) { return k; };

  // A + B -> C with solvent S
  DissolvedReaction reaction{ { { "DROP", rate } },
                              { a, b },  // 2 reactants
                              { c },     // 1 product
                              solvent,
                              phase };

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  std::string k_param = "DROP." + phase.name_ + "." + reaction.uuid_ + ".k";

  std::unordered_map<std::string, std::size_t> state_parameter_indices;
  state_parameter_indices[k_param] = 0;

  std::unordered_map<std::string, std::size_t> state_variable_indices;
  state_variable_indices["DROP.AQUEOUS.A"] = 0;
  state_variable_indices["DROP.AQUEOUS.B"] = 1;
  state_variable_indices["DROP.AQUEOUS.C"] = 2;
  state_variable_indices["DROP.AQUEOUS.S"] = 3;

  auto forcing_func =
      reaction.ForcingFunction<MatrixPolicy>(phase_prefixes, state_parameter_indices, state_variable_indices);

  MatrixPolicy state_parameters(1, 1);
  state_parameters[0][0] = k;

  MatrixPolicy state_variables(1, 4);
  state_variables[0][0] = 0.001;  // [A]
  state_variables[0][1] = 0.002;  // [B]
  state_variables[0][2] = 0.0;    // [C]
  state_variables[0][3] = 50.0;   // [S]

  MatrixPolicy forcing_terms(1, 4, 0.0);

  forcing_func(state_parameters, state_variables, forcing_terms);

  // rate = k / [S]^(2-1) * [A] * [B] = k / [S] * [A] * [B]
  double expected_rate = k / 50.0 * 0.001 * 0.002;

  EXPECT_NEAR(forcing_terms[0][0], -expected_rate, 1e-15);  // A (reactant)
  EXPECT_NEAR(forcing_terms[0][1], -expected_rate, 1e-15);  // B (reactant)
  EXPECT_NEAR(forcing_terms[0][2], expected_rate, 1e-15);   // C (product)
  EXPECT_EQ(forcing_terms[0][3], 0.0);                      // S (solvent only)
}

TEST(DissolvedReaction, ForcingFunctionMultipleProducts)
{
  using MatrixPolicy = micm::VectorMatrix<double>;

  auto a = micm::Species{ "A" };
  auto b = micm::Species{ "B" };
  auto c = micm::Species{ "C" };
  auto solvent = micm::Species{ "S" };

  auto phase = micm::Phase{ "AQUEOUS", { { a }, { b }, { c }, { solvent } } };

  double k = 2.0;
  auto rate = [k](const micm::Conditions& conditions) { return k; };

  // A -> B + C with solvent S
  DissolvedReaction reaction{ { { "DROP", rate } },
                              { a },     // 1 reactant
                              { b, c },  // 2 products
                              solvent,
                              phase };

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  std::string k_param = "DROP." + phase.name_ + "." + reaction.uuid_ + ".k";

  std::unordered_map<std::string, std::size_t> state_parameter_indices;
  state_parameter_indices[k_param] = 0;

  std::unordered_map<std::string, std::size_t> state_variable_indices;
  state_variable_indices["DROP.AQUEOUS.A"] = 0;
  state_variable_indices["DROP.AQUEOUS.B"] = 1;
  state_variable_indices["DROP.AQUEOUS.C"] = 2;
  state_variable_indices["DROP.AQUEOUS.S"] = 3;

  auto forcing_func =
      reaction.ForcingFunction<MatrixPolicy>(phase_prefixes, state_parameter_indices, state_variable_indices);

  MatrixPolicy state_parameters(1, 1);
  state_parameters[0][0] = k;

  MatrixPolicy state_variables(1, 4);
  state_variables[0][0] = 3.0;   // [A]
  state_variables[0][1] = 1.0;   // [B]
  state_variables[0][2] = 0.5;   // [C]
  state_variables[0][3] = 40.0;  // [S]

  MatrixPolicy forcing_terms(1, 4, 0.0);

  forcing_func(state_parameters, state_variables, forcing_terms);

  // rate = k / [S]^(1-1) * [A] = k * [A] = 2.0 * 3.0 = 6.0
  double expected_rate = k * 3.0;

  EXPECT_NEAR(forcing_terms[0][0], -expected_rate, 1e-10);  // A (reactant): -rate
  EXPECT_NEAR(forcing_terms[0][1], expected_rate, 1e-10);   // B (product): +rate
  EXPECT_NEAR(forcing_terms[0][2], expected_rate, 1e-10);   // C (product): +rate
  EXPECT_EQ(forcing_terms[0][3], 0.0);                      // S (solvent only)
}

TEST(DissolvedReaction, ForcingFunctionMultipleCells)
{
  using MatrixPolicy = micm::VectorMatrix<double>;

  auto a = micm::Species{ "A" };
  auto b = micm::Species{ "B" };
  auto solvent = micm::Species{ "S" };

  auto phase = micm::Phase{ "AQUEOUS", { { a }, { b }, { solvent } } };

  double k = 0.1;
  auto rate = [k](const micm::Conditions& conditions) { return k; };

  DissolvedReaction reaction{ { { "MODE1", rate } }, { a }, { b }, solvent, phase };

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("MODE1");

  std::string k_param = "MODE1." + phase.name_ + "." + reaction.uuid_ + ".k";

  std::unordered_map<std::string, std::size_t> state_parameter_indices;
  state_parameter_indices[k_param] = 0;

  std::unordered_map<std::string, std::size_t> state_variable_indices;
  state_variable_indices["MODE1.AQUEOUS.A"] = 0;
  state_variable_indices["MODE1.AQUEOUS.B"] = 1;
  state_variable_indices["MODE1.AQUEOUS.S"] = 2;

  auto forcing_func =
      reaction.ForcingFunction<MatrixPolicy>(phase_prefixes, state_parameter_indices, state_variable_indices);

  MatrixPolicy state_parameters(3, 1);
  MatrixPolicy state_variables(3, 3);
  MatrixPolicy forcing_terms(3, 3, 0.0);

  for (std::size_t i = 0; i < 3; ++i)
  {
    state_parameters[i][0] = k;
    state_variables[i][0] = 1.0 + i * 0.5;  // Different [A]
    state_variables[i][1] = 0.0;            // [B]
    state_variables[i][2] = 55.0;           // [S]
  }

  forcing_func(state_parameters, state_variables, forcing_terms);

  for (std::size_t i = 0; i < 3; ++i)
  {
    double a_conc = 1.0 + i * 0.5;
    double expected_rate = k * a_conc;

    EXPECT_NEAR(forcing_terms[i][0], -expected_rate, 1e-10);
    EXPECT_NEAR(forcing_terms[i][1], expected_rate, 1e-10);
    EXPECT_EQ(forcing_terms[i][2], 0.0);
  }
}

TEST(DissolvedReaction, ForcingFunctionMultiplePhaseInstances)
{
  using MatrixPolicy = micm::VectorMatrix<double>;

  auto a = micm::Species{ "A" };
  auto b = micm::Species{ "B" };
  auto solvent = micm::Species{ "S" };

  auto phase = micm::Phase{ "AQUEOUS", { { a }, { b }, { solvent } } };

  double k = 0.1;
  auto rate = [k](const micm::Conditions& conditions) { return k; };

  DissolvedReaction reaction{ { { "SMALL_DROP", rate }, { "LARGE_DROP", rate } }, { a }, { b }, solvent, phase };

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("SMALL_DROP");
  phase_prefixes["AQUEOUS"].insert("LARGE_DROP");

  std::string k_small = "SMALL_DROP." + phase.name_ + "." + reaction.uuid_ + ".k";
  std::string k_large = "LARGE_DROP." + phase.name_ + "." + reaction.uuid_ + ".k";

  std::unordered_map<std::string, std::size_t> state_parameter_indices;
  state_parameter_indices[k_small] = 0;
  state_parameter_indices[k_large] = 1;

  std::unordered_map<std::string, std::size_t> state_variable_indices;
  state_variable_indices["SMALL_DROP.AQUEOUS.A"] = 0;
  state_variable_indices["SMALL_DROP.AQUEOUS.B"] = 1;
  state_variable_indices["SMALL_DROP.AQUEOUS.S"] = 2;
  state_variable_indices["LARGE_DROP.AQUEOUS.A"] = 3;
  state_variable_indices["LARGE_DROP.AQUEOUS.B"] = 4;
  state_variable_indices["LARGE_DROP.AQUEOUS.S"] = 5;

  auto forcing_func =
      reaction.ForcingFunction<MatrixPolicy>(phase_prefixes, state_parameter_indices, state_variable_indices);

  MatrixPolicy state_parameters(1, 2);
  state_parameters[0][0] = k;
  state_parameters[0][1] = k;

  MatrixPolicy state_variables(1, 6);
  // Small drop
  state_variables[0][0] = 2.0;   // [A]
  state_variables[0][1] = 0.0;   // [B]
  state_variables[0][2] = 50.0;  // [S]
  // Large drop
  state_variables[0][3] = 5.0;   // [A]
  state_variables[0][4] = 0.0;   // [B]
  state_variables[0][5] = 60.0;  // [S]

  MatrixPolicy forcing_terms(1, 6, 0.0);

  forcing_func(state_parameters, state_variables, forcing_terms);

  // Small drop: rate = k * [A] = 0.1 * 2.0 = 0.2
  double rate_small = k * 2.0;
  EXPECT_NEAR(forcing_terms[0][0], -rate_small, 1e-10);
  EXPECT_NEAR(forcing_terms[0][1], rate_small, 1e-10);
  EXPECT_EQ(forcing_terms[0][2], 0.0);

  // Large drop: rate = k * [A] = 0.1 * 5.0 = 0.5
  double rate_large = k * 5.0;
  EXPECT_NEAR(forcing_terms[0][3], -rate_large, 1e-10);
  EXPECT_NEAR(forcing_terms[0][4], rate_large, 1e-10);
  EXPECT_EQ(forcing_terms[0][5], 0.0);
}

// ============================================================================
// JacobianFunction Tests
// ============================================================================

TEST(DissolvedReaction, JacobianFunctionBasicPartials)
{
  using MatrixPolicy = micm::Matrix<double>;
  using SparseMatrixPolicy = micm::SparseMatrix<double, micm::SparseMatrixStandardOrderingCompressedSparseRow>;

  auto a = micm::Species{ "A" };
  auto b = micm::Species{ "B" };
  auto solvent = micm::Species{ "S" };

  auto phase = micm::Phase{ "AQUEOUS", { { a }, { b }, { solvent } } };

  double k = 0.1;
  auto rate = [k](const micm::Conditions& conditions) { return k; };

  // A -> B with solvent S
  DissolvedReaction reaction{ { { "MODE1", rate } }, { a }, { b }, solvent, phase };

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("MODE1");

  std::string k_param = "MODE1." + phase.name_ + "." + reaction.uuid_ + ".k";

  std::unordered_map<std::string, std::size_t> state_parameter_indices;
  state_parameter_indices[k_param] = 0;

  std::unordered_map<std::string, std::size_t> state_variable_indices;
  state_variable_indices["MODE1.AQUEOUS.A"] = 0;
  state_variable_indices["MODE1.AQUEOUS.B"] = 1;
  state_variable_indices["MODE1.AQUEOUS.S"] = 2;

  auto jacobian_elements = reaction.NonZeroJacobianElements(phase_prefixes, state_variable_indices);
  auto jacobian_builder = SparseMatrixPolicy::Create(3).SetNumberOfBlocks(1);
  for (const auto& elem : jacobian_elements)
  {
    jacobian_builder.WithElement(elem.first, elem.second);
  }
  SparseMatrixPolicy jacobian(jacobian_builder);
  jacobian.Fill(0.0);

  auto jacobian_func = reaction.JacobianFunction<MatrixPolicy, SparseMatrixPolicy>(
      phase_prefixes, state_parameter_indices, state_variable_indices, jacobian);

  MatrixPolicy state_parameters(1, 1);
  state_parameters[0][0] = k;

  MatrixPolicy state_variables(1, 3);
  state_variables[0][0] = 2.0;   // [A]
  state_variables[0][1] = 0.5;   // [B]
  state_variables[0][2] = 55.0;  // [S]

  jacobian_func(state_parameters, state_variables, jacobian);

  // For A -> B (1 reactant, n_r = 1):
  // rate = k / [S]^0 * [A] = k * [A]
  // d(rate)/d[A] = k
  // d(rate)/d[S] = k * (1-1) / [S]^1 * [A] = 0

  // -J[A,A] = +d(rate)/d[A] = +k (stores negative Jacobian)
  EXPECT_NEAR(jacobian[0][0][0], k, 1e-9);
  // -J[B,A] = -d(rate)/d[A] = -k
  EXPECT_NEAR(jacobian[0][1][0], -k, 1e-9);
  // -J[A,S] = +d(rate)/d[S] ~ 0 (small damping residual)
  EXPECT_NEAR(jacobian[0][0][2], 0.0, 1e-9);
  // -J[B,S] = -d(rate)/d[S] ~ 0 (small damping residual)
  EXPECT_NEAR(jacobian[0][1][2], 0.0, 1e-9);
}

TEST(DissolvedReaction, JacobianFunctionMultipleReactants)
{
  using MatrixPolicy = micm::Matrix<double>;
  using SparseMatrixPolicy = micm::SparseMatrix<double, micm::SparseMatrixStandardOrderingCompressedSparseRow>;

  auto a = micm::Species{ "A" };
  auto b = micm::Species{ "B" };
  auto c = micm::Species{ "C" };
  auto solvent = micm::Species{ "S" };

  auto phase = micm::Phase{ "AQUEOUS", { { a }, { b }, { c }, { solvent } } };

  double k = 1.0;
  auto rate = [k](const micm::Conditions& conditions) { return k; };

  // A + B -> C with solvent S
  DissolvedReaction reaction{ { { "DROP", rate } },
                              { a, b },  // 2 reactants
                              { c },     // 1 product
                              solvent,
                              phase };

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  std::string k_param = "DROP." + phase.name_ + "." + reaction.uuid_ + ".k";

  std::unordered_map<std::string, std::size_t> state_parameter_indices;
  state_parameter_indices[k_param] = 0;

  std::unordered_map<std::string, std::size_t> state_variable_indices;
  state_variable_indices["DROP.AQUEOUS.A"] = 0;
  state_variable_indices["DROP.AQUEOUS.B"] = 1;
  state_variable_indices["DROP.AQUEOUS.C"] = 2;
  state_variable_indices["DROP.AQUEOUS.S"] = 3;

  auto jacobian_elements = reaction.NonZeroJacobianElements(phase_prefixes, state_variable_indices);
  auto jacobian_builder = SparseMatrixPolicy::Create(4).SetNumberOfBlocks(1);
  for (const auto& elem : jacobian_elements)
  {
    jacobian_builder.WithElement(elem.first, elem.second);
  }
  SparseMatrixPolicy jacobian(jacobian_builder);
  jacobian.Fill(0.0);

  auto jacobian_func = reaction.JacobianFunction<MatrixPolicy, SparseMatrixPolicy>(
      phase_prefixes, state_parameter_indices, state_variable_indices, jacobian);

  MatrixPolicy state_parameters(1, 1);
  state_parameters[0][0] = k;

  MatrixPolicy state_variables(1, 4);
  state_variables[0][0] = 0.001;  // [A]
  state_variables[0][1] = 0.002;  // [B]
  state_variables[0][2] = 0.0;    // [C]
  state_variables[0][3] = 50.0;   // [S]

  jacobian_func(state_parameters, state_variables, jacobian);

  // rate = k / [S] * [A] * [B]
  // d(rate)/d[A] = k / [S] * [B] = 1.0 / 50.0 * 0.002 = 4e-5
  double d_rate_d_a = k / 50.0 * 0.002;
  // d(rate)/d[B] = k / [S] * [A] = 1.0 / 50.0 * 0.001 = 2e-5
  double d_rate_d_b = k / 50.0 * 0.001;
  // d(rate)/d[S] = k * (1-2) / [S]^2 * [A] * [B] = -1.0 / 2500 * 0.001 * 0.002
  double d_rate_d_s = k * (1 - 2) / std::pow(50.0, 2) * 0.001 * 0.002;

  // Stores -J convention: reactant rows get +, product rows get -
  // -J[A,A]
  EXPECT_NEAR(jacobian[0][0][0], d_rate_d_a, 1e-15);
  // -J[A,B]
  EXPECT_NEAR(jacobian[0][0][1], d_rate_d_b, 1e-15);
  // -J[B,A]
  EXPECT_NEAR(jacobian[0][1][0], d_rate_d_a, 1e-15);
  // -J[B,B]
  EXPECT_NEAR(jacobian[0][1][1], d_rate_d_b, 1e-15);
  // -J[C,A]
  EXPECT_NEAR(jacobian[0][2][0], -d_rate_d_a, 1e-15);
  // -J[C,B]
  EXPECT_NEAR(jacobian[0][2][1], -d_rate_d_b, 1e-15);
  // -J[A,S]
  EXPECT_NEAR(jacobian[0][0][3], d_rate_d_s, 1e-15);
  // -J[B,S]
  EXPECT_NEAR(jacobian[0][1][3], d_rate_d_s, 1e-15);
  // -J[C,S]
  EXPECT_NEAR(jacobian[0][2][3], -d_rate_d_s, 1e-15);
}

TEST(DissolvedReaction, JacobianFunctionSolventPartial)
{
  using MatrixPolicy = micm::Matrix<double>;
  using SparseMatrixPolicy = micm::SparseMatrix<double, micm::SparseMatrixStandardOrderingCompressedSparseRow>;

  auto a = micm::Species{ "A" };
  auto b = micm::Species{ "B" };
  auto c = micm::Species{ "C" };

  auto phase = micm::Phase{ "AQUEOUS", { { a }, { b }, { c } } };

  double k = 2.0;
  auto rate = [k](const micm::Conditions& conditions) { return k; };

  // A + C -> B with solvent C (solvent is a reactant)
  // n_r = 2, rate = k / [C]^1 * [A] * [C] = k * [A]
  DissolvedReaction reaction{ { { "MODE1", rate } },
                              { a, c },  // 2 reactants (C is also solvent)
                              { b },     // 1 product
                              c,         // solvent = reactant C
                              phase };

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("MODE1");

  std::string k_param = "MODE1." + phase.name_ + "." + reaction.uuid_ + ".k";

  std::unordered_map<std::string, std::size_t> state_parameter_indices;
  state_parameter_indices[k_param] = 0;

  std::unordered_map<std::string, std::size_t> state_variable_indices;
  state_variable_indices["MODE1.AQUEOUS.A"] = 0;
  state_variable_indices["MODE1.AQUEOUS.C"] = 1;
  state_variable_indices["MODE1.AQUEOUS.B"] = 2;

  auto jacobian_elements = reaction.NonZeroJacobianElements(phase_prefixes, state_variable_indices);
  auto jacobian_builder = SparseMatrixPolicy::Create(3).SetNumberOfBlocks(1);
  for (const auto& elem : jacobian_elements)
  {
    jacobian_builder.WithElement(elem.first, elem.second);
  }
  SparseMatrixPolicy jacobian(jacobian_builder);
  jacobian.Fill(0.0);

  auto jacobian_func = reaction.JacobianFunction<MatrixPolicy, SparseMatrixPolicy>(
      phase_prefixes, state_parameter_indices, state_variable_indices, jacobian);

  MatrixPolicy state_parameters(1, 1);
  state_parameters[0][0] = k;

  MatrixPolicy state_variables(1, 3);
  state_variables[0][0] = 0.5;   // [A]
  state_variables[0][1] = 40.0;  // [C] (solvent and reactant)
  state_variables[0][2] = 0.0;   // [B]

  jacobian_func(state_parameters, state_variables, jacobian);

  // rate = k / [C]^(2-1) * [A] * [C] = k * [A]
  // d(rate)/d[A] = k / [C] * [C] = k
  double d_rate_d_a = k / 40.0 * 40.0;
  // d(rate)/d[C] (as reactant) = k / [C] * [A]
  double d_rate_d_c_reactant = k / 40.0 * 0.5;
  // d(rate)/d[C] (as solvent) = k * (1-2) / [C]^2 * [A] * [C] = -k * [A] / [C]
  double d_rate_d_c_solvent = k * (1 - 2) / std::pow(40.0, 2) * 0.5 * 40.0;

  // Stores -J convention
  // -J[A,A] = +d(rate)/d[A]
  EXPECT_NEAR(jacobian[0][0][0], d_rate_d_a, 1e-10);
  // -J[B,A] = -d(rate)/d[A]
  EXPECT_NEAR(jacobian[0][2][0], -d_rate_d_a, 1e-10);
  // -J[A,C] = +d(rate)/d[C]_reactant + d(rate)/d[C]_solvent
  EXPECT_NEAR(jacobian[0][0][1], d_rate_d_c_reactant + d_rate_d_c_solvent, 1e-10);
  // -J[B,C] = -d(rate)/d[C]_reactant - d(rate)/d[C]_solvent
  EXPECT_NEAR(jacobian[0][2][1], -d_rate_d_c_reactant - d_rate_d_c_solvent, 1e-10);
}

// ============================================================================
// FD Jacobian Tests
// ============================================================================

// FD check: A -> B (unimolecular, 1 cell)
TEST(DissolvedReaction, JacobianFDUnimolecular)
{
  auto a = micm::Species{ "A" };
  auto b = micm::Species{ "B" };
  auto s = micm::Species{ "S" };
  auto phase = micm::Phase{ "AQUEOUS", { { a }, { b }, { s } } };

  double k = 0.1;
  DissolvedReaction reaction{ { { "DROP", [k](const micm::Conditions&) { return k; } } }, { a }, { b }, s, phase };

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  std::string k_param = "DROP." + phase.name_ + "." + reaction.uuid_ + ".k";
  std::unordered_map<std::string, std::size_t> spi{ { k_param, 0 } };
  std::unordered_map<std::string, std::size_t> svi{ { "DROP.AQUEOUS.A", 0 },
                                                    { "DROP.AQUEOUS.B", 1 },
                                                    { "DROP.AQUEOUS.S", 2 } };

  MatrixPolicy params(1, 1);
  params[0][0] = k;
  MatrixPolicy vars(1, 3);
  vars[0][0] = 2.0;
  vars[0][1] = 0.5;
  vars[0][2] = 55.0;

  CheckFiniteDifferenceJacobian(reaction, phase_prefixes, spi, svi, params, vars);
}

// FD check: A + B -> C (bimolecular, solvent normalization active)
TEST(DissolvedReaction, JacobianFDBimolecular)
{
  auto a = micm::Species{ "A" };
  auto b = micm::Species{ "B" };
  auto c = micm::Species{ "C" };
  auto s = micm::Species{ "S" };
  auto phase = micm::Phase{ "AQUEOUS", { { a }, { b }, { c }, { s } } };

  double k = 1.0;
  DissolvedReaction reaction{ { { "DROP", [k](const micm::Conditions&) { return k; } } }, { a, b }, { c }, s, phase };

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  std::string k_param = "DROP." + phase.name_ + "." + reaction.uuid_ + ".k";
  std::unordered_map<std::string, std::size_t> spi{ { k_param, 0 } };
  std::unordered_map<std::string, std::size_t> svi{
    { "DROP.AQUEOUS.A", 0 }, { "DROP.AQUEOUS.B", 1 }, { "DROP.AQUEOUS.C", 2 }, { "DROP.AQUEOUS.S", 3 }
  };

  MatrixPolicy params(1, 1);
  params[0][0] = k;
  MatrixPolicy vars(1, 4);
  vars[0][0] = 0.001;
  vars[0][1] = 0.002;
  vars[0][2] = 0.0;
  vars[0][3] = 50.0;

  CheckFiniteDifferenceJacobian(reaction, phase_prefixes, spi, svi, params, vars);
}

// FD check: multiple grid cells
TEST(DissolvedReaction, JacobianFDMultiCell)
{
  auto a = micm::Species{ "A" };
  auto b = micm::Species{ "B" };
  auto s = micm::Species{ "S" };
  auto phase = micm::Phase{ "AQUEOUS", { { a }, { b }, { s } } };

  double k = 0.5;
  DissolvedReaction reaction{ { { "DROP", [k](const micm::Conditions&) { return k; } } }, { a }, { b }, s, phase };

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  std::string k_param = "DROP." + phase.name_ + "." + reaction.uuid_ + ".k";
  std::unordered_map<std::string, std::size_t> spi{ { k_param, 0 } };
  std::unordered_map<std::string, std::size_t> svi{ { "DROP.AQUEOUS.A", 0 },
                                                    { "DROP.AQUEOUS.B", 1 },
                                                    { "DROP.AQUEOUS.S", 2 } };

  std::size_t n_cells = 4;
  MatrixPolicy params(n_cells, 1);
  MatrixPolicy vars(n_cells, 3);
  double gas_concs[] = { 1.0e-3, 5.0e-3, 1.0e-2, 5.0e-2 };
  double solvent_concs[] = { 55.0, 50.0, 45.0, 40.0 };
  for (std::size_t i = 0; i < n_cells; ++i)
  {
    params[i][0] = k;
    vars[i][0] = gas_concs[i];
    vars[i][1] = 0.0;
    vars[i][2] = solvent_concs[i];
  }

  CheckFiniteDifferenceJacobian(reaction, phase_prefixes, spi, svi, params, vars);
}

// FD check: two phase instances (SMALL_DROP + LARGE_DROP)
TEST(DissolvedReaction, JacobianFDMultiPhaseInstance)
{
  auto a = micm::Species{ "A" };
  auto b = micm::Species{ "B" };
  auto s = micm::Species{ "S" };
  auto phase = micm::Phase{ "AQUEOUS", { { a }, { b }, { s } } };

  double k = 0.2;
  DissolvedReaction reaction{ { { "SMALL", [k](const micm::Conditions&) { return k; } },
                                { "LARGE", [k](const micm::Conditions&) { return k; } } },
                              { a },
                              { b },
                              s,
                              phase };

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("SMALL");
  phase_prefixes["AQUEOUS"].insert("LARGE");

  std::string k_small = "SMALL." + phase.name_ + "." + reaction.uuid_ + ".k";
  std::string k_large = "LARGE." + phase.name_ + "." + reaction.uuid_ + ".k";
  std::unordered_map<std::string, std::size_t> spi{ { k_small, 0 }, { k_large, 1 } };
  std::unordered_map<std::string, std::size_t> svi{ { "SMALL.AQUEOUS.A", 0 }, { "SMALL.AQUEOUS.B", 1 },
                                                    { "SMALL.AQUEOUS.S", 2 }, { "LARGE.AQUEOUS.A", 3 },
                                                    { "LARGE.AQUEOUS.B", 4 }, { "LARGE.AQUEOUS.S", 5 } };

  MatrixPolicy params(1, 2);
  params[0][0] = k;
  params[0][1] = k;
  MatrixPolicy vars(1, 6);
  vars[0][0] = 2.0e-3;
  vars[0][1] = 0.0;
  vars[0][2] = 55.0;
  vars[0][3] = 5.0e-3;
  vars[0][4] = 0.0;
  vars[0][5] = 50.0;

  CheckFiniteDifferenceJacobian(reaction, phase_prefixes, spi, svi, params, vars);
}

// FD check: solvent also appears as a reactant (n_r=2, [S] normalization non-trivial)
TEST(DissolvedReaction, JacobianFDSolventIsReactant)
{
  auto a = micm::Species{ "A" };
  auto c = micm::Species{ "C" };  // solvent and reactant
  auto b = micm::Species{ "B" };
  auto phase = micm::Phase{ "AQUEOUS", { { a }, { b }, { c } } };

  double k = 2.0;
  DissolvedReaction reaction{ { { "DROP", [k](const micm::Conditions&) { return k; } } }, { a, c }, { b }, c, phase };

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  std::string k_param = "DROP." + phase.name_ + "." + reaction.uuid_ + ".k";
  std::unordered_map<std::string, std::size_t> spi{ { k_param, 0 } };
  std::unordered_map<std::string, std::size_t> svi{ { "DROP.AQUEOUS.A", 0 },
                                                    { "DROP.AQUEOUS.C", 1 },
                                                    { "DROP.AQUEOUS.B", 2 } };

  MatrixPolicy params(1, 1);
  params[0][0] = k;
  MatrixPolicy vars(1, 3);
  vars[0][0] = 0.5;
  vars[0][1] = 40.0;
  vars[0][2] = 0.0;

  CheckFiniteDifferenceJacobian(reaction, phase_prefixes, spi, svi, params, vars);
}

// ============================================================================
// Solvent-floor singularity-protection tests
// ============================================================================

// When [S] → 0, rate should go to zero (not infinity/NaN)
TEST(DissolvedReaction, ForcingFunctionSolventFloorZeroSolvent)
{
  using VM = micm::VectorMatrix<double>;

  auto a = micm::Species{ "A" };
  auto b = micm::Species{ "B" };
  auto c = micm::Species{ "C" };
  auto s = micm::Species{ "S" };
  auto phase = micm::Phase{ "AQUEOUS", { { a }, { b }, { c }, { s } } };

  double k = 1.0;
  DissolvedReaction reaction{ { { "DROP", [k](const micm::Conditions&) { return k; } } }, { a, b }, { c }, s, phase };

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  std::string k_param = "DROP." + phase.name_ + "." + reaction.uuid_ + ".k";
  std::unordered_map<std::string, std::size_t> spi{ { k_param, 0 } };
  std::unordered_map<std::string, std::size_t> svi{
    { "DROP.AQUEOUS.A", 0 }, { "DROP.AQUEOUS.B", 1 }, { "DROP.AQUEOUS.C", 2 }, { "DROP.AQUEOUS.S", 3 }
  };

  auto ff = reaction.ForcingFunction<VM>(phase_prefixes, spi, svi);

  VM params(1, 1);
  params[0][0] = k;
  VM vars(1, 4);
  vars[0][0] = 0.001;
  vars[0][1] = 0.002;
  vars[0][2] = 0.0;
  vars[0][3] = 0.0;  // [S] = 0
  VM forcing(1, 4, 0.0);

  ff(params, vars, forcing);

  // rate = k * 0 / (0+eps)^2 * [A]*[B] → 0
  EXPECT_NEAR(forcing[0][0], 0.0, 1e-10) << "Forcing should be zero when [S]=0";
  EXPECT_NEAR(forcing[0][1], 0.0, 1e-10);
  EXPECT_NEAR(forcing[0][2], 0.0, 1e-10);
  EXPECT_FALSE(std::isnan(forcing[0][0])) << "Should not produce NaN";
  EXPECT_FALSE(std::isinf(forcing[0][0])) << "Should not produce Inf";
}

// Jacobian should also be finite and well-defined at [S]=0
TEST(DissolvedReaction, JacobianFunctionSolventFloorZeroSolvent)
{
  auto a = micm::Species{ "A" };
  auto b = micm::Species{ "B" };
  auto c = micm::Species{ "C" };
  auto s = micm::Species{ "S" };
  auto phase = micm::Phase{ "AQUEOUS", { { a }, { b }, { c }, { s } } };

  double k = 1.0;
  DissolvedReaction reaction{ [k](const micm::Conditions&) { return k; }, { a, b }, { c }, s, phase };

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  std::string k_param = phase.name_ + "." + reaction.uuid_ + ".k";
  std::unordered_map<std::string, std::size_t> spi{ { k_param, 0 } };
  std::unordered_map<std::string, std::size_t> svi{
    { "DROP.AQUEOUS.A", 0 }, { "DROP.AQUEOUS.B", 1 }, { "DROP.AQUEOUS.C", 2 }, { "DROP.AQUEOUS.S", 3 }
  };

  auto jac_elements = reaction.NonZeroJacobianElements(phase_prefixes, svi);
  auto jac_builder = SparseMatrixPolicy::Create(4).SetNumberOfBlocks(1);
  for (const auto& elem : jac_elements)
    jac_builder.WithElement(elem.first, elem.second);
  SparseMatrixPolicy jacobian(jac_builder);
  jacobian.Fill(0.0);

  auto jac_func = reaction.JacobianFunction<MatrixPolicy, SparseMatrixPolicy>(phase_prefixes, spi, svi, jacobian);

  MatrixPolicy params(1, 1);
  params[0][0] = k;
  MatrixPolicy vars(1, 4);
  vars[0][0] = 0.001;
  vars[0][1] = 0.002;
  vars[0][2] = 0.0;
  vars[0][3] = 0.0;  // [S] = 0
  jac_func(params, vars, jacobian);

  // When [S] = 0 the rate r = k*[S]/([S]+δ)^n_r * ∏Ri = 0, so
  // ∂r/∂[Ri] = k*[S]/([S]+δ)^n_r * ∏_{j≠i}Rj → 0 naturally (reactant/product columns).
  // ∂r/∂[S]  = k*(δ+(1-n_r)*[S])/([S]+δ)^{n_r+1} * ∏Ri → k/δ^n_r * ∏Ri (finite, but large).
  // Verify: no NaN/Inf anywhere; reactant/product columns near zero; solvent column finite.
  std::size_t solvent_col = svi.at("DROP.AQUEOUS.S");
  const double eps = 1.0e-20;
  const std::size_t n_r = 2;  // reactants = {A, B}
  const double expected_dr_dS =
      1.0 * (eps + (1.0 - n_r) * 0.0) / std::pow(eps, n_r + 1) * vars[0][0] * vars[0][1];  // k * (δ/δ^{n_r+1}) * [A]*[B]
  for (std::size_t i = 0; i < 4; ++i)
    for (std::size_t j = 0; j < 4; ++j)
    {
      try
      {
        double val = jacobian[0][i][j];
        EXPECT_FALSE(std::isnan(val)) << "NaN at [" << i << "," << j << "]";
        EXPECT_FALSE(std::isinf(val)) << "Inf at [" << i << "," << j << "]";
        if (j == solvent_col)
        {
          // ∂r/∂[S] is large but finite; sign flips between reactant and product rows
          double sign = (i == 2) ? -1.0 : 1.0;  // product C is row 2
          EXPECT_NEAR(val, sign * expected_dr_dS, std::abs(expected_dr_dS) * 1e-10)
              << "Unexpected solvent-column Jacobian at [" << i << "," << j << "]";
        }
        else
        {
          EXPECT_NEAR(val, 0.0, 1e-5) << "Expected near-zero reactant/product-column Jacobian when [S]=0 at [" << i << ","
                                      << j << "]";
        }
      }
      catch (...)
      { /* sparse zero — expected */
      }
    }
}

// ============================================================================
// Rate-cap (min_halflife_) tests
// ============================================================================

// In the uncapped regime (rate << r_max), capped == uncapped
TEST(DissolvedReaction, ForcingFunctionCappedUncappedRegime)
{
  using VM = micm::VectorMatrix<double>;

  auto a = micm::Species{ "A" };
  auto b = micm::Species{ "B" };
  auto s = micm::Species{ "S" };
  auto phase = micm::Phase{ "AQUEOUS", { { a }, { b }, { s } } };

  double k = 1.0e-6;    // very slow
  double t_half = 1.0;  // short half-life → r_max = [A] / t_half = large compared to r

  DissolvedReaction uncapped{ [k](const micm::Conditions&) { return k; }, { a }, { b }, s, phase };
  DissolvedReaction capped{ [k](const micm::Conditions&) { return k; }, { a }, { b }, s, phase, 1.0e-20, t_half };

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");
  std::string k_param = phase.name_ + "." + uncapped.uuid_ + ".k";
  std::string k_param_c = phase.name_ + "." + capped.uuid_ + ".k";
  std::unordered_map<std::string, std::size_t> spi_u{ { k_param, 0 } };
  std::unordered_map<std::string, std::size_t> spi_c{ { k_param_c, 0 } };
  std::unordered_map<std::string, std::size_t> svi{ { "DROP.AQUEOUS.A", 0 },
                                                    { "DROP.AQUEOUS.B", 1 },
                                                    { "DROP.AQUEOUS.S", 2 } };

  auto ff_u = uncapped.ForcingFunction<VM>(phase_prefixes, spi_u, svi);
  auto ff_c = capped.ForcingFunction<VM>(phase_prefixes, spi_c, svi);

  VM params(1, 1);
  params[0][0] = k;
  VM vars(1, 3);
  vars[0][0] = 1.0;
  vars[0][1] = 0.0;
  vars[0][2] = 55.0;  // [A]=1

  // r = k * [A] = 1e-6; r_max = [A] / t_half = 1.0 → r << r_max
  VM forcing_u(1, 3, 0.0);
  VM forcing_c(1, 3, 0.0);
  ff_u(params, vars, forcing_u);
  ff_c(params, vars, forcing_c);

  EXPECT_NEAR(forcing_c[0][0], forcing_u[0][0], std::abs(forcing_u[0][0]) * 1e-10)
      << "Capped should equal uncapped when rate << r_max";
  EXPECT_NEAR(forcing_c[0][1], forcing_u[0][1], std::abs(forcing_u[0][1]) * 1e-10);
}

// In the saturated regime (rate >> r_max), |forcing| ≈ r_max = [A] / t_half
TEST(DissolvedReaction, ForcingFunctionCappedSaturatedRegime)
{
  using VM = micm::VectorMatrix<double>;

  auto a = micm::Species{ "A" };
  auto b = micm::Species{ "B" };
  auto s = micm::Species{ "S" };
  auto phase = micm::Phase{ "AQUEOUS", { { a }, { b }, { s } } };

  double k = 1.0e6;        // very fast reaction
  double t_half = 1000.0;  // long half-life → r_max = [A] / t_half = small

  DissolvedReaction reaction{ [k](const micm::Conditions&) { return k; }, { a }, { b }, s, phase, 1.0e-20, t_half };

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");
  std::string k_param = phase.name_ + "." + reaction.uuid_ + ".k";
  std::unordered_map<std::string, std::size_t> spi{ { k_param, 0 } };
  std::unordered_map<std::string, std::size_t> svi{ { "DROP.AQUEOUS.A", 0 },
                                                    { "DROP.AQUEOUS.B", 1 },
                                                    { "DROP.AQUEOUS.S", 2 } };

  auto ff = reaction.ForcingFunction<VM>(phase_prefixes, spi, svi);

  VM params(1, 1);
  params[0][0] = k;
  VM vars(1, 3);
  double A = 1.0;
  vars[0][0] = A;
  vars[0][1] = 0.0;
  vars[0][2] = 55.0;

  // r = k * [A] = 1e6; r_max = [A] / t_half = 1e-3 → rate >> r_max
  // expected capped rate ≈ r_max = A / t_half (tanh(u) ≈ 1)
  VM forcing(1, 3, 0.0);
  ff(params, vars, forcing);

  double r_max = A / t_half;
  EXPECT_NEAR(forcing[0][0], -r_max, r_max * 1e-6) << "Capped forcing should approach -r_max in saturated regime";
  EXPECT_NEAR(forcing[0][1], r_max, r_max * 1e-6);
}

// FD check: capped Jacobian (uncapped regime) — should match uncapped analytical
TEST(DissolvedReaction, JacobianFDCappedUncappedRegime)
{
  auto a = micm::Species{ "A" };
  auto b = micm::Species{ "B" };
  auto s = micm::Species{ "S" };
  auto phase = micm::Phase{ "AQUEOUS", { { a }, { b }, { s } } };

  double k = 1.0e-5;    // slow reaction, uncapped
  double t_half = 0.1;  // r_max = [A]/t_half >> r

  DissolvedReaction reaction{ [k](const micm::Conditions&) { return k; }, { a }, { b }, s, phase, 1.0e-20, t_half };

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");
  std::string k_param = phase.name_ + "." + reaction.uuid_ + ".k";
  std::unordered_map<std::string, std::size_t> spi{ { k_param, 0 } };
  std::unordered_map<std::string, std::size_t> svi{ { "DROP.AQUEOUS.A", 0 },
                                                    { "DROP.AQUEOUS.B", 1 },
                                                    { "DROP.AQUEOUS.S", 2 } };

  MatrixPolicy params(1, 1);
  params[0][0] = k;
  MatrixPolicy vars(1, 3);
  vars[0][0] = 1.0;
  vars[0][1] = 0.0;
  vars[0][2] = 55.0;

  CheckFiniteDifferenceJacobian(reaction, phase_prefixes, spi, svi, params, vars);
}

// FD check: capped Jacobian (saturated regime) — sech² ≈ 0, correction term dominates
TEST(DissolvedReaction, JacobianFDCappedSaturatedRegime)
{
  auto a = micm::Species{ "A" };
  auto b = micm::Species{ "B" };
  auto s = micm::Species{ "S" };
  auto phase = micm::Phase{ "AQUEOUS", { { a }, { b }, { s } } };

  double k = 1.0e5;  // fast reaction, deeply capped
  double t_half = 1000.0;

  DissolvedReaction reaction{ [k](const micm::Conditions&) { return k; }, { a }, { b }, s, phase, 1.0e-20, t_half };

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");
  std::string k_param = phase.name_ + "." + reaction.uuid_ + ".k";
  std::unordered_map<std::string, std::size_t> spi{ { k_param, 0 } };
  std::unordered_map<std::string, std::size_t> svi{ { "DROP.AQUEOUS.A", 0 },
                                                    { "DROP.AQUEOUS.B", 1 },
                                                    { "DROP.AQUEOUS.S", 2 } };

  MatrixPolicy params(1, 1);
  params[0][0] = k;
  MatrixPolicy vars(1, 3);
  vars[0][0] = 1.0;
  vars[0][1] = 0.0;
  vars[0][2] = 55.0;

  // In the saturated regime the Jacobian is dominated by the correction term.
  // FD uses h = max(|y|*1e-7, 1e-7), which is small enough that tanh is still
  // smooth. Tolerance is relaxed slightly because tanh curvature gives 2nd-order
  // FD errors proportional to h² * d³r/dy³.
  CheckFiniteDifferenceJacobian(reaction, phase_prefixes, spi, svi, params, vars);
}

// FD check: capped bimolecular — both reactant partial derivatives exercised
TEST(DissolvedReaction, JacobianFDCappedBimolecular)
{
  auto a = micm::Species{ "A" };
  auto b = micm::Species{ "B" };
  auto c = micm::Species{ "C" };
  auto s = micm::Species{ "S" };
  auto phase = micm::Phase{ "AQUEOUS", { { a }, { b }, { c }, { s } } };

  double k = 50.0;       // moderately fast
  double t_half = 10.0;  // r_max = min(A,B) / t_half

  DissolvedReaction reaction{ [k](const micm::Conditions&) { return k; }, { a, b }, { c }, s, phase, 1.0e-20, t_half };

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");
  std::string k_param = phase.name_ + "." + reaction.uuid_ + ".k";
  std::unordered_map<std::string, std::size_t> spi{ { k_param, 0 } };
  std::unordered_map<std::string, std::size_t> svi{
    { "DROP.AQUEOUS.A", 0 }, { "DROP.AQUEOUS.B", 1 }, { "DROP.AQUEOUS.C", 2 }, { "DROP.AQUEOUS.S", 3 }
  };

  MatrixPolicy params(1, 1);
  params[0][0] = k;
  MatrixPolicy vars(1, 4);
  // r = k/[S]*[A]*[B] = 50/50 * 0.01 * 0.02 = 2e-4; r_max = min(A,B)/t_half = 1e-3 → mild cap
  vars[0][0] = 0.01;
  vars[0][1] = 0.02;
  vars[0][2] = 0.0;
  vars[0][3] = 50.0;

  CheckFiniteDifferenceJacobian(reaction, phase_prefixes, spi, svi, params, vars);
}
