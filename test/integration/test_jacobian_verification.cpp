// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Finite-difference Jacobian verification tests for all MIAM processes and constraints.
// Uses MICM's FiniteDifferenceJacobian / CompareJacobianToFiniteDifference utilities to
// verify that analytical Jacobian implementations match numerical approximations.

#include <miam/miam.hpp>
#include <miam/processes/constants/equilibrium_constant.hpp>
#include <miam/processes/constants/henrys_law_constant.hpp>
#include <micm/CPU.hpp>
#include <micm/util/jacobian_verification.hpp>

#include <gtest/gtest.h>

#include <cmath>

using namespace micm;
using namespace miam;

using DenseMatrix = micm::Matrix<double>;
using SparseMatrixFD = micm::SparseMatrix<double, micm::SparseMatrixStandardOrdering>;

// ═══════════════════════════════════════════════════════════════════════════════
// Helper: build state index maps from a Model, allocating variable and parameter
// indices in a deterministic order.
// ═══════════════════════════════════════════════════════════════════════════════
namespace
{
  struct IndexMaps
  {
    std::unordered_map<std::string, std::size_t> variable_indices;
    std::unordered_map<std::string, std::size_t> parameter_indices;
    std::size_t num_variables;
    std::size_t num_parameters;
  };

  IndexMaps BuildIndexMaps(const Model& model)
  {
    IndexMaps result;
    auto var_names = model.StateVariableNames();
    // Include species used by processes/constraints that aren't part of any
    // representation (e.g. gas-phase species used by HenryLawPhaseTransfer)
    auto species_used = model.SpeciesUsed();
    for (const auto& name : species_used)
      var_names.insert(name);
    auto param_names = model.StateParameterNames();

    std::size_t idx = 0;
    for (const auto& name : var_names)
      result.variable_indices[name] = idx++;
    result.num_variables = idx;

    idx = 0;
    for (const auto& name : param_names)
      result.parameter_indices[name] = idx++;
    result.num_parameters = idx;

    return result;
  }

  /// Helper to verify process Jacobian (forcing-based) for a given model
  void VerifyProcessJacobian(
      const Model& model,
      const IndexMaps& maps,
      const DenseMatrix& variables,
      const DenseMatrix& parameters,
      const std::vector<micm::Conditions>& conditions,
      double atol = 0,
      double rtol = 0)
  {
    const std::size_t num_blocks = variables.NumRows();
    const std::size_t num_species = maps.num_variables;

    // Update temperature-dependent parameters
    auto update_fn = model.UpdateStateParametersFunction<DenseMatrix>(maps.parameter_indices);
    DenseMatrix params_copy(parameters);
    update_fn(conditions, params_copy);

    // Build sparse Jacobian with declared sparsity
    auto nz_elements = model.NonZeroJacobianElements(maps.variable_indices);
    auto builder = SparseMatrixFD::Create(num_species).SetNumberOfBlocks(num_blocks).InitialValue(0.0);
    for (const auto& elem : nz_elements)
      builder = builder.WithElement(elem.first, elem.second);
    SparseMatrixFD analytical_jac(builder);

    // Compute analytical Jacobian
    auto jacobian_fn =
        model.JacobianFunction<DenseMatrix, SparseMatrixFD>(maps.parameter_indices, maps.variable_indices, analytical_jac);
    jacobian_fn(params_copy, variables, analytical_jac);

    // Compute finite-difference Jacobian
    auto forcing_fn = model.ForcingFunction<DenseMatrix>(maps.parameter_indices, maps.variable_indices);
    auto fd_wrapper = [&](const DenseMatrix& vars, DenseMatrix& forcing)
    { forcing_fn(params_copy, vars, forcing); };

    auto fd_jac = micm::FiniteDifferenceJacobian<DenseMatrix>(fd_wrapper, variables, num_species);

    // Compare
    auto comparison = (atol > 0 && rtol > 0)
        ? micm::CompareJacobianToFiniteDifference<DenseMatrix, SparseMatrixFD>(
              analytical_jac, fd_jac, num_species, atol, rtol)
        : micm::CompareJacobianToFiniteDifference<DenseMatrix, SparseMatrixFD>(
              analytical_jac, fd_jac, num_species);

    EXPECT_TRUE(comparison.passed) << "Process Jacobian mismatch: block=" << comparison.worst_block
                                   << " row=" << comparison.worst_row << " col=" << comparison.worst_col
                                   << " analytical=" << comparison.worst_analytical << " fd=" << comparison.worst_fd;

    // Check sparsity completeness
    auto sparsity =
        micm::CheckJacobianSparsityCompleteness<DenseMatrix, SparseMatrixFD>(analytical_jac, fd_jac, num_species);

    EXPECT_TRUE(sparsity.passed) << "Missing sparsity at block=" << sparsity.worst_block
                                 << " row=" << sparsity.worst_row << " col=" << sparsity.worst_col
                                 << " fd_value=" << sparsity.worst_fd;
  }

  /// Helper to verify constraint Jacobian (residual-based) for a given model
  void VerifyConstraintJacobian(
      const Model& model,
      const IndexMaps& maps,
      const DenseMatrix& variables,
      const DenseMatrix& parameters,
      const std::vector<micm::Conditions>& conditions,
      double atol = 0,
      double rtol = 0)
  {
    const std::size_t num_blocks = variables.NumRows();
    const std::size_t num_species = maps.num_variables;

    // Update temperature-dependent constraint parameters
    auto update_fn = model.UpdateStateParametersFunction<DenseMatrix>(maps.parameter_indices);
    DenseMatrix params_copy(parameters);
    update_fn(conditions, params_copy);
    auto constraint_update_fn = model.ConstraintUpdateStateParametersFunction<DenseMatrix>(maps.parameter_indices);
    constraint_update_fn(conditions, params_copy);

    // Build sparse Jacobian with declared sparsity
    auto nz_elements = model.NonZeroConstraintJacobianElements(maps.variable_indices);
    auto builder = SparseMatrixFD::Create(num_species).SetNumberOfBlocks(num_blocks).InitialValue(0.0);
    for (const auto& elem : nz_elements)
      builder = builder.WithElement(elem.first, elem.second);
    SparseMatrixFD analytical_jac(builder);

    // Compute analytical Jacobian
    auto jac_fn =
        model.ConstraintJacobianFunction<DenseMatrix, SparseMatrixFD>(maps.parameter_indices, maps.variable_indices, analytical_jac);
    jac_fn(variables, DenseMatrix(variables.NumRows(), std::max(maps.num_parameters, std::size_t(1)), 0.0), analytical_jac);

    // Compute finite-difference Jacobian from residual function
    auto residual_fn = model.ConstraintResidualFunction<DenseMatrix>(maps.parameter_indices, maps.variable_indices);
    auto fd_wrapper = [&](const DenseMatrix& vars, DenseMatrix& forcing)
    { residual_fn(vars, DenseMatrix(vars.NumRows(), std::max(maps.num_parameters, std::size_t(1)), 0.0), forcing); };

    auto fd_jac = micm::FiniteDifferenceJacobian<DenseMatrix>(fd_wrapper, variables, num_species);

    // Compare
    auto comparison = (atol > 0 && rtol > 0)
        ? micm::CompareJacobianToFiniteDifference<DenseMatrix, SparseMatrixFD>(
              analytical_jac, fd_jac, num_species, atol, rtol)
        : micm::CompareJacobianToFiniteDifference<DenseMatrix, SparseMatrixFD>(
              analytical_jac, fd_jac, num_species);

    EXPECT_TRUE(comparison.passed) << "Constraint Jacobian mismatch: block=" << comparison.worst_block
                                   << " row=" << comparison.worst_row << " col=" << comparison.worst_col
                                   << " analytical=" << comparison.worst_analytical << " fd=" << comparison.worst_fd;

    // Check sparsity completeness
    auto sparsity =
        micm::CheckJacobianSparsityCompleteness<DenseMatrix, SparseMatrixFD>(analytical_jac, fd_jac, num_species);

    EXPECT_TRUE(sparsity.passed) << "Missing sparsity at block=" << sparsity.worst_block
                                 << " row=" << sparsity.worst_row << " col=" << sparsity.worst_col
                                 << " fd_value=" << sparsity.worst_fd;
  }
}  // namespace

// ═══════════════════════════════════════════════════════════════════════════════
// Process Jacobian Tests
// ═══════════════════════════════════════════════════════════════════════════════

/// @brief DissolvedReaction process Jacobian: A → B with solvent C
TEST(JacobianVerification, DissolvedReactionProcess)
{
  auto A = Species{ "A" };
  auto B = Species{ "B" };
  auto C = Species{ "C" };

  auto aqueous_phase = Phase{ "AQUEOUS", { { A }, { B }, { C } } };
  auto droplet = representation::UniformSection{ "DROPLET", { aqueous_phase } };

  double k = 0.1;
  auto rate = [k](const Conditions&) { return k; };
  auto reaction = process::DissolvedReaction{ rate, { A }, { B }, C, aqueous_phase };

  auto model = Model{ .name_ = "AEROSOL", .representations_ = { droplet } };
  model.AddProcesses({ reaction });

  auto maps = BuildIndexMaps(model);

  // Two blocks with different concentrations
  DenseMatrix variables(2, maps.num_variables, 0.0);
  variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.A")] = 0.8;
  variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.B")] = 0.1;
  variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.C")] = 1.0e-4;
  variables[1][maps.variable_indices.at("DROPLET.AQUEOUS.A")] = 0.3;
  variables[1][maps.variable_indices.at("DROPLET.AQUEOUS.B")] = 0.5;
  variables[1][maps.variable_indices.at("DROPLET.AQUEOUS.C")] = 2.0e-4;

  DenseMatrix parameters(2, std::max(maps.num_parameters, std::size_t(1)), 0.0);
  std::vector<Conditions> conditions(2);
  conditions[0].temperature_ = 298.15;
  conditions[0].pressure_ = 101325.0;
  conditions[1].temperature_ = 310.0;
  conditions[1].pressure_ = 101325.0;

  VerifyProcessJacobian(model, maps, variables, parameters, conditions);
}

/// @brief DissolvedReversibleReaction process Jacobian: A ⇌ B with solvent C
TEST(JacobianVerification, DissolvedReversibleReactionProcess)
{
  auto A = Species{ "A" };
  auto B = Species{ "B" };
  auto C = Species{ "C" };

  auto aqueous_phase = Phase{ "AQUEOUS", { { A }, { B }, { C } } };
  auto droplet = representation::UniformSection{ "DROPLET", { aqueous_phase } };

  double k_f = 0.1, k_r = 0.05;
  auto forward_rate = [k_f](const Conditions&) { return k_f; };
  auto reverse_rate = [k_r](const Conditions&) { return k_r; };
  auto reaction = process::DissolvedReversibleReaction{
    forward_rate, reverse_rate, { A }, { B }, C, aqueous_phase
  };

  auto model = Model{ .name_ = "AEROSOL", .representations_ = { droplet } };
  model.AddProcesses({ reaction });

  auto maps = BuildIndexMaps(model);

  DenseMatrix variables(2, maps.num_variables, 0.0);
  variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.A")] = 0.6;
  variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.B")] = 0.3;
  variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.C")] = 1.0e-4;
  variables[1][maps.variable_indices.at("DROPLET.AQUEOUS.A")] = 0.2;
  variables[1][maps.variable_indices.at("DROPLET.AQUEOUS.B")] = 0.7;
  variables[1][maps.variable_indices.at("DROPLET.AQUEOUS.C")] = 3.0e-4;

  DenseMatrix parameters(2, std::max(maps.num_parameters, std::size_t(1)), 0.0);
  std::vector<Conditions> conditions(2);
  conditions[0].temperature_ = 298.15;
  conditions[0].pressure_ = 101325.0;
  conditions[1].temperature_ = 280.0;
  conditions[1].pressure_ = 101325.0;

  VerifyProcessJacobian(model, maps, variables, parameters, conditions);
}

/// @brief HenryLawPhaseTransfer process Jacobian with SingleMomentMode
TEST(JacobianVerification, HenryLawPhaseTransferProcess)
{
  double Mw_gas = 0.044;
  double Mw_solvent = 0.018;
  double rho_solvent = 1000.0;
  double D_g = 1.5e-5;
  double alpha = 0.05;
  double HLC_val = 3.4e-2;

  auto A_g = Species{ "A_g", { { "molecular weight [kg mol-1]", Mw_gas } } };
  auto A_aq = Species{ "A_aq", { { "molecular weight [kg mol-1]", Mw_gas }, { "density [kg m-3]", 1800.0 } } };
  auto H2O = Species{ "H2O", { { "molecular weight [kg mol-1]", Mw_solvent }, { "density [kg m-3]", rho_solvent } } };

  Phase gas_phase{ "GAS", { { A_g } } };
  Phase aqueous_phase{ "AQUEOUS", { { A_aq }, { H2O } } };

  auto droplet = representation::SingleMomentMode{ "DROPLET", { aqueous_phase }, 5.0e-6, 1.2 };

  auto transfer = process::HenryLawPhaseTransferBuilder()
      .SetCondensedPhase(aqueous_phase)
      .SetGasSpecies(A_g)
      .SetCondensedSpecies(A_aq)
      .SetSolvent(H2O)
      .SetHenrysLawConstant(process::constant::HenrysLawConstant(
          process::constant::HenrysLawConstantParameters{ .HLC_ref_ = HLC_val }))
      .SetDiffusionCoefficient(D_g)
      .SetAccommodationCoefficient(alpha)
      .Build();

  auto model = Model{ .name_ = "AEROSOL", .representations_ = { droplet } };
  model.AddProcesses({ transfer });

  auto maps = BuildIndexMaps(model);

  DenseMatrix variables(2, maps.num_variables, 0.0);
  variables[0][maps.variable_indices.at("A_g")] = 1.0e-3;
  variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.A_aq")] = 1.0e-5;
  variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.H2O")] = 55.5;
  variables[1][maps.variable_indices.at("A_g")] = 5.0e-4;
  variables[1][maps.variable_indices.at("DROPLET.AQUEOUS.A_aq")] = 2.0e-4;
  variables[1][maps.variable_indices.at("DROPLET.AQUEOUS.H2O")] = 55.5;

  DenseMatrix parameters(2, std::max(maps.num_parameters, std::size_t(1)), 0.0);
  std::vector<Conditions> conditions(2);
  conditions[0].temperature_ = 298.15;
  conditions[0].pressure_ = 101325.0;
  conditions[1].temperature_ = 280.0;
  conditions[1].pressure_ = 101325.0;

  // Use slightly relaxed tolerances due to the multi-step rate computation
  VerifyProcessJacobian(model, maps, variables, parameters, conditions, 1.0e-4, 1.0e-3);
}

/// @brief HenryLawPhaseTransfer with TwoMomentMode (multi-instance)
TEST(JacobianVerification, HenryLawPhaseTransferTwoMomentMode)
{
  double Mw_gas = 0.044;
  double Mw_solvent = 0.018;
  double rho_solvent = 1000.0;
  double D_g = 1.5e-5;
  double alpha = 0.05;
  double HLC_val = 3.4e-2;

  auto A_g = Species{ "A_g", { { "molecular weight [kg mol-1]", Mw_gas } } };
  auto A_aq = Species{ "A_aq", { { "molecular weight [kg mol-1]", Mw_gas }, { "density [kg m-3]", 1800.0 } } };
  auto H2O = Species{ "H2O", { { "molecular weight [kg mol-1]", Mw_solvent }, { "density [kg m-3]", rho_solvent } } };

  Phase gas_phase{ "GAS", { { A_g } } };
  Phase aqueous_phase{ "AQUEOUS", { { A_aq }, { H2O } } };

  auto droplet = representation::TwoMomentMode{
    "DROPLET", { aqueous_phase }, 1.2
  };

  auto transfer = process::HenryLawPhaseTransferBuilder()
      .SetCondensedPhase(aqueous_phase)
      .SetGasSpecies(A_g)
      .SetCondensedSpecies(A_aq)
      .SetSolvent(H2O)
      .SetHenrysLawConstant(process::constant::HenrysLawConstant(
          process::constant::HenrysLawConstantParameters{ .HLC_ref_ = HLC_val }))
      .SetDiffusionCoefficient(D_g)
      .SetAccommodationCoefficient(alpha)
      .Build();

  auto model = Model{ .name_ = "AEROSOL", .representations_ = { droplet } };
  model.AddProcesses({ transfer });

  auto maps = BuildIndexMaps(model);

  DenseMatrix variables(1, maps.num_variables, 0.0);
  variables[0][maps.variable_indices.at("A_g")] = 1.0e-3;
  variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.A_aq")] = 1.0e-5;
  variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.H2O")] = 55.5;
  // TwoMomentMode also has a number concentration variable
  for (const auto& [name, idx] : maps.variable_indices)
  {
    if (name.find("NumberConcentration") != std::string::npos)
      variables[0][idx] = 1.0e7;
  }

  DenseMatrix parameters(1, std::max(maps.num_parameters, std::size_t(1)), 0.0);
  std::vector<Conditions> conditions(1);
  conditions[0].temperature_ = 298.15;
  conditions[0].pressure_ = 101325.0;

  VerifyProcessJacobian(model, maps, variables, parameters, conditions, 1.0e-4, 1.0e-3);
}

/// @brief Multiple processes in a single model: DissolvedReaction + DissolvedReversibleReaction
TEST(JacobianVerification, MultipleProcessesCombined)
{
  auto A = Species{ "A" };
  auto B = Species{ "B" };
  auto C = Species{ "C" };
  auto D = Species{ "D" };
  auto S = Species{ "S" };

  auto aqueous_phase = Phase{ "AQUEOUS", { { A }, { B }, { C }, { D }, { S } } };
  auto droplet = representation::UniformSection{ "DROPLET", { aqueous_phase } };

  // A → B (irreversible)
  auto reaction1 = process::DissolvedReaction{
    [](const Conditions&) { return 0.1; }, { A }, { B }, S, aqueous_phase
  };

  // C ⇌ D (reversible)
  auto reaction2 = process::DissolvedReversibleReaction{
    [](const Conditions&) { return 0.2; },
    [](const Conditions&) { return 0.05; },
    { C }, { D }, S, aqueous_phase
  };

  auto model = Model{ .name_ = "AEROSOL", .representations_ = { droplet } };
  model.AddProcesses(reaction1, reaction2);

  auto maps = BuildIndexMaps(model);

  DenseMatrix variables(2, maps.num_variables, 0.0);
  variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.A")] = 0.5;
  variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.B")] = 0.3;
  variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.C")] = 0.4;
  variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.D")] = 0.2;
  variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.S")] = 1.0e-4;
  variables[1][maps.variable_indices.at("DROPLET.AQUEOUS.A")] = 0.1;
  variables[1][maps.variable_indices.at("DROPLET.AQUEOUS.B")] = 0.8;
  variables[1][maps.variable_indices.at("DROPLET.AQUEOUS.C")] = 0.9;
  variables[1][maps.variable_indices.at("DROPLET.AQUEOUS.D")] = 0.1;
  variables[1][maps.variable_indices.at("DROPLET.AQUEOUS.S")] = 2.0e-4;

  DenseMatrix parameters(2, std::max(maps.num_parameters, std::size_t(1)), 0.0);
  std::vector<Conditions> conditions(2);
  conditions[0].temperature_ = 298.15;
  conditions[0].pressure_ = 101325.0;
  conditions[1].temperature_ = 273.15;
  conditions[1].pressure_ = 101325.0;

  VerifyProcessJacobian(model, maps, variables, parameters, conditions);
}

// ═══════════════════════════════════════════════════════════════════════════════
// Constraint Jacobian Tests
// ═══════════════════════════════════════════════════════════════════════════════

/// @brief DissolvedEquilibriumConstraint: B ⇌ C (K_eq), C algebraic
TEST(JacobianVerification, DissolvedEquilibriumConstraint)
{
  auto A = Species{ "A" };
  auto B = Species{ "B" };
  auto C = Species{ "C" };
  auto S = Species{ "S" };

  auto aqueous_phase = Phase{ "AQUEOUS", { { A }, { B }, { C }, { S } } };
  auto droplet = representation::UniformSection{ "DROPLET", { aqueous_phase } };

  double K_eq = 2.0;
  auto equil = constraint::DissolvedEquilibriumConstraintBuilder()
      .SetPhase(aqueous_phase)
      .SetReactants({ B })
      .SetProducts({ C })
      .SetAlgebraicSpecies(C)
      .SetSolvent(S)
      .SetEquilibriumConstant(process::constant::EquilibriumConstant(
          process::constant::EquilibriumConstantParameters{ .A_ = K_eq }))
      .Build();

  auto model = Model{ .name_ = "AEROSOL", .representations_ = { droplet } };
  model.AddConstraints({ equil });

  auto maps = BuildIndexMaps(model);

  DenseMatrix variables(2, maps.num_variables, 0.0);
  variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.A")] = 0.3;
  variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.B")] = 0.4;
  variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.C")] = 0.2;
  variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.S")] = 1.0e-4;
  variables[1][maps.variable_indices.at("DROPLET.AQUEOUS.A")] = 0.1;
  variables[1][maps.variable_indices.at("DROPLET.AQUEOUS.B")] = 0.6;
  variables[1][maps.variable_indices.at("DROPLET.AQUEOUS.C")] = 0.8;
  variables[1][maps.variable_indices.at("DROPLET.AQUEOUS.S")] = 2.0e-4;

  DenseMatrix parameters(2, std::max(maps.num_parameters, std::size_t(1)), 0.0);
  std::vector<Conditions> conditions(2);
  conditions[0].temperature_ = 298.15;
  conditions[0].pressure_ = 101325.0;
  conditions[1].temperature_ = 310.0;
  conditions[1].pressure_ = 101325.0;

  VerifyConstraintJacobian(model, maps, variables, parameters, conditions);
}

/// @brief LinearConstraint: [A] + [B] + [C] = total, B algebraic
TEST(JacobianVerification, LinearConstraint)
{
  auto A = Species{ "A" };
  auto B = Species{ "B" };
  auto C = Species{ "C" };

  auto aqueous_phase = Phase{ "AQUEOUS", { { A }, { B }, { C } } };
  auto droplet = representation::UniformSection{ "DROPLET", { aqueous_phase } };

  double total = 1.0;
  auto mass_cons = constraint::LinearConstraintBuilder()
      .SetAlgebraicSpecies(aqueous_phase, B)
      .AddTerm(aqueous_phase, A, 1.0)
      .AddTerm(aqueous_phase, B, 1.0)
      .AddTerm(aqueous_phase, C, 1.0)
      .SetConstant(total)
      .Build();

  auto model = Model{ .name_ = "AEROSOL", .representations_ = { droplet } };
  model.AddConstraints({ mass_cons });

  auto maps = BuildIndexMaps(model);

  DenseMatrix variables(2, maps.num_variables, 0.0);
  variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.A")] = 0.3;
  variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.B")] = 0.5;
  variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.C")] = 0.2;
  variables[1][maps.variable_indices.at("DROPLET.AQUEOUS.A")] = 0.1;
  variables[1][maps.variable_indices.at("DROPLET.AQUEOUS.B")] = 0.7;
  variables[1][maps.variable_indices.at("DROPLET.AQUEOUS.C")] = 0.2;

  DenseMatrix parameters(2, std::max(maps.num_parameters, std::size_t(1)), 0.0);
  std::vector<Conditions> conditions(2);
  conditions[0].temperature_ = 298.15;
  conditions[0].pressure_ = 101325.0;
  conditions[1].temperature_ = 298.15;
  conditions[1].pressure_ = 101325.0;

  VerifyConstraintJacobian(model, maps, variables, parameters, conditions);
}

/// @brief HenryLawEquilibriumConstraint: A_g ⇌ A_aq, A_aq algebraic
TEST(JacobianVerification, HenryLawEquilibriumConstraint)
{
  double Mw_solvent = 0.018;
  double rho_solvent = 1000.0;
  double HLC = 4.0e-4;

  auto A_g = Species{ "A_g" };
  auto A_aq = Species{ "A_aq" };
  auto H2O = Species{ "H2O",
      { { "molecular weight [kg mol-1]", Mw_solvent },
        { "density [kg m-3]", rho_solvent } } };

  Phase gas_phase{ "GAS", { { A_g } } };
  Phase aqueous_phase{ "AQUEOUS", { { A_aq }, { H2O } } };

  auto droplet = representation::UniformSection{ "DROPLET", { aqueous_phase } };

  auto hl_constraint = constraint::HenryLawEquilibriumConstraintBuilder()
      .SetGasSpecies(A_g)
      .SetCondensedSpecies(A_aq)
      .SetSolvent(H2O)
      .SetCondensedPhase(aqueous_phase)
      .SetHenryLawConstant(process::constant::HenrysLawConstant(
          process::constant::HenrysLawConstantParameters{ .HLC_ref_ = HLC }))
      .SetMwSolvent(Mw_solvent)
      .SetRhoSolvent(rho_solvent)
      .Build();

  auto model = Model{ .name_ = "AEROSOL", .representations_ = { droplet } };
  model.AddConstraints({ hl_constraint });

  auto maps = BuildIndexMaps(model);

  DenseMatrix variables(2, maps.num_variables, 0.0);
  variables[0][maps.variable_indices.at("A_g")] = 1.0e-3;
  variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.A_aq")] = 5.0e-4;
  variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.H2O")] = 55555.0;
  variables[1][maps.variable_indices.at("A_g")] = 2.0e-3;
  variables[1][maps.variable_indices.at("DROPLET.AQUEOUS.A_aq")] = 1.0e-3;
  variables[1][maps.variable_indices.at("DROPLET.AQUEOUS.H2O")] = 55555.0;

  DenseMatrix parameters(2, std::max(maps.num_parameters, std::size_t(1)), 0.0);
  std::vector<Conditions> conditions(2);
  conditions[0].temperature_ = 298.15;
  conditions[0].pressure_ = 101325.0;
  conditions[1].temperature_ = 280.0;
  conditions[1].pressure_ = 101325.0;

  VerifyConstraintJacobian(model, maps, variables, parameters, conditions);
}

// ═══════════════════════════════════════════════════════════════════════════════
// Combined Process + Constraint Jacobian Tests
// ═══════════════════════════════════════════════════════════════════════════════

/// @brief DissolvedReaction + DissolvedEquilibriumConstraint + LinearConstraint
///        A → B (kinetic), B ⇌ C (equilibrium), [A]+[B]+[C] = total (conservation)
TEST(JacobianVerification, ProcessAndConstraintsCombined)
{
  auto A = Species{ "A" };
  auto B = Species{ "B" };
  auto C = Species{ "C" };
  auto S = Species{ "S" };

  auto aqueous_phase = Phase{ "AQUEOUS", { { A }, { B }, { C }, { S } } };
  auto droplet = representation::UniformSection{ "DROPLET", { aqueous_phase } };

  double k = 0.1;
  double K_eq = 2.0;
  double total = 1.0;

  auto reaction = process::DissolvedReaction{
    [k](const Conditions&) { return k; }, { A }, { B }, S, aqueous_phase
  };

  auto equil = constraint::DissolvedEquilibriumConstraintBuilder()
      .SetPhase(aqueous_phase)
      .SetReactants({ B })
      .SetProducts({ C })
      .SetAlgebraicSpecies(C)
      .SetSolvent(S)
      .SetEquilibriumConstant(process::constant::EquilibriumConstant(
          process::constant::EquilibriumConstantParameters{ .A_ = K_eq }))
      .Build();

  auto mass_cons = constraint::LinearConstraintBuilder()
      .SetAlgebraicSpecies(aqueous_phase, B)
      .AddTerm(aqueous_phase, A, 1.0)
      .AddTerm(aqueous_phase, B, 1.0)
      .AddTerm(aqueous_phase, C, 1.0)
      .SetConstant(total)
      .Build();

  auto model = Model{ .name_ = "AEROSOL", .representations_ = { droplet } };
  model.AddProcesses({ reaction });
  model.AddConstraints(equil, mass_cons);

  auto maps = BuildIndexMaps(model);

  DenseMatrix variables(2, maps.num_variables, 0.0);
  variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.A")] = 0.5;
  variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.B")] = 0.3;
  variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.C")] = 0.2;
  variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.S")] = 1.0e-4;
  variables[1][maps.variable_indices.at("DROPLET.AQUEOUS.A")] = 0.1;
  variables[1][maps.variable_indices.at("DROPLET.AQUEOUS.B")] = 0.4;
  variables[1][maps.variable_indices.at("DROPLET.AQUEOUS.C")] = 0.5;
  variables[1][maps.variable_indices.at("DROPLET.AQUEOUS.S")] = 2.0e-4;

  DenseMatrix parameters(2, std::max(maps.num_parameters, std::size_t(1)), 0.0);
  std::vector<Conditions> conditions(2);
  conditions[0].temperature_ = 298.15;
  conditions[0].pressure_ = 101325.0;
  conditions[1].temperature_ = 310.0;
  conditions[1].pressure_ = 101325.0;

  // Verify both process and constraint Jacobians
  VerifyProcessJacobian(model, maps, variables, parameters, conditions);
  VerifyConstraintJacobian(model, maps, variables, parameters, conditions);
}

/// @brief HenryLawEquilibriumConstraint + LinearConstraint with gas-phase species
TEST(JacobianVerification, HenryLawEquilibriumWithConservation)
{
  double Mw_solvent = 0.018;
  double rho_solvent = 1000.0;
  double HLC = 4.0e-4;

  auto Precursor = Species{ "Precursor" };
  auto A_g = Species{ "A_g" };
  auto A_aq = Species{ "A_aq" };
  auto H2O = Species{ "H2O",
      { { "molecular weight [kg mol-1]", Mw_solvent },
        { "density [kg m-3]", rho_solvent } } };

  Phase gas_phase{ "GAS", { { Precursor }, { A_g } } };
  Phase aqueous_phase{ "AQUEOUS", { { A_aq }, { H2O } } };

  auto droplet = representation::UniformSection{ "DROPLET", { aqueous_phase } };

  auto hl_constraint = constraint::HenryLawEquilibriumConstraintBuilder()
      .SetGasSpecies(A_g)
      .SetCondensedSpecies(A_aq)
      .SetSolvent(H2O)
      .SetCondensedPhase(aqueous_phase)
      .SetHenryLawConstant(process::constant::HenrysLawConstant(
          process::constant::HenrysLawConstantParameters{ .HLC_ref_ = HLC }))
      .SetMwSolvent(Mw_solvent)
      .SetRhoSolvent(rho_solvent)
      .Build();

  double total = 1.0;
  auto mass_cons = constraint::LinearConstraintBuilder()
      .SetAlgebraicSpecies(gas_phase, A_g)
      .AddTerm(gas_phase, Precursor, 1.0)
      .AddTerm(gas_phase, A_g, 1.0)
      .AddTerm(aqueous_phase, A_aq, 1.0)
      .SetConstant(total)
      .Build();

  auto model = Model{ .name_ = "AEROSOL", .representations_ = { droplet } };
  model.AddConstraints(hl_constraint, mass_cons);

  auto maps = BuildIndexMaps(model);

  DenseMatrix variables(2, maps.num_variables, 0.0);
  variables[0][maps.variable_indices.at("Precursor")] = 0.8;
  variables[0][maps.variable_indices.at("A_g")] = 0.1;
  variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.A_aq")] = 0.1;
  variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.H2O")] = 55555.0;
  variables[1][maps.variable_indices.at("Precursor")] = 0.3;
  variables[1][maps.variable_indices.at("A_g")] = 0.4;
  variables[1][maps.variable_indices.at("DROPLET.AQUEOUS.A_aq")] = 0.3;
  variables[1][maps.variable_indices.at("DROPLET.AQUEOUS.H2O")] = 55555.0;

  DenseMatrix parameters(2, std::max(maps.num_parameters, std::size_t(1)), 0.0);
  std::vector<Conditions> conditions(2);
  conditions[0].temperature_ = 298.15;
  conditions[0].pressure_ = 101325.0;
  conditions[1].temperature_ = 280.0;
  conditions[1].pressure_ = 101325.0;

  VerifyConstraintJacobian(model, maps, variables, parameters, conditions);
}
