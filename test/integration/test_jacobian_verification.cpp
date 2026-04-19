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
  variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.H2O")] = 0.017;  // mol/m³ air
  variables[1][maps.variable_indices.at("A_g")] = 5.0e-4;
  variables[1][maps.variable_indices.at("DROPLET.AQUEOUS.A_aq")] = 2.0e-4;
  variables[1][maps.variable_indices.at("DROPLET.AQUEOUS.H2O")] = 0.017;  // mol/m³ air

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
  variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.H2O")] = 0.017;  // mol/m³ air
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
  variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.H2O")] = 0.017;  // mol/m³ air
  variables[1][maps.variable_indices.at("A_g")] = 2.0e-3;
  variables[1][maps.variable_indices.at("DROPLET.AQUEOUS.A_aq")] = 1.0e-3;
  variables[1][maps.variable_indices.at("DROPLET.AQUEOUS.H2O")] = 0.017;  // mol/m³ air

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
  variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.H2O")] = 0.017;  // mol/m³ air
  variables[1][maps.variable_indices.at("Precursor")] = 0.3;
  variables[1][maps.variable_indices.at("A_g")] = 0.4;
  variables[1][maps.variable_indices.at("DROPLET.AQUEOUS.A_aq")] = 0.3;
  variables[1][maps.variable_indices.at("DROPLET.AQUEOUS.H2O")] = 0.017;  // mol/m³ air

  DenseMatrix parameters(2, std::max(maps.num_parameters, std::size_t(1)), 0.0);
  std::vector<Conditions> conditions(2);
  conditions[0].temperature_ = 298.15;
  conditions[0].pressure_ = 101325.0;
  conditions[1].temperature_ = 280.0;
  conditions[1].pressure_ = 101325.0;

  VerifyConstraintJacobian(model, maps, variables, parameters, conditions);
}

// ═══════════════════════════════════════════════════════════════════════════════
// Solvent Damping Range Tests
// ═══════════════════════════════════════════════════════════════════════════════

/// @brief Verify DissolvedReaction Jacobian across the damping range (solvent from 1e-4 down to 0)
TEST(JacobianVerification, DissolvedReactionDampingRange)
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

  // FD verification at solvent levels well above eps (1e-10) where central differences are accurate
  for (double sol : { 55.0, 1.0, 1.0e-2, 1.0e-4 })
  {
    SCOPED_TRACE("solvent = " + std::to_string(sol));

    DenseMatrix variables(1, maps.num_variables, 0.0);
    variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.A")] = 0.5;
    variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.B")] = 0.3;
    variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.C")] = sol;

    DenseMatrix parameters(1, std::max(maps.num_parameters, std::size_t(1)), 0.0);
    std::vector<Conditions> conditions(1);
    conditions[0].temperature_ = 298.15;
    conditions[0].pressure_ = 101325.0;

    VerifyProcessJacobian(model, maps, variables, parameters, conditions);
  }

  // At extreme low solvent (near/below eps), FD can't resolve the steep damping gradient.
  // Verify finiteness instead.
  auto forcing_fn = model.ForcingFunction<DenseMatrix>(maps.parameter_indices, maps.variable_indices);
  auto jac_nz = model.NonZeroJacobianElements(maps.variable_indices);
  auto jac_builder = SparseMatrixFD::Create(maps.num_variables).SetNumberOfBlocks(1).InitialValue(0.0);
  for (const auto& elem : jac_nz)
    jac_builder = jac_builder.WithElement(elem.first, elem.second);
  SparseMatrixFD jac(jac_builder);
  auto jac_fn = model.JacobianFunction<DenseMatrix, SparseMatrixFD>(
      maps.parameter_indices, maps.variable_indices, jac);

  for (double sol : { 1.0e-10, 1.0e-15, 0.0 })
  {
    SCOPED_TRACE("solvent (finiteness) = " + std::to_string(sol));

    DenseMatrix variables(1, maps.num_variables, 0.0);
    variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.A")] = 0.5;
    variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.B")] = 0.3;
    variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.C")] = sol;

    DenseMatrix parameters(1, std::max(maps.num_parameters, std::size_t(1)), 0.0);
    std::vector<Conditions> conditions(1);
    conditions[0].temperature_ = 298.15;
    conditions[0].pressure_ = 101325.0;

    auto update_fn = model.UpdateStateParametersFunction<DenseMatrix>(maps.parameter_indices);
    update_fn(conditions, parameters);

    DenseMatrix forcing(1, maps.num_variables, 0.0);
    forcing_fn(parameters, variables, forcing);
    for (std::size_t j = 0; j < maps.num_variables; ++j)
      EXPECT_TRUE(std::isfinite(forcing[0][j])) << "forcing[" << j << "] is not finite at sol=" << sol;

    jac.Fill(0.0);
    jac_fn(parameters, variables, jac);
    for (const auto& v : jac.AsVector())
      EXPECT_TRUE(std::isfinite(v)) << "Jacobian element is not finite at sol=" << sol;
  }
}

/// @brief Verify DissolvedReversibleReaction Jacobian across the damping range
TEST(JacobianVerification, DissolvedReversibleReactionDampingRange)
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

  // FD verification at solvent levels well above eps (1e-10)
  for (double sol : { 55.0, 1.0, 1.0e-2, 1.0e-4 })
  {
    SCOPED_TRACE("solvent = " + std::to_string(sol));

    DenseMatrix variables(1, maps.num_variables, 0.0);
    variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.A")] = 0.6;
    variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.B")] = 0.3;
    variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.C")] = sol;

    DenseMatrix parameters(1, std::max(maps.num_parameters, std::size_t(1)), 0.0);
    std::vector<Conditions> conditions(1);
    conditions[0].temperature_ = 298.15;
    conditions[0].pressure_ = 101325.0;

    VerifyProcessJacobian(model, maps, variables, parameters, conditions);
  }

  // At extreme low solvent, verify finiteness only
  auto forcing_fn = model.ForcingFunction<DenseMatrix>(maps.parameter_indices, maps.variable_indices);
  auto jac_nz = model.NonZeroJacobianElements(maps.variable_indices);
  auto jac_builder = SparseMatrixFD::Create(maps.num_variables).SetNumberOfBlocks(1).InitialValue(0.0);
  for (const auto& elem : jac_nz)
    jac_builder = jac_builder.WithElement(elem.first, elem.second);
  SparseMatrixFD jac(jac_builder);
  auto jac_fn = model.JacobianFunction<DenseMatrix, SparseMatrixFD>(
      maps.parameter_indices, maps.variable_indices, jac);

  for (double sol : { 1.0e-10, 1.0e-15, 0.0 })
  {
    SCOPED_TRACE("solvent (finiteness) = " + std::to_string(sol));

    DenseMatrix variables(1, maps.num_variables, 0.0);
    variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.A")] = 0.6;
    variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.B")] = 0.3;
    variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.C")] = sol;

    DenseMatrix parameters(1, std::max(maps.num_parameters, std::size_t(1)), 0.0);
    std::vector<Conditions> conditions(1);
    conditions[0].temperature_ = 298.15;
    conditions[0].pressure_ = 101325.0;

    auto update_fn = model.UpdateStateParametersFunction<DenseMatrix>(maps.parameter_indices);
    update_fn(conditions, parameters);

    DenseMatrix forcing(1, maps.num_variables, 0.0);
    forcing_fn(parameters, variables, forcing);
    for (std::size_t j = 0; j < maps.num_variables; ++j)
      EXPECT_TRUE(std::isfinite(forcing[0][j])) << "forcing[" << j << "] is not finite at sol=" << sol;

    jac.Fill(0.0);
    jac_fn(parameters, variables, jac);
    for (const auto& v : jac.AsVector())
      EXPECT_TRUE(std::isfinite(v)) << "Jacobian element is not finite at sol=" << sol;
  }
}

/// @brief Verify DissolvedEquilibriumConstraint Jacobian across the damping range
TEST(JacobianVerification, DissolvedEquilibriumConstraintDampingRange)
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

  // FD verification at solvent levels well above eps (1e-10)
  for (double sol : { 55.0, 1.0, 1.0e-2, 1.0e-4 })
  {
    SCOPED_TRACE("solvent = " + std::to_string(sol));

    DenseMatrix variables(1, maps.num_variables, 0.0);
    variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.A")] = 0.3;
    variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.B")] = 0.4;
    variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.C")] = 0.2;
    variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.S")] = sol;

    DenseMatrix parameters(1, std::max(maps.num_parameters, std::size_t(1)), 0.0);
    std::vector<Conditions> conditions(1);
    conditions[0].temperature_ = 298.15;
    conditions[0].pressure_ = 101325.0;

    VerifyConstraintJacobian(model, maps, variables, parameters, conditions);
  }

  // At extreme low solvent, verify finiteness only
  auto residual_fn = model.ConstraintResidualFunction<DenseMatrix>(maps.parameter_indices, maps.variable_indices);
  auto jac_nz = model.NonZeroConstraintJacobianElements(maps.variable_indices);
  auto jac_builder = SparseMatrixFD::Create(maps.num_variables).SetNumberOfBlocks(1).InitialValue(0.0);
  for (const auto& elem : jac_nz)
    jac_builder = jac_builder.WithElement(elem.first, elem.second);
  SparseMatrixFD jac(jac_builder);
  auto jac_fn = model.ConstraintJacobianFunction<DenseMatrix, SparseMatrixFD>(
      maps.parameter_indices, maps.variable_indices, jac);

  for (double sol : { 1.0e-10, 1.0e-15, 0.0 })
  {
    SCOPED_TRACE("solvent (finiteness) = " + std::to_string(sol));

    DenseMatrix variables(1, maps.num_variables, 0.0);
    variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.A")] = 0.3;
    variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.B")] = 0.4;
    variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.C")] = 0.2;
    variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.S")] = sol;

    DenseMatrix parameters(1, std::max(maps.num_parameters, std::size_t(1)), 0.0);
    std::vector<Conditions> conditions(1);
    conditions[0].temperature_ = 298.15;
    conditions[0].pressure_ = 101325.0;

    auto update_fn = model.ConstraintUpdateStateParametersFunction<DenseMatrix>(maps.parameter_indices);
    update_fn(conditions, parameters);

    DenseMatrix residual(1, maps.num_variables, 0.0);
    residual_fn(variables, parameters, residual);
    for (std::size_t j = 0; j < maps.num_variables; ++j)
      EXPECT_TRUE(std::isfinite(residual[0][j])) << "residual[" << j << "] is not finite at sol=" << sol;

    jac.Fill(0.0);
    jac_fn(variables, parameters, jac);
    for (const auto& v : jac.AsVector())
      EXPECT_TRUE(std::isfinite(v)) << "Constraint Jacobian element is not finite at sol=" << sol;
  }
}

/// @brief Verify combined process+constraint Jacobian at zero solvent
TEST(JacobianVerification, CombinedProcessAndConstraintZeroSolvent)
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

  // Test at zero solvent — should NOT produce NaN anymore
  DenseMatrix variables(1, maps.num_variables, 0.0);
  variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.A")] = 0.5;
  variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.B")] = 0.3;
  variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.C")] = 0.2;
  variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.S")] = 0.0;

  DenseMatrix parameters(1, std::max(maps.num_parameters, std::size_t(1)), 0.0);
  std::vector<Conditions> conditions(1);
  conditions[0].temperature_ = 298.15;
  conditions[0].pressure_ = 101325.0;

  // At zero solvent, FD can't resolve the steep damping gradient.
  // Verify finiteness of forcing/Jacobian instead of FD accuracy.
  auto update_fn = model.UpdateStateParametersFunction<DenseMatrix>(maps.parameter_indices);
  update_fn(conditions, parameters);
  auto cons_update_fn = model.ConstraintUpdateStateParametersFunction<DenseMatrix>(maps.parameter_indices);
  cons_update_fn(conditions, parameters);

  // Process forcing/Jacobian finiteness
  auto forcing_fn = model.ForcingFunction<DenseMatrix>(maps.parameter_indices, maps.variable_indices);
  DenseMatrix forcing(1, maps.num_variables, 0.0);
  forcing_fn(parameters, variables, forcing);
  for (std::size_t j = 0; j < maps.num_variables; ++j)
    EXPECT_TRUE(std::isfinite(forcing[0][j])) << "Process forcing[" << j << "] is not finite at sol=0";

  auto proc_nz = model.NonZeroJacobianElements(maps.variable_indices);
  auto proc_jac_builder = SparseMatrixFD::Create(maps.num_variables).SetNumberOfBlocks(1).InitialValue(0.0);
  for (const auto& elem : proc_nz)
    proc_jac_builder = proc_jac_builder.WithElement(elem.first, elem.second);
  SparseMatrixFD proc_jac(proc_jac_builder);
  auto proc_jac_fn = model.JacobianFunction<DenseMatrix, SparseMatrixFD>(
      maps.parameter_indices, maps.variable_indices, proc_jac);
  proc_jac_fn(parameters, variables, proc_jac);
  for (const auto& v : proc_jac.AsVector())
    EXPECT_TRUE(std::isfinite(v)) << "Process Jacobian element is not finite at sol=0";

  // Constraint residual/Jacobian finiteness
  auto residual_fn = model.ConstraintResidualFunction<DenseMatrix>(maps.parameter_indices, maps.variable_indices);
  DenseMatrix residual(1, maps.num_variables, 0.0);
  residual_fn(variables, parameters, residual);
  for (std::size_t j = 0; j < maps.num_variables; ++j)
    EXPECT_TRUE(std::isfinite(residual[0][j])) << "Constraint residual[" << j << "] is not finite at sol=0";

  auto cons_nz = model.NonZeroConstraintJacobianElements(maps.variable_indices);
  auto cons_jac_builder = SparseMatrixFD::Create(maps.num_variables).SetNumberOfBlocks(1).InitialValue(0.0);
  for (const auto& elem : cons_nz)
    cons_jac_builder = cons_jac_builder.WithElement(elem.first, elem.second);
  SparseMatrixFD cons_jac(cons_jac_builder);
  auto cons_jac_fn = model.ConstraintJacobianFunction<DenseMatrix, SparseMatrixFD>(
      maps.parameter_indices, maps.variable_indices, cons_jac);
  cons_jac_fn(variables, parameters, cons_jac);
  for (const auto& v : cons_jac.AsVector())
    EXPECT_TRUE(std::isfinite(v)) << "Constraint Jacobian element is not finite at sol=0";
}

// ═══════════════════════════════════════════════════════════════════════════════
// Rate-Capped Dissolved Reaction Jacobian Tests
// ═══════════════════════════════════════════════════════════════════════════════

/// @brief DissolvedReaction with max_halflife: single reactant, sweep rate regimes
TEST(JacobianVerification, DissolvedReactionCappedSingleReactant)
{
  auto A = Species{ "A" };
  auto B = Species{ "B" };
  auto C = Species{ "C" };

  auto aqueous_phase = Phase{ "AQUEOUS", { { A }, { B }, { C } } };
  auto droplet = representation::UniformSection{ "DROPLET", { aqueous_phase } };

  double t_half = 1.0;
  auto reaction = process::DissolvedReactionBuilder()
      .SetPhase(aqueous_phase)
      .SetReactants({ A })
      .SetProducts({ B })
      .SetSolvent(C)
      .SetRateConstant([](const Conditions&) { return 0.5; })
      .SetMaxHalflife(t_half)
      .Build();

  auto model = Model{ .name_ = "AEROSOL", .representations_ = { droplet } };
  model.AddProcesses({ reaction });

  auto maps = BuildIndexMaps(model);

  // Test across rate regimes:
  //   Low  k[A]: r << r_max  (uncapped, tanh(x)≈x)
  //   Mid  k[A]: r ≈  r_max  (transition region)
  //   High k[A]: r >> r_max  (heavily capped, tanh(x)≈1)
  for (double conc_A : { 0.01, 0.5, 1.0, 5.0, 100.0 })
  {
    SCOPED_TRACE("A = " + std::to_string(conc_A));

    DenseMatrix variables(1, maps.num_variables, 0.0);
    variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.A")] = conc_A;
    variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.B")] = 0.1;
    variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.C")] = 1.0;

    DenseMatrix parameters(1, std::max(maps.num_parameters, std::size_t(1)), 0.0);
    std::vector<Conditions> conditions(1);
    conditions[0].temperature_ = 298.15;
    conditions[0].pressure_ = 101325.0;

    VerifyProcessJacobian(model, maps, variables, parameters, conditions);
  }
}

/// @brief DissolvedReaction with max_halflife: two reactants, sweep rate regimes
TEST(JacobianVerification, DissolvedReactionCappedTwoReactants)
{
  auto A = Species{ "A" };
  auto B = Species{ "B" };
  auto P = Species{ "P" };
  auto S = Species{ "S" };

  auto aqueous_phase = Phase{ "AQUEOUS", { { A }, { B }, { P }, { S } } };
  auto droplet = representation::UniformSection{ "DROPLET", { aqueous_phase } };

  double t_half = 0.5;
  auto reaction = process::DissolvedReactionBuilder()
      .SetPhase(aqueous_phase)
      .SetReactants({ A, B })
      .SetProducts({ P })
      .SetSolvent(S)
      .SetRateConstant([](const Conditions&) { return 1.0; })
      .SetMaxHalflife(t_half)
      .Build();

  auto model = Model{ .name_ = "AEROSOL", .representations_ = { droplet } };
  model.AddProcesses({ reaction });

  auto maps = BuildIndexMaps(model);

  // Sweep both reactant concentrations to test soft-min transitions
  struct TestPoint { double A; double B; };
  for (auto [cA, cB] : std::vector<TestPoint>{{ 0.01, 0.01 }, { 0.5, 0.5 }, { 10.0, 0.01 }, { 0.01, 10.0 }, { 5.0, 5.0 }})
  {
    SCOPED_TRACE("A=" + std::to_string(cA) + " B=" + std::to_string(cB));

    DenseMatrix variables(1, maps.num_variables, 0.0);
    variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.A")] = cA;
    variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.B")] = cB;
    variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.P")] = 0.1;
    variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.S")] = 1.0;

    DenseMatrix parameters(1, std::max(maps.num_parameters, std::size_t(1)), 0.0);
    std::vector<Conditions> conditions(1);
    conditions[0].temperature_ = 298.15;
    conditions[0].pressure_ = 101325.0;

    VerifyProcessJacobian(model, maps, variables, parameters, conditions);
  }
}

/// @brief DissolvedReaction with max_halflife: vary solvent across damping range
TEST(JacobianVerification, DissolvedReactionCappedSolventRange)
{
  auto A = Species{ "A" };
  auto B = Species{ "B" };
  auto C = Species{ "C" };

  auto aqueous_phase = Phase{ "AQUEOUS", { { A }, { B }, { C } } };
  auto droplet = representation::UniformSection{ "DROPLET", { aqueous_phase } };

  auto reaction = process::DissolvedReactionBuilder()
      .SetPhase(aqueous_phase)
      .SetReactants({ A })
      .SetProducts({ B })
      .SetSolvent(C)
      .SetRateConstant([](const Conditions&) { return 1.0; })
      .SetMaxHalflife(1.0)
      .Build();

  auto model = Model{ .name_ = "AEROSOL", .representations_ = { droplet } };
  model.AddProcesses({ reaction });

  auto maps = BuildIndexMaps(model);

  // FD verification at solvent levels above eps where central differences work
  for (double sol : { 55.0, 1.0, 1.0e-2, 1.0e-4 })
  {
    SCOPED_TRACE("solvent = " + std::to_string(sol));

    DenseMatrix variables(1, maps.num_variables, 0.0);
    variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.A")] = 0.5;
    variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.B")] = 0.3;
    variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.C")] = sol;

    DenseMatrix parameters(1, std::max(maps.num_parameters, std::size_t(1)), 0.0);
    std::vector<Conditions> conditions(1);
    conditions[0].temperature_ = 298.15;
    conditions[0].pressure_ = 101325.0;

    VerifyProcessJacobian(model, maps, variables, parameters, conditions);
  }

  // At extreme low solvent, verify finiteness
  auto forcing_fn = model.ForcingFunction<DenseMatrix>(maps.parameter_indices, maps.variable_indices);
  auto jac_nz = model.NonZeroJacobianElements(maps.variable_indices);
  auto jac_builder = SparseMatrixFD::Create(maps.num_variables).SetNumberOfBlocks(1).InitialValue(0.0);
  for (const auto& elem : jac_nz)
    jac_builder = jac_builder.WithElement(elem.first, elem.second);
  SparseMatrixFD jac(jac_builder);
  auto jac_fn = model.JacobianFunction<DenseMatrix, SparseMatrixFD>(
      maps.parameter_indices, maps.variable_indices, jac);

  for (double sol : { 1.0e-10, 1.0e-15, 0.0 })
  {
    SCOPED_TRACE("solvent (finiteness) = " + std::to_string(sol));

    DenseMatrix variables(1, maps.num_variables, 0.0);
    variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.A")] = 0.5;
    variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.B")] = 0.3;
    variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.C")] = sol;

    DenseMatrix parameters(1, std::max(maps.num_parameters, std::size_t(1)), 0.0);
    std::vector<Conditions> conditions(1);
    conditions[0].temperature_ = 298.15;
    conditions[0].pressure_ = 101325.0;

    auto update_fn = model.UpdateStateParametersFunction<DenseMatrix>(maps.parameter_indices);
    update_fn(conditions, parameters);

    DenseMatrix forcing(1, maps.num_variables, 0.0);
    forcing_fn(parameters, variables, forcing);
    for (std::size_t j = 0; j < maps.num_variables; ++j)
      EXPECT_TRUE(std::isfinite(forcing[0][j])) << "forcing[" << j << "] is not finite at sol=" << sol;

    jac.Fill(0.0);
    jac_fn(parameters, variables, jac);
    for (const auto& v : jac.AsVector())
      EXPECT_TRUE(std::isfinite(v)) << "Jacobian element is not finite at sol=" << sol;
  }
}

/// @brief DissolvedReaction with max_halflife: multi-block (2 grid cells)
TEST(JacobianVerification, DissolvedReactionCappedMultiBlock)
{
  auto A = Species{ "A" };
  auto B = Species{ "B" };
  auto C = Species{ "C" };

  auto aqueous_phase = Phase{ "AQUEOUS", { { A }, { B }, { C } } };
  auto droplet = representation::UniformSection{ "DROPLET", { aqueous_phase } };

  auto reaction = process::DissolvedReactionBuilder()
      .SetPhase(aqueous_phase)
      .SetReactants({ A })
      .SetProducts({ B })
      .SetSolvent(C)
      .SetRateConstant([](const Conditions&) { return 2.0; })
      .SetMaxHalflife(0.1)
      .Build();

  auto model = Model{ .name_ = "AEROSOL", .representations_ = { droplet } };
  model.AddProcesses({ reaction });

  auto maps = BuildIndexMaps(model);

  // Two blocks: one uncapped regime, one heavily capped
  DenseMatrix variables(2, maps.num_variables, 0.0);
  variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.A")] = 0.001;  // low rate, uncapped
  variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.B")] = 0.1;
  variables[0][maps.variable_indices.at("DROPLET.AQUEOUS.C")] = 1.0;
  variables[1][maps.variable_indices.at("DROPLET.AQUEOUS.A")] = 50.0;   // high rate, capped
  variables[1][maps.variable_indices.at("DROPLET.AQUEOUS.B")] = 0.1;
  variables[1][maps.variable_indices.at("DROPLET.AQUEOUS.C")] = 1.0;

  DenseMatrix parameters(2, std::max(maps.num_parameters, std::size_t(1)), 0.0);
  std::vector<Conditions> conditions(2);
  conditions[0].temperature_ = 298.15;
  conditions[0].pressure_ = 101325.0;
  conditions[1].temperature_ = 298.15;
  conditions[1].pressure_ = 101325.0;

  VerifyProcessJacobian(model, maps, variables, parameters, conditions);
}

// ═══════════════════════════════════════════════════════════════════════════════
// Rate-Capped → Uncapped Convergence Tests
// ═══════════════════════════════════════════════════════════════════════════════

/// @brief Helper: build a capped and uncapped model from the same reaction definition,
///        evaluate forcing and Jacobian at the same state, and verify convergence.
///
///        When the natural half-life t_nat = C_min / r  >>  max_halflife, the argument
///        u = r / r_max  is small and tanh(u) ≈ u − u³/3. The relative difference
///        between capped and uncapped rates is therefore ≈ u²/3.
///
///        For the Jacobian, the same Taylor argument gives sech²(u) ≈ 1 − 2u²/3 and
///        the correction term scales as u³, so element-wise relative error is also O(u²).
namespace
{
  struct ConvergenceResult
  {
    double max_forcing_rel_error;
    double max_jacobian_rel_error;
  };

  ConvergenceResult CompareCappedToUncapped(
      double k,
      double max_halflife,
      const std::vector<micm::Species>& reactants,
      const std::vector<micm::Species>& products,
      const micm::Species& solvent,
      const micm::Phase& phase,
      const representation::UniformSection& droplet,
      const std::unordered_map<std::string, double>& concentrations)
  {
    auto rate_fn = [k](const Conditions&) { return k; };

    // Build uncapped model
    auto rxn_uncapped = process::DissolvedReactionBuilder()
        .SetPhase(phase)
        .SetReactants(reactants)
        .SetProducts(products)
        .SetSolvent(solvent)
        .SetRateConstant(rate_fn)
        .Build();
    auto model_uncapped = Model{ .name_ = "AEROSOL", .representations_ = { droplet } };
    model_uncapped.AddProcesses({ rxn_uncapped });

    // Build capped model
    auto rxn_capped = process::DissolvedReactionBuilder()
        .SetPhase(phase)
        .SetReactants(reactants)
        .SetProducts(products)
        .SetSolvent(solvent)
        .SetRateConstant(rate_fn)
        .SetMaxHalflife(max_halflife)
        .Build();
    auto model_capped = Model{ .name_ = "AEROSOL", .representations_ = { droplet } };
    model_capped.AddProcesses({ rxn_capped });

    auto maps_u = BuildIndexMaps(model_uncapped);
    auto maps_c = BuildIndexMaps(model_capped);

    // Verify index maps are identical (same species, same ordering)
    EXPECT_EQ(maps_u.num_variables, maps_c.num_variables);
    EXPECT_EQ(maps_u.variable_indices, maps_c.variable_indices);

    const std::size_t n = maps_u.num_variables;

    DenseMatrix variables(1, n, 0.0);
    for (const auto& [name, val] : concentrations)
      variables[0][maps_u.variable_indices.at(name)] = val;

    DenseMatrix params_u(1, std::max(maps_u.num_parameters, std::size_t(1)), 0.0);
    DenseMatrix params_c(1, std::max(maps_c.num_parameters, std::size_t(1)), 0.0);
    std::vector<Conditions> conditions(1);
    conditions[0].temperature_ = 298.15;
    conditions[0].pressure_ = 101325.0;

    model_uncapped.UpdateStateParametersFunction<DenseMatrix>(maps_u.parameter_indices)(conditions, params_u);
    model_capped.UpdateStateParametersFunction<DenseMatrix>(maps_c.parameter_indices)(conditions, params_c);

    // Compare forcing
    DenseMatrix forcing_u(1, n, 0.0), forcing_c(1, n, 0.0);
    model_uncapped.ForcingFunction<DenseMatrix>(maps_u.parameter_indices, maps_u.variable_indices)(params_u, variables, forcing_u);
    model_capped.ForcingFunction<DenseMatrix>(maps_c.parameter_indices, maps_c.variable_indices)(params_c, variables, forcing_c);

    double max_forcing_rel = 0.0;
    for (std::size_t j = 0; j < n; ++j)
    {
      double ref = std::abs(forcing_u[0][j]);
      if (ref > 0.0)
        max_forcing_rel = std::max(max_forcing_rel, std::abs(forcing_c[0][j] - forcing_u[0][j]) / ref);
    }

    // Compare Jacobian
    auto nz = model_uncapped.NonZeroJacobianElements(maps_u.variable_indices);
    auto builder = SparseMatrixFD::Create(n).SetNumberOfBlocks(1).InitialValue(0.0);
    for (const auto& elem : nz)
      builder = builder.WithElement(elem.first, elem.second);
    SparseMatrixFD jac_u(builder), jac_c(builder);

    model_uncapped.JacobianFunction<DenseMatrix, SparseMatrixFD>(
        maps_u.parameter_indices, maps_u.variable_indices, jac_u)(params_u, variables, jac_u);
    model_capped.JacobianFunction<DenseMatrix, SparseMatrixFD>(
        maps_c.parameter_indices, maps_c.variable_indices, jac_c)(params_c, variables, jac_c);

    double max_jac_rel = 0.0;
    const auto& vec_u = jac_u.AsVector();
    const auto& vec_c = jac_c.AsVector();
    for (std::size_t i = 0; i < vec_u.size(); ++i)
    {
      double ref = std::abs(vec_u[i]);
      if (ref > 0.0)
        max_jac_rel = std::max(max_jac_rel, std::abs(vec_c[i] - vec_u[i]) / ref);
    }

    return { max_forcing_rel, max_jac_rel };
  }
}  // namespace

/// @brief Single-reactant convergence: verify O(u²) approach to uncapped as u = r/r_max → 0.
///        k=0.1, [A]=0.5, [C]=1.0  →  r = 0.05.  r_max = [A]/t_half = 0.5/t_half.
///        u = r/r_max = 0.05 * t_half / 0.5 = 0.1 * t_half.
///        For small u, relative error ≈ u²/3.
///        We decrease t_half so that u shrinks and verify monotonic convergence.
TEST(JacobianVerification, CappedConvergesToUncappedSingleReactant)
{
  auto A = Species{ "A" };
  auto B = Species{ "B" };
  auto C = Species{ "C" };
  auto phase = Phase{ "AQUEOUS", { { A }, { B }, { C } } };
  auto droplet = representation::UniformSection{ "DROPLET", { phase } };

  // u = 0.1 * t_half, so t_half ∈ {1, 0.1, 0.01} gives u ∈ {0.1, 0.01, 0.001}
  double prev_forcing_err = 1.0;
  double prev_jac_err = 1.0;
  for (double t_half : { 1.0, 0.1, 0.01 })
  {
    SCOPED_TRACE("max_halflife = " + std::to_string(t_half));
    double u = 0.1 * t_half;
    // Forcing: r_c = r_max*tanh(u) ≈ r*(1 − u²/3), so rel error ≈ u²/3
    double expected_forcing_rel = u * u / 3.0;
    // Jacobian: worst element is ∂r_c/∂S = sech²(u)*∂r/∂S, where |sech²(u)−1| ≈ u²
    double expected_jac_rel = u * u;

    auto result = CompareCappedToUncapped(
        0.1, t_half, { A }, { B }, C, phase, droplet,
        { { "DROPLET.AQUEOUS.A", 0.5 }, { "DROPLET.AQUEOUS.B", 0.1 }, { "DROPLET.AQUEOUS.C", 1.0 } });

    // Verify the error is bounded by the Taylor prediction (with 2x margin for higher-order terms)
    EXPECT_LT(result.max_forcing_rel_error, 2.0 * expected_forcing_rel)
        << "Forcing error " << result.max_forcing_rel_error << " exceeds 2*u^2/3 = " << 2.0 * expected_forcing_rel;
    EXPECT_LT(result.max_jacobian_rel_error, 2.0 * expected_jac_rel)
        << "Jacobian error " << result.max_jacobian_rel_error << " exceeds 2*u^2 = " << 2.0 * expected_jac_rel;

    // Verify monotonic convergence (errors shrink as t_half shrinks → u shrinks)
    EXPECT_LT(result.max_forcing_rel_error, prev_forcing_err);
    EXPECT_LT(result.max_jacobian_rel_error, prev_jac_err);
    prev_forcing_err = result.max_forcing_rel_error;
    prev_jac_err = result.max_jacobian_rel_error;
  }

  // At the smallest t_half (0.01), u = 0.001, expected rel ≈ u² = 1e-6
  EXPECT_LT(prev_forcing_err, 1.0e-6);
  EXPECT_LT(prev_jac_err, 1.0e-5);
}

/// @brief Two-reactant convergence with asymmetric concentrations.
///        [A]=10.0, [B]=0.01  →  C_min ≈ 0.01 (limited by B).
///        k=0.001  →  r = 0.001 * 10.0 * 0.01 = 1e-4.
///        u = r * t_half / C_min = 1e-4 * t_half / 0.01 = 0.01 * t_half.
///        Decrease t_half to verify convergence and final precision.
TEST(JacobianVerification, CappedConvergesToUncappedTwoReactants)
{
  auto A = Species{ "A" };
  auto B = Species{ "B" };
  auto P = Species{ "P" };
  auto S = Species{ "S" };
  auto phase = Phase{ "AQUEOUS", { { A }, { B }, { P }, { S } } };
  auto droplet = representation::UniformSection{ "DROPLET", { phase } };

  double k = 0.001;
  // u = 0.01 * t_half
  double prev_forcing_err = 1.0;
  double prev_jac_err = 1.0;
  for (double t_half : { 10.0, 1.0, 0.1, 0.01 })
  {
    SCOPED_TRACE("max_halflife = " + std::to_string(t_half));

    auto result = CompareCappedToUncapped(
        k, t_half, { A, B }, { P }, S, phase, droplet,
        { { "DROPLET.AQUEOUS.A", 10.0 },
          { "DROPLET.AQUEOUS.B", 0.01 },
          { "DROPLET.AQUEOUS.P", 0.1 },
          { "DROPLET.AQUEOUS.S", 1.0 } });

    // Verify monotonic convergence as t_half shrinks
    EXPECT_LT(result.max_forcing_rel_error, prev_forcing_err);
    EXPECT_LT(result.max_jacobian_rel_error, prev_jac_err);
    prev_forcing_err = result.max_forcing_rel_error;
    prev_jac_err = result.max_jacobian_rel_error;
  }

  // At t_half=0.01, u ≈ 1e-4, so relative error ≈ u²/3 ≈ 3e-9
  EXPECT_LT(prev_forcing_err, 1.0e-6);
  EXPECT_LT(prev_jac_err, 1.0e-6);
}
