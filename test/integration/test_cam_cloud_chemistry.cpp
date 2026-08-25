// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// CPU driver for the CAM cloud aqueous-phase sulfate-chemistry integration
// tests. Each TEST() below forwards to a templated policy function in
// cam_cloud_chemistry_policy.hpp with a `CpuSolverBuilder`. The FD Jacobian
// verification (Step 5) uses fixed CPU matrix types and lives here since it
// exercises Model's original std::function-returning API directly rather
// than driving a solver.

#include "cam_cloud_chemistry_policy.hpp"

#include <miam/miam.hpp>

#include <micm/CPU.hpp>
#include <micm/util/jacobian_verification.hpp>

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

namespace policy = miam_test_cam_cloud_chemistry;

using CpuBuilder = micm::CpuSolverBuilder<micm::RosenbrockSolverParameters>;

namespace
{
  CpuBuilder StandardBuilder()
  {
    return CpuBuilder(micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters());
  }

  // Step 1c needs tighter constraint-initialisation and more max steps
  CpuBuilder Step1cBuilder()
  {
    auto params = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
    params.constraint_init_max_iterations_ = 200;
    params.constraint_init_tolerance_ = 1e-20;
    params.max_number_of_steps_ = 1500;
    return CpuBuilder(params);
  }
}  // namespace

TEST(CamCloudChemistry, Step1_SingleHLC)
{
  policy::Step1_SingleHLC(StandardBuilder());
}

TEST(CamCloudChemistry, Step1b_KwOnly)
{
  policy::Step1b_KwOnly(StandardBuilder());
}

TEST(CamCloudChemistry, Step1c_KwNaiveIC)
{
  policy::Step1c_KwNaiveIC(Step1cBuilder());
}

TEST(CamCloudChemistry, Step2_HLC_Plus_Dissociation)
{
  policy::Step2_HLC_Plus_Dissociation(StandardBuilder());
}

TEST(CamCloudChemistry, Step3_FullEquilibrium)
{
  policy::Step3_FullEquilibrium(StandardBuilder());
}

TEST(CamCloudChemistry, Step3b_NaiveInitialConditions)
{
  policy::Step3b_NaiveInitialConditions(StandardBuilder());
}

TEST(CamCloudChemistry, Step4_FullSystemWithKinetics)
{
  policy::Step4_FullSystemWithKinetics(StandardBuilder());
}

TEST(CamCloudChemistry, Step4b_NaiveInitialConditions)
{
  policy::Step4b_NaiveInitialConditions(StandardBuilder());
}

// ════════════════════════════════════════════════════════════════════════
// TEST 5: FD Jacobian verification against the analytical Model Jacobian.
// CPU-only — exercises Model's std::function-returning API to build a
// finite-difference reference and does not run a solver.
// ════════════════════════════════════════════════════════════════════════
TEST(CamCloudChemistry, Step5_JacobianVerification)
{
  using namespace micm;
  using namespace miam;

  using DenseMatrix = policy::DenseMatrix;
  constexpr double M_ATM_TO_MOL_M3_PA = policy::M_ATM_TO_MOL_M3_PA;
  constexpr double c_H2O_M = policy::c_H2O_M;
  constexpr double C_H2O = policy::C_H2O;

  auto so2_g = Species{ "SO2" };
  auto h2o2_g = Species{ "H2O2" };
  auto o3_g = Species{ "O3" };
  auto so2_aq = Species{ "SO2_aq" };
  auto h2o2_aq = Species{ "H2O2_aq" };
  auto o3_aq = Species{ "O3_aq" };
  auto hp = Species{ "Hp" };
  auto ohm = Species{ "OHm" };
  auto hso3m = Species{ "HSO3m" };
  auto so3mm = Species{ "SO3mm" };
  auto so4mm = Species{ "SO4mm" };
  auto so2oohm = Species{ "SO2OOHm" };
  auto h2o = Species{ "H2O", { { "molecular weight [kg mol-1]", 0.018 }, { "density [kg m-3]", 1000.0 } } };

  Phase gas_phase{ "GAS", { so2_g, h2o2_g, o3_g } };
  Phase aqueous_phase{ "AQUEOUS", { h2o, so2_aq, h2o2_aq, o3_aq, hp, ohm, hso3m, so3mm, so4mm, so2oohm } };
  auto cloud = UniformSection{ "CLOUD", { aqueous_phase } };

  auto hl_so2 = HenrysLawEquilibriumConstraintBuilder()
                    .SetGasSpecies(so2_g)
                    .SetCondensedSpecies(so2_aq)
                    .SetSolvent(h2o)
                    .SetCondensedPhase(aqueous_phase)
                    .SetHenrysLawConstant(HenrysLawConstant({ .HLC_ref_ = 1.23 * M_ATM_TO_MOL_M3_PA, .C_ = 3120.0 }))
                    .Build();

  auto hl_h2o2 = HenrysLawEquilibriumConstraintBuilder()
                     .SetGasSpecies(h2o2_g)
                     .SetCondensedSpecies(h2o2_aq)
                     .SetSolvent(h2o)
                     .SetCondensedPhase(aqueous_phase)
                     .SetHenrysLawConstant(HenrysLawConstant({ .HLC_ref_ = 7.4e4 * M_ATM_TO_MOL_M3_PA, .C_ = 6621.0 }))
                     .Build();

  auto hl_o3 = HenrysLawEquilibriumConstraintBuilder()
                   .SetGasSpecies(o3_g)
                   .SetCondensedSpecies(o3_aq)
                   .SetSolvent(h2o)
                   .SetCondensedPhase(aqueous_phase)
                   .SetHenrysLawConstant(HenrysLawConstant({ .HLC_ref_ = 1.15e-2 * M_ATM_TO_MOL_M3_PA, .C_ = 2560.0 }))
                   .Build();

  auto eq_kw = DissolvedEquilibriumConstraintBuilder()
                   .SetPhase(aqueous_phase)
                   .SetReactants({ h2o })
                   .SetProducts({ hp, ohm })
                   .SetAlgebraicSpecies(ohm)
                   .SetSolvent(h2o)
                   .SetEquilibriumConstant(EquilibriumConstant({ .A_ = 1.0e-14 / (c_H2O_M * c_H2O_M), .C_ = 6710.0 }))
                   .Build();

  auto eq_ka1 = DissolvedEquilibriumConstraintBuilder()
                    .SetPhase(aqueous_phase)
                    .SetReactants({ so2_aq })
                    .SetProducts({ hso3m, hp })
                    .SetAlgebraicSpecies(hso3m)
                    .SetSolvent(h2o)
                    .SetEquilibriumConstant(EquilibriumConstant({ .A_ = 1.7e-2 / c_H2O_M, .C_ = 2090.0 }))
                    .Build();

  auto eq_ka2 = DissolvedEquilibriumConstraintBuilder()
                    .SetPhase(aqueous_phase)
                    .SetReactants({ hso3m })
                    .SetProducts({ so3mm, hp })
                    .SetAlgebraicSpecies(so3mm)
                    .SetSolvent(h2o)
                    .SetEquilibriumConstant(EquilibriumConstant({ .A_ = 6.0e-8 / c_H2O_M, .C_ = 1120.0 }))
                    .Build();

  double total_S = 3.01e-8 + 1.0;
  auto mass_S = LinearConstraintBuilder()
                    .SetAlgebraicSpecies(gas_phase, so2_g)
                    .AddTerm(gas_phase, so2_g, 1.0)
                    .AddTerm(aqueous_phase, so2_aq, 1.0)
                    .AddTerm(aqueous_phase, hso3m, 1.0)
                    .AddTerm(aqueous_phase, so3mm, 1.0)
                    .AddTerm(aqueous_phase, so4mm, 1.0)
                    .AddTerm(aqueous_phase, so2oohm, 1.0)
                    .DiagnoseConstantFromState()
                    .Build();

  auto mass_H2O2 = LinearConstraintBuilder()
                       .SetAlgebraicSpecies(gas_phase, h2o2_g)
                       .AddTerm(gas_phase, h2o2_g, 1.0)
                       .AddTerm(aqueous_phase, h2o2_aq, 1.0)
                       .DiagnoseConstantFromState()
                       .Build();

  auto mass_O3 = LinearConstraintBuilder()
                     .SetAlgebraicSpecies(gas_phase, o3_g)
                     .AddTerm(gas_phase, o3_g, 1.0)
                     .AddTerm(aqueous_phase, o3_aq, 1.0)
                     .DiagnoseConstantFromState()
                     .Build();

  auto charge = LinearConstraintBuilder()
                    .SetAlgebraicSpecies(aqueous_phase, hp)
                    .AddTerm(aqueous_phase, hp, 1.0)
                    .AddTerm(aqueous_phase, ohm, -1.0)
                    .AddTerm(aqueous_phase, hso3m, -1.0)
                    .AddTerm(aqueous_phase, so3mm, -2.0)
                    .AddTerm(aqueous_phase, so4mm, -2.0)
                    .AddTerm(aqueous_phase, so2oohm, -1.0)
                    .SetConstant(0.0)
                    .Build();

  auto rxn1a = DissolvedReversibleReactionBuilder()
                   .SetPhase(aqueous_phase)
                   .SetReactants({ hso3m, h2o2_aq })
                   .SetProducts({ so2oohm, h2o })
                   .SetSolvent(h2o)
                   .SetForwardRateConstant(EquilibriumConstant({ .A_ = c_H2O_M * (7.45e7 / 13.0), .C_ = 4430.0 }))
                   .SetEquilibriumConstant(EquilibriumConstant({ .A_ = 1725.0 }))
                   .Build();

  auto rxn1b = DissolvedReactionBuilder()
                   .SetPhase(aqueous_phase)
                   .SetReactants({ so2oohm, hp })
                   .SetProducts({ so4mm })
                   .SetSolvent(h2o)
                   .SetRateConstant(
                       [](const Conditions& c) -> double
                       { return c_H2O_M * 2.4e6 * std::exp(-4430.0 * (1.0 / c.temperature_ - 1.0 / 298.0)); })
                   .Build();

  auto rxn2 = DissolvedReactionBuilder()
                  .SetPhase(aqueous_phase)
                  .SetReactants({ hso3m, o3_aq })
                  .SetProducts({ so4mm, hp })
                  .SetSolvent(h2o)
                  .SetRateConstant(
                      [](const Conditions& c) -> double
                      { return c_H2O_M * 3.75e5 * std::exp(-5530.0 * (1.0 / c.temperature_ - 1.0 / 298.0)); })
                  .Build();

  auto rxn3 = DissolvedReactionBuilder()
                  .SetPhase(aqueous_phase)
                  .SetReactants({ so3mm, o3_aq })
                  .SetProducts({ so4mm })
                  .SetSolvent(h2o)
                  .SetRateConstant(
                      [](const Conditions& c) -> double
                      { return c_H2O_M * 1.59e9 * std::exp(-5280.0 * (1.0 / c.temperature_ - 1.0 / 298.0)); })
                  .Build();

  auto model = Model{ .name_ = "CLOUD", .representations_ = { cloud } };
  model.AddProcesses(rxn1a, rxn1b, rxn2, rxn3);
  model.AddConstraints(hl_so2, hl_h2o2, hl_o3, eq_kw, eq_ka1, eq_ka2, mass_S, mass_H2O2, mass_O3, charge);

  auto maps = policy::BuildIndexMaps(model);

  DenseMatrix variables(2, maps.num_variables, 0.0);
  auto set_var = [&](int block, const std::string& name, double val)
  { variables[block][maps.variable_indices.at(name)] = val; };

  set_var(0, "SO2", 2.5e-8);
  set_var(0, "H2O2", 2.0e-8);
  set_var(0, "O3", 1.5e-6);
  set_var(0, "CLOUD.AQUEOUS.H2O", C_H2O);
  set_var(0, "CLOUD.AQUEOUS.SO2_aq", 1.0e-6);
  set_var(0, "CLOUD.AQUEOUS.H2O2_aq", 0.1);
  set_var(0, "CLOUD.AQUEOUS.O3_aq", 5.0e-7);
  set_var(0, "CLOUD.AQUEOUS.Hp", 5.0e-3);
  set_var(0, "CLOUD.AQUEOUS.OHm", 2.0e-6);
  set_var(0, "CLOUD.AQUEOUS.HSO3m", 3.0e-3);
  set_var(0, "CLOUD.AQUEOUS.SO3mm", 3.0e-5);
  set_var(0, "CLOUD.AQUEOUS.SO4mm", 1.0);
  set_var(0, "CLOUD.AQUEOUS.SO2OOHm", 1.0e-6);

  set_var(1, "SO2", 1.0e-8);
  set_var(1, "H2O2", 1.0e-8);
  set_var(1, "O3", 1.4e-6);
  set_var(1, "CLOUD.AQUEOUS.H2O", C_H2O);
  set_var(1, "CLOUD.AQUEOUS.SO2_aq", 5.0e-7);
  set_var(1, "CLOUD.AQUEOUS.H2O2_aq", 0.05);
  set_var(1, "CLOUD.AQUEOUS.O3_aq", 3.0e-7);
  set_var(1, "CLOUD.AQUEOUS.Hp", 8.0e-3);
  set_var(1, "CLOUD.AQUEOUS.OHm", 1.0e-6);
  set_var(1, "CLOUD.AQUEOUS.HSO3m", 1.5e-3);
  set_var(1, "CLOUD.AQUEOUS.SO3mm", 1.0e-5);
  set_var(1, "CLOUD.AQUEOUS.SO4mm", 1.5);
  set_var(1, "CLOUD.AQUEOUS.SO2OOHm", 5.0e-7);

  DenseMatrix parameters(2, std::max(maps.num_parameters, std::size_t(1)), 0.0);
  policy::Vector<Conditions> conditions(2);
  conditions[0].temperature_ = 280.0;
  conditions[0].pressure_ = 70000.0;
  conditions[1].temperature_ = 290.0;
  conditions[1].pressure_ = 80000.0;

  std::cout << "\n=== Step 5: FD Jacobian Verification ===" << std::endl;

  std::cout << "Checking process Jacobian..." << std::endl;
  policy::VerifyProcessJacobian(model, maps, variables, parameters, conditions, 1.0e-4, 1.0e-4);

  std::cout << "Checking constraint Jacobian..." << std::endl;
  policy::VerifyConstraintJacobian(model, maps, variables, parameters, conditions);

  std::cout << "=== Step 5 PASSED ===" << std::endl;
}
