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
  HenryLawEquilibriumConstraint constraint(hlc, A_g, A_aq, h2o, aqueous_phase, Mw_water, rho_water);

  std::map<std::string, std::set<std::string>> phase_prefixes;
  phase_prefixes["AQUEOUS"].insert("DROP");

  auto names = constraint.ConstraintAlgebraicVariableNames(phase_prefixes);
  EXPECT_EQ(names.size(), 1);
  EXPECT_TRUE(names.count("DROP.AQUEOUS.A_aq"));
}

TEST(HenryLawEquilibriumConstraint, AlgebraicVariableNamesMultiplePrefixes)
{
  auto hlc = [](const micm::Conditions&) { return 5.0e3; };
  HenryLawEquilibriumConstraint constraint(hlc, A_g, A_aq, h2o, aqueous_phase, Mw_water, rho_water);

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
  HenryLawEquilibriumConstraint constraint(hlc, A_g, A_aq, h2o, aqueous_phase, Mw_water, rho_water);

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
  HenryLawEquilibriumConstraint constraint(hlc, A_g, A_aq, h2o, aqueous_phase, Mw_water, rho_water);

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
  HenryLawEquilibriumConstraint constraint(hlc, A_g, A_aq, h2o, aqueous_phase, Mw_water, rho_water);

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
  HenryLawEquilibriumConstraint constraint(hlc, A_g, A_aq, h2o, aqueous_phase, Mw_water, rho_water);

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
  auto residual_fn = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, state_indices);

  double gas_conc = 1.0e-6;
  double aq_conc = 0.5;
  double solvent_conc = 300.0;

  DMP state_variables{ 1, 3, 0.0 };
  state_variables[0][0] = gas_conc;
  state_variables[0][1] = aq_conc;
  state_variables[0][2] = solvent_conc;

  DMP residual{ 1, 3, 0.0 };
  residual_fn(state_variables, residual);

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
  HenryLawEquilibriumConstraint constraint(hlc, A_g, A_aq, h2o, aqueous_phase, Mw_water, rho_water);

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
  auto residual_fn = constraint.ConstraintResidualFunction<DMP>(phase_prefixes, state_indices);

  DMP state_variables{ 1, 5, 0.0 };
  state_variables[0][0] = 1.0e-5;  // [A_g]
  state_variables[0][1] = 0.3;     // LARGE [A_aq]
  state_variables[0][2] = 400.0;   // LARGE [H2O]
  state_variables[0][3] = 0.1;     // SMALL [A_aq]
  state_variables[0][4] = 200.0;   // SMALL [H2O]

  DMP residual{ 1, 5, 0.0 };
  residual_fn(state_variables, residual);

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
  HenryLawEquilibriumConstraint constraint(hlc, A_g, A_aq, h2o, aqueous_phase, Mw_water, rho_water);

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
  HenryLawEquilibriumConstraint constraint(hlc, A_g, A_aq, h2o, aqueous_phase, Mw_water, rho_water);

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
