// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <miam/process_set.hpp>
#include <miam/processes/dissolved_reversible_reaction.hpp>
#include <micm/system/conditions.hpp>
#include <micm/system/phase.hpp>
#include <micm/system/species.hpp>
#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>

#include <gtest/gtest.h>

using namespace miam;

using MatrixPolicy = micm::Matrix<double>;
using SparseMatrixPolicy = micm::SparseMatrix<double, micm::SparseMatrixStandardOrderingCompressedSparseRow>;
using ProcessSet = MiamProcessSet<MatrixPolicy, SparseMatrixPolicy>;

namespace
{
  // Helper to create a simple A <-> B reaction for tests
  struct TestFixture
  {
    micm::Species a{ "A" };
    micm::Species b{ "B" };
    micm::Species solvent{ "SOLVENT" };
    micm::Phase aqueous_phase{ "AQUEOUS", { { a }, { b }, { solvent } } };
    double k_forward = 0.1;
    double k_reverse = 0.05;

    process::DissolvedReversibleReaction MakeReaction() const
    {
      return process::DissolvedReversibleReaction{
        [kf = k_forward](const micm::Conditions&) { return kf; },
        [kr = k_reverse](const micm::Conditions&) { return kr; },
        { a },
        { b },
        solvent,
        aqueous_phase
      };
    }

    ProcessSet::PhaseMap MakePhaseMap() const
    {
      ProcessSet::PhaseMap phase_prefixes;
      phase_prefixes["AQUEOUS"].insert("DROP");
      return phase_prefixes;
    }

    ProcessSet::IndexMap MakeVariableIndices() const
    {
      return {
          { "DROP.AQUEOUS.A", 0 },
          { "DROP.AQUEOUS.B", 1 },
          { "DROP.AQUEOUS.SOLVENT", 2 }
      };
    }
  };
}  // namespace

TEST(MiamProcessSet, ConstructFromDissolvedReversibleReaction)
{
    TestFixture fix;
    // Should compile and not throw
    ProcessSet ps(fix.MakeReaction());

    // All function members should be assigned
    EXPECT_TRUE(ps.process_parameter_names_);
    EXPECT_TRUE(ps.species_used_);
    EXPECT_TRUE(ps.required_aerosol_properties_);
    EXPECT_TRUE(ps.non_zero_jacobian_elements_);
    EXPECT_TRUE(ps.update_state_parameters_function_);
    EXPECT_TRUE(ps.get_forcing_function_);
    EXPECT_TRUE(ps.get_jacobian_function_);
}

TEST(MiamProcessSet, ProcessParameterNames)
{
    TestFixture fix;
    ProcessSet ps(fix.MakeReaction());
    auto phase_prefixes = fix.MakePhaseMap();

    auto names = ps.process_parameter_names_(phase_prefixes);

    // Should have 2 parameter names (k_forward and k_reverse) for the single phase instance
    EXPECT_EQ(names.size(), 2);
    // Names follow the pattern AQUEOUS.<uuid>.k_forward / k_reverse — just check count
    for (const auto& name : names)
    {
        EXPECT_TRUE(name.find("AQUEOUS.") == 0);
    }
}

TEST(MiamProcessSet, SpeciesUsed)
{
    TestFixture fix;
    ProcessSet ps(fix.MakeReaction());
    auto phase_prefixes = fix.MakePhaseMap();

    auto species = ps.species_used_(phase_prefixes);

    EXPECT_EQ(species.size(), 3);
    EXPECT_TRUE(species.count("DROP.AQUEOUS.A"));
    EXPECT_TRUE(species.count("DROP.AQUEOUS.B"));
    EXPECT_TRUE(species.count("DROP.AQUEOUS.SOLVENT"));
}

TEST(MiamProcessSet, RequiredAerosolProperties)
{
    TestFixture fix;
    ProcessSet ps(fix.MakeReaction());

    auto required = ps.required_aerosol_properties_();

    // DissolvedReversibleReaction needs no aerosol properties
    EXPECT_TRUE(required.empty());
}

TEST(MiamProcessSet, NonZeroJacobianElements)
{
    TestFixture fix;
    ProcessSet ps(fix.MakeReaction());
    auto phase_prefixes = fix.MakePhaseMap();
    auto var_indices = fix.MakeVariableIndices();
    ProcessSet::ProviderMap providers;  // empty — not used by this process

    auto elements = ps.non_zero_jacobian_elements_(phase_prefixes, var_indices, providers);

    // A <-> B with solvent: the reaction touches all 3 species
    // All pairs of (reactant, product, solvent) × (reactant, product, solvent) should appear
    EXPECT_FALSE(elements.empty());
    // At minimum: d[A]/d[A], d[A]/d[B], d[B]/d[A], d[B]/d[B], plus solvent couplings
    EXPECT_TRUE(elements.count({ 0, 0 }));  // d[A]/d[A]
    EXPECT_TRUE(elements.count({ 1, 0 }));  // d[B]/d[A]
    EXPECT_TRUE(elements.count({ 0, 1 }));  // d[A]/d[B]
    EXPECT_TRUE(elements.count({ 1, 1 }));  // d[B]/d[B]
    EXPECT_TRUE(elements.count({ 0, 2 }));  // d[A]/d[SOLVENT]
    EXPECT_TRUE(elements.count({ 1, 2 }));  // d[B]/d[SOLVENT]
}

TEST(MiamProcessSet, UpdateStateParametersFunction)
{
    TestFixture fix;
    auto reaction = fix.MakeReaction();
    auto phase_prefixes = fix.MakePhaseMap();

    // Build parameter names so we can construct the index map
    auto param_names = reaction.ProcessParameterNames(phase_prefixes);
    ProcessSet::IndexMap param_indices;
    std::size_t idx = 0;
    for (const auto& name : param_names)
    {
        param_indices[name] = idx++;
    }

    ProcessSet ps(std::move(reaction));
    auto update_fn = ps.update_state_parameters_function_(phase_prefixes, param_indices);

    EXPECT_TRUE(update_fn);

    // 2 grid cells, 2 parameters
    MatrixPolicy state_parameters(2, param_names.size(), 0.0);
    std::vector<micm::Conditions> conditions(2);
    conditions[0].temperature_ = 298.15;
    conditions[1].temperature_ = 310.0;

    update_fn(conditions, state_parameters);

    // Rate constants should be written (0.1 and 0.05 for both cells)
    // Parameters are not necessarily in k_forward, k_reverse order, but values should be 0.1 and 0.05
    for (std::size_t i_cell = 0; i_cell < 2; ++i_cell)
    {
        double val0 = state_parameters[i_cell][0];
        double val1 = state_parameters[i_cell][1];
        EXPECT_TRUE((val0 == 0.1 && val1 == 0.05) || (val0 == 0.05 && val1 == 0.1));
    }
}

TEST(MiamProcessSet, ForcingFunction)
{
    TestFixture fix;
    auto reaction = fix.MakeReaction();
    auto phase_prefixes = fix.MakePhaseMap();

    auto param_names = reaction.ProcessParameterNames(phase_prefixes);
    ProcessSet::IndexMap param_indices;
    std::size_t idx = 0;
    for (const auto& name : param_names)
    {
        param_indices[name] = idx++;
    }
    auto var_indices = fix.MakeVariableIndices();

    // Wrap the same reaction so UUIDs match the parameter index map
    ProcessSet ps(reaction);
    ProcessSet::ProviderMap providers;

    auto forcing_fn = ps.get_forcing_function_(phase_prefixes, param_indices, var_indices, std::move(providers));
    EXPECT_TRUE(forcing_fn);

    // Set up state: 1 grid cell, 2 params, 3 variables
    MatrixPolicy state_parameters(1, param_names.size(), 0.0);
    MatrixPolicy state_variables(1, 3, 0.0);
    MatrixPolicy forcing(1, 3, 0.0);

    // Fill parameters using the same wrapped process
    auto update_fn = ps.update_state_parameters_function_(phase_prefixes, param_indices);
    std::vector<micm::Conditions> conditions(1);
    conditions[0].temperature_ = 298.15;
    update_fn(conditions, state_parameters);

    // Set concentration: A=1.0, B=0.0, SOLVENT=0.017
    state_variables[0][0] = 1.0;
    state_variables[0][1] = 0.0;
    state_variables[0][2] = 0.017;

    forcing_fn(state_parameters, state_variables, forcing);

    // For A <-> B with 1 reactant and 1 product:
    //   forward = k_f * [A] / [SOLVENT]^(n_reactants-1) = k_f * [A]
    //   reverse = k_r * [B] / [SOLVENT]^(n_products-1)  = k_r * [B]
    //   f_A = -forward + reverse, f_B = forward - reverse
    double forward_val = fix.k_forward * 1.0;
    double reverse_val = fix.k_reverse * 0.0;
    double expected_f_A = -forward_val + reverse_val;
    double expected_f_B = forward_val - reverse_val;

    EXPECT_NEAR(forcing[0][0], expected_f_A, 1e-10);
    EXPECT_NEAR(forcing[0][1], expected_f_B, 1e-10);
}

TEST(MiamProcessSet, JacobianFunction)
{
    TestFixture fix;
    auto reaction = fix.MakeReaction();
    auto phase_prefixes = fix.MakePhaseMap();

    auto param_names = reaction.ProcessParameterNames(phase_prefixes);
    ProcessSet::IndexMap param_indices;
    std::size_t idx = 0;
    for (const auto& name : param_names)
    {
        param_indices[name] = idx++;
    }
    auto var_indices = fix.MakeVariableIndices();

    // Build the sparse Jacobian structure
    auto nz_elements = reaction.NonZeroJacobianElements(phase_prefixes, var_indices);
    auto jacobian_builder = SparseMatrixPolicy::Create(3).SetNumberOfBlocks(1);
    for (const auto& [row, col] : nz_elements)
    {
        jacobian_builder = jacobian_builder.WithElement(row, col);
    }
    SparseMatrixPolicy jacobian(jacobian_builder);
    jacobian.Fill(0.0);

    // Wrap the same reaction so UUIDs match
    ProcessSet ps(reaction);
    ProcessSet::ProviderMap providers;

    auto jacobian_fn =
        ps.get_jacobian_function_(phase_prefixes, param_indices, var_indices, jacobian, std::move(providers));
    EXPECT_TRUE(jacobian_fn);

    // Set up state: 1 grid cell
    MatrixPolicy state_parameters(1, param_names.size(), 0.0);
    MatrixPolicy state_variables(1, 3, 0.0);

    // Fill parameters using the same wrapped process
    auto update_fn = ps.update_state_parameters_function_(phase_prefixes, param_indices);
    std::vector<micm::Conditions> conditions(1);
    conditions[0].temperature_ = 298.15;
    update_fn(conditions, state_parameters);

    state_variables[0][0] = 1.0;   // A
    state_variables[0][1] = 0.5;   // B
    state_variables[0][2] = 0.017;  // SOLVENT

    jacobian_fn(state_parameters, state_variables, jacobian);

    // Verify non-zero Jacobian entries were written
    bool has_nonzero = false;
    for (const auto& [row, col] : nz_elements)
    {
        double val = jacobian[0][row][col];
        if (val != 0.0)
        {
            has_nonzero = true;
            break;
        }
    }
    EXPECT_TRUE(has_nonzero);
}

TEST(MiamProcessSet, MoveConstruction)
{
    TestFixture fix;
    ProcessSet ps1(fix.MakeReaction());

    // Move construct
    ProcessSet ps2(std::move(ps1));

    auto phase_prefixes = fix.MakePhaseMap();
    auto species = ps2.species_used_(phase_prefixes);
    EXPECT_EQ(species.size(), 3);
}

TEST(MiamProcessSet, StorableInVector)
{
    TestFixture fix;
    std::vector<ProcessSet> processes;
    processes.push_back(ProcessSet(fix.MakeReaction()));
    processes.push_back(ProcessSet(fix.MakeReaction()));

    auto phase_prefixes = fix.MakePhaseMap();
    for (const auto& ps : processes)
    {
        auto species = ps.species_used_(phase_prefixes);
        EXPECT_EQ(species.size(), 3);
    }
}
