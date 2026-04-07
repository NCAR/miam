// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/aerosol_property.hpp>

#include <micm/system/phase.hpp>
#include <micm/util/matrix.hpp>

#include <cmath>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <gtest/gtest.h>

constexpr std::size_t minStateSize = 3;
const std::string test_model_name = "TEST_MODEL";

std::vector<micm::Phase> getTestPhases()
{
    micm::Phase phase1{ "PHASE1", { { micm::Species{ "SPECIES_A" } }, { micm::Species{ "SPECIES_B" } } } };
    micm::Phase phase2{ "PHASE2", { { micm::Species{ "SPECIES_C" } } } };
    return { phase1, phase2 };
}

template<typename ModelPolicy>
std::tuple<std::size_t, std::size_t> testStateSize()
{
    const auto phases = getTestPhases();
    auto model = ModelPolicy{ test_model_name, phases };
    EXPECT_GE(std::get<0>(model.StateSize()), minStateSize);
    EXPECT_GE(std::get<1>(model.StateSize()), 0);
    return model.StateSize();
}

template<typename ModelPolicy>
std::set<std::string> testStateVariableNames()
{
    const auto phases = getTestPhases();
    auto model = ModelPolicy{ test_model_name, phases };
    const auto names = model.StateVariableNames();
    EXPECT_EQ(names.size(), std::get<0>(model.StateSize()));

    return names;
}

template<typename ModelPolicy>
std::set<std::string> testStateParameterNames()
{
    const auto phases = getTestPhases();
    auto model = ModelPolicy{ test_model_name, phases };
    const auto names = model.StateParameterNames();
    EXPECT_EQ(names.size(), std::get<1>(model.StateSize()));

    return names;
}

template<typename ModelPolicy>
std::map<std::string, double> testDefaultParameters()
{
    const auto phases = getTestPhases();
    auto model = ModelPolicy{ test_model_name, phases };
    auto params = model.DefaultParameters();
    auto param_names = model.StateParameterNames();
    
    // Make sure the number of default parameters matches the number of state parameters
    EXPECT_EQ(params.size(), param_names.size());
    
    // Make sure all state parameter names have default values
    for (const auto& name : param_names)
    {
        EXPECT_TRUE(params.find(name) != params.end()); // name should exist in params
    }
    return params;
}

// --- Phases whose species carry molecular weight and density ---

inline std::vector<micm::Phase> getPhysicalTestPhases()
{
    micm::Species species_a{ "SPECIES_A", { { "molecular weight [kg mol-1]", 0.1 },
                                             { "density [kg m-3]", 1000.0 } } };
    micm::Species species_b{ "SPECIES_B", { { "molecular weight [kg mol-1]", 0.2 },
                                             { "density [kg m-3]", 2000.0 } } };
    micm::Species species_c{ "SPECIES_C", { { "molecular weight [kg mol-1]", 0.3 },
                                             { "density [kg m-3]", 1500.0 } } };
    micm::Phase phase1{ "PHASE1", { species_a, species_b } };
    micm::Phase phase2{ "PHASE2", { species_c } };
    return { phase1, phase2 };
}

inline std::vector<micm::Phase> getSinglePhysicalTestPhase()
{
    micm::Species species_a{ "SPECIES_A", { { "molecular weight [kg mol-1]", 0.1 },
                                             { "density [kg m-3]", 1000.0 } } };
    micm::Phase phase1{ "PHASE1", { species_a } };
    return { phase1 };
}

// --- Helpers for provider tests ---

inline void setTestParameters(
    micm::Matrix<double>& params,
    const std::unordered_map<std::string, std::size_t>& param_indices)
{
    for (const auto& [name, idx] : param_indices)
    {
        if (name.find("GEOMETRIC_MEAN_RADIUS") != std::string::npos)
            params[0][idx] = 1.0e-7;
        else if (name.find("GEOMETRIC_STANDARD_DEVIATION") != std::string::npos)
            params[0][idx] = 1.5;
        else if (name.find("MIN_RADIUS") != std::string::npos)
            params[0][idx] = 5.0e-8;
        else if (name.find("MAX_RADIUS") != std::string::npos)
            params[0][idx] = 1.5e-7;
    }
}

inline void setTestVariables(
    micm::Matrix<double>& vars,
    const std::unordered_map<std::string, std::size_t>& var_indices)
{
    vars[0][var_indices.at(test_model_name + ".PHASE1.SPECIES_A")] = 100.0;
    vars[0][var_indices.at(test_model_name + ".PHASE1.SPECIES_B")] = 200.0;
    auto c_it = var_indices.find(test_model_name + ".PHASE2.SPECIES_C");
    if (c_it != var_indices.end())
        vars[0][c_it->second] = 300.0;
    auto nc_it = var_indices.find(test_model_name + ".NUMBER_CONCENTRATION");
    if (nc_it != var_indices.end())
        vars[0][nc_it->second] = 1.0e9;
}

inline void buildIndexMaps(
    const std::set<std::string>& var_names,
    const std::set<std::string>& param_names,
    std::unordered_map<std::string, std::size_t>& var_indices,
    std::unordered_map<std::string, std::size_t>& param_indices)
{
    std::size_t idx = 0;
    for (const auto& n : var_names)
        var_indices[n] = idx++;
    idx = 0;
    for (const auto& n : param_names)
        param_indices[n] = idx++;
}

// Check analytic partials against central finite differences
inline void checkFiniteDifferences(
    const miam::AerosolPropertyProvider<micm::Matrix<double>>& provider,
    micm::Matrix<double>& params,
    micm::Matrix<double>& vars)
{
    using Mat = micm::Matrix<double>;
    std::size_t n_deps = provider.dependent_variable_indices.size();
    Mat result(1, 1, 0.0);
    Mat partials(1, n_deps, 0.0);
    provider.ComputeValueAndDerivatives(params, vars, result, partials);

    for (std::size_t k = 0; k < n_deps; ++k)
    {
        std::size_t var_idx = provider.dependent_variable_indices[k];
        double orig = vars[0][var_idx];
        double h = std::max(std::abs(orig) * 1.0e-5, 1.0e-8);

        Mat r_plus(1, 1, 0.0);
        Mat r_minus(1, 1, 0.0);
        vars[0][var_idx] = orig + h;
        provider.ComputeValue(params, vars, r_plus);
        vars[0][var_idx] = orig - h;
        provider.ComputeValue(params, vars, r_minus);
        vars[0][var_idx] = orig;

        double fd = (r_plus[0][0] - r_minus[0][0]) / (2.0 * h);
        double tol = std::max(std::abs(fd) * 1.0e-4, 1.0e-12);
        EXPECT_NEAR(partials[0][k], fd, tol)
            << "partial[" << k << "] (var idx " << var_idx << "): analytic=" << partials[0][k] << " fd=" << fd;
    }
}

// --- Phase structure tests ---

template<typename ModelPolicy>
void testPhaseStatePrefixes()
{
    const auto phases = getPhysicalTestPhases();
    auto model = ModelPolicy{ test_model_name, phases };
    auto prefixes = model.PhaseStatePrefixes();

    EXPECT_EQ(prefixes.size(), 2u);
    ASSERT_TRUE(prefixes.find("PHASE1") != prefixes.end());
    ASSERT_TRUE(prefixes.find("PHASE2") != prefixes.end());
    EXPECT_EQ(prefixes["PHASE1"].size(), 1u);
    EXPECT_EQ(prefixes["PHASE2"].size(), 1u);
    EXPECT_TRUE(prefixes["PHASE1"].count(test_model_name));
    EXPECT_TRUE(prefixes["PHASE2"].count(test_model_name));
}

// --- PhaseVolumeFraction: single phase ---

template<typename ModelPolicy>
void testPhaseVolumeFractionSinglePhase()
{
    using Mat = micm::Matrix<double>;
    const auto phases = getSinglePhysicalTestPhase();
    auto model = ModelPolicy{ test_model_name, phases };

    auto var_names = model.StateVariableNames();
    auto param_names = model.StateParameterNames();
    std::unordered_map<std::string, std::size_t> var_indices, param_indices;
    buildIndexMaps(var_names, param_names, var_indices, param_indices);

    auto provider = model.template GetPropertyProvider<Mat>(
        miam::AerosolProperty::PhaseVolumeFraction, param_indices, var_indices, "PHASE1");

    // Single phase → no dependent variables, always 1.0
    EXPECT_TRUE(provider.dependent_variable_indices.empty());
    EXPECT_TRUE(provider.ComputeValue);
    EXPECT_TRUE(provider.ComputeValueAndDerivatives);

    Mat params(1, param_names.size(), 0.0);
    Mat vars(1, var_names.size(), 0.0);
    setTestParameters(params, param_indices);
    vars[0][var_indices.at(test_model_name + ".PHASE1.SPECIES_A")] = 42.0;
    auto nc_it = var_indices.find(test_model_name + ".NUMBER_CONCENTRATION");
    if (nc_it != var_indices.end())
        vars[0][nc_it->second] = 1.0e9;

    Mat result(1, 1, 0.0);
    provider.ComputeValue(params, vars, result);
    EXPECT_DOUBLE_EQ(result[0][0], 1.0);

    Mat partials(1, 0, 0.0);
    provider.ComputeValueAndDerivatives(params, vars, result, partials);
    EXPECT_DOUBLE_EQ(result[0][0], 1.0);
}

// --- PhaseVolumeFraction: multi-phase values and derivatives ---

template<typename ModelPolicy>
void testPhaseVolumeFractionMultiPhase()
{
    using Mat = micm::Matrix<double>;
    const auto phases = getPhysicalTestPhases();
    auto model = ModelPolicy{ test_model_name, phases };

    auto var_names = model.StateVariableNames();
    auto param_names = model.StateParameterNames();
    std::unordered_map<std::string, std::size_t> var_indices, param_indices;
    buildIndexMaps(var_names, param_names, var_indices, param_indices);

    auto phi1_provider = model.template GetPropertyProvider<Mat>(
        miam::AerosolProperty::PhaseVolumeFraction, param_indices, var_indices, "PHASE1");
    auto phi2_provider = model.template GetPropertyProvider<Mat>(
        miam::AerosolProperty::PhaseVolumeFraction, param_indices, var_indices, "PHASE2");

    EXPECT_FALSE(phi1_provider.dependent_variable_indices.empty());
    EXPECT_FALSE(phi2_provider.dependent_variable_indices.empty());
    EXPECT_TRUE(phi1_provider.ComputeValue);
    EXPECT_TRUE(phi2_provider.ComputeValue);
    EXPECT_TRUE(phi1_provider.ComputeValueAndDerivatives);
    EXPECT_TRUE(phi2_provider.ComputeValueAndDerivatives);

    Mat params(1, param_names.size(), 0.0);
    Mat vars(1, var_names.size(), 0.0);
    setTestParameters(params, param_indices);
    setTestVariables(vars, var_indices);

    // mw/rho: A = 0.1/1000 = 1e-4, B = 0.2/2000 = 1e-4, C = 0.3/1500 = 2e-4
    // V_phase1 = 100*1e-4 + 200*1e-4 = 0.03
    // V_phase2 = 300*2e-4 = 0.06
    // V_total  = 0.09
    // phi_1 = 1/3,  phi_2 = 2/3
    Mat result1(1, 1, 0.0);
    Mat result2(1, 1, 0.0);

    phi1_provider.ComputeValue(params, vars, result1);
    phi2_provider.ComputeValue(params, vars, result2);

    EXPECT_NEAR(result1[0][0], 1.0 / 3.0, 1.0e-12);
    EXPECT_NEAR(result2[0][0], 2.0 / 3.0, 1.0e-12);
    EXPECT_NEAR(result1[0][0] + result2[0][0], 1.0, 1.0e-12);

    // Verify analytic partials via finite differences
    checkFiniteDifferences(phi1_provider, params, vars);
    checkFiniteDifferences(phi2_provider, params, vars);
}

// --- EffectiveRadius: positive value and derivatives ---

template<typename ModelPolicy>
void testEffectiveRadiusProvider()
{
    using Mat = micm::Matrix<double>;
    const auto phases = getPhysicalTestPhases();
    auto model = ModelPolicy{ test_model_name, phases };

    auto var_names = model.StateVariableNames();
    auto param_names = model.StateParameterNames();
    std::unordered_map<std::string, std::size_t> var_indices, param_indices;
    buildIndexMaps(var_names, param_names, var_indices, param_indices);

    auto provider = model.template GetPropertyProvider<Mat>(
        miam::AerosolProperty::EffectiveRadius, param_indices, var_indices);

    EXPECT_TRUE(provider.ComputeValue);
    EXPECT_TRUE(provider.ComputeValueAndDerivatives);

    Mat params(1, param_names.size(), 0.0);
    Mat vars(1, var_names.size(), 0.0);
    setTestParameters(params, param_indices);
    setTestVariables(vars, var_indices);

    Mat result(1, 1, 0.0);
    provider.ComputeValue(params, vars, result);
    EXPECT_GT(result[0][0], 0.0);

    checkFiniteDifferences(provider, params, vars);
}

// --- NumberConcentration: positive value and derivatives ---

template<typename ModelPolicy>
void testNumberConcentrationProvider()
{
    using Mat = micm::Matrix<double>;
    const auto phases = getPhysicalTestPhases();
    auto model = ModelPolicy{ test_model_name, phases };

    auto var_names = model.StateVariableNames();
    auto param_names = model.StateParameterNames();
    std::unordered_map<std::string, std::size_t> var_indices, param_indices;
    buildIndexMaps(var_names, param_names, var_indices, param_indices);

    auto provider = model.template GetPropertyProvider<Mat>(
        miam::AerosolProperty::NumberConcentration, param_indices, var_indices);

    EXPECT_TRUE(provider.ComputeValue);
    EXPECT_TRUE(provider.ComputeValueAndDerivatives);
    EXPECT_FALSE(provider.dependent_variable_indices.empty());

    Mat params(1, param_names.size(), 0.0);
    Mat vars(1, var_names.size(), 0.0);
    setTestParameters(params, param_indices);
    setTestVariables(vars, var_indices);

    Mat result(1, 1, 0.0);
    provider.ComputeValue(params, vars, result);
    EXPECT_GT(result[0][0], 0.0);

    checkFiniteDifferences(provider, params, vars);
}