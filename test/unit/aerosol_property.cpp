// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <miam/aerosol_property.hpp>
#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>

using namespace miam;

TEST(AerosolProperty, EnumValues)
{
    // Verify distinct enum values
    EXPECT_NE(
        static_cast<int>(AerosolProperty::EffectiveRadius),
        static_cast<int>(AerosolProperty::NumberConcentration));
    EXPECT_NE(
        static_cast<int>(AerosolProperty::EffectiveRadius),
        static_cast<int>(AerosolProperty::PhaseVolumeFraction));
    EXPECT_NE(
        static_cast<int>(AerosolProperty::NumberConcentration),
        static_cast<int>(AerosolProperty::PhaseVolumeFraction));
}

TEST(AerosolProperty, EnumUsableAsMapKey)
{
    std::map<AerosolProperty, std::string> property_names;
    property_names[AerosolProperty::EffectiveRadius] = "r_eff";
    property_names[AerosolProperty::NumberConcentration] = "N";
    property_names[AerosolProperty::PhaseVolumeFraction] = "phi";

    EXPECT_EQ(property_names.size(), 3);
    EXPECT_EQ(property_names[AerosolProperty::EffectiveRadius], "r_eff");
    EXPECT_EQ(property_names[AerosolProperty::NumberConcentration], "N");
    EXPECT_EQ(property_names[AerosolProperty::PhaseVolumeFraction], "phi");
}

TEST(AerosolPropertyProvider, DefaultConstruction)
{
    using MatrixPolicy = micm::VectorMatrix<double>;
    AerosolPropertyProvider<MatrixPolicy> provider;

    // Default-constructed provider should have empty dependent variable indices
    EXPECT_TRUE(provider.dependent_variable_indices.empty());
    // Function members are empty (not assigned)
    EXPECT_FALSE(provider.ComputeValue);
    EXPECT_FALSE(provider.ComputeValueAndDerivatives);
}

TEST(AerosolPropertyProvider, ComputeValueAssignment)
{
    using MatrixPolicy = micm::VectorMatrix<double>;

    AerosolPropertyProvider<MatrixPolicy> provider;
    provider.dependent_variable_indices = { 0, 2 };

    bool was_called = false;
    provider.ComputeValue = [&was_called](
                                const MatrixPolicy& /* params */,
                                const MatrixPolicy& /* vars */,
                                MatrixPolicy& /* result */)
    { was_called = true; };

    MatrixPolicy params(1, 1, 0.0);
    MatrixPolicy vars(1, 1, 0.0);
    MatrixPolicy result(1, 1, 0.0);

    provider.ComputeValue(params, vars, result);
    EXPECT_TRUE(was_called);
    EXPECT_EQ(provider.dependent_variable_indices.size(), 2);
    EXPECT_EQ(provider.dependent_variable_indices[0], 0);
    EXPECT_EQ(provider.dependent_variable_indices[1], 2);
}

TEST(AerosolPropertyProvider, ComputeValueAndDerivativesAssignment)
{
    using MatrixPolicy = micm::VectorMatrix<double>;

    AerosolPropertyProvider<MatrixPolicy> provider;
    provider.dependent_variable_indices = { 1 };

    bool was_called = false;
    provider.ComputeValueAndDerivatives =
        [&was_called](
            const MatrixPolicy& /* params */,
            const MatrixPolicy& /* vars */,
            MatrixPolicy& /* result */,
            MatrixPolicy& /* partials */)
    { was_called = true; };

    MatrixPolicy params(1, 1, 0.0);
    MatrixPolicy vars(1, 1, 0.0);
    MatrixPolicy result(1, 1, 0.0);
    MatrixPolicy partials(1, 1, 0.0);

    provider.ComputeValueAndDerivatives(params, vars, result, partials);
    EXPECT_TRUE(was_called);
}

TEST(AerosolPropertyProvider, StorableInProviderMap)
{
    using MatrixPolicy = micm::VectorMatrix<double>;
    using ProviderMap = std::map<std::string, std::vector<AerosolPropertyProvider<MatrixPolicy>>>;

    AerosolPropertyProvider<MatrixPolicy> provider1;
    provider1.dependent_variable_indices = { 0 };
    provider1.ComputeValue = [](const MatrixPolicy&, const MatrixPolicy&, MatrixPolicy&) {};

    AerosolPropertyProvider<MatrixPolicy> provider2;
    provider2.dependent_variable_indices = { 1, 2 };
    provider2.ComputeValue = [](const MatrixPolicy&, const MatrixPolicy&, MatrixPolicy&) {};

    ProviderMap providers;
    providers["DROPLET"].push_back(provider1);
    providers["DROPLET"].push_back(provider2);

    EXPECT_EQ(providers.size(), 1);
    EXPECT_EQ(providers["DROPLET"].size(), 2);
    EXPECT_EQ(providers["DROPLET"][0].dependent_variable_indices.size(), 1);
    EXPECT_EQ(providers["DROPLET"][1].dependent_variable_indices.size(), 2);
}
