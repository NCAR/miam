// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <miam/representations/aerosol_property.hpp>
#include <miam/representations/aerosol_property_descriptor.hpp>

#include <micm/util/vector_matrix.hpp>

#include <gtest/gtest.h>

#include <map>

using namespace miam;

TEST(AerosolProperty, EnumValues)
{
  EXPECT_NE(static_cast<int>(AerosolProperty::EffectiveRadius), static_cast<int>(AerosolProperty::NumberConcentration));
  EXPECT_NE(static_cast<int>(AerosolProperty::EffectiveRadius), static_cast<int>(AerosolProperty::PhaseVolumeFraction));
  EXPECT_NE(static_cast<int>(AerosolProperty::NumberConcentration), static_cast<int>(AerosolProperty::PhaseVolumeFraction));
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

TEST(AerosolPropertyDescriptor, ConstructAndQuery)
{
  using MatrixPolicy = micm::VectorMatrix<double>;
  SingleMomentModeEffectiveRadiusDescriptor<MatrixPolicy> effective_radius{ 0, 1 };
  EXPECT_TRUE(effective_radius.DependentVariableIndices().empty());

  SingleMomentModeNumberConcentrationDescriptor<MatrixPolicy> number_concentration{ 0, 1, { 2, 3 }, { 1.0e-5, 2.0e-5 } };
  ASSERT_EQ(number_concentration.DependentVariableIndices().size(), 2u);
  EXPECT_EQ(number_concentration.DependentVariableIndices()[0], 2u);
  EXPECT_EQ(number_concentration.DependentVariableIndices()[1], 3u);

  PhaseVolumeFractionDescriptor<MatrixPolicy> phase_volume_fraction{ { 4, 5 }, { 1.0, 2.0 }, 1, { 4, 5 } };
  ASSERT_EQ(phase_volume_fraction.DependentVariableIndices().size(), 2u);
}

TEST(AerosolPropertyDescriptor, StorableInDescriptorMap)
{
  using MatrixPolicy = micm::VectorMatrix<double>;
  using Descriptor = AerosolPropertyDescriptor<MatrixPolicy>;
  using DescriptorMap = std::map<std::string, std::map<AerosolProperty, Descriptor>>;

  DescriptorMap descriptors;
  descriptors["DROPLET"][AerosolProperty::EffectiveRadius] =
      SingleMomentModeEffectiveRadiusDescriptor<MatrixPolicy>{ 0, 1 };
  descriptors["DROPLET"][AerosolProperty::NumberConcentration] =
      SingleMomentModeNumberConcentrationDescriptor<MatrixPolicy>{ 0, 1, { 2, 3 }, { 1.0e-5, 2.0e-5 } };

  EXPECT_EQ(descriptors.size(), 1u);
  EXPECT_EQ(descriptors["DROPLET"].size(), 2u);
  EXPECT_EQ(DependentVariableIndices(descriptors["DROPLET"].at(AerosolProperty::NumberConcentration)).size(), 2u);
}
