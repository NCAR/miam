#include <miam/processes/constants/henrys_law_constant.hpp>

#include <micm/system/conditions.hpp>

#include <cmath>

#include <gtest/gtest.h>

using namespace miam::process::constant;

TEST(HenrysLawConstant, DefaultParameters)
{
  HenrysLawConstant hlc;
  EXPECT_DOUBLE_EQ(hlc.parameters_.HLC_ref_, 1.0);
  EXPECT_DOUBLE_EQ(hlc.parameters_.C_, 0.0);
  EXPECT_DOUBLE_EQ(hlc.parameters_.T0_, 298.15);
}

TEST(HenrysLawConstant, CustomParameters)
{
  HenrysLawConstantParameters params{ .HLC_ref_ = 3.4e-2, .C_ = 2500.0, .T0_ = 298.15 };
  HenrysLawConstant hlc(params);
  EXPECT_DOUBLE_EQ(hlc.parameters_.HLC_ref_, 3.4e-2);
  EXPECT_DOUBLE_EQ(hlc.parameters_.C_, 2500.0);
  EXPECT_DOUBLE_EQ(hlc.parameters_.T0_, 298.15);
}

TEST(HenrysLawConstant, AtReferenceTemperature)
{
  // At reference temperature, exp term should be 1, so HLC = HLC_ref
  HenrysLawConstantParameters params{ .HLC_ref_ = 3.4e-2, .C_ = 2500.0, .T0_ = 298.15 };
  HenrysLawConstant hlc(params);

  EXPECT_DOUBLE_EQ(hlc.Calculate(298.15), 3.4e-2);

  micm::Conditions conditions;
  conditions.temperature_ = 298.15;
  EXPECT_DOUBLE_EQ(hlc.Calculate(conditions), 3.4e-2);
}

TEST(HenrysLawConstant, ZeroTemperatureDependence)
{
  // When C = 0, HLC is constant regardless of temperature
  HenrysLawConstantParameters params{ .HLC_ref_ = 5.0e-3, .C_ = 0.0, .T0_ = 298.15 };
  HenrysLawConstant hlc(params);

  EXPECT_DOUBLE_EQ(hlc.Calculate(200.0), 5.0e-3);
  EXPECT_DOUBLE_EQ(hlc.Calculate(298.15), 5.0e-3);
  EXPECT_DOUBLE_EQ(hlc.Calculate(400.0), 5.0e-3);
}

TEST(HenrysLawConstant, TemperatureDependence)
{
  // HLC(T) = HLC_ref * exp(C * (1/T - 1/T0))
  // For SO2: HLC_ref = 1.23 mol/(m³·Pa), C = 3120 K, T0 = 298.15 K
  HenrysLawConstantParameters params{ .HLC_ref_ = 1.23, .C_ = 3120.0, .T0_ = 298.15 };
  HenrysLawConstant hlc(params);

  double T = 280.0;
  double expected = 1.23 * std::exp(3120.0 * (1.0 / T - 1.0 / 298.15));
  EXPECT_DOUBLE_EQ(hlc.Calculate(T), expected);

  // At lower temperature (positive C), HLC should increase (more soluble)
  EXPECT_GT(hlc.Calculate(280.0), hlc.Calculate(298.15));

  // At higher temperature (positive C), HLC should decrease (less soluble)
  EXPECT_LT(hlc.Calculate(320.0), hlc.Calculate(298.15));
}

TEST(HenrysLawConstant, ConditionsOverload)
{
  HenrysLawConstantParameters params{ .HLC_ref_ = 0.8, .C_ = 2000.0, .T0_ = 298.15 };
  HenrysLawConstant hlc(params);

  micm::Conditions conditions;
  conditions.temperature_ = 310.0;

  double expected = 0.8 * std::exp(2000.0 * (1.0 / 310.0 - 1.0 / 298.15));
  EXPECT_DOUBLE_EQ(hlc.Calculate(conditions), expected);
  EXPECT_DOUBLE_EQ(hlc.Calculate(conditions), hlc.Calculate(310.0));
}

TEST(HenrysLawConstant, NegativeTemperatureDependence)
{
  // Some species have negative C (solubility increases with temperature)
  HenrysLawConstantParameters params{ .HLC_ref_ = 2.0, .C_ = -1500.0, .T0_ = 298.15 };
  HenrysLawConstant hlc(params);

  // At lower temperature, HLC should decrease (negative C)
  EXPECT_LT(hlc.Calculate(280.0), hlc.Calculate(298.15));

  // At higher temperature, HLC should increase (negative C)
  EXPECT_GT(hlc.Calculate(320.0), hlc.Calculate(298.15));
}

TEST(HenrysLawConstant, KnownValues)
{
  // Verify against hand-calculated value
  // HLC_ref = 1.0, C = 1000.0, T0 = 300.0, T = 250.0
  // 1/T - 1/T0 = 1/250 - 1/300 = 0.004 - 0.003333... = 0.000666...
  // exp(1000 * 0.000666...) = exp(0.666...) ≈ 1.94773
  HenrysLawConstantParameters params{ .HLC_ref_ = 1.0, .C_ = 1000.0, .T0_ = 300.0 };
  HenrysLawConstant hlc(params);

  double result = hlc.Calculate(250.0);
  double expected = std::exp(1000.0 * (1.0 / 250.0 - 1.0 / 300.0));
  EXPECT_DOUBLE_EQ(result, expected);
  EXPECT_NEAR(result, 1.9477, 0.001);
}
