// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Kokkos driver for the CAM cloud aqueous-phase sulfate-chemistry integration
// tests. Each TEST() forwards to a templated policy function in
// cam_cloud_chemistry_policy.hpp with a `KokkosSolverBuilder`. This exercises
// MIAM's Model interface end-to-end with Kokkos-backed dense and sparse
// matrix classes.

#include "../cam_cloud_chemistry_policy.hpp"

#include <miam/miam.hpp>

#include <micm/Kokkos.hpp>

#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>

namespace policy = miam_test_cam_cloud_chemistry;

using KokkosBuilder = micm::KokkosSolverBuilder<micm::RosenbrockSolverParameters>;

namespace
{
  KokkosBuilder StandardBuilder()
  {
    return KokkosBuilder(micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters());
  }

  KokkosBuilder Step1cBuilder()
  {
    auto params = micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters();
    params.constraint_init_max_iterations_ = 200;
    params.constraint_init_tolerance_ = 1e-20;
    params.max_number_of_steps_ = 1500;
    return KokkosBuilder(params);
  }
}  // namespace

TEST(KokkosCamCloudChemistry, Step1_SingleHLC)
{
  policy::Step1_SingleHLC(StandardBuilder());
}

TEST(KokkosCamCloudChemistry, Step1b_KwOnly)
{
  policy::Step1b_KwOnly(StandardBuilder());
}

TEST(KokkosCamCloudChemistry, Step1c_KwNaiveIC)
{
  policy::Step1c_KwNaiveIC(Step1cBuilder());
}

TEST(KokkosCamCloudChemistry, Step2_HLC_Plus_Dissociation)
{
  policy::Step2_HLC_Plus_Dissociation(StandardBuilder());
}

TEST(KokkosCamCloudChemistry, Step3_FullEquilibrium)
{
  policy::Step3_FullEquilibrium(StandardBuilder());
}

TEST(KokkosCamCloudChemistry, Step3b_NaiveInitialConditions)
{
  policy::Step3b_NaiveInitialConditions(StandardBuilder());
}

TEST(KokkosCamCloudChemistry, Step4_FullSystemWithKinetics)
{
  policy::Step4_FullSystemWithKinetics(StandardBuilder());
}

TEST(KokkosCamCloudChemistry, Step4b_NaiveInitialConditions)
{
  policy::Step4b_NaiveInitialConditions(StandardBuilder());
}

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  Kokkos::initialize(argc, argv);
  int result = RUN_ALL_TESTS();
  Kokkos::finalize();
  return result;
}
