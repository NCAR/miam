// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include "cam_cloud_chemistry_policy.hpp"

#include <miam/miam.hpp>

#include <micm/CPU.hpp>

#include <gtest/gtest.h>

namespace policy = miam_test_cam_cloud_chemistry;

namespace
{
  template<std::size_t L>
  using VectorBuilder = micm::CpuSolverBuilder<
      micm::RosenbrockSolverParameters,
      micm::VectorMatrix<double, L>,
      micm::SparseMatrix<double, micm::SparseMatrixVectorOrdering<L>>>;

  template<std::size_t L>
  VectorBuilder<L> MakeVectorBuilder()
  {
    return VectorBuilder<L>(micm::RosenbrockSolverParameters::FourStageDifferentialAlgebraicRosenbrockParameters());
  }
}  // namespace

TEST(VectorMatrixModel, ConstraintsOnly_L1)
{
  policy::Step1b_KwOnly(MakeVectorBuilder<1>());
}

TEST(VectorMatrixModel, ConstraintsOnly_L4)
{
  policy::Step1b_KwOnly(MakeVectorBuilder<4>());
}

TEST(VectorMatrixModel, ProcessesAndConstraints_L1)
{
  policy::Step4_FullSystemWithKinetics(MakeVectorBuilder<1>());
}

TEST(VectorMatrixModel, ProcessesAndConstraints_L4)
{
  policy::Step4_FullSystemWithKinetics(MakeVectorBuilder<4>());
}
