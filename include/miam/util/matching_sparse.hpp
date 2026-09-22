// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <micm/util/matrix.hpp>
#include <micm/util/sparse_matrix.hpp>
#include <micm/util/sparse_matrix_standard_ordering_compressed_sparse_row.hpp>

#ifdef MICM_USE_KOKKOS
  #include <micm/kokkos/util/kokkos_dense_matrix.hpp>
  #include <micm/kokkos/util/kokkos_sparse_matrix.hpp>
  #include <micm/util/sparse_matrix_vector_ordering_compressed_sparse_row.hpp>
#endif

namespace miam::detail
{
  /// @brief Maps a `DenseMatrixPolicy` to the `SparseMatrixPolicy` MIAM pairs with it.
  /// @details MICM's ExternalModel dispatcher hands MIAM only DP in `AddForcingTerms<DP>`,
  ///          but Sets are templated on `<DP, SP>`. This trait recovers SP from DP for the
  ///          combinations MIAM supports (CPU + Kokkos-serial/CUDA).
  template<class DenseMatrixPolicy>
  struct MatchingSparse;

  template<class T>
  struct MatchingSparse<micm::Matrix<T>>
  {
    using type = micm::SparseMatrix<T, micm::SparseMatrixStandardOrderingCompressedSparseRow>;
  };

#ifdef MICM_USE_KOKKOS
  template<class T, micm::Index L>
  struct MatchingSparse<micm::KokkosDenseMatrix<T, L>>
  {
    using type = micm::KokkosSparseMatrix<T, micm::SparseMatrixVectorOrderingCompressedSparseRow<L>>;
  };
#endif

  template<class DenseMatrixPolicy>
  using MatchingSparseT = typename MatchingSparse<DenseMatrixPolicy>::type;

  /// @brief Inverse of `MatchingSparse`. `Model::FinalizeConstraintSetup<SP>` receives only SP;
  ///        this trait recovers the DP MIAM expects to pair with it, so the SP-typed collection
  ///        can be templated as `<DP, SP>` at Finalize time.
  template<class SparseMatrixPolicy>
  struct MatchingDense;

  template<class T, class OrderingPolicy>
  struct MatchingDense<micm::SparseMatrix<T, OrderingPolicy>>
  {
    using type = micm::Matrix<T>;
  };

#ifdef MICM_USE_KOKKOS
  template<class T, micm::Index L>
  struct MatchingDense<micm::KokkosSparseMatrix<T, micm::SparseMatrixVectorOrderingCompressedSparseRow<L>>>
  {
    using type = micm::KokkosDenseMatrix<T, L>;
  };
#endif

  template<class SparseMatrixPolicy>
  using MatchingDenseT = typename MatchingDense<SparseMatrixPolicy>::type;
}  // namespace miam::detail
