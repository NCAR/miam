// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#pragma once

#include <miam/representations/aerosol_property.hpp>
#include <miam/representations/phase_volume_fraction_descriptor.hpp>
#include <miam/representations/single_moment_mode.hpp>
#include <miam/representations/two_moment_mode.hpp>
#include <miam/representations/uniform_section.hpp>

#include <concepts>
#include <utility>
#include <variant>
#include <vector>

namespace miam
{
  /// @brief Concrete descriptors expose these three methods (checked structurally).
  /// @details Each per-representation descriptor type in this file's transitive includes satisfies
  ///          the concept, as does each host- or test-side synthetic descriptor.
  template<class Descriptor, class DenseMatrixPolicy>
  concept AerosolPropertyDescriptorLike =
      requires(const Descriptor& d, const DenseMatrixPolicy& p, DenseMatrixPolicy& r, DenseMatrixPolicy& jac) {
        { d.Evaluate(p, p, r) };
        { d.EvaluateAndDerivatives(p, p, r, jac) };
        { d.DependentVariableIndices() } -> std::convertible_to<const std::vector<std::size_t>&>;
      };

  /// @brief Closed set of production descriptor types accepted by `Model::BuildDescriptors`.
  /// @details Tests may pass their own concrete descriptor type or a different variant; the
  ///          overloads below accept both shapes.
  template<typename DenseMatrixPolicy>
  using AerosolPropertyDescriptor = std::variant<
      SingleMomentModeEffectiveRadiusDescriptor<DenseMatrixPolicy>,
      SingleMomentModeNumberConcentrationDescriptor<DenseMatrixPolicy>,
      TwoMomentModeEffectiveRadiusDescriptor<DenseMatrixPolicy>,
      TwoMomentModeNumberConcentrationDescriptor<DenseMatrixPolicy>,
      UniformSectionEffectiveRadiusDescriptor<DenseMatrixPolicy>,
      UniformSectionNumberConcentrationDescriptor<DenseMatrixPolicy>,
      PhaseVolumeFractionDescriptor<DenseMatrixPolicy>>;

  /// @brief Widen a rep-local descriptor variant into `AerosolPropertyDescriptor<Policy>`.
  /// @details Consumed by `Model::BuildDescriptors`; each representation's `GetPropertyDescriptor`
  ///          returns its own three-alternative variant which is coerced here into the wider
  ///          production variant that consumer methods take.
  template<typename DenseMatrixPolicy, class... Ds>
  AerosolPropertyDescriptor<DenseMatrixPolicy> WidenDescriptorVariant(std::variant<Ds...> narrow)
  {
    return std::visit(
        [](auto&& d) -> AerosolPropertyDescriptor<DenseMatrixPolicy> { return { std::move(d) }; }, std::move(narrow));
  }

  template<class Descriptor, class DenseMatrixPolicy>
    requires AerosolPropertyDescriptorLike<Descriptor, DenseMatrixPolicy>
  void EvaluateAerosolProperty(
      const Descriptor& descriptor,
      const DenseMatrixPolicy& state_parameters,
      const DenseMatrixPolicy& state_variables,
      DenseMatrixPolicy& result)
  {
    descriptor.Evaluate(state_parameters, state_variables, result);
  }

  template<class Descriptor, class DenseMatrixPolicy>
    requires AerosolPropertyDescriptorLike<Descriptor, DenseMatrixPolicy>
  void EvaluateAerosolPropertyAndDerivatives(
      const Descriptor& descriptor,
      const DenseMatrixPolicy& state_parameters,
      const DenseMatrixPolicy& state_variables,
      DenseMatrixPolicy& result,
      DenseMatrixPolicy& partials)
  {
    descriptor.EvaluateAndDerivatives(state_parameters, state_variables, result, partials);
  }

  template<class Descriptor, class DenseMatrixPolicy>
    requires AerosolPropertyDescriptorLike<Descriptor, DenseMatrixPolicy>
  const std::vector<std::size_t>& DependentVariableIndices(const Descriptor& descriptor)
  {
    return descriptor.DependentVariableIndices();
  }

  template<class... Ds, class DenseMatrixPolicy>
  void EvaluateAerosolProperty(
      const std::variant<Ds...>& descriptor,
      const DenseMatrixPolicy& state_parameters,
      const DenseMatrixPolicy& state_variables,
      DenseMatrixPolicy& result)
  {
    std::visit([&](const auto& d) { d.Evaluate(state_parameters, state_variables, result); }, descriptor);
  }

  template<class... Ds, class DenseMatrixPolicy>
  void EvaluateAerosolPropertyAndDerivatives(
      const std::variant<Ds...>& descriptor,
      const DenseMatrixPolicy& state_parameters,
      const DenseMatrixPolicy& state_variables,
      DenseMatrixPolicy& result,
      DenseMatrixPolicy& partials)
  {
    std::visit(
        [&](const auto& d) { d.EvaluateAndDerivatives(state_parameters, state_variables, result, partials); },
        descriptor);
  }

  template<class... Ds>
  const std::vector<std::size_t>& DependentVariableIndices(const std::variant<Ds...>& descriptor)
  {
    return std::visit(
        [](const auto& d) -> const std::vector<std::size_t>& { return d.DependentVariableIndices(); }, descriptor);
  }
}  // namespace miam
