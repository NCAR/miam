# Kokkos-Compatible `Model` — Design & Implementation Plan

## Overview

Rework MIAM's `Model`, per-process, and per-constraint classes so that a single
`Model` instance dispatches into MICM's solver hot loop without any
`std::function` on the solve path, and works uniformly with `micm::Matrix`,
`micm::VectorMatrix`, and `micm::KokkosDenseMatrix` / `micm::KokkosSparseMatrix`
backings. The design mirrors MICM's own `micm::ProcessSet` and the
`test/integration/stub_aerosol_*` external-model examples (MICM PR
[NCAR/micm#1080](https://github.com/NCAR/micm/pull/1080)).

The end state:

- Inside `Solver::Solve()`, execution is entirely on device. `AddForcingTerms`,
  `SubtractJacobianTerms`, `AddConstraintResidual`, and
  `SubtractConstraintJacobian` never touch host state.
- Inside `Solver::UpdateStateParameters()`, host work is a single loop over
  expression buckets, followed by a single `Kokkos::parallel_for` per external
  model.
- No `std::function` anywhere in the process, constraint, aerosol-property, or
  `Model` code paths.

## Diagnosis

The current `Model` (added in the "new-MICM-external-model-API" PR) caches
`std::vector<std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)>>`
inside seven `mutable std::any` members. Each stored `std::function` wraps a
lambda produced by `DenseMatrixPolicy::Function(...)`. Inside those lambdas the
per-process code captures:

- `this` by pointer.
- `StateVariableIndices` and `JacobianIndices`, both of which hold a
  `micm::Matrix<std::size_t>` — a non-trivially-copyable owning matrix.
- `rate_constant_`, a `std::function<double(const micm::Conditions&)>`, invoked
  from inside another lambda.

For `micm::Matrix` this works because `DenseMatrixPolicy::Function(...)` runs
the lambda on the host, so every capture stays on the host side. For
`micm::KokkosDenseMatrix`, `Function(...)` launches a
`Kokkos::parallel_for(TeamPolicy)` and the lambda body executes on device. At
that point every one of the captures above is either undefined-behaviour
device code (`std::function`, `micm::Matrix` heap owner) or a host pointer
masquerading as device state (`this`).

The failure mode observed in `test_kokkos_cam_cloud_chemistry` (introduced in
the "kokkos-cam-cloud-chemistry-test" PR) is exactly this: everything compiles,
everything links, all eight steps return `SolverState != Converged` at
`t = 0` because the solve path is reading garbage from device memory that was
supposed to hold `this`-relative offsets.

## Target architecture

The refactor is broken into four independently-shippable layers.

### Layer A — POD rate/equilibrium/HLC "expression" classes

Each temperature-dependent constant becomes its own POD class with a
`MICM_DEVICE_FUNCTION Calculate(const micm::Conditions&) const` method. Every
field of every expression class is `micm::Real` or `micm::Index`, so the class
is trivially copyable and captures cleanly into a `MICM_LAMBDA`.

```cpp
namespace miam
{
  class ArrheniusExpression
  {
   public:
    micm::ArrheniusRateConstantParameters params_;
    MICM_DEVICE_FUNCTION micm::Real Calculate(const micm::Conditions& conditions) const
    { return micm::CalculateArrhenius(params_, conditions.temperature_, conditions.pressure_); }
  };

  class VantHoffExpression
  {
   public:
    micm::Real A_ = 1.0, C_ = 0.0, T0_ = 298.15;
    MICM_DEVICE_FUNCTION micm::Real Calculate(const micm::Conditions& conditions) const;
  };

  class HenrysLawExpression
  {
   public:
    micm::Real HLC_ref_ = 0.0, C_ = 0.0, T0_ = 298.15;
    MICM_DEVICE_FUNCTION micm::Real Calculate(const micm::Conditions& conditions) const;
  };

  class UserDefinedConstantExpression
  {
   public:
    micm::Real value_ = 0.0;
    MICM_DEVICE_FUNCTION micm::Real Calculate(const micm::Conditions&) const { return value_; }
  };
}
```

The list of allowed alternatives per site is fixed at compile time in a
`std::variant`, one variant per constant kind:

```cpp
namespace miam
{
  using RateConstantExpression        = std::variant<ArrheniusExpression, UserDefinedConstantExpression>;
  using EquilibriumConstantExpression = std::variant<VantHoffExpression, UserDefinedConstantExpression>;
  using HenrysLawConstantExpression   = std::variant<HenrysLawExpression, UserDefinedConstantExpression>;
}
```

The config-level variant is **host-only**. It is populated by the builders on
the host, visited exactly once at `Model::FinalizeProcessSetup` /
`Model::FinalizeConstraintSetup` time to append into the constants buckets
(Layer D), and never inspected at solve time.

**Every `std::function<double(const micm::Conditions&)>` on a MIAM config
class is deleted.** The six affected members are:

| Config class | Member |
|-|-|
| `DissolvedReaction` | `rate_constant_` |
| `DissolvedReversibleReaction` | `forward_rate_constant_`, `reverse_rate_constant_`, `equilibrium_constant_` |
| `HenrysLawPhaseTransfer` | `henrys_law_constant_` |
| `DissolvedEquilibriumConstraint` | `equilibrium_constant_` |
| `HenrysLawEquilibriumConstraint` | `henrys_law_constant_` |

Every `SetRateConstant(std::function<double(...)>)`,
`SetEquilibriumConstant(std::function<double(...)>)`, and
`SetHenrysLawConstant(std::function<double(...)>)` builder overload is deleted
outright. The new overloads accept an expression variant (or one of its
alternatives) directly. Builders that previously accepted concrete parameter
types such as `micm::ArrheniusRateConstantParameters` keep working, since
those overloads simply wrap the parameters into the corresponding expression
class.

`DissolvedReversibleReactionBuilder::Build()`'s current "derive the missing
rate constant" logic — which today wraps two `std::function`s into a third —
moves into `Update`-time evaluation: when both `equilibrium_constant_` and
`reverse_rate_constant_` are set, the missing `forward_rate_constant_` is
computed inside the `MICM_LAMBDA` as `K_eq * k_r`, with no closure wrapping.

### Layer B — POD aerosol property descriptors

`include/miam/representations/aerosol_property.hpp` currently defines an
`AerosolPropertyProvider<DenseMatrixPolicy>` struct with two
`std::function` fields (`ComputeValue`, `ComputeValueAndDerivatives`). Each
representation (`SingleMomentMode`, `TwoMomentMode`, `UniformSection`) fills
those `std::function`s by wrapping `DenseMatrixPolicy::Function(...)`. The
only consumer, `HenrysLawPhaseTransfer`, calls
`provider.ComputeValueAndDerivatives(...)` from inside its own
`ForcingFunction` lambda — a `std::function` call from inside a lambda that
is itself captured in another `std::function`.

`AerosolPropertyProvider` is deleted. It is replaced by a tagged POD
descriptor plus a pair of `MICM_DEVICE_FUNCTION` free evaluators that
`HenrysLawPhaseTransfer` calls inline from its own `MICM_LAMBDA`:

```cpp
namespace miam
{
  enum class RepresentationKind : std::uint8_t
  { SingleMomentMode, TwoMomentMode, UniformSection };

  struct AerosolPropertyDescriptor
  {
    AerosolProperty     property_;
    RepresentationKind  representation_;
    // Column indices into state_parameters / state_variables used by the
    // evaluator (unused slots are `-1`).
    std::array<micm::Index, 4> parameter_indices_{};
    std::array<micm::Index, 4> variable_indices_{};
    std::array<micm::Index, 4> dependent_variable_indices_{};
    std::uint8_t number_of_dependent_variables_ = 0;
  };

  MICM_DEVICE_FUNCTION inline void EvaluateAerosolProperty(
      const AerosolPropertyDescriptor& descriptor,
      const auto& parameters_view,
      const auto& variables_view,
      auto& value_out);

  MICM_DEVICE_FUNCTION inline void EvaluateAerosolPropertyAndDerivatives(
      const AerosolPropertyDescriptor& descriptor,
      const auto& parameters_view,
      const auto& variables_view,
      auto& value_out,
      auto& partials_out);
}
```

Each representation's `GetPropertyProvider<DenseMatrixPolicy>(...)` is
renamed to `GetPropertyDescriptor(...)` and returns a POD
`AerosolPropertyDescriptor` — no template parameter, no `std::function`, no
allocation. `Model::BuildProviders<DenseMatrixPolicy>` disappears; the
equivalent host-side descriptor gathering happens once at Finalize and is
handed to each process's Set constructor.

### Layer C — Per-process and per-constraint "Set" objects

Mirror `micm::ProcessSet<DenseMatrixPolicy, SparseMatrixPolicy>` and MICM's
`stub_aerosol_1` external model. For each of the three process types and
three constraint types, introduce a `<Config>Set<DenseMatrixPolicy,
SparseMatrixPolicy>` template class that owns policy-typed index caches:

```cpp
namespace miam
{
  template<class DenseMatrixPolicy, class SparseMatrixPolicy>
  class DissolvedReactionSet
  {
   private:
    template<class U>
    using Vector = typename SparseMatrixPolicy::template VectorType<U>;

    // Flattened per-phase caches (num_phases * per-phase-stride).
    Vector<micm::Index> reactant_indices_;
    Vector<micm::Index> product_indices_;
    Vector<micm::Index> solvent_indices_;
    Vector<micm::Index> jacobian_flat_ids_;

    micm::Index k_state_parameter_index_ = -1;
    micm::Index num_phases_ = 0;
    micm::Index num_reactants_ = 0;
    micm::Index num_products_ = 0;
    micm::Real  solvent_floor_ = 1.0e-20;
    micm::Real  min_halflife_ = 0.0;

   public:
    DissolvedReactionSet(
        const DissolvedReaction& config,
        const std::map<std::string, std::set<std::string>>& phase_prefixes,
        const std::unordered_map<std::string, micm::Index>& state_parameter_indices,
        const std::unordered_map<std::string, micm::Index>& state_variable_indices,
        const std::unordered_map<std::string, std::map<AerosolProperty, AerosolPropertyDescriptor>>& providers,
        const SparseMatrixPolicy& jacobian);

    void AddForcingTerms(
        const DenseMatrixPolicy& state_parameters,
        const DenseMatrixPolicy& state_variables,
        DenseMatrixPolicy& forcing) const;

    void SubtractJacobianTerms(
        const DenseMatrixPolicy& state_parameters,
        const DenseMatrixPolicy& state_variables,
        SparseMatrixPolicy& jacobian) const;
  };
}
```

`Vector<U>` is the same alias `micm::ProcessSet` uses. On CPU it resolves to
`micm::PaddedVector<U, 1>`; on Kokkos it resolves to
`micm::KokkosPaddedVector<U, L>`. Both expose `.CopyToDevice()`,
`.CopyToHost()`, and a device-safe `.GetView()`.

The constructor does every string-map lookup, calls `jacobian.VectorIndex(...)`
for every non-zero Jacobian element, appends into the four `Vector<micm::Index>`
members, and calls `.CopyToDevice()` on each — exactly the pattern of
`micm::ProcessSet::ProcessSet(...)` + `SetJacobianFlatIds(...)`.

`AddForcingTerms` iterates phases on the **host** (matching `stub_aerosol_1`)
and per phase issues one `DenseMatrixPolicy::Function(...)` call. Only
trivially copyable data (a couple of `micm::Index` values and
`Vector<micm::Index>::GetView()` handles) crosses into the `MICM_LAMBDA`:

```cpp
template<class DenseMatrixPolicy, class SparseMatrixPolicy>
void DissolvedReactionSet<DenseMatrixPolicy, SparseMatrixPolicy>::AddForcingTerms(
    const DenseMatrixPolicy& state_parameters,
    const DenseMatrixPolicy& state_variables,
    DenseMatrixPolicy& forcing) const
{
  const auto reactant_view = reactant_indices_.GetView();
  const auto product_view  = product_indices_.GetView();
  const auto solvent_view  = solvent_indices_.GetView();
  const micm::Index num_reactants        = num_reactants_;
  const micm::Index num_products         = num_products_;
  const micm::Real  solvent_floor        = solvent_floor_;
  const micm::Index k_parameter_index    = k_state_parameter_index_;

  for (micm::Index phase = 0; phase < num_phases_; ++phase)
  {
    DenseMatrixPolicy::Function(
        MICM_LAMBDA(
            const typename DenseMatrixPolicy::ViewType&      forcing_view,
            const typename DenseMatrixPolicy::ConstViewType& state_view,
            const typename DenseMatrixPolicy::ConstViewType& params_view)
        {
          auto rate = forcing_view.GetRowVariable();
          forcing_view.ForEachRow(
              [=](micm::Real& out, const micm::Real& k, const micm::Real& solvent)
              { out = k * solvent / std::pow(solvent + solvent_floor, num_reactants); },
              rate,
              params_view.GetConstColumnView(k_parameter_index),
              state_view.GetConstColumnView(solvent_view[phase]));
          for (micm::Index r = 0; r < num_reactants; ++r)
            forcing_view.ForEachRow(
                [](micm::Real& out, const micm::Real& reactant) { out *= reactant; },
                rate,
                state_view.GetConstColumnView(reactant_view[phase * num_reactants + r]));
          for (micm::Index r = 0; r < num_reactants; ++r)
            state_view.ForEachRow(
                [](const micm::Real& rate, micm::Real& forcing) { forcing -= rate; },
                rate,
                forcing_view.GetColumnView(reactant_view[phase * num_reactants + r]));
          for (micm::Index p = 0; p < num_products; ++p)
            state_view.ForEachRow(
                [](const micm::Real& rate, micm::Real& forcing) { forcing += rate; },
                rate,
                forcing_view.GetColumnView(product_view[phase * num_products + p]));
        },
        forcing, state_variables, state_parameters)(forcing, state_variables, state_parameters);
  }
}
```

`SubtractJacobianTerms` reads pre-computed flat IDs from
`jacobian_flat_ids_.GetView()` and calls `jacobian_view.GetBlockView(flat_id)`.

`DissolvedReversibleReactionSet` and `HenrysLawPhaseTransferSet` follow the
same shape with their own `Vector<micm::Index>` caches; the constraint Sets
(`DissolvedEquilibriumConstraintSet`, `HenrysLawEquilibriumConstraintSet`,
`LinearConstraintSet`) each expose `AddConstraintResidual` and
`SubtractConstraintJacobian` instead of `AddForcingTerms` and
`SubtractJacobianTerms`, matching MICM's
`stub_aerosol_with_constraints::AddConstraintResidual` and
`SubtractConstraintJacobian`.

`HenrysLawPhaseTransferSet` additionally captures its
`AerosolPropertyDescriptor` by value into `MICM_LAMBDA` and calls
`EvaluateAerosolPropertyAndDerivatives(...)` inline from the kernel body.

**Note:** No Set class has an `UpdateStateParameters` method. That path is
handled by the shared constants buckets in Layer D.

### Layer D — `Model` collapse plus `ConstantsBucket`

The variant of process/constraint config objects on `Model` stays as-is (it
is matrix-policy-independent and defines the user-facing API). The heavy
state lives in per-`(DenseMatrixPolicy, SparseMatrixPolicy)` collections
built lazily on first solve-time call:

```cpp
namespace miam
{
  template<class DenseMatrixPolicy, class SparseMatrixPolicy>
  struct ProcessSetCollection
  {
    std::vector<DissolvedReactionSet<DenseMatrixPolicy, SparseMatrixPolicy>>            dissolved_reactions_;
    std::vector<DissolvedReversibleReactionSet<DenseMatrixPolicy, SparseMatrixPolicy>>  reversible_reactions_;
    std::vector<HenrysLawPhaseTransferSet<DenseMatrixPolicy, SparseMatrixPolicy>>       phase_transfers_;
  };

  template<class DenseMatrixPolicy, class SparseMatrixPolicy>
  struct ConstraintSetCollection
  {
    std::vector<DissolvedEquilibriumConstraintSet<DenseMatrixPolicy, SparseMatrixPolicy>> dissolved_equilibrium_constraints_;
    std::vector<HenrysLawEquilibriumConstraintSet<DenseMatrixPolicy, SparseMatrixPolicy>> henrys_law_equilibrium_constraints_;
    std::vector<LinearConstraintSet<DenseMatrixPolicy, SparseMatrixPolicy>>               linear_constraints_;
  };
}
```

`Model::FinalizeProcessSetup<SparseMatrixPolicy>(...)` stashes
`state_parameter_indices`, `state_variable_indices`, the sparse pattern (as
`std::any` holding the `SparseMatrixPolicy` itself), and resets any stale
D-dependent cache. `Model::FinalizeConstraintSetup<SparseMatrixPolicy>(...)`
does the same. `Model::AddForcingTerms<DenseMatrixPolicy>(...)`,
`SubtractJacobianTerms<DenseMatrixPolicy, SparseMatrixPolicy>(...)`,
`UpdateStateParameters<DenseMatrixPolicy>(...)`, and their three constraint
counterparts all delegate to the two collections, constructing them on first
call from the stashed inputs.

#### `ConstantsBucket`

To keep `UpdateStateParameters` on-device, add a **bucketed** constants
store alongside each `ProcessSetCollection` and `ConstraintSetCollection`.
One `Vector<>` per concrete expression class, all destination columns in
the state-parameters matrix collected into per-bucket `Vector<micm::Index>`s:

```cpp
namespace miam
{
  template<class DenseMatrixPolicy, class SparseMatrixPolicy, class... ExpressionTypes>
  class ConstantsBucket
  {
   private:
    template<class U>
    using Vector = typename SparseMatrixPolicy::template VectorType<U>;

    std::tuple<Vector<ExpressionTypes>...> expressions_;
    std::tuple<Vector<micm::Index>,
               /* ... one per expression ... */
               Vector<micm::Index>>       destination_columns_;

   public:
    // Host-only, called by Model::FinalizeProcessSetup after visiting each
    // config-level `std::variant<...>` on the host.
    template<class Expression>
    void Append(const Expression& expression, micm::Index destination_column);

    void CopyToDevice();

    // Called on-device: one `Kokkos::parallel_for` per external model.
    void UpdateStateParameters(
        const auto& conditions,
        DenseMatrixPolicy& state_parameters) const;
  };
}
```

`Append` is called on the host at Finalize time. It is the **only** place
that ever touches the config-level `std::variant<...>`, and it happens
exactly once per `Model` build.

`UpdateStateParameters` issues **one** `Kokkos::parallel_for` (via one
`DenseMatrixPolicy::Function(...)`) and, inside the single `MICM_LAMBDA`,
runs a compile-time-unrolled fold over the buckets:

```cpp
template<class DenseMatrixPolicy, class SparseMatrixPolicy, class... ExpressionTypes>
void ConstantsBucket<DenseMatrixPolicy, SparseMatrixPolicy, ExpressionTypes...>::UpdateStateParameters(
    const auto& conditions, DenseMatrixPolicy& state_parameters) const
{
  DenseMatrixPolicy::Function(
      MICM_LAMBDA(
          const typename DenseMatrixPolicy::ViewType&                          params_view,
          const typename decltype(conditions)::ConstViewType&                  conditions_view)
      {
        [&]<std::size_t... Is>(std::index_sequence<Is...>)
        {
          ( /* one bucket iteration per expression type */
            [&]
            {
              const auto expressions_view = std::get<Is>(expressions_).GetView();
              const auto columns_view     = std::get<Is>(destination_columns_).GetView();
              const micm::Index count     = static_cast<micm::Index>(expressions_view.size());
              for (micm::Index i = 0; i < count; ++i)
              {
                const auto expression         = expressions_view[i];
                const auto destination_column = columns_view[i];
                params_view.ForEachRow(
                    [=](micm::Real& out, const micm::Conditions& c)
                    { out = expression.Calculate(c); },
                    params_view.GetColumnView(destination_column),
                    conditions_view);
              }
            }()
          , ...);
        }(std::index_sequence_for<ExpressionTypes...>{});
      },
      state_parameters, conditions)(state_parameters, conditions);
}
```

Every capture is trivially copyable: `expression` is one alternative POD,
`destination_column` is a `micm::Index`, `count` is a `micm::Index`. No
`std::visit`, no `std::function`, no host allocation on the solve path.

Adding a new expression class (e.g. `TroeExpression`) is a small, targeted
change: write the POD class + its `Calculate` method, add it to the
relevant `std::variant<...>` alias, and add it to the `ExpressionTypes...`
pack of the corresponding `ConstantsBucket` instantiation.

#### `Model` public methods after the collapse

```cpp
class Model
{
 public:
  template<class SparseMatrixPolicy>
  void FinalizeProcessSetup(
      const std::unordered_map<std::string, micm::Index>& state_parameter_indices,
      const std::unordered_map<std::string, micm::Index>& state_variable_indices,
      const SparseMatrixPolicy& jacobian);

  template<class SparseMatrixPolicy>
  void FinalizeConstraintSetup(
      const std::unordered_map<std::string, micm::Index>& state_parameter_indices,
      const std::unordered_map<std::string, micm::Index>& state_variable_indices,
      const SparseMatrixPolicy& jacobian);

  template<class DenseMatrixPolicy>
  void UpdateStateParameters(
      const typename DenseMatrixPolicy::template VectorType<micm::Conditions>& conditions,
      DenseMatrixPolicy& state_parameters) const;

  template<class DenseMatrixPolicy>
  void AddForcingTerms(
      const DenseMatrixPolicy& state_parameters,
      const DenseMatrixPolicy& state_variables,
      DenseMatrixPolicy& forcing) const;

  template<class DenseMatrixPolicy, class SparseMatrixPolicy>
  void SubtractJacobianTerms(
      const DenseMatrixPolicy& state_parameters,
      const DenseMatrixPolicy& state_variables,
      SparseMatrixPolicy& jacobian) const;

  template<class DenseMatrixPolicy>
  void UpdateConstraintStateParameters(
      const typename DenseMatrixPolicy::template VectorType<micm::Conditions>& conditions,
      DenseMatrixPolicy& state_parameters) const;

  template<class DenseMatrixPolicy>
  void InitializeConstraintParameters(
      const DenseMatrixPolicy& state_variables,
      DenseMatrixPolicy& state_parameters) const;

  template<class DenseMatrixPolicy>
  void AddConstraintResidual(
      const DenseMatrixPolicy& state_parameters,
      const DenseMatrixPolicy& state_variables,
      DenseMatrixPolicy& forcing) const;

  template<class DenseMatrixPolicy, class SparseMatrixPolicy>
  void SubtractConstraintJacobian(
      const DenseMatrixPolicy& state_parameters,
      const DenseMatrixPolicy& state_variables,
      SparseMatrixPolicy& jacobian) const;
};
```

**Deleted by this collapse:**

- All seven `mutable std::any` caches on `Model` (`cached_process_update_fns_`,
  `cached_forcing_fns_`, `cached_jacobian_fns_`, plus four constraint mates).
- `Model::BuildProviders<DenseMatrixPolicy>`.
- Every `Model` factory that returned `std::function`:
  `UpdateStateParametersFunction`, `ForcingFunction`, `JacobianFunction`,
  `ConstraintUpdateStateParametersFunction`,
  `InitializeConstraintParametersFunction`, `ConstraintResidualFunction`,
  `ConstraintJacobianFunction`.
- Every per-process and per-constraint `std::function`-returning factory
  (`UpdateStateParametersFunction`, `ForcingFunction`, `JacobianFunction` on
  the three process configs; the constraint-side equivalents on the three
  constraint configs).

## Invariants that hold by design

1. **During `Solver::Solve()` (Rosenbrock hot loop): zero host residency.** All
   work is Kokkos kernels launched by `AddForcingTerms<DenseMatrixPolicy>` and
   `SubtractJacobianTerms<DenseMatrixPolicy, SparseMatrixPolicy>`. Every input
   is either a POD `micm::Index` / `micm::Real` captured by value or a
   `Vector<T>::GetView()` handle already resident on device.
2. **During `Solver::UpdateStateParameters()`: one Kokkos kernel launch per
   external model.** No `std::visit` at solve time, no `std::function`, no
   host allocation.
3. **Config-level `std::variant<...>` is host-only.** It is visited exactly
   once, at `Model::FinalizeProcessSetup` / `FinalizeConstraintSetup` time,
   and used only to populate the constants buckets. It never crosses to
   device.

## Per-file change list

| Layer | Files under `include/miam/` |
|-|-|
| A | `processes/constants/rate_expression.hpp` (new); the six config `.hpp` + `_builder.hpp` pairs — drop `std::function` members and `std::function`-taking builder overloads outright |
| B | `representations/aerosol_property.hpp` (rewrite — POD descriptor); `representations/aerosol_property_evaluator.hpp` (new); `representations/single_moment_mode.hpp`, `two_moment_mode.hpp`, `uniform_section.hpp` (rename `GetPropertyProvider` → `GetPropertyDescriptor`) |
| C | Six new headers: `processes/dissolved_reaction_set.hpp`, `processes/dissolved_reversible_reaction_set.hpp`, `processes/henrys_law_phase_transfer_set.hpp`, `constraints/dissolved_equilibrium_constraint_set.hpp`, `constraints/henrys_law_equilibrium_constraint_set.hpp`, `constraints/linear_constraint_set.hpp`. The corresponding six config headers keep their build-time query methods and public data; every `*Function<...>()` factory is deleted. |
| D | `model/model.hpp` — delete all `std::function`-returning factories and the seven `mutable std::any` caches; add `ProcessSetCollection<DenseMatrixPolicy, SparseMatrixPolicy>`, `ConstraintSetCollection<DenseMatrixPolicy, SparseMatrixPolicy>`, and the `ConstantsBucket` companions; add the lazy-init trigger inside `AddForcingTerms<DenseMatrixPolicy>` / `SubtractJacobianTerms<DenseMatrixPolicy, SparseMatrixPolicy>` / their constraint mates. |

## Test migration

Every test that currently calls a `Model::*Function<DenseMatrixPolicy>()` or
`Model::*Function<DenseMatrixPolicy, SparseMatrixPolicy>()` migrates to the
direct-dispatch API. MICM's `micm::FiniteDifferenceJacobian`,
`micm::CompareJacobianToFiniteDifference`, and
`micm::CheckJacobianSparsityCompleteness` (in
`micm/util/jacobian_verification.hpp`) are matrix-policy-templated, accept
any callable, and internally invoke `.CopyToHost()` on their perturbation
matrices, so they work uniformly on CPU and Kokkos backings.

The migration pattern is:

```cpp
// Before:
auto forcing_fn = model.ForcingFunction<DenseMatrix>(parameter_map, variable_map);
forcing_fn(parameters, variables, forcing);

// After:
model.FinalizeProcessSetup(parameter_map, variable_map, jacobian);
model.AddForcingTerms(parameters, variables, forcing);
```

And for the FD wrappers:

```cpp
auto forcing_wrapper = [&](const DenseMatrixPolicy& variables, DenseMatrixPolicy& forcing)
{ model.AddForcingTerms(parameters, variables, forcing); };
auto finite_difference_jacobian =
    micm::FiniteDifferenceJacobian<DenseMatrixPolicy>(forcing_wrapper, variables, number_of_species);
```

Files migrated:

- `test/integration/test_jacobian_verification.cpp` (39 hits of the old API).
- `test/integration/test_cam_cloud_chemistry.cpp` — Step 5 (`Step5_JacobianVerification`)
  moves out of the CPU driver and into `cam_cloud_chemistry_policy.hpp` as a
  templated policy function alongside the other eight steps. Both drivers
  (CPU and Kokkos) then exercise all nine steps.
- `test/unit/processes/dissolved_reaction.cpp`,
  `test/unit/processes/dissolved_reversible_reaction.cpp`,
  `test/unit/processes/henrys_law_phase_transfer.cpp` — each unit test does
  a `FinalizeProcessSetup` on the raw process (with a hand-built
  `SparseMatrixFD`) or on a small `Model` containing just that process, then
  calls `AddForcingTerms` / `SubtractJacobianTerms` /
  `UpdateStateParameters` directly.
- `test/unit/constraints/dissolved_equilibrium_constraint.cpp`,
  `test/unit/constraints/henrys_law_equilibrium_constraint.cpp`,
  `test/unit/constraints/linear_constraint.cpp` — same pattern for
  constraints.

## Staging plan

Five small, reviewable PRs. Each ends with all CPU tests green, and each
expands the set of `test_kokkos_cam_cloud_chemistry` steps that pass.

| PR | Content | Kokkos-serial state at end |
|-|-|-|
| 1 | Layer A. Add `RateExpression` classes + variant aliases. Update all six builders' constant-setters to require an expression variant (or a concrete alternative). `std::function` overloads deleted. Ad-hoc lambda constants in tests rewritten to `UserDefinedConstantExpression{ .value_ = k }` or `ArrheniusExpression{ .params_ = ... }`. | 0/9 — plumbing only |
| 2 | Layer B. POD `AerosolPropertyDescriptor` + evaluators. Three representations updated. `HenrysLawPhaseTransfer` consumer swap. | 0/9 |
| 3 | Layer C for `DissolvedReactionSet`. Build the Set class + `.CopyToDevice()`. `Model` routes `DissolvedReaction` through the new Set; leaves the other two process types on the existing paths behind a variant + `if constexpr`. Delete all OLD `DissolvedReaction::*Function<...>()` factories. Migrate `DissolvedReaction` unit tests. | 0/9 — no test exercises only `DissolvedReaction` |
| 4 | Layer C for `DissolvedReversibleReactionSet` and `HenrysLawPhaseTransferSet`. Delete their OLD `*Function<...>()` factories. Delete `Model`'s process-side `mutable std::any` glue. Migrate their unit tests. | Processes are Kokkos-safe; steps that use only constraints still fail. |
| 5 | Layer C for the three constraint Sets + Layer D (`Model` collapse: delete the seven `mutable std::any` caches, replace with `ProcessSetCollection` / `ConstraintSetCollection` + `ConstantsBucket`). Delete every remaining `*Function<...>()` factory on both `Model` and every config class. Migrate `test_jacobian_verification`, remaining unit tests, and Step 5 of the CAM cloud policy. | **All 9/9 `KokkosCamCloudChemistry.*` steps pass on Kokkos-serial.** |

## Decisions locked in

1. **Constants storage.** Each rate/equilibrium/HLC expression is its own POD
   class with a `MICM_DEVICE_FUNCTION Calculate(const micm::Conditions&) const`
   method. Concrete alternatives are grouped in a `std::variant<...>` alias
   per constant kind. Config-side variants are host-only and never touched
   at solve time.
2. **All `std::function` usage removed outright**, including the builder
   overloads. No deprecation window.
3. **FD Jacobian tests migrate to MICM's current API** and become
   matrix-policy-templated; they run against both CPU and Kokkos backends
   using the same policy header.
4. **Cache storage uses `SparseMatrixPolicy::template VectorType<U>`** as
   the `Vector<U>` alias — the same pattern `micm::ProcessSet` uses. Provides
   device-mirrored views and `.CopyTo{Device,Host}()` at no extra cost.
5. **`UpdateStateParameters` is fully device-resident** via a bucketed
   `ConstantsBucket<DenseMatrixPolicy, SparseMatrixPolicy, ExpressionTypes...>`,
   populated once at Finalize and iterated inside a single `MICM_LAMBDA`.
