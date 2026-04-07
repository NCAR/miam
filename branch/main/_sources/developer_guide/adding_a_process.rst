=================
Adding a Process
=================

This guide walks through the steps required to add a new process type to
MIAM.

Process Interface
=================

Every process type must be a struct or class that provides the following
methods. MIAM uses ``std::variant`` (not virtual inheritance), so the
interface is enforced by the ``Model::ForEachProcess`` visitor pattern.

Required Methods
----------------

1. **CopyWithNewUuid()**

   .. code-block:: c++

      MyProcess CopyWithNewUuid() const;

   Return a copy of the process with a fresh UUID. Called by
   ``Model::AddProcesses`` to ensure each registered process instance has
   a unique identifier.

2. **ProcessParameterNames(phase_prefixes)**

   .. code-block:: c++

      std::set<std::string> ProcessParameterNames(
          const std::map<std::string, std::set<std::string>>& phase_prefixes) const;

   Return unique names for all state parameters this process needs (e.g.,
   rate constants). The ``phase_prefixes`` map provides the set of
   mode/section prefixes that own each phase name.

3. **SpeciesUsed(phase_prefixes)**

   .. code-block:: c++

      std::set<std::string> SpeciesUsed(
          const std::map<std::string, std::set<std::string>>& phase_prefixes) const;

   Return the set of fully-qualified state variable names that this
   process reads or writes.

4. **RequiredAerosolProperties()**

   .. code-block:: c++

      std::map<std::string, std::vector<AerosolProperty>> RequiredAerosolProperties() const;

   Return a map from phase name to the list of aerosol properties this
   process needs. Return an empty map if none are required.

5. **NonZeroJacobianElements(phase_prefixes, state_indices)**

   .. code-block:: c++

      std::set<std::pair<std::size_t, std::size_t>> NonZeroJacobianElements(
          const std::map<std::string, std::set<std::string>>& phase_prefixes,
          const std::unordered_map<std::string, std::size_t>& state_indices) const;

   Return all (row, col) pairs where this process contributes non-zero
   Jacobian entries. Include both direct and indirect (through provider)
   entries.

6. **UpdateStateParametersFunction<DenseMatrixPolicy>(...)**

   .. code-block:: c++

      template<typename DenseMatrixPolicy>
      std::function<void(const std::vector<micm::Conditions>&, DenseMatrixPolicy&)>
      UpdateStateParametersFunction(
          const std::map<std::string, std::set<std::string>>& phase_prefixes,
          const std::unordered_map<std::string, std::size_t>& param_indices) const;

   Return a closure that evaluates temperature-dependent constants and
   writes them into the state parameter matrix.

7. **ForcingFunction<DenseMatrixPolicy>(...)**

   .. code-block:: c++

      template<typename DenseMatrixPolicy>
      std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)>
      ForcingFunction(
          const std::map<std::string, std::set<std::string>>& phase_prefixes,
          const std::unordered_map<std::string, std::size_t>& param_indices,
          const std::unordered_map<std::string, std::size_t>& var_indices,
          const ProviderMap& providers) const;

   Return a closure that adds this process's contribution to forcing
   terms. The closure signature is ``(params, vars, forcing) -> void``.

8. **JacobianFunction<DenseMatrixPolicy, SparseMatrixPolicy>(...)**

   .. code-block:: c++

      template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
      std::function<void(const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)>
      JacobianFunction(
          const std::map<std::string, std::set<std::string>>& phase_prefixes,
          const std::unordered_map<std::string, std::size_t>& param_indices,
          const std::unordered_map<std::string, std::size_t>& var_indices,
          const SparseMatrixPolicy& jacobian,
          const ProviderMap& providers) const;

   Return a closure that adds Jacobian contributions. **Remember**: MICM
   stores **−J**, so negate all entries relative to the mathematical
   derivative.

Registering the New Type
========================

1. Add the header to [include/miam/process.hpp](include/miam/process.hpp).

2. Add the type to the ``ProcessVariant`` in
   [include/miam/model.hpp](include/miam/model.hpp):

   .. code-block:: c++

      using ProcessVariant = std::variant<
          process::DissolvedReversibleReaction,
          process::HenryLawPhaseTransfer,
          process::MyNewProcess        // ← add here
      >;

That's it — the ``std::visit`` in ``ForEachProcess`` will automatically
dispatch to the new type.

Sign Convention Checklist
=========================

- Forcing: compute the mathematical rate and write it directly (positive
  for production, negative for consumption).
- Jacobian: compute the mathematical ∂f/∂x, then **negate** before
  writing into the sparse Jacobian matrix.
- Verify with a unit test that compares analytical Jacobian entries
  against finite-difference approximations.
