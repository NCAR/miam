=========================
Adding a Representation
=========================

This guide covers the interface a new representation type must satisfy
and how to register it with the Model.

Representation Interface
========================

Like processes, representations use ``std::variant`` dispatch. A new
representation must provide the following methods:

Required Methods
----------------

1. **StateSize()**

   .. code-block:: c++

      std::tuple<std::size_t, std::size_t> StateSize() const;

   Return ``(num_variables, num_parameters)``—the count of state variables
   and state parameters introduced by this representation.

2. **StateVariableNames()**

   .. code-block:: c++

      std::set<std::string> StateVariableNames() const;

   Return fully-qualified names for all state variables
   (e.g., ``"MODE.PHASE.SPECIES"``).

3. **StateParameterNames()**

   .. code-block:: c++

      std::set<std::string> StateParameterNames() const;

   Return names for all state parameters (e.g., ``"MODE.geometric_mean_radius"``).

4. **PhaseStatePrefixes()**

   .. code-block:: c++

      std::map<std::string, std::set<std::string>> PhaseStatePrefixes() const;

   Return a map from phase name to the set of mode/section prefixes.
   Processes use this to expand reactions across all instances.

5. **GetPropertyProvider<DenseMatrixPolicy>(...)**

   .. code-block:: c++

      template<typename DenseMatrixPolicy>
      AerosolPropertyProvider<DenseMatrixPolicy> GetPropertyProvider(
          AerosolProperty property,
          const std::unordered_map<std::string, std::size_t>& param_indices,
          const std::unordered_map<std::string, std::size_t>& var_indices,
          const std::string& phase_name) const;

   Create a provider for the requested property. The provider must:

   - Populate ``dependent_variable_indices`` with the state variable
     indices that the property depends on.
   - Implement ``ComputeValue`` for the forcing function path.
   - Implement ``ComputeValueAndDerivatives`` for the Jacobian path,
     writing partial derivatives into the partials matrix.

6. **SetDefaultParameters(state)** (optional but recommended)

   Write default parameter values (GMD, GSD, section bounds) into a
   solver state. Convenience for users setting up initial conditions.

Provider Contract
=================

The ``AerosolPropertyProvider`` struct is the bridge between
representations and processes:

- ``dependent_variable_indices``: fixed at creation time, used by the
  process to determine Jacobian sparsity and map partials to state
  variable indices.
- ``ComputeValue``: called on the forcing path (no derivatives needed).
- ``ComputeValueAndDerivatives``: called on the Jacobian path. Column
  :math:`k` of the partials matrix corresponds to
  :math:`\partial P / \partial y_k` where :math:`y_k` is the state
  variable at ``dependent_variable_indices[k]``.

Registering the New Type
========================

1. Add the header to [include/miam/representation.hpp](include/miam/representation.hpp).

2. Add the type to the ``RepresentationVariant`` in
   [include/miam/model.hpp](include/miam/model.hpp):

   .. code-block:: c++

      using RepresentationVariant = std::variant<
          representation::SingleMomentMode,
          representation::TwoMomentMode,
          representation::UniformSection,
          representation::MyNewRepresentation   // ← add here
      >;

Testing Providers
=================

Verify each provider with:

1. **Value test**: set known state values, call ``ComputeValue``, compare
   against an analytical expectation.
2. **Derivative test**: call ``ComputeValueAndDerivatives``, compare
   partial derivatives against finite-difference approximations:

   .. code-block:: c++

      double dP_dy = (P(y + h) - P(y - h)) / (2 * h);

3. **Sparsity test**: verify that ``dependent_variable_indices`` is
   consistent with which partials are non-zero.
