# Henry's Law Phase Transfer — Design & Implementation Plan

## Overview

Introduce a Henry's Law phase transfer process to MIAM, enabling gas-particle
partitioning between a single gas-phase species and one or more condensed-phase
instances. Follows the pattern established by `DissolvedReversibleReaction`
(PR #15) and the CAMP `rxn_HL_phase_transfer` treatment.

## Physics

For gas species A partitioning into condensed phase(s):

```
A(gas) ⇌ A(condensed)
```

### Concentration units

All state variables are stored in `mol m⁻³ air` (moles per cubic meter of air).
Henry's Law, however, relates gas-phase concentration to **dissolved**
concentration (`mol m⁻³ solvent` — moles per cubic meter of solvent). To
convert, we need the **condensed-phase solvent species** (typically water) and
its **volume fraction**:

```
f_v = [solvent] · solvent_molecular_weight / solvent_density
```

where `[solvent]` is the solvent concentration (`mol m⁻³ air`), `solvent_molecular_weight` is
the solvent molar mass (`kg mol⁻¹`), and `solvent_density` is the solvent density
(`kg m⁻³`). The dissolved concentration of species A is then:

```
[A]_dissolved = [A]_aq / f_v     (mol m⁻³ solvent)
```

### Mass transfer rates

**Condensation rate** (per condensed-phase instance):

```
k_cond = 4π · r_eff · N · D_g · f(Kn)
```

where `r_eff` = effective radius, `N` = number concentration, `D_g` = gas-phase
diffusion coefficient, and `f(Kn)` is the regime correction factor.

**Continuum regime**: `f(Kn) = 1`

**Transition regime** (Fuchs-Sutugin):

```
f(Kn) = (1 + Kn) / (1 + 2·Kn·(1 + Kn) / α)
```

where `Kn = λ / r_eff`, `λ = 3·D_g / c̄`, `c̄ = √(8RT / (π·Mw))`,
`α` = mass accommodation coefficient.

### Phase volume fraction scaling

All three representations assume spherical particles, and the exposed surface
area of a phase within a particle is proportional to its volume fraction in the
mode/section. When a mode hosts multiple phases, the effective condensation rate
for a particular phase is scaled by that phase’s volume fraction:

```
φ_p = V_phase / V_total
```

where `V_phase = Σ_i ( [species_{p,i}] · Mw_{p,i} / ρ_{p,i} )` is the total
volume of species in phase `p`, and `V_total` is the total volume across
all phases in the mode/section (see § Species property requirements below).

For example, if a mode has aqueous volume 2 and organic volume 5, then
φ_aq = 2/7 and a gas species partitioning into the aqueous phase sees only
2/7 of the total particle surface area.

The **effective condensation rate** for phase `p` is:

```
k_cond_eff = φ_p · k_cond
k_evap_eff = k_cond_eff / (HLC · R · T) = φ_p · k_evap
```

**Evaporation rate constant** (from Henry's Law):

```
k_evap = k_cond / (HLC · R · T)
```

where `HLC = H(T)` is the temperature-dependent Henry's Law constant
(`mol m⁻³ Pa⁻¹`), recalculated whenever conditions change and stored as a
state parameter (see § Process setup phase below).

### Net rate of change

Using the dissolved concentration for the equilibrium comparison,
and scaling by the phase volume fraction φ_p:

```
d[A]_gas/dt  = -φ_p · k_cond · [A]_gas + φ_p · k_evap · [A]_aq / f_v
d[A]_aq/dt   = +φ_p · k_cond · [A]_gas - φ_p · k_evap · [A]_aq / f_v
d[solvent]/dt = 0   (solvent not consumed by this process)
```

Note: the forcing terms for gas and aq are in `mol m⁻³ air` per unit time.
The `1/f_v` factor converts from the air-volume basis back to the dissolved
basis as needed by Henry's Law equilibrium. For a mode with a single phase,
φ_p = 1 and the expressions reduce to the un-scaled form.

## Mathematical Reference

This section collects every equation used in the Henry's Law phase transfer
rate calculation, its partial derivatives with respect to state variables,
and the Jacobian sign convention used by the MICM solver.

### Constants and notation

| Symbol | Description | Units |
|--------|-------------|-------|
| R      | Universal gas constant (8.314462618) | J mol⁻¹ K⁻¹ |
| T      | Temperature | K |
| D_g    | Gas-phase diffusion coefficient | m² s⁻¹ |
| α      | Mass accommodation coefficient | dimensionless (0–1) |
| gas_molecular_weight | Molecular weight of gas species | kg mol⁻¹ |
| solvent_molecular_weight | Molecular weight of solvent | kg mol⁻¹ |
| solvent_density | Density of solvent | kg m⁻³ |
| [A]_gas | Gas-phase concentration | mol m⁻³ air |
| [A]_aq | Aqueous-phase concentration | mol m⁻³ air |
| [solvent] | Solvent concentration | mol m⁻³ air |

### 1. Henry's Law Constant (temperature-dependent)

```
HLC(T) = HLC_ref · exp( C · (1/T − 1/T₀) )
```

| Variable | Description | Units |
|----------|-------------|-------|
| HLC_ref  | Reference Henry's Law constant at T₀ | mol m⁻³ Pa⁻¹ |
| C        | Temperature dependence parameter | K |
| T₀       | Reference temperature (default 298.15) | K |
| HLC(T)   | Henry's Law constant at temperature T | mol m⁻³ Pa⁻¹ |

HLC is evaluated once per time step in `UpdateStateParametersFunction` and
stored as a state parameter — it is not differentiated with respect to state
variables.

### 2. Mean molecular speed

```
c̄ = √(8 R T / (π gas_molecular_weight))
```

Units: m s⁻¹. Used to derive the mean free path.

### 3. Mean free path

```
λ = 3 D_g / c̄
```

Units: m.

### 4. Knudsen number

```
Kn = λ / r_eff
```

Dimensionless. Characterizes the gas-particle interaction regime
(Kn ≪ 1 → continuum, Kn ≫ 1 → free-molecular).

### 5. Fuchs-Sutugin transition-regime correction

```
f(Kn) = (1 + Kn) / (1 + 2 Kn (1 + Kn) / α)
```

Dimensionless. Interpolates between the continuum (f → 1) and
free-molecular (f → α / (2 Kn)) limits.

**Partial derivative with respect to Kn** (needed for the chain rule
through r_eff):

```
Let denom = 1 + 2 Kn (1 + Kn) / α

df/dKn = [α − 2 − 4 Kn − 2 Kn²] / [α · denom²]
```

### 6. Condensation rate

```
k_c = 4π · r_eff · N · D_g · f(Kn)
```

Units: s⁻¹. The first-order rate constant for gas-to-condensed transfer.

**Partial derivative with respect to r_eff** (chain rule through Kn):

```
dKn/dr_eff = −Kn / r_eff

dk_c/dr_eff = 4π N D_g · [ f(Kn) + r_eff · df/dKn · dKn/dr_eff ]
            = 4π N D_g · [ f(Kn) − Kn · df/dKn ]
```

Units: s⁻¹ m⁻¹.

**Partial derivative with respect to N** (linear dependence):

```
dk_c/dN = k_c / N
```

Units: s⁻¹ # ⁻¹ m³. Because k_c is linear in N, the derivative is
simply the per-particle condensation rate.

### 7. Evaporation rate

```
k_e = k_c / (HLC · R · T)
```

Units: s⁻¹. Derived from Henry's Law equilibrium — at equilibrium the
condensation and evaporation fluxes balance.

**Partial derivatives** — obtained by dividing the k_c partial by the
same constant factor:

```
dk_e/dr_eff = dk_c/dr_eff / (HLC · R · T)
dk_e/dN     = dk_c/dN     / (HLC · R · T) = k_e / N
```

### 8. Solvent volume fraction

```
f_v = [solvent] · solvent_molecular_weight / solvent_density
```

Dimensionless. Converts the dissolved-phase concentration basis
(`mol m⁻³ solvent`) back to the air-volume basis (`mol m⁻³ air`).

### 9. Net transfer rate (forcing function)

For a single condensed-phase instance with phase volume fraction φ_p:

```
R_net = φ_p · k_c · [A]_gas − φ_p · k_e · [A]_aq / f_v
```

Units: mol m⁻³ s⁻¹. The ODE right-hand sides are:

```
d[A]_gas / dt = −R_net      (gas is consumed)
d[A]_aq  / dt = +R_net      (condensed species is produced)
d[solvent] / dt = 0          (solvent unchanged by this process)
```

When multiple condensed-phase instances exist (e.g., different aerosol
modes containing the same phase), each contributes its own R_net with its
own r_eff, N, and φ_p; the contributions are summed into the same gas
species forcing term.

### 10. Jacobian entries

The MICM Rosenbrock solver stores **−J** (negative Jacobian). All entries
produced by the Jacobian function are therefore negated relative to the
mathematical derivatives. The table below shows the mathematical derivative
and the value actually stored.

#### 10a. Direct entries (w.r.t. state variables [A]_gas, [A]_aq, [solvent])

| Entry | Mathematical J | Stored −J |
|-------|---------------|-----------|
| J[gas, gas]     | −φ_p · k_c               | +φ_p · k_c               |
| J[gas, aq]      | +φ_p · k_e / f_v         | −φ_p · k_e / f_v         |
| J[gas, solvent] | −φ_p · k_e · [A]_aq / (f_v · [solvent]) | +φ_p · k_e · [A]_aq / (f_v · [solvent]) |
| J[aq, gas]      | +φ_p · k_c               | −φ_p · k_c               |
| J[aq, aq]       | −φ_p · k_e / f_v         | +φ_p · k_e / f_v         |
| J[aq, solvent]  | +φ_p · k_e · [A]_aq / (f_v · [solvent]) | −φ_p · k_e · [A]_aq / (f_v · [solvent]) |

**Derivations (mathematical J, before negation):**

```
J[gas,gas] = ∂(−R_net)/∂[A]_gas = −φ_p · k_c

J[gas,aq]  = ∂(−R_net)/∂[A]_aq  = +φ_p · k_e / f_v

J[gas,solvent] = ∂(−R_net)/∂[solvent]
    R_net contains the term  −φ_p · k_e · [A]_aq / f_v
    where f_v = [solvent] · Mw_s/ρ_s, so:
    ∂/∂[solvent](−(−φ_p k_e [A]_aq / f_v))
      = −φ_p · k_e · [A]_aq / (f_v · [solvent])

J[aq,x] = −J[gas,x]   for all x  (mass conservation)
```

Note that J[gas,x] = −J[aq,x] for every column x — this antisymmetry is
a direct consequence of mass conservation (the gas loss equals the
condensed-phase gain) and holds regardless of the sign convention.

#### 10b. Indirect entries through aerosol properties

When an aerosol property (r_eff, N, or φ_p) depends on a state variable
y_j, the Jacobian gains additional entries via the chain rule.

**Through r_eff** (applies when r_eff depends on species concentrations,
e.g., TwoMomentMode):

```
R_net = φ_p · (k_c · [A]_gas − k_e · [A]_aq / f_v)

∂R_net/∂y_j via r_eff
    = φ_p · (dk_c/dr_eff · [A]_gas − dk_e/dr_eff · [A]_aq / f_v) · ∂r_eff/∂y_j

Stored −J[gas, y_j] += +φ_p · (dk_c/dr · [A]_gas − dk_e/dr · [A]_aq / f_v) · ∂r/∂y_j
Stored −J[aq,  y_j] += −φ_p · (dk_c/dr · [A]_gas − dk_e/dr · [A]_aq / f_v) · ∂r/∂y_j
```

**Through N** (applies when N depends on species concentrations, e.g.,
SingleMomentMode and UniformSection, or on the number concentration
state variable in TwoMomentMode):

```
∂R_net/∂y_j via N
    = φ_p · (dk_c/dN · [A]_gas − dk_e/dN · [A]_aq / f_v) · ∂N/∂y_j

Stored −J[gas, y_j] += +φ_p · (dk_c/dN · [A]_gas − dk_e/dN · [A]_aq / f_v) · ∂N/∂y_j
Stored −J[aq,  y_j] += −φ_p · (dk_c/dN · [A]_gas − dk_e/dN · [A]_aq / f_v) · ∂N/∂y_j
```

**Through φ_p** (applies when the mode contains multiple phases and φ_p
depends on species concentrations):

```
Let R = k_c · [A]_gas − k_e · [A]_aq / f_v   (un-scaled net rate)
Then R_net = φ_p · R

∂R_net/∂y_j via φ_p = R · ∂φ_p/∂y_j

Stored −J[gas, y_j] += +R · ∂φ_p/∂y_j
Stored −J[aq,  y_j] += −R · ∂φ_p/∂y_j
```

The unscaled rate R is computed once and reused for all φ_p-dependent
variables. The ∂φ_p/∂y_j values come from the phase volume fraction
provider (see § Aerosol property partial derivatives below).

#### 10c. Jacobian sparsity

The set of nonzero Jacobian elements is determined at setup time from the
union of:

- **Direct**: (gas, gas), (gas, aq), (gas, solvent), (aq, gas), (aq, aq),
  (aq, solvent) — 6 entries per instance.
- **Via r_eff**: (gas, y_j), (aq, y_j) for each y_j in
  `r_eff_provider.dependent_variable_indices` — 2 entries per dependency.
- **Via N**: (gas, y_j), (aq, y_j) for each y_j in
  `N_provider.dependent_variable_indices` — 2 entries per dependency.
- **Via φ_p**: (gas, y_j), (aq, y_j) for each y_j in
  `phi_provider.dependent_variable_indices` — 2 entries per dependency.

The total number of nonzero entries depends on the representation type.
For SingleMomentMode with a single phase the provider dependencies are:

| Property | Dependencies | Entries |
|----------|-------------|---------|
| r_eff | none (parameterized) | 0 |
| N | all species in mode | 2 × n_species |
| φ_p | none (single phase → φ = 1) | 0 |

For TwoMomentMode with two phases, each containing m_p and m_q species
respectively, plus a number concentration variable:

| Property | Dependencies | Entries |
|----------|-------------|---------|
| r_eff | all species + N_var | 2 × (m_p + m_q + 1) |
| N | N_var | 2 × 1 |
| φ_p | all species | 2 × (m_p + m_q) |

### 11. Aerosol property partial derivatives

These derivatives are computed by representation-specific providers and
consumed by the process's Jacobian function via the chain rule
(§ 10b above).

#### Effective radius (r_eff)

**SingleMomentMode** — parameterized, no state variable dependencies:

```
r_eff = GMD · exp(2.5 · ln²(GSD))

∂r_eff/∂y_j = 0   for all state variables y_j
```

**TwoMomentMode** — depends on total volume and number concentration:

```
V_total = Σ_p Σ_i [species_{p,i}] · Mw_{p,i} / ρ_{p,i}
r_mean  = (3 V_total / (4π N))^(1/3)
r_eff   = r_mean · exp(2.5 · ln²(GSD))

∂r_eff/∂[species_{p,i}] = r_eff · (Mw_{p,i} / ρ_{p,i}) / (3 V_total)
∂r_eff/∂N               = −r_eff / (3 N)
```

**UniformSection** — parameterized, no state variable dependencies:

```
r_eff = (r_min + r_max) / 2

∂r_eff/∂y_j = 0   for all state variables y_j
```

#### Number concentration (N)

**SingleMomentMode** — derived from total volume and fixed single-particle
volume:

```
V_single = (4/3)π · GMD³ · exp(4.5 · ln²(GSD))
V_total  = Σ_p Σ_i [species_{p,i}] · Mw_{p,i} / ρ_{p,i}
N        = V_total / V_single

∂N/∂[species_{p,i}] = (Mw_{p,i} / ρ_{p,i}) / V_single
```

**TwoMomentMode** — prognostic state variable:

```
N = [N_var]      (a directly tracked state variable)

∂N/∂[N_var] = 1
```

**UniformSection** — derived from total volume and fixed single-particle
volume:

```
V_single = (4/3)π · r_eff³
V_total  = Σ_p Σ_i [species_{p,i}] · Mw_{p,i} / ρ_{p,i}
N        = V_total / V_single

∂N/∂[species_{p,i}] = (Mw_{p,i} / ρ_{p,i}) / V_single
```

#### Phase volume fraction (φ_p)

All representation types use the same formula. When a mode or section
contains only a single phase, φ_p = 1 and all partials are zero.

For multi-phase modes/sections:

```
V_phase = Σ_i [species_{p,i}] · Mw_{p,i} / ρ_{p,i}    (target phase only)
V_total = Σ_q Σ_i [species_{q,i}] · Mw_{q,i} / ρ_{q,i}  (all phases)
φ_p     = V_phase / V_total

Same-phase species:
  ∂φ_p/∂[species_{p,i}] = (Mw_{p,i} / ρ_{p,i}) · (1 − φ_p) / V_total

Other-phase species:
  ∂φ_p/∂[species_{q,i}] = −(Mw_{q,i} / ρ_{q,i}) · φ_p / V_total
```

### 12. Summary of the complete derivative chain

For a single condensed-phase instance i, the full Jacobian contribution to
the gas species equation is (mathematical J, before MICM negation):

```
∂(d[A]_gas/dt) / ∂y_j =

  DIRECT (y_j ∈ {[A]_gas, [A]_aq, [solvent]}):
    see § 10a table

  + VIA r_eff:
    −φ_p · (dk_c/dr · [A]_gas − dk_e/dr · [A]_aq / f_v) · ∂r/∂y_j

  + VIA N:
    −φ_p · (dk_c/dN · [A]_gas − dk_e/dN · [A]_aq / f_v) · ∂N/∂y_j

  + VIA φ_p:
    −R · ∂φ_p/∂y_j     where R = k_c · [A]_gas − k_e · [A]_aq / f_v
```

The aqueous species row is always the negative of the gas row
(J[aq, y_j] = −J[gas, y_j]), so only one set of chain-rule products needs
to be computed. All stored values are then negated for the MICM −J
convention.

## Key Differences from DissolvedReversibleReaction

| Aspect | Dissolved Reversible Reaction | Henry's Law Phase Transfer |
|--------|-------------------------------|----------------------------|
| Phases involved | Single condensed phase | Gas phase + 1+ condensed phases |
| Gas-phase species | None | Single gas-phase species |
| Solvent species | Implicit (not tracked) | Explicit condensed-phase solvent per instance (needed for `f_v`) |
| Rate dependence | Rate constants from `Conditions` only | Rate depends on particle properties (r_eff, N) from representation state |
| Phase instances | Same reaction replicated per instance | Each instance has distinct r_eff/N/solvent, affecting rate per-instance |
| Condensed species | Multiple reactants/products | Single condensed solute species + solvent species per instance |

**Note on modes/sections and phases:** Each representation instance (mode or
section) can host one or more phases — e.g., an accumulation mode might contain
both an aqueous phase and an organic phase. All species across all phases in a
mode/section contribute to the total particle volume used to derive size
properties like effective radius and number concentration. Because the particles
are assumed spherical, each phase’s exposed surface area is proportional to its
volume fraction φ_p = V_phase / V_total within the mode/section (see § Phase
volume fraction scaling).

## Common Process Interface

All MIAM process types implement a common interface so that the Model can
store, iterate, and dispatch to them generically. This makes it possible to
add new process types without modifying Model.

### Required methods

Every MIAM process must implement these methods:

| Method | Signature | Purpose |
|--------|-----------|--------|
| `ProcessParameterNames` | `(phase_prefixes) → set<string>` | Names of parameters this process contributes to the state |
| `SpeciesUsed` | `(phase_prefixes) → set<string>` | State variable names this process reads/writes |
| `RequiredAerosolProperties` | `() → map<string, vector<AerosolProperty>>` | Properties needed from representations, keyed by phase name |
| `NonZeroJacobianElements` | `(phase_prefixes, var_indices, providers) → set<pair>` | Jacobian sparsity (may depend on provider dependencies) |
| `UpdateStateParametersFunction` | `(phase_prefixes, param_indices) → function` | Updates rate parameters from `Conditions` |
| `ForcingFunction` | `(phase_prefixes, param_indices, var_indices, providers) → function` | Builds the forcing lambda |
| `JacobianFunction` | `(phase_prefixes, param_indices, var_indices, jacobian, providers) → function` | Builds the Jacobian lambda |

`providers` is a `map<string, vector<AerosolPropertyProvider>>` — one
vector of providers per phase prefix. Processes that don't need aerosol
properties (e.g., `DissolvedReversibleReaction`) return an empty map from
`RequiredAerosolProperties()` and ignore the `providers` argument.

### Type erasure

Model stores processes in a single collection using a type-erased wrapper.
This follows the same pattern MICM uses for `ExternalModelProcessSet`:

```cpp
// Type-erased wrapper around any MIAM process
template<typename DenseMatrixPolicy, typename SparseMatrixPolicy>
struct MiamProcessSet
{
  using PhaseMap = std::map<std::string, std::set<std::string>>;
  using IndexMap = std::unordered_map<std::string, std::size_t>;
  using ProviderMap = std::map<std::string,
      std::vector<AerosolPropertyProvider<DenseMatrixPolicy>>>;

  std::function<std::set<std::string>(const PhaseMap&)>
      process_parameter_names_;
  std::function<std::set<std::string>(const PhaseMap&)>
      species_used_;
  std::function<std::map<std::string, std::vector<AerosolProperty>>()>
      required_aerosol_properties_;
  std::function<std::set<std::pair<std::size_t, std::size_t>>(
      const PhaseMap&, const IndexMap&, const ProviderMap&)>
      non_zero_jacobian_elements_;
  std::function<std::function<void(
      const std::vector<micm::Conditions>&, DenseMatrixPolicy&)>(
      const PhaseMap&, const IndexMap&)>
      update_state_parameters_function_;
  std::function<std::function<void(
      const DenseMatrixPolicy&, const DenseMatrixPolicy&, DenseMatrixPolicy&)>(
      const PhaseMap&, const IndexMap&, const IndexMap&, ProviderMap)>
      get_forcing_function_;
  std::function<std::function<void(
      const DenseMatrixPolicy&, const DenseMatrixPolicy&, SparseMatrixPolicy&)>(
      const PhaseMap&, const IndexMap&, const IndexMap&,
      const SparseMatrixPolicy&, ProviderMap)>
      get_jacobian_function_;

  template<typename ProcessType>
  MiamProcessSet(ProcessType&& process)
  {
    auto shared = std::make_shared<std::decay_t<ProcessType>>(
        std::forward<ProcessType>(process));

    process_parameter_names_ = [shared](const PhaseMap& pp)
    { return shared->ProcessParameterNames(pp); };

    species_used_ = [shared](const PhaseMap& pp)
    { return shared->SpeciesUsed(pp); };

    required_aerosol_properties_ = [shared]()
    { return shared->RequiredAerosolProperties(); };

    non_zero_jacobian_elements_ =
        [shared](const PhaseMap& pp, const IndexMap& vi,
                 const ProviderMap& prov)
    { return shared->NonZeroJacobianElements(pp, vi, prov); };

    update_state_parameters_function_ =
        [shared](const PhaseMap& pp, const IndexMap& pi)
    {
      return shared->template
          UpdateStateParametersFunction<DenseMatrixPolicy>(pp, pi);
    };

    get_forcing_function_ =
        [shared](const PhaseMap& pp, const IndexMap& pi,
                 const IndexMap& vi, ProviderMap prov)
    {
      return shared->template ForcingFunction<DenseMatrixPolicy>(
          pp, pi, vi, std::move(prov));
    };

    get_jacobian_function_ =
        [shared](const PhaseMap& pp, const IndexMap& pi,
                 const IndexMap& vi, const SparseMatrixPolicy& jac,
                 ProviderMap prov)
    {
      return shared->template
          JacobianFunction<DenseMatrixPolicy, SparseMatrixPolicy>(
              pp, pi, vi, jac, std::move(prov));
    };
  }
};
```

### Adapting existing processes

`DissolvedReversibleReaction` currently doesn't take a `providers` argument.
To satisfy the common interface:

```cpp
// Returns empty — no aerosol properties needed
std::map<std::string, std::vector<AerosolProperty>>
RequiredAerosolProperties() const { return {}; }

// Existing methods gain a providers parameter they ignore:
template<typename DenseMatrixPolicy>
auto ForcingFunction(
    const PhaseMap& phase_prefixes,
    const IndexMap& param_indices,
    const IndexMap& var_indices,
    const ProviderMap& /* providers */) const
{
  // Delegate to existing implementation (unchanged)
  return ForcingFunctionImpl<DenseMatrixPolicy>(
      phase_prefixes, param_indices, var_indices);
}
```

## Aerosol Property Provider Design

### Requirements

Processes like Henry's Law phase transfer need aerosol/cloud particle properties
(effective radius, number concentration) that come from the representations.
These properties:

1. May depend on state parameters (e.g., GMD, GSD) or state variables (e.g.,
   number concentration), depending on the representation type.
2. Need to provide partial derivatives w.r.t. state variables for Jacobian
   computation.
3. Must be callable within the solver hot loop with **zero memory allocations**.
4. Should be **one provider per property** for scalability—new properties can be
   added without changing existing providers.
5. Must work across all representation types (SingleMomentMode, TwoMomentMode,
   UniformSection) with representation-specific implementations.

### Pseudocode: Aerosol Property API from a Process's Perspective

The following shows how a generic process would use the aerosol property
provider API. This is written from the process's perspective—only the
interface it needs is shown; internal representation details are hidden.

The pseudocode follows the patterns established in `DissolvedReversibleReaction`:
- `DenseMatrixPolicy::Function()` wraps the outer lambda, receives GroupView
  objects (opaque views over vectorized row groups).
- `ForEachRow()` iterates across grid cells within a group, receiving element
  references from ColumnViews and RowVariables.
- `ForEachBlock()` does the same for sparse Jacobian blocks.
- `GetRowVariable()` creates a stack-allocated per-group temporary.
- `GetColumnView(idx)` / `GetConstColumnView(idx)` reference state data columns.
- No `matrix[i][j]` indexing — all access through the vectorized API.

```cpp
// ====================================================================
//  TYPES (defined once in the aerosol property header)
// ====================================================================

enum class AerosolProperty {
    EffectiveRadius,       // [m]
    NumberConcentration,   // [# m^-3]
    PhaseVolumeFraction    // [dimensionless, 0–1]
    // future: Density, SurfaceArea, ...
};

// A provider for a single aerosol property, created at setup time by a
// representation instance. Captures all needed parameter/variable
// column indices internally.
//
// Operates on ForEachRow-compatible column views — no per-cell indexing.
// Partial derivatives are written into columns of a pre-allocated
// DenseMatrixPolicy, one column per dependent variable.
template<typename DenseMatrixPolicy>
struct AerosolPropertyProvider {
    // State variable indices that this property has non-zero partial
    // derivatives with respect to. Fixed at creation time. Used by
    // processes to:
    //   1. Determine Jacobian sparsity (NonZeroJacobianElements)
    //   2. Know how many columns the partials matrix needs
    //   3. Map partials columns back to state variable indices
    std::vector<std::size_t> dependent_variable_indices;

    // Compute the property value for all grid cells in the current group.
    // Called inside a ForEachRow loop — receives column views and writes
    // into a RowVariable.
    //
    //   result:      mutable RowVariable to write the property value into
    //   params_view: const GroupView of state parameters
    //   vars_view:   const GroupView of state variables
    //
    // The provider internally calls params_view.GetConstColumnView(...)
    // and vars_view.GetConstColumnView(...) using its captured indices.
    std::function<void(
        auto&& params_view,
        auto&& vars_view,
        auto&& result          // RowVariable — written to as output
    )> ComputeValue;

    // Compute the property value AND partial derivatives for all grid
    // cells in the current group. Called inside a ForEachRow loop.
    //
    //   result:           mutable RowVariable for the property value
    //   partials_matrix:  mutable GroupView of the partials DenseMatrix
    //                     (num_cells × num_dependent_variables)
    //                     Column k corresponds to d(property)/d(var[dependent_variable_indices[k]])
    //   params_view:      const GroupView of state parameters
    //   vars_view:        const GroupView of state variables
    std::function<void(
        auto&& params_view,
        auto&& vars_view,
        auto&& result,             // RowVariable — value output
        auto&& partials_matrix     // mutable GroupView — partials output
    )> ComputeValueAndDerivatives;
};


// ====================================================================
//  REPRESENTATION INTERFACE (each representation type implements this)
// ====================================================================

// Each representation provides:
template<typename DenseMatrixPolicy>
AerosolPropertyProvider<DenseMatrixPolicy> GetPropertyProvider(
    AerosolProperty property,
    const auto& state_parameter_indices,   // unordered_map<string, size_t>
    const auto& state_variable_indices     // unordered_map<string, size_t>
) const;
// Throws std::runtime_error if the representation can't provide the
// requested property.


// ====================================================================
//  EXAMPLE PROVIDER IMPLEMENTATIONS (inside representations)
// ====================================================================

// --- SingleMomentMode: EffectiveRadius from GMD, GSD parameters ---
// r_eff = GMD * exp(2.5 * ln²(GSD))
// No state variable dependencies → dependent_variable_indices is empty,
// partials matrix has zero columns.

provider.dependent_variable_indices = {};

provider.ComputeValue = [gmd_param_idx, gsd_param_idx]
    (auto&& params_view, auto&& vars_view, auto&& result)
{
    params_view.ForEachRow(
        [](const double& gmd, const double& gsd, double& r_eff) {
            double ln_gsd = std::log(gsd);
            r_eff = gmd * std::exp(2.5 * ln_gsd * ln_gsd);
        },
        params_view.GetConstColumnView(gmd_param_idx),
        params_view.GetConstColumnView(gsd_param_idx),
        result);
};

provider.ComputeValueAndDerivatives = [gmd_param_idx, gsd_param_idx]
    (auto&& params_view, auto&& vars_view, auto&& result, auto&& partials)
{
    // Same computation; partials matrix has 0 columns, nothing to write
    params_view.ForEachRow(
        [](const double& gmd, const double& gsd, double& r_eff) {
            double ln_gsd = std::log(gsd);
            r_eff = gmd * std::exp(2.5 * ln_gsd * ln_gsd);
        },
        params_view.GetConstColumnView(gmd_param_idx),
        params_view.GetConstColumnView(gsd_param_idx),
        result);
};


// --- SingleMomentMode: NumberConcentration from total volume ---
// V_total = Σ_p Σ_i [species_{p,i}] · Mw_{p,i} / ρ_{p,i}
//   (sum over all phases p in the mode, all species i in each phase)
// V_single = (4/3)π · GMD³ · exp(4.5 · ln²(GSD))
// N = V_total / V_single
// ∂N/∂[species_{p,i}] = (Mw_{p,i} / ρ_{p,i}) / V_single
//
// dependent_variable_indices = { all species variable indices across all phases }
// Mw_over_rho[k] = Mw_k / ρ_k (pre-computed at setup)

provider.dependent_variable_indices = species_var_indices;  // one per species

provider.ComputeValue = [gmd_param_idx, gsd_param_idx,
                         species_var_indices, mw_over_rho]
    (auto&& params_view, auto&& vars_view, auto&& result)
{
    // Zero result before accumulating
    params_view.ForEachRow(
        [](double& N) { N = 0.0; },
        result);

    // Accumulate V_total = Σ [species_i] * Mw_i/ρ_i
    for (size_t k = 0; k < species_var_indices.size(); ++k) {
        params_view.ForEachRow(
            [mwr = mw_over_rho[k]](const double& conc, double& V) {
                V += conc * mwr;
            },
            vars_view.GetConstColumnView(species_var_indices[k]),
            result);  // result temporarily holds V_total
    }

    // N = V_total / V_single
    params_view.ForEachRow(
        [](const double& gmd, const double& gsd, double& N) {
            double ln_gsd = std::log(gsd);
            double V_single = (4.0/3.0) * M_PI * gmd * gmd * gmd
                             * std::exp(4.5 * ln_gsd * ln_gsd);
            N = N / V_single;  // N held V_total from accumulation
        },
        params_view.GetConstColumnView(gmd_param_idx),
        params_view.GetConstColumnView(gsd_param_idx),
        result);
};

provider.ComputeValueAndDerivatives = [gmd_param_idx, gsd_param_idx,
                                       species_var_indices, mw_over_rho]
    (auto&& params_view, auto&& vars_view, auto&& result, auto&& partials)
{
    // Compute V_total (same as ComputeValue)
    params_view.ForEachRow(
        [](double& N) { N = 0.0; },
        result);
    for (size_t k = 0; k < species_var_indices.size(); ++k) {
        params_view.ForEachRow(
            [mwr = mw_over_rho[k]](const double& conc, double& V) {
                V += conc * mwr;
            },
            vars_view.GetConstColumnView(species_var_indices[k]),
            result);
    }

    // N = V_total / V_single, ∂N/∂[species_k] = (Mw_k/ρ_k) / V_single
    params_view.ForEachRow(
        [](const double& gmd, const double& gsd, double& N) {
            double ln_gsd = std::log(gsd);
            double V_single = (4.0/3.0) * M_PI * gmd * gmd * gmd
                             * std::exp(4.5 * ln_gsd * ln_gsd);
            N = N / V_single;
        },
        params_view.GetConstColumnView(gmd_param_idx),
        params_view.GetConstColumnView(gsd_param_idx),
        result);

    // Partials: ∂N/∂[species_k] = (Mw_k/ρ_k) / V_single
    for (size_t k = 0; k < species_var_indices.size(); ++k) {
        params_view.ForEachRow(
            [mwr = mw_over_rho[k]](const double& gmd, const double& gsd,
                                   double& dN_dspecies) {
                double ln_gsd = std::log(gsd);
                double V_single = (4.0/3.0) * M_PI * gmd * gmd * gmd
                                 * std::exp(4.5 * ln_gsd * ln_gsd);
                dN_dspecies = mwr / V_single;
            },
            params_view.GetConstColumnView(gmd_param_idx),
            params_view.GetConstColumnView(gsd_param_idx),
            partials.GetColumnView(k));
    }
};


// --- TwoMomentMode: NumberConcentration from state variable ---
// N = state_variables[NUMBER_CONCENTRATION]
// d(N)/d(N) = 1, so partials matrix has 1 column.

provider.dependent_variable_indices = { nc_var_idx };

provider.ComputeValue = [nc_var_idx]
    (auto&& params_view, auto&& vars_view, auto&& result)
{
    vars_view.ForEachRow(
        [](const double& nc, double& N) { N = nc; },
        vars_view.GetConstColumnView(nc_var_idx),
        result);
};

provider.ComputeValueAndDerivatives = [nc_var_idx]
    (auto&& params_view, auto&& vars_view, auto&& result, auto&& partials)
{
    vars_view.ForEachRow(
        [](const double& nc, double& N, double& dN_dN) {
            N = nc;
            dN_dN = 1.0;
        },
        vars_view.GetConstColumnView(nc_var_idx),
        result,
        partials.GetColumnView(0));  // column 0 = d(N)/d(nc_var)
};


// --- TwoMomentMode: EffectiveRadius from V_total, N, GSD ---
// V_total = Σ_p Σ_i [species_{p,i}] · Mw_{p,i} / ρ_{p,i}
//   (sum over all phases p in the mode, all species i in each phase)
// V_mean = V_total / N
// r_mean = (3 · V_mean / (4π))^(1/3)
// r_eff = r_mean · exp(2.5 · ln²(GSD))
//
// dependent_variable_indices = { all species indices across all phases..., nc_var_idx }
// Partials: ∂r_eff/∂[species_{p,i}] = r_eff / (3·V_total) · (Mw_{p,i}/ρ_{p,i})
//           ∂r_eff/∂N = -r_eff / (3·N)

provider.dependent_variable_indices = species_and_N_indices;  // species first, N last

provider.ComputeValue = [gsd_param_idx, species_var_indices,
                         nc_var_idx, mw_over_rho]
    (auto&& params_view, auto&& vars_view, auto&& result)
{
    // Accumulate V_total
    params_view.ForEachRow(
        [](double& r) { r = 0.0; },
        result);
    for (size_t k = 0; k < species_var_indices.size(); ++k) {
        params_view.ForEachRow(
            [mwr = mw_over_rho[k]](const double& conc, double& V) {
                V += conc * mwr;
            },
            vars_view.GetConstColumnView(species_var_indices[k]),
            result);
    }

    // r_eff = (3·V_total/(4π·N))^(1/3) · exp(2.5·ln²(GSD))
    params_view.ForEachRow(
        [](const double& gsd, const double& nc, double& r_eff) {
            double ln_gsd = std::log(gsd);
            double V_mean = r_eff / nc;  // r_eff held V_total
            double r_mean = std::cbrt(3.0 * V_mean / (4.0 * M_PI));
            r_eff = r_mean * std::exp(2.5 * ln_gsd * ln_gsd);
        },
        params_view.GetConstColumnView(gsd_param_idx),
        vars_view.GetConstColumnView(nc_var_idx),
        result);
};

provider.ComputeValueAndDerivatives = [gsd_param_idx, species_var_indices,
                                       nc_var_idx, mw_over_rho]
    (auto&& params_view, auto&& vars_view, auto&& result, auto&& partials)
{
    // ComputeValue (same as above)
    params_view.ForEachRow([](double& r) { r = 0.0; }, result);
    for (size_t k = 0; k < species_var_indices.size(); ++k) {
        params_view.ForEachRow(
            [mwr = mw_over_rho[k]](const double& conc, double& V) {
                V += conc * mwr;
            },
            vars_view.GetConstColumnView(species_var_indices[k]),
            result);
    }
    // Store V_total before overwriting result with r_eff
    // (need V_total for partials — use a separate temporary or
    //  recompute; shown here with the V_total still in result)
    // After computing r_eff, partials use the formula directly:
    //   ∂r_eff/∂[species_k] = r_eff / (3·V_total) · (Mw_k/ρ_k)
    //   ∂r_eff/∂N = -r_eff / (3·N)

    // Compute r_eff and write partials for species in one pass
    // (V_total is in result at this point)
    for (size_t k = 0; k < species_var_indices.size(); ++k) {
        params_view.ForEachRow(
            [mwr = mw_over_rho[k]](const double& gsd, const double& nc,
                                   const double& V_total,
                                   double& dr_dspecies) {
                double ln_gsd = std::log(gsd);
                double V_mean = V_total / nc;
                double r_mean = std::cbrt(3.0 * V_mean / (4.0 * M_PI));
                double r_eff = r_mean * std::exp(2.5 * ln_gsd * ln_gsd);
                dr_dspecies = r_eff * mwr / (3.0 * V_total);
            },
            params_view.GetConstColumnView(gsd_param_idx),
            vars_view.GetConstColumnView(nc_var_idx),
            result,
            partials.GetColumnView(k));
    }

    // Partial w.r.t. N (last column in partials)
    size_t N_col = species_var_indices.size();  // N is last in dependent list
    params_view.ForEachRow(
        [](const double& gsd, const double& nc, const double& V_total,
           double& dr_dN) {
            double ln_gsd = std::log(gsd);
            double V_mean = V_total / nc;
            double r_mean = std::cbrt(3.0 * V_mean / (4.0 * M_PI));
            double r_eff = r_mean * std::exp(2.5 * ln_gsd * ln_gsd);
            dr_dN = -r_eff / (3.0 * nc);
        },
        params_view.GetConstColumnView(gsd_param_idx),
        vars_view.GetConstColumnView(nc_var_idx),
        result,
        partials.GetColumnView(N_col));

    // Now overwrite result with r_eff
    params_view.ForEachRow(
        [](const double& gsd, const double& nc, double& r_eff) {
            double ln_gsd = std::log(gsd);
            double V_mean = r_eff / nc;  // r_eff still holds V_total
            double r_mean = std::cbrt(3.0 * V_mean / (4.0 * M_PI));
            r_eff = r_mean * std::exp(2.5 * ln_gsd * ln_gsd);
        },
        params_view.GetConstColumnView(gsd_param_idx),
        vars_view.GetConstColumnView(nc_var_idx),
        result);
};


// --- UniformSection: EffectiveRadius from min/max radius parameters ---
// r_eff = (r_min + r_max) / 2
// No state variable dependencies → partials matrix has zero columns.

provider.dependent_variable_indices = {};

provider.ComputeValue = [rmin_param_idx, rmax_param_idx]
    (auto&& params_view, auto&& vars_view, auto&& result)
{
    params_view.ForEachRow(
        [](const double& r_min, const double& r_max, double& r_eff) {
            r_eff = 0.5 * (r_min + r_max);
        },
        params_view.GetConstColumnView(rmin_param_idx),
        params_view.GetConstColumnView(rmax_param_idx),
        result);
};

provider.ComputeValueAndDerivatives = [rmin_param_idx, rmax_param_idx]
    (auto&& params_view, auto&& vars_view, auto&& result, auto&& partials)
{
    // Same computation; partials matrix has 0 columns, nothing to write
    params_view.ForEachRow(
        [](const double& r_min, const double& r_max, double& r_eff) {
            r_eff = 0.5 * (r_min + r_max);
        },
        params_view.GetConstColumnView(rmin_param_idx),
        params_view.GetConstColumnView(rmax_param_idx),
        result);
};


// --- UniformSection: NumberConcentration from total volume ---
// V_total = Σ_p Σ_i [species_{p,i}] · Mw_{p,i} / ρ_{p,i}
// V_single = (4/3)π · r_eff³
// N = V_total / V_single
// ∂N/∂[species_{p,i}] = (Mw_{p,i} / ρ_{p,i}) / V_single
//
// dependent_variable_indices = { all species variable indices across all phases }

provider.dependent_variable_indices = species_var_indices;

provider.ComputeValue = [rmin_param_idx, rmax_param_idx,
                         species_var_indices, mw_over_rho]
    (auto&& params_view, auto&& vars_view, auto&& result)
{
    // Accumulate V_total
    params_view.ForEachRow(
        [](double& N) { N = 0.0; },
        result);
    for (size_t k = 0; k < species_var_indices.size(); ++k) {
        params_view.ForEachRow(
            [mwr = mw_over_rho[k]](const double& conc, double& V) {
                V += conc * mwr;
            },
            vars_view.GetConstColumnView(species_var_indices[k]),
            result);
    }

    // N = V_total / V_single
    params_view.ForEachRow(
        [](const double& r_min, const double& r_max, double& N) {
            double r_eff = 0.5 * (r_min + r_max);
            double V_single = (4.0/3.0) * M_PI * r_eff * r_eff * r_eff;
            N = N / V_single;  // N held V_total
        },
        params_view.GetConstColumnView(rmin_param_idx),
        params_view.GetConstColumnView(rmax_param_idx),
        result);
};

provider.ComputeValueAndDerivatives = [rmin_param_idx, rmax_param_idx,
                                       species_var_indices, mw_over_rho]
    (auto&& params_view, auto&& vars_view, auto&& result, auto&& partials)
{
    // Compute V_total (same as ComputeValue)
    params_view.ForEachRow(
        [](double& N) { N = 0.0; },
        result);
    for (size_t k = 0; k < species_var_indices.size(); ++k) {
        params_view.ForEachRow(
            [mwr = mw_over_rho[k]](const double& conc, double& V) {
                V += conc * mwr;
            },
            vars_view.GetConstColumnView(species_var_indices[k]),
            result);
    }

    // N = V_total / V_single
    params_view.ForEachRow(
        [](const double& r_min, const double& r_max, double& N) {
            double r_eff = 0.5 * (r_min + r_max);
            double V_single = (4.0/3.0) * M_PI * r_eff * r_eff * r_eff;
            N = N / V_single;
        },
        params_view.GetConstColumnView(rmin_param_idx),
        params_view.GetConstColumnView(rmax_param_idx),
        result);

    // Partials: ∂N/∂[species_k] = (Mw_k/ρ_k) / V_single
    for (size_t k = 0; k < species_var_indices.size(); ++k) {
        params_view.ForEachRow(
            [mwr = mw_over_rho[k]](const double& r_min, const double& r_max,
                                   double& dN_dspecies) {
                double r_eff = 0.5 * (r_min + r_max);
                double V_single = (4.0/3.0) * M_PI * r_eff * r_eff * r_eff;
                dN_dspecies = mwr / V_single;
            },
            params_view.GetConstColumnView(rmin_param_idx),
            params_view.GetConstColumnView(rmax_param_idx),
            partials.GetColumnView(k));
    }
};


// --- All representations: PhaseVolumeFraction ---
// V_phase = Σ_i ( [species_{p,i}] · Mw_{p,i} / ρ_{p,i} )  (target phase only)
// V_total = Σ_q Σ_i ( [species_{q,i}] · Mw_{q,i} / ρ_{q,i} )  (all phases)
// φ_p = V_phase / V_total
//
// ∂φ_p/∂[species_{p,i}] = (Mw_{p,i}/ρ_{p,i}) · (1 - φ_p) / V_total  (same phase)
// ∂φ_p/∂[species_{q,i}] = -(Mw_{q,i}/ρ_{q,i}) · φ_p / V_total        (other phase)
//
// dependent_variable_indices = { all species across all phases }
// phase_species_count = number of species in the target phase (first entries)
// mw_over_rho[k] = Mw_k / ρ_k for species k (all phases, target phase first)
//
// Note: when a mode has a single phase, φ_p = 1 and all partials are zero.
// The provider can short-circuit in that case.

provider.dependent_variable_indices = all_species_var_indices;
// Ordered: target phase species first (indices 0..phase_species_count-1),
// then other phases' species.

provider.ComputeValue = [all_species_var_indices, mw_over_rho,
                         phase_species_count]
    (auto&& params_view, auto&& vars_view, auto&& result)
{
    // result will hold V_total; use a temporary for V_phase
    auto V_phase = result;  // RowVariable — stack allocated

    params_view.ForEachRow(
        [](double& vt, double& vp) { vt = 0.0; vp = 0.0; },
        result, V_phase);

    for (size_t k = 0; k < all_species_var_indices.size(); ++k) {
        if (k < phase_species_count) {
            // Target phase: contributes to both V_phase and V_total
            params_view.ForEachRow(
                [mwr = mw_over_rho[k]](const double& conc,
                                       double& vt, double& vp) {
                    double vol = conc * mwr;
                    vt += vol;
                    vp += vol;
                },
                vars_view.GetConstColumnView(all_species_var_indices[k]),
                result, V_phase);
        } else {
            // Other phases: contribute to V_total only
            params_view.ForEachRow(
                [mwr = mw_over_rho[k]](const double& conc, double& vt) {
                    vt += conc * mwr;
                },
                vars_view.GetConstColumnView(all_species_var_indices[k]),
                result);
        }
    }

    // φ_p = V_phase / V_total
    params_view.ForEachRow(
        [](double& phi, const double& vp) {
            phi = (phi > 0.0) ? vp / phi : 1.0;  // phi held V_total
        },
        result, V_phase);
};

provider.ComputeValueAndDerivatives = [all_species_var_indices, mw_over_rho,
                                       phase_species_count]
    (auto&& params_view, auto&& vars_view, auto&& result, auto&& partials)
{
    auto V_phase = result;  // RowVariable — stack allocated
    auto V_total = result;  // will reuse result for this

    // Compute V_total and V_phase (same as ComputeValue)
    params_view.ForEachRow(
        [](double& vt, double& vp) { vt = 0.0; vp = 0.0; },
        V_total, V_phase);

    for (size_t k = 0; k < all_species_var_indices.size(); ++k) {
        if (k < phase_species_count) {
            params_view.ForEachRow(
                [mwr = mw_over_rho[k]](const double& conc,
                                       double& vt, double& vp) {
                    double vol = conc * mwr;
                    vt += vol;
                    vp += vol;
                },
                vars_view.GetConstColumnView(all_species_var_indices[k]),
                V_total, V_phase);
        } else {
            params_view.ForEachRow(
                [mwr = mw_over_rho[k]](const double& conc, double& vt) {
                    vt += conc * mwr;
                },
                vars_view.GetConstColumnView(all_species_var_indices[k]),
                V_total);
        }
    }

    // φ_p = V_phase / V_total
    params_view.ForEachRow(
        [](double& phi, const double& vp, const double& vt) {
            phi = (vt > 0.0) ? vp / vt : 1.0;
        },
        result, V_phase, V_total);

    // Partials
    for (size_t k = 0; k < all_species_var_indices.size(); ++k) {
        if (k < phase_species_count) {
            // Same phase: ∂φ_p/∂[species_{p,i}] = (Mw/ρ)·(1-φ_p)/V_total
            params_view.ForEachRow(
                [mwr = mw_over_rho[k]](const double& phi, const double& vt,
                                       double& dphi) {
                    dphi = mwr * (1.0 - phi) / vt;
                },
                result, V_total,
                partials.GetColumnView(k));
        } else {
            // Other phase: ∂φ_p/∂[species_{q,i}] = -(Mw/ρ)·φ_p/V_total
            params_view.ForEachRow(
                [mwr = mw_over_rho[k]](const double& phi, const double& vt,
                                       double& dphi) {
                    dphi = -mwr * phi / vt;
                },
                result, V_total,
                partials.GetColumnView(k));
        }
    }
};
//  PROCESS PARAMETER NAMES — declares state parameters that this
//  process contributes (one HLC per phase instance)
// ====================================================================

// Following the DissolvedReversibleReaction pattern: each temperature-
// dependent constant gets a unique state parameter name built from
// prefix + phase name + uuid + constant name.

std::set<std::string> ProcessParameterNames(
    const std::map<std::string, std::set<std::string>>& phase_prefixes) const
{
    std::set<std::string> names;
    auto it = phase_prefixes.find(condensed_phase_.name_);
    if (it != phase_prefixes.end()) {
        for (const auto& prefix : it->second) {
            // e.g., "SMALL_DROP.AQUEOUS.abc-123.hlc"
            names.insert(prefix + "." + condensed_phase_.name_
                         + "." + uuid_ + ".hlc");
            names.insert(prefix + "." + condensed_phase_.name_
                         + "." + uuid_ + ".temperature");
        }
    }
    return names;
}


// ====================================================================
//  UPDATE STATE PARAMETERS — recalculates HLC(T) from Conditions
//  Called whenever Conditions change (before each time step).
// ====================================================================

template<typename DenseMatrixPolicy>
auto UpdateStateParametersFunction(
    const std::map<std::string, std::set<std::string>>& phase_prefixes,
    const auto& state_parameter_indices) const
{
    // Collect per-instance HLC and temperature parameter indices
    std::vector<size_t> hlc_indices;
    std::vector<size_t> temp_indices;
    auto it = phase_prefixes.find(condensed_phase_.name_);
    if (it != phase_prefixes.end()) {
        for (const auto& prefix : it->second) {
            hlc_indices.push_back(state_parameter_indices.at(
                prefix + "." + condensed_phase_.name_
                + "." + uuid_ + ".hlc"));
            temp_indices.push_back(state_parameter_indices.at(
                prefix + "." + condensed_phase_.name_
                + "." + uuid_ + ".temperature"));
        }
    }

    DenseMatrixPolicy dummy{ 1, state_parameter_indices.size(), 0.0 };
    std::vector<micm::Conditions> dummy_conditions;

    return DenseMatrixPolicy::Function(
        [this, hlc_indices, temp_indices](auto&& conditions, auto&& params)
        {
            for (size_t i = 0; i < hlc_indices.size(); ++i) {
                params.ForEachRow(
                    [&](const micm::Conditions& cond, double& hlc, double& T)
                    {
                        hlc = henry_law_constant_(cond);
                        T = cond.temperature_;
                    },
                    conditions,
                    params.GetColumnView(hlc_indices[i]),
                    params.GetColumnView(temp_indices[i]));
            }
        },
        dummy_conditions, dummy);
}


// ====================================================================
//  PROCESS SETUP PHASE — called once when building solver functions
//  (inside the process's ForcingFunction() / JacobianFunction())
// ====================================================================

// The process receives a ProviderMap from the Model layer through the
// common interface:
//   ProviderMap = map<string, vector<AerosolPropertyProvider>>
//   key = phase prefix (e.g., "SMALL_DROP"), value = ordered providers
//   matching the order returned by RequiredAerosolProperties()
//
// The process unpacks this into PerInstanceData at function-build time
// (once), then the hot-path lambda captures only the flat vector.

struct PerInstanceData {
    size_t aq_species_idx;       // state variable index for condensed solute
    size_t solvent_species_idx;  // state variable index for condensed solvent
    size_t hlc_param_idx;        // state parameter index for HLC(T) [mol m⁻³ Pa⁻¹]
    size_t temperature_param_idx; // state parameter index for T [K]
    double solvent_molecular_weight;           // solvent molar mass [kg mol⁻¹]
    double solvent_density;          // solvent density [kg m⁻³]
    AerosolPropertyProvider<DenseMatrixPolicy> r_eff_provider;
    AerosolPropertyProvider<DenseMatrixPolicy> N_provider;
    AerosolPropertyProvider<DenseMatrixPolicy> phi_provider;  // phase volume fraction
    CondensationRateProvider cond_rate_provider;
};

// Built from the ProviderMap inside ForcingFunction() / JacobianFunction():
std::vector<PerInstanceData> instances;
for (const auto& [prefix, provider_vec] : providers) {
    PerInstanceData inst;
    inst.aq_species_idx = var_indices.at(
        prefix + "." + condensed_phase_.name_ + "." + condensed_species_.name_);
    inst.solvent_species_idx = var_indices.at(
        prefix + "." + condensed_phase_.name_ + "." + solvent_.name_);
    inst.hlc_param_idx = param_indices.at(
        prefix + "." + condensed_phase_.name_ + "." + uuid_ + ".hlc");
    inst.temperature_param_idx = param_indices.at(
        prefix + "." + condensed_phase_.name_ + "." + uuid_ + ".temperature");
    inst.solvent_molecular_weight = solvent_molecular_weight_;
    inst.solvent_density = solvent_density_;
    // Providers arrive in RequiredAerosolProperties() order:
    //   [0] = EffectiveRadius, [1] = NumberConcentration,
    //   [2] = PhaseVolumeFraction
    inst.r_eff_provider = provider_vec[0];
    inst.N_provider     = provider_vec[1];
    inst.phi_provider   = provider_vec[2];
    inst.cond_rate_provider = MakeCondensationRateProvider(D_g_, alpha_, gas_molecular_weight_);
    instances.push_back(std::move(inst));
}

// --- Determine Jacobian sparsity from providers ---
std::set<std::pair<size_t, size_t>> jacobian_nonzeros;

for (const auto& inst : instances) {
    // Direct: gas ↔ aq (solute)
    jacobian_nonzeros.insert({gas_idx, gas_idx});
    jacobian_nonzeros.insert({gas_idx, inst.aq_species_idx});
    jacobian_nonzeros.insert({inst.aq_species_idx, gas_idx});
    jacobian_nonzeros.insert({inst.aq_species_idx, inst.aq_species_idx});

    // Direct: rate depends on solvent through f_v = [solvent] · Mw/ρ
    jacobian_nonzeros.insert({gas_idx, inst.solvent_species_idx});
    jacobian_nonzeros.insert({inst.aq_species_idx, inst.solvent_species_idx});

    // Indirect: through aerosol properties
    for (size_t var_j : inst.r_eff_provider.dependent_variable_indices) {
        jacobian_nonzeros.insert({gas_idx, var_j});
        jacobian_nonzeros.insert({inst.aq_species_idx, var_j});
    }
    for (size_t var_j : inst.N_provider.dependent_variable_indices) {
        jacobian_nonzeros.insert({gas_idx, var_j});
        jacobian_nonzeros.insert({inst.aq_species_idx, var_j});
    }
    for (size_t var_j : inst.phi_provider.dependent_variable_indices) {
        jacobian_nonzeros.insert({gas_idx, var_j});
        jacobian_nonzeros.insert({inst.aq_species_idx, var_j});
    }
}

// --- Pre-allocate working DenseMatrixPolicy instances (ONCE) ---
// These are captured by the returned lambda and reused every call.
// Dimensions: (num_cells × num_dependent_variables_for_that_property)
// The number of columns matches dependent_variable_indices.size().

size_t max_r_deps = 0, max_N_deps = 0, max_phi_deps = 0;
for (const auto& inst : instances) {
    max_r_deps = std::max(max_r_deps,
                          inst.r_eff_provider.dependent_variable_indices.size());
    max_N_deps = std::max(max_N_deps,
                          inst.N_provider.dependent_variable_indices.size());
    max_phi_deps = std::max(max_phi_deps,
                            inst.phi_provider.dependent_variable_indices.size());
}
// Use dummy 1-row matrices at setup; real num_rows determined at invocation
DenseMatrixPolicy r_eff_partials{ 1, max_r_deps, 0.0 };
DenseMatrixPolicy N_partials{ 1, max_N_deps, 0.0 };
DenseMatrixPolicy phi_partials{ 1, max_phi_deps, 0.0 };


// ====================================================================
//  FORCING FUNCTION (HOT PATH)
//  Uses DenseMatrixPolicy::Function() — operates on GroupView objects.
//  Zero allocations. All providers/indices captured at setup time.
// ====================================================================

DenseMatrixPolicy dummy_params{ 1, state_parameter_indices.size(), 0.0 };
DenseMatrixPolicy dummy_vars{ 1, state_variable_indices.size(), 0.0 };

return DenseMatrixPolicy::Function(
    [instances, gas_idx]
    (auto&& state_parameters, auto&& state_variables, auto&& forcing_terms)
    {
        // Stack-allocated per-group temporaries (no heap allocation)
        auto r_eff = forcing_terms.GetRowVariable();
        auto N     = forcing_terms.GetRowVariable();
        auto phi   = forcing_terms.GetRowVariable();  // phase volume fraction
        auto k_cond_var = forcing_terms.GetRowVariable();
        auto k_evap_var = forcing_terms.GetRowVariable();
        auto net_rate   = forcing_terms.GetRowVariable();

        for (const auto& inst : instances) {

            // --- Compute aerosol properties for all cells in group ---
            inst.r_eff_provider.ComputeValue(
                state_parameters, state_variables, r_eff);
            inst.N_provider.ComputeValue(
                state_parameters, state_variables, N);
            inst.phi_provider.ComputeValue(
                state_parameters, state_variables, phi);

            // --- Compute condensation/evaporation rates ---
            // k_cond via provider (encapsulates Fuchs-Sutugin regime)
            inst.cond_rate_provider.ComputeValue(r_eff, N, k_cond_var);

            // net = φ_p · k_cond · [gas] - φ_p · k_evap · [aq] / f_v
            // where k_evap = k_cond / (HLC·R·T), f_v = [solvent]·Mw/ρ
            // HLC is pre-computed in UpdateStateParametersFunction
            double molar_volume = inst.solvent_molecular_weight / inst.solvent_density;

            state_parameters.ForEachRow(
                [&](const double& kc, const double& phi_p,
                    const double& hlc, const double& T,
                    const double& gas_conc, const double& aq_conc,
                    const double& solvent_conc,
                    double& k_evap, double& net)
                {
                    double kc_eff = phi_p * kc;
                    k_evap = kc_eff / (hlc * R_gas * T);
                    double f_v = solvent_conc * molar_volume;
                    net = kc_eff * gas_conc - k_evap * aq_conc / f_v;
                },
                k_cond_var,
                phi,
                state_parameters.GetConstColumnView(inst.hlc_param_idx),
                state_parameters.GetConstColumnView(inst.temperature_param_idx),
                state_variables.GetConstColumnView(gas_idx),
                state_variables.GetConstColumnView(inst.aq_species_idx),
                state_variables.GetConstColumnView(inst.solvent_species_idx),
                k_evap_var,
                net_rate);

            // --- Apply forcing: gas species loses, aq species gains ---
            // d[gas]/dt -= net,  d[aq]/dt += net
            // (solvent is NOT consumed by this process)
            state_variables.ForEachRow(
                [](const double& net, double& gas_forcing) {
                    gas_forcing -= net;
                },
                net_rate,
                forcing_terms.GetColumnView(gas_idx));

            state_variables.ForEachRow(
                [](const double& net, double& aq_forcing) {
                    aq_forcing += net;
                },
                net_rate,
                forcing_terms.GetColumnView(inst.aq_species_idx));
        }
    },
    dummy_params, dummy_vars, dummy_vars);


// ====================================================================
//  JACOBIAN FUNCTION (HOT PATH)
//  Uses SparseMatrixPolicy::Function() — operates on block views.
//  Partial derivatives stored in DenseMatrixPolicy instances
//  pre-allocated at setup time (r_eff_partials, N_partials, phi_partials).
// ====================================================================

// net = φ_p · k_cond · [gas] - φ_p · k_evap · [aq] / f_v
// where f_v = [solvent] · solvent_molecular_weight / solvent_density
//       φ_p = V_phase / V_total  (phase volume fraction)
//       k_evap = k_cond / (HLC · R · T)
//       HLC = pre-computed Henry's Law constant (state parameter)
//
// Let R = k_cond · [gas] - k_evap · [aq] / f_v  (un-scaled net rate)
// Then net = φ_p · R
//
// Derivatives of net rate (for one condensed-phase instance):
//   d(net)/d([gas])     = φ_p · k_cond
//   d(net)/d([aq])      = -φ_p · k_evap / f_v
//   d(net)/d([solvent]) = φ_p · k_evap · [aq] / (f_v · [solvent])
//   d(net)/d(var_j) via r_eff:
//       = φ_p · (dk_cond/d(r_eff)·[gas] - dk_evap/d(r_eff)·[aq]/f_v) · dr/dj
//   d(net)/d(var_j) via N:
//       = φ_p · (dk_cond/d(N)·[gas] - dk_evap/d(N)·[aq]/f_v) · dN/dj
//   d(net)/d(var_j) via φ_p:
//       = R · dφ_p/d(var_j)

return SparseMatrixPolicy::Function(
    [instances, gas_idx, jacobian_indices,
     r_eff_partials = std::move(r_eff_partials),
     N_partials = std::move(N_partials),
     phi_partials = std::move(phi_partials)]
    (auto&& state_parameters, auto&& state_variables,
     auto&& jacobian_values) mutable
    {
        auto r_eff    = jacobian_values.GetBlockVariable();
        auto N        = jacobian_values.GetBlockVariable();
        auto phi      = jacobian_values.GetBlockVariable();  // phase volume fraction
        auto k_cond   = jacobian_values.GetBlockVariable();
        auto k_evap   = jacobian_values.GetBlockVariable();
        auto f_v_var  = jacobian_values.GetBlockVariable();  // solvent volume fraction
        auto net_unscaled = jacobian_values.GetBlockVariable();  // R = k_cond·[gas] - k_evap·[aq]/f_v
        auto dk_cond_prop = jacobian_values.GetBlockVariable();
        auto dk_evap_prop = jacobian_values.GetBlockVariable();

        auto jac_id = jacobian_indices.indices_.AsVector().begin();

        for (const auto& inst : instances) {
            double molar_volume = inst.solvent_molecular_weight / inst.solvent_density;

            // --- Compute property values ---
            inst.r_eff_provider.ComputeValue(
                state_parameters, state_variables, r_eff);
            inst.N_provider.ComputeValue(
                state_parameters, state_variables, N);
            inst.phi_provider.ComputeValue(
                state_parameters, state_variables, phi);

            // --- Compute rate constants, volume fraction, and un-scaled rate ---
            inst.cond_rate_provider.ComputeValue(r_eff, N, k_cond);

            jacobian_values.ForEachBlock(
                [&](const double& kc, const double& hlc, const double& T,
                    const double& solvent_conc,
                    const double& gas_conc, const double& aq_conc,
                    double& ke, double& fv, double& R)
                {
                    ke = kc / (hlc * R_gas * T);
                    fv = solvent_conc * molar_volume;
                    R = kc * gas_conc - ke * aq_conc / fv;
                },
                k_cond,
                state_parameters.GetConstColumnView(inst.hlc_param_idx),
                state_parameters.GetConstColumnView(inst.temperature_param_idx),
                state_variables.GetConstColumnView(inst.solvent_species_idx),
                state_variables.GetConstColumnView(gas_idx),
                state_variables.GetConstColumnView(inst.aq_species_idx),
                k_evap, f_v_var, net_unscaled);

            // --- Direct Jacobian entries (scaled by φ_p) ---
            // J[gas, gas] -= φ_p · k_cond
            jacobian_values.ForEachBlock(
                [](const double& phi_p, const double& kc, double& jac) {
                    jac -= phi_p * kc;
                },
                phi, k_cond,
                jacobian_values.GetBlockView(*jac_id++));

            // J[gas, aq] += φ_p · k_evap / f_v
            jacobian_values.ForEachBlock(
                [](const double& phi_p, const double& ke, const double& fv,
                   double& jac) {
                    jac += phi_p * ke / fv;
                },
                phi, k_evap, f_v_var,
                jacobian_values.GetBlockView(*jac_id++));

            // J[gas, solvent] -= φ_p · k_evap · [aq] / (f_v · [solvent])
            jacobian_values.ForEachBlock(
                [](const double& phi_p, const double& ke, const double& fv,
                   const double& aq_conc, const double& solvent_conc,
                   double& jac)
                {
                    jac -= phi_p * ke * aq_conc / (fv * solvent_conc);
                },
                phi, k_evap, f_v_var,
                state_variables.GetConstColumnView(inst.aq_species_idx),
                state_variables.GetConstColumnView(inst.solvent_species_idx),
                jacobian_values.GetBlockView(*jac_id++));

            // J[aq, gas] += φ_p · k_cond
            jacobian_values.ForEachBlock(
                [](const double& phi_p, const double& kc, double& jac) {
                    jac += phi_p * kc;
                },
                phi, k_cond,
                jacobian_values.GetBlockView(*jac_id++));

            // J[aq, aq] -= φ_p · k_evap / f_v
            jacobian_values.ForEachBlock(
                [](const double& phi_p, const double& ke, const double& fv,
                   double& jac) {
                    jac -= phi_p * ke / fv;
                },
                phi, k_evap, f_v_var,
                jacobian_values.GetBlockView(*jac_id++));

            // J[aq, solvent] += φ_p · k_evap · [aq] / (f_v · [solvent])
            jacobian_values.ForEachBlock(
                [](const double& phi_p, const double& ke, const double& fv,
                   const double& aq_conc, const double& solvent_conc,
                   double& jac)
                {
                    jac += phi_p * ke * aq_conc / (fv * solvent_conc);
                },
                phi, k_evap, f_v_var,
                state_variables.GetConstColumnView(inst.aq_species_idx),
                state_variables.GetConstColumnView(inst.solvent_species_idx),
                jacobian_values.GetBlockView(*jac_id++));

            // --- Indirect Jacobian entries through r_eff (scaled by φ_p) ---
            if (!inst.r_eff_provider.dependent_variable_indices.empty()) {
                inst.r_eff_provider.ComputeValueAndDerivatives(
                    state_parameters, state_variables,
                    r_eff, r_eff_partials);

                // dk_cond_prop receives dk_cond/dr_eff (4th arg)
                // dk_evap_prop receives dk_cond/dN (5th arg, unused here)
                inst.cond_rate_provider.ComputeValueAndDerivatives(
                    r_eff, N, k_cond, dk_cond_prop, dk_evap_prop);

                // Convert dk_cond/dr to dk_evap/dr using HLC and T from state params
                jacobian_values.ForEachBlock(
                    [&](const double& dk_dr, const double& hlc, const double& T,
                        double& dk_evap_dr)
                    {
                        dk_evap_dr = dk_dr / (hlc * R_gas * T);
                    },
                    dk_cond_prop,
                    state_parameters.GetConstColumnView(inst.hlc_param_idx),
                    state_parameters.GetConstColumnView(inst.temperature_param_idx),
                    dk_evap_prop);

                for (size_t k = 0; k < inst.r_eff_provider.dependent_variable_indices.size(); ++k) {
                    // Chain rule: d(net)/d(var_j) = φ_p · (dk_cond/dr·[gas]
                    //                              - dk_evap/dr·[aq]/f_v) · dr/dj

                    // J[gas, var_j]
                    jacobian_values.ForEachBlock(
                        [](const double& phi_p,
                           const double& dk_dr, const double& dk_evap_dr,
                           const double& dr_dvar,
                           const double& gas_conc, const double& aq_conc,
                           const double& fv,
                           double& jac)
                        {
                            double dk_cond_dj = dk_dr * dr_dvar;
                            double dk_evap_dj = dk_evap_dr * dr_dvar;
                            jac -= phi_p * (dk_cond_dj * gas_conc
                                           - dk_evap_dj * aq_conc / fv);
                        },
                        phi,
                        dk_cond_prop,
                        dk_evap_prop,
                        r_eff_partials.GetConstColumnView(k),
                        state_variables.GetConstColumnView(gas_idx),
                        state_variables.GetConstColumnView(inst.aq_species_idx),
                        f_v_var,
                        jacobian_values.GetBlockView(*jac_id++));

                    // J[aq, var_j]
                    jacobian_values.ForEachBlock(
                        [](const double& phi_p,
                           const double& dk_dr, const double& dk_evap_dr,
                           const double& dr_dvar,
                           const double& gas_conc, const double& aq_conc,
                           const double& fv,
                           double& jac)
                        {
                            double dk_cond_dj = dk_dr * dr_dvar;
                            double dk_evap_dj = dk_evap_dr * dr_dvar;
                            jac += phi_p * (dk_cond_dj * gas_conc
                                           - dk_evap_dj * aq_conc / fv);
                        },
                        phi,
                        dk_cond_prop,
                        dk_evap_prop,
                        r_eff_partials.GetConstColumnView(k),
                        state_variables.GetConstColumnView(gas_idx),
                        state_variables.GetConstColumnView(inst.aq_species_idx),
                        f_v_var,
                        jacobian_values.GetBlockView(*jac_id++));
                }
            }

            // --- Indirect Jacobian entries through N (scaled by φ_p) ---
            if (!inst.N_provider.dependent_variable_indices.empty()) {
                inst.N_provider.ComputeValueAndDerivatives(
                    state_parameters, state_variables,
                    N, N_partials);

                // Same arg order as r_eff section:
                // dk_cond_prop receives dk_cond/dr_eff (4th arg, unused here)
                // dk_evap_prop receives dk_cond/dN (5th arg)
                inst.cond_rate_provider.ComputeValueAndDerivatives(
                    r_eff, N, k_cond, dk_cond_prop, dk_evap_prop);

                // Convert dk_cond/dN to dk_evap/dN using HLC and T from state params
                // dk_evap_prop already holds dk_cond/dN; overwrite dk_cond_prop
                // with dk_evap/dN so we can use dk_evap_prop and dk_cond_prop
                // as dk_cond/dN and dk_evap/dN respectively below.
                jacobian_values.ForEachBlock(
                    [&](const double& dk_dN, const double& hlc, const double& T,
                        double& dk_evap_dN)
                    {
                        dk_evap_dN = dk_dN / (hlc * R_gas * T);
                    },
                    dk_evap_prop,
                    state_parameters.GetConstColumnView(inst.hlc_param_idx),
                    state_parameters.GetConstColumnView(inst.temperature_param_idx),
                    dk_cond_prop);

                for (size_t k = 0; k < inst.N_provider.dependent_variable_indices.size(); ++k) {
                    // J[gas, var_j]
                    jacobian_values.ForEachBlock(
                        [](const double& phi_p,
                           const double& dk_dN, const double& dk_evap_dN,
                           const double& dN_dvar,
                           const double& gas_conc, const double& aq_conc,
                           const double& fv,
                           double& jac)
                        {
                            double dk_cond_dj = dk_dN * dN_dvar;
                            double dk_evap_dj = dk_evap_dN * dN_dvar;
                            jac -= phi_p * (dk_cond_dj * gas_conc
                                           - dk_evap_dj * aq_conc / fv);
                        },
                        phi,
                        dk_evap_prop,   // holds dk_cond/dN
                        dk_cond_prop,   // holds dk_evap/dN (computed above)
                        N_partials.GetConstColumnView(k),
                        state_variables.GetConstColumnView(gas_idx),
                        state_variables.GetConstColumnView(inst.aq_species_idx),
                        f_v_var,
                        jacobian_values.GetBlockView(*jac_id++));

                    // J[aq, var_j]
                    jacobian_values.ForEachBlock(
                        [](const double& phi_p,
                           const double& dk_dN, const double& dk_evap_dN,
                           const double& dN_dvar,
                           const double& gas_conc, const double& aq_conc,
                           const double& fv,
                           double& jac)
                        {
                            double dk_cond_dj = dk_dN * dN_dvar;
                            double dk_evap_dj = dk_evap_dN * dN_dvar;
                            jac += phi_p * (dk_cond_dj * gas_conc
                                           - dk_evap_dj * aq_conc / fv);
                        },
                        phi,
                        dk_evap_prop,   // holds dk_cond/dN
                        dk_cond_prop,   // holds dk_evap/dN
                        N_partials.GetConstColumnView(k),
                        state_variables.GetConstColumnView(gas_idx),
                        state_variables.GetConstColumnView(inst.aq_species_idx),
                        f_v_var,
                        jacobian_values.GetBlockView(*jac_id++));
                }
            }

            // --- Indirect Jacobian entries through φ_p ---
            // net = φ_p · R where R = k_cond·[gas] - k_evap·[aq]/f_v
            // d(net)/d(var_j) via φ_p = R · dφ_p/d(var_j)
            if (!inst.phi_provider.dependent_variable_indices.empty()) {
                inst.phi_provider.ComputeValueAndDerivatives(
                    state_parameters, state_variables,
                    phi, phi_partials);

                for (size_t k = 0; k < inst.phi_provider.dependent_variable_indices.size(); ++k) {
                    // J[gas, var_j] -= R · dφ_p/d(var_j)
                    jacobian_values.ForEachBlock(
                        [](const double& R, const double& dphi_dvar,
                           double& jac)
                        {
                            jac -= R * dphi_dvar;
                        },
                        net_unscaled,
                        phi_partials.GetConstColumnView(k),
                        jacobian_values.GetBlockView(*jac_id++));

                    // J[aq, var_j] += R · dφ_p/d(var_j)
                    jacobian_values.ForEachBlock(
                        [](const double& R, const double& dphi_dvar,
                           double& jac)
                        {
                            jac += R * dphi_dvar;
                        },
                        net_unscaled,
                        phi_partials.GetConstColumnView(k),
                        jacobian_values.GetBlockView(*jac_id++));
                }
            }
        }
    },
    dummy_params, dummy_vars, jacobian);
```

### Performance characteristics

- **Setup time** (called once): Index lookups, provider creation,
  `DenseMatrixPolicy` allocation for partials all happen here.
  `std::map` lookups for prefix→representation mapping are acceptable.
- **Hot path** (forcing/Jacobian, called every solver iteration):
  - Zero memory allocations — all working storage pre-allocated.
  - Vectorized across grid cells: `ForEachRow()` / `ForEachBlock()` process
    L cells per group (matching `VectorMatrix` group size).
  - `std::function` overhead: one indirect call per group per property
    (amortized over L cells).
  - Partial derivatives stored in `DenseMatrixPolicy` columns — same
    vectorized access pattern, no raw pointer arithmetic.
  - `RowVariable` / `BlockVariable` temporaries are stack-allocated per group.
  - Per-instance loop is inside the `Function()` lambda — all instances
    share the same group iteration without redundant setup.

### Species property requirements

Aerosol property calculations require per-species physical constants that are
stored as properties on the `micm::Species` class via its `properties_double_`
map (accessed with `species.GetProperty<double>(key)`).

| Property | Key | Units | Used for |
|----------|-----|-------|----------|
| Density | `"density"` | `kg m⁻³` | Converting species concentration to volume |
| Molecular weight | `"molecular weight"` | `kg mol⁻¹` | Converting species concentration to volume |

The total volume of a mode/section is computed from the species concentrations
across all phase instances it contains. A single mode or section may include
one or more phases (e.g., an accumulation mode with both aqueous and organic
phases). Every species in every phase contributes to the total particle
volume:

```
V_total = Σ_p Σ_i ( [species_{p,i}] · Mw_{p,i} / ρ_{p,i} )   [m³ m⁻³ air]
```

where `p` indexes phases within the mode/section, `i` indexes species within
phase `p`, `[species_{p,i}]` is in `mol m⁻³ air`, `Mw_{p,i}` is molecular
weight (`kg mol⁻¹`), and `ρ_{p,i}` is density (`kg m⁻³`).

These properties must be set on species before they are used in representations
that compute volume-dependent aerosol properties:

```cpp
auto h2o = Species{ "H2O", {{ "density", 1000.0 }, { "molecular weight", 0.018 }} };
auto nacl = Species{ "NaCl", {{ "density", 2165.0 }, { "molecular weight", 0.058 }} };
```

### Property availability by representation type

| Property             | SingleMomentMode              | TwoMomentMode                         | UniformSection                 |
|----------------------|-------------------------------|---------------------------------------|--------------------------------|
| EffectiveRadius      | From GMD, GSD (params); no state var deps | From V_total, N, GSD; depends on all species (all phases) + N | From min/max radius (params); no state var deps |
| NumberConcentration   | From V_total, GMD, GSD; depends on all species (all phases) | Direct state variable (∂N/∂N = 1)  | From V_total, r_min, r_max; depends on all species (all phases) |
| PhaseVolumeFraction   | V_phase / V_total; depends on all species (all phases) | V_phase / V_total; depends on all species (all phases) | V_phase / V_total; depends on all species (all phases) |

**Formulas:**

*SingleMomentMode — EffectiveRadius (parameterized):*
```
r_eff = GMD · exp(2.5 · ln²(GSD))
∂r_eff/∂[species_i] = 0  (no state variable dependencies)
```

*SingleMomentMode — NumberConcentration (from total volume):*
```
V_total = Σ_p Σ_i ( [species_{p,i}] · Mw_{p,i} / ρ_{p,i} )   (sum over all phases p, species i)
V_single = (4/3)π · GMD³ · exp(4.5 · ln²(GSD))
N = V_total / V_single
∂N/∂[species_{p,i}] = (Mw_{p,i} / ρ_{p,i}) / V_single
```

*TwoMomentMode — EffectiveRadius (from volume and number):*
```
V_total = Σ_p Σ_i ( [species_{p,i}] · Mw_{p,i} / ρ_{p,i} )   (sum over all phases p, species i)
V_mean = V_total / N
r_mean = (3 · V_mean / (4π))^(1/3)
r_eff = r_mean · exp(2.5 · ln²(GSD))
∂r_eff/∂[species_{p,i}] = r_eff / (3 · V_total) · (Mw_{p,i} / ρ_{p,i})
∂r_eff/∂N = -r_eff / (3 · N)
```

*TwoMomentMode — NumberConcentration (state variable):*
```
N = state_variables[NUMBER_CONCENTRATION]
∂N/∂N = 1
```

*UniformSection — EffectiveRadius (parameterized):*
```
r_eff = (r_min + r_max) / 2   (or geometric mean, depending on convention)
∂r_eff/∂[species_i] = 0
```

*UniformSection — NumberConcentration (from total volume):*
```
V_total = Σ_p Σ_i ( [species_{p,i}] · Mw_{p,i} / ρ_{p,i} )   (sum over all phases p, species i)
V_single = (4/3)π · r_eff³
N = V_total / V_single
∂N/∂[species_{p,i}] = (Mw_{p,i} / ρ_{p,i}) / V_single
```

*All representations — PhaseVolumeFraction (from species volumes):*
```
V_phase = Σ_i ( [species_{p,i}] · Mw_{p,i} / ρ_{p,i} )   (species in target phase p only)
V_total = Σ_q Σ_i ( [species_{q,i}] · Mw_{q,i} / ρ_{q,i} )   (all phases)
φ_p = V_phase / V_total

∂φ_p/∂[species_{p,i}] = (Mw_{p,i}/ρ_{p,i}) · (1 - φ_p) / V_total   (same phase)
∂φ_p/∂[species_{q,i}] = -(Mw_{q,i}/ρ_{q,i}) · φ_p / V_total          (other phase q ≠ p)
```

Note: when a mode has only a single phase, φ_p = 1, all partials are zero,
and the provider can short-circuit.

### Open questions

1. **Model-level routing** *(resolved)*: How should the Model route a
   `(prefix, property)` request to the correct representation?

   **Resolution: Model builds providers generically; processes receive them
   through the common interface.**

   No MICM changes are required. MICM's `ExternalModelProcessSet` already calls
   `Model.ForcingFunction(state_parameter_indices, state_variable_indices)` with
   **unified** global index maps containing both gas-phase species and all
   external model variables. The Model handles routing internally.

   **Workflow (generic for all process types):**

   1. Model stores all processes in a single `processes_` vector using the
      type-erased `MiamProcessSet` wrapper (see Common Process Interface above).
   2. When MICM calls `Model.ForcingFunction()` (or `JacobianFunction()`),
      the Model iterates `processes_` uniformly. For each process:
      - Calls `required_aerosol_properties_()` to discover what properties
        the process needs (may be empty).
      - For each requested `(phase_name, property)`, iterates
        `representations_` via `std::visit`, finds the representation owning
        each prefix, and calls `GetPropertyProvider()` to build providers.
      - Bundles providers into a `ProviderMap` (keyed by phase prefix).
      - Passes the `ProviderMap` to the process's `get_forcing_function_()`
        (or `get_jacobian_function_()`).
   3. The process builds its lambda using the providers. Processes that
      returned empty from `RequiredAerosolProperties()` receive an empty map
      and ignore it.

   This runs at function-build time (once), not on the hot path. The
   `std::visit` over representations and string lookups are acceptable here.

   **Model-level pseudocode:**

   ```cpp
   // Model now stores a single process collection:
   //   std::vector<process::DissolvedReversibleReaction> dissolved_reactions_;
   // is replaced by the type-erased collection built at AddProcesses() time.

   template<typename DenseMatrixPolicy>
   auto ForcingFunction(
       const std::unordered_map<std::string, std::size_t>& state_parameter_indices,
       const std::unordered_map<std::string, std::size_t>& state_variable_indices) const
   {
     auto phase_prefixes = CollectPhaseStatePrefixes();
     std::vector<std::function<void(const DenseMatrixPolicy&,
                                    const DenseMatrixPolicy&,
                                    DenseMatrixPolicy&)>> forcing_functions;

     // --- Single loop over all process types ---
     for (const auto& process : processes_)
     {
       // Build providers for whatever this process needs
       auto providers = BuildProviders<DenseMatrixPolicy>(
           process.required_aerosol_properties_(),
           phase_prefixes,
           state_parameter_indices, state_variable_indices);

       forcing_functions.push_back(
           process.get_forcing_function_(
               phase_prefixes, state_parameter_indices,
               state_variable_indices, std::move(providers)));
     }

     return [forcing_functions](const DenseMatrixPolicy& params,
                                const DenseMatrixPolicy& vars,
                                DenseMatrixPolicy& forcing) {
       for (const auto& fn : forcing_functions)
         fn(params, vars, forcing);
     };
   }
   ```

   **`BuildProviders` helper (private Model method):**

   ```cpp
   // Builds aerosol property providers for the properties a process requests.
   // Returns an empty map if the process doesn't need any.
   template<typename DenseMatrixPolicy>
   ProviderMap<DenseMatrixPolicy> BuildProviders(
       const std::map<std::string, std::vector<AerosolProperty>>& required,
       const std::map<std::string, std::set<std::string>>& phase_prefixes,
       const std::unordered_map<std::string, std::size_t>& param_indices,
       const std::unordered_map<std::string, std::size_t>& var_indices) const
   {
     ProviderMap<DenseMatrixPolicy> result;
     if (required.empty()) return result;

     for (const auto& [phase_name, properties] : required)
     {
       auto pp_it = phase_prefixes.find(phase_name);
       if (pp_it == phase_prefixes.end())
         throw std::runtime_error("Phase not found: " + phase_name);

       for (const auto& prefix : pp_it->second)
       {
         // Find the representation that owns this prefix
         for (const auto& repr : representations_)
         {
           std::visit(
               [&](const auto& r)
               {
                 auto repr_prefixes = r.PhaseStatePrefixes();
                 auto phase_it = repr_prefixes.find(phase_name);
                 if (phase_it != repr_prefixes.end() &&
                     phase_it->second.count(prefix))
                 {
                   for (const auto& prop : properties)
                   {
                     result[prefix].push_back(
                         r.template GetPropertyProvider<DenseMatrixPolicy>(
                             prop, param_indices, var_indices));
                   }
                 }
               },
               repr);
         }
       }
     }
     return result;
   }
   ```

   **Why no MICM changes:**

   - The `ExternalModelProcessSet` interface (`ForcingFunction`, `JacobianFunction`,
     etc.) already receives unified state index maps covering all gas-phase and
     external model variables.
   - The gas species index (e.g., `"CO2"`) is in the unified
     `state_variable_indices` map because MICM's `SetSystem()` registers both gas
     and external model variables.
   - All routing complexity is inside `Model`, invisible to MICM.

   **Implications for implementation steps:**

   No separate MICM PR is needed. The common process interface (Step 0) and
   Model integration (Step 5) cover this:
   - Step 0: Define `MiamProcessSet` and the common interface; refactor
     `DissolvedReversibleReaction` to satisfy it.
   - Step 4–5: `HenryLawPhaseTransfer` implements the same interface;
     Model's single loop handles both process types identically.

2. **Provider GroupView compatibility** *(resolved)*: The providers receive
   `auto&&` view objects and call `ForEachRow()` / `GetConstColumnView()` on
   them. Can these be called from inside a `SparseMatrixPolicy::Function()`
   lambda (for the Jacobian) where the dense state data arrives as a different
   view type?

   **Resolution: Yes — dense matrices can be passed to
   `SparseMatrixPolicy::Function()`.**

   MICM's `SparseMatrixPolicy::Function()` accepts mixed argument types: sparse
   matrices, dense matrices, and vectors. Inside the lambda, dense column views
   (`dense.GetConstColumnView()`) can be used alongside sparse block views
   (`sparse.GetConstBlockView()`) in the same `ForEachBlock()` call. The only
   constraints are:
   - Row/block counts must match across all arguments.
   - Ordering policy `L` values must match (validated at function-creation time).

   This is confirmed by MICM's `testSparseAndDenseMatrixFunction` and
   `testSparseAndVectorMatrixFunction` tests
   ([test_sparse_matrix_policy.hpp](https://github.com/NCAR/micm/blob/bd284ea053c10c6e961aca2da86cb4dbf2e00aaf/test/unit/util/test_sparse_matrix_policy.hpp#L1350-L1442)).

   For the Jacobian function, this means the provider lambdas can use
   `GetConstColumnView()` on the dense state parameter/variable views passed
   to `SparseMatrixPolicy::Function()` without any re-wrapping. No special
   handling is needed.

3. **Future properties** *(resolved)*: Properties like `Density` or
   `SurfaceArea` may depend on multiple state variables with more complex
   derivative expressions.

   **Resolution: The provider design already handles multi-variable
   dependencies.** `EffectiveRadius` (TwoMomentMode) depends on all species
   concentrations plus number concentration. `NumberConcentration`
   (SingleMomentMode, UniformSection) depends on all species concentrations.
   The `DenseMatrixPolicy` partials matrix scales naturally — one column per
   dependent variable. The `ForEachRow` / `ForEachBlock` arg-passing pattern
   works by looping over species indices at setup time, not by passing
   variable numbers of args in a single call. See the updated provider
   implementations above for concrete examples.

   Species-level physical constants (`"density"` and `"molecular weight"`)
   are obtained from `micm::Species::GetProperty<double>()` at provider
   creation time and captured as `mw_over_rho` vectors in the provider
   lambdas — no runtime lookups.

4. **Transition regime derivatives** *(resolved)*: `dk_cond/dr_eff` in the
   transition regime involves `df/dKn · dKn/dr`. Should this be the process's
   responsibility or should the mass transfer utility provide a combined
   function?

   **Resolution: Use a `CondensationRateProvider` — the same provider pattern
   as aerosol properties.** The mass transfer utility provides a factory
   function that returns a provider struct with `ComputeValue` and
   `ComputeValueAndDerivatives` lambdas:

   ```cpp
   struct CondensationRateProvider {
       // Compute k_cond from r_eff, N, and physical parameters.
       // result receives k_cond.
       std::function<void(/* r_eff view, N view, result view */)> ComputeValue;

       // Compute k_cond and partial derivatives dk_cond/dr_eff, dk_cond/dN
       // in a single pass, sharing intermediates (Kn, f(Kn), df/dKn).
       // result receives k_cond; dk_dr receives dk_cond/dr_eff;
       // dk_dN receives dk_cond/dN.
       std::function<void(/* r_eff view, N view,
                             result view, dk_dr view, dk_dN view */)>
           ComputeValueAndDerivatives;
   };
   ```

   **Factory function** (regime-specific logic baked in at setup):

   ```cpp
   // include/miam/util/condensation_rate.hpp
   CondensationRateProvider MakeCondensationRateProvider(
       double D_g, double alpha, double gas_molecular_weight);
   ```

   The returned lambdas capture the physical constants and encapsulate the
   full Fuchs-Sutugin transition-regime calculation:

   ```
   c̄  = √(8RT / (πMw))
   λ  = 3·D_g / c̄
   Kn = λ / r_eff
   f  = (1 + Kn) / (1 + 2·Kn·(1 + Kn) / α)
   k_cond = 4π · r_eff · N · D_g · f

   Derivatives (shared intermediates — computed once):
     df/dKn  = [α - 2 - 4·Kn - 2·Kn²] / [α·(1 + 2·Kn·(1+Kn)/α)²]
     dKn/dr  = -Kn / r_eff

     dk_cond/dr_eff = 4π·N·D_g · (f + r_eff · df/dKn · dKn/dr)
     dk_cond/dN     = k_cond / N   (linear in N)
   ```

   **Benefits:**
   - Forcing function calls `ComputeValue` — no derivative overhead
   - Jacobian calls `ComputeValueAndDerivatives` — intermediates computed
     once, both partials extracted in the same pass
   - Regime-specific math is self-contained in the utility, not spread
     across the process
   - Same structural pattern as `AerosolPropertyProvider`, so the process
     code treats both uniformly

---

## Implementation Steps

| PR  | Steps | Description |
|-----|-------|-------------|
| PR A | 0     | Common process interface + refactor DissolvedReversibleReaction |
| PR B | 1     | Aerosol property providers on representations |
| PR C | 2–3   | Henry's Law constant + mass transfer utilities |
| PR D | 4     | Core HenryLawPhaseTransfer process |
| PR E | 5–6   | Builder, model integration, integration tests |
| PR F | 7     | Documentation |

### Step 0 — Common Process Interface

Define the `MiamProcessSet` type-erased wrapper and the common method
interface (see Common Process Interface section above). Refactor Model to
store a single `std::vector<MiamProcessSet>` and use one loop for all
process methods (`StateSize`, `ProcessParameterNames`, `SpeciesUsed`,
`NonZeroJacobianElements`, `ForcingFunction`, `JacobianFunction`). Add the
generic `BuildProviders()` helper to Model.

Refactor `DissolvedReversibleReaction` to satisfy the interface:
- Add `RequiredAerosolProperties()` returning `{}`.
- Add thin overloads of `ForcingFunction`, `JacobianFunction`, and
  `NonZeroJacobianElements` that accept a `ProviderMap` (ignored) and
  delegate to the existing implementations.

Existing tests must continue to pass unchanged after the refactor.

### Step 1 — Aerosol Property Providers

Add the `AerosolProperty` enum, `AerosolPropertyProvider` struct, and
`GetPropertyProvider()` method to each representation type. Providers use
`ForEachRow()` with `RowVariable` outputs and `GetConstColumnView()` inputs.
Partial derivatives written into `DenseMatrixPolicy` column views. Unit
test each provider's `ComputeValue()` and `ComputeValueAndDerivatives()`
methods.

### Step 2 — Henry's Law Constant

`include/miam/processes/constants/henrys_law_constant.hpp` — Temperature-dependent
Henry's Law constant: `HLC(T) = HLC_ref · exp(C · (1/T - 1/T_ref))`, where
`HLC_ref` is the Henry's Law constant at reference temperature `T_ref` and `C`
is a parameter defining the temperature dependence.
The process stores a `henry_law_constant_` callable (same pattern as
`DissolvedReversibleReaction`'s `forward_rate_constant_`), which is invoked by
`UpdateStateParametersFunction` each time conditions change. The computed value
is written to a state parameter column (`hlc`) for each phase instance.

### Step 3 — Mass Transfer Rate Utilities

`include/miam/util/condensation_rate.hpp` — `CondensationRateProvider` struct
with `ComputeValue` and `ComputeValueAndDerivatives` lambdas, plus factory
function `MakeCondensationRateProvider(D_g, alpha, gas_molecular_weight)`. Encapsulates
mean molecular speed, mean free path, Knudsen number, Fuchs-Sutugin correction,
and all derivatives (`dk_cond/dr_eff`, `dk_cond/dN`) — same provider pattern as
aerosol properties. Helper functions (mean free path, Fuchs-Sutugin factor) are
internal implementation details, not part of the public interface.

### Step 4 — HenryLawPhaseTransfer Process

Core process struct implementing the common interface:
`RequiredAerosolProperties()` returns `{phase_name → {EffectiveRadius,
NumberConcentration, PhaseVolumeFraction}}`, `ProcessParameterNames()`, `SpeciesUsed()`,
`NonZeroJacobianElements()`, `UpdateStateParametersFunction()`,
`ForcingFunction()`, `JacobianFunction()`. Receives aerosol property
providers from Model through the `ProviderMap` argument — same path as
every other process.

### Step 5 — Builder + Model Integration

Fluent builder API. `Model::AddProcesses()` wraps the process in a
`MiamProcessSet` and appends to the unified `processes_` vector. Model's
single loop in `ForcingFunction()` / `JacobianFunction()` calls
`BuildProviders()` and dispatches to the process — no process-type-specific
code in Model.

### Step 6 — Integration Tests

1. Simple 1-species, 1-instance with analytical solution comparison
2. Multi-instance with mass conservation verification
3. Temperature-dependent Henry's Law constant verification
4. Continuum vs. transition regime rate comparison

### Step 7 — Documentation

Update docs/ and README with new process type and usage examples.
