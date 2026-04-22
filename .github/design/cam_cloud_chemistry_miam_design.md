# CAM Cloud Chemistry → MIAM Design Document

## 1. Overview

This document describes the cloud aqueous chemistry system currently implemented in
CAM (`mo_setsox.F90` / `cloud_aqueous_chemistry.F90`), maps it onto MIAM's process
and constraint framework, identifies gaps, and provides a draft MIAM configuration
script.

**CAM's approach**: Iterative bisection for pH (charge balance), followed by
forward-Euler oxidation reactions, all within a single operator-split timestep.

**MIAM's approach**: DAE (differential-algebraic equation) system solved by a
Rosenbrock method. Equilibrium partitioning and acid-base dissociation are algebraic
constraints; oxidation reactions are ODE terms. The solver enforces all constraints
simultaneously.

---

## 2. Chemical System in CAM Cloud Chemistry

### 2.1 Gas-Phase Species

| Species | Role | Notes |
|---------|------|-------|
| SO₂(g) | Primary reactant | Dissolves, oxidized to S(VI) |
| H₂O₂(g) | Oxidant | Dissolves, reacts with S(IV) |
| O₃(g) | Oxidant | Dissolves, reacts with S(IV) |
| HO₂(g) | Radical | Dissolves, self-reacts → H₂O₂(aq) |
| NH₃(g) | Base | Dissolves, protonates |
| HNO₃(g) | Acid | Dissolves, speciation affects pH |
| CO₂(g) | Acid | Dissolves, speciation affects pH |

### 2.2 Aqueous-Phase Species

| Species | Formula | Role |
|---------|---------|------|
| SO₂·H₂O | SO₂(aq) | Dissolved SO₂ |
| HSO₃⁻ | HSO₃⁻ | Bisulfite (reacts with H₂O₂ and O₃) |
| SO₃²⁻ | SO₃²⁻ | Sulfite (reacts with O₃) |
| SO₄²⁻ | SO₄²⁻ | Sulfate product |
| H₂O₂(aq) | H₂O₂(aq) | Dissolved H₂O₂ |
| HO₂⁻ | HO₂⁻ | Hydroperoxide anion |
| O₃(aq) | O₃(aq) | Dissolved O₃ |
| HO₂(aq) | HO₂(aq) | Dissolved HO₂ |
| O₂⁻ | O₂⁻ | Superoxide radical anion |
| H⁺ | H⁺ | Hydronium |
| OH⁻ | OH⁻ | Hydroxide |
| NH₃·H₂O | NH₃(aq) | Dissolved ammonia |
| NH₄⁺ | NH₄⁺ | Ammonium |
| HNO₃(aq) | HNO₃(aq) | Dissolved nitric acid |
| NO₃⁻ | NO₃⁻ | Nitrate |
| CO₂·H₂O | CO₂(aq) | Dissolved CO₂ |
| HCO₃⁻ | HCO₃⁻ | Bicarbonate |
| H₂O | H₂O(l) | Solvent |

### 2.3 Henry's Law Partitioning Parameters

All Henry's Law constants follow the van 't Hoff form:

$$H(T) = H_{ref} \cdot \exp\!\left(-\Delta H / R \cdot (1/T - 1/T_{ref})\right)$$

In MIAM's notation: $\text{HLC}(T) = \text{HLC}_{ref} \cdot \exp\!\left(C \cdot (1/T - 1/T_0)\right)$

where $C = \Delta H / R$ (enthalpy of dissolution divided by gas constant). CAM
uses $T_{ref} = 298$ K.

**Important**: MIAM expects HLC in **mol m⁻³ Pa⁻¹**. CAM uses **M atm⁻¹**
(= mol L⁻¹ atm⁻¹). The conversion is:

$$\text{HLC}_{\text{MIAM}} = \text{HLC}_{\text{CAM}} \times \frac{1000}{101325} \approx 9.8692 \times 10^{-3} \times \text{HLC}_{\text{CAM}}$$

| Species | Type | $H_{ref}$ (M atm⁻¹) | $\Delta H/R$ (K) | Reference |
|---------|------|----------------------|-------------------|-----------|
| HNO₃ | Monoprotic acid | 2.1×10⁵ | 8700 | Sander (2015) / Schwartz & White (1981) |
| SO₂ | Diprotic acid | 1.23 | 3120 | Chameides (1984) / NBS (1965) |
| NH₃ | Base | 58.0 | 4085 | Chameides (1984) / NBS (1965) |
| CO₂ | Diprotic acid | 3.1×10⁻² | 2423 | Chameides (1984) / NBS (1965) |
| H₂O₂ | Monoprotic acid | 7.4×10⁴ | 6621 | (source uncertain) |
| O₃ | Neutral | 1.15×10⁻² | 2560 | Chameides (1984) / NBS (1965) |
| HO₂ | Special | 9.0×10³ | 0 | Chameides (1984); Schwartz (1984) |

### 2.4 Acid Dissociation / Protonation Equilibria

All equilibrium constants follow:

$$K_{eq}(T) = A \cdot \exp\!\left(C \cdot (1/T_0 - 1/T)\right)$$

which matches MIAM's `EquilibriumConstant` directly.

| # | Equilibrium | $A$ | $C$ (K) | $T_0$ (K) | Units of $A$ |
|---|-------------|-----|---------|------------|-------------|
| E1 | H₂O ⇌ H⁺ + OH⁻ | 1.0×10⁻¹⁴ | 0 | 298.15 | mol² L⁻² |
| E2 | SO₂·H₂O ⇌ HSO₃⁻ + H⁺ | 1.7×10⁻² | 2090 | 298.15 | mol L⁻¹ |
| E3 | HSO₃⁻ ⇌ SO₃²⁻ + H⁺ | 6.0×10⁻⁸ | 1120 | 298.15 | mol L⁻¹ |
| E4 | HNO₃(aq) ⇌ NO₃⁻ + H⁺ | 15.4 | 0 | 298.15 | mol L⁻¹ |
| E5 | NH₃·H₂O + H⁺ ⇌ NH₄⁺ + H₂O | * | * | 298.15 | See note |
| E6 | CO₂·H₂O ⇌ HCO₃⁻ + H⁺ | 4.3×10⁻⁷ | −913 | 298.15 | mol L⁻¹ |
| E7 | H₂O₂(aq) ⇌ HO₂⁻ + H⁺ | 2.2×10⁻¹² | −3730 | 298.15 | mol L⁻¹ |
| E8 | HO₂(aq) ⇌ O₂⁻ + H⁺ | 2.05×10⁻⁵ | 0 | 298.15 | mol L⁻¹ |

**Note on E5 (NH₃)**: CAM treats NH₃ as a base. The effective equilibrium is:
NH₃(aq) + H₂O ⇌ NH₄⁺ + OH⁻, with $K_b$ = 1.7×10⁻⁵, $C$ = −4325 K.
This is combined with water dissociation to give the protonation equilibrium.
The effective protonation constant is: $K_{prot} = K_b / K_w$.

### 2.5 Kinetic Reactions (SO₂ Oxidation)

#### R1: S(IV) + H₂O₂ → S(VI)

$$k_{H_2O_2} = 7.45 \times 10^{7} \cdot \exp\!\left(-4430 \cdot (1/T - 1/298)\right) \cdot \frac{[\text{H}^+]}{1 + 13 \cdot [\text{H}^+]}$$

The reactive species is HSO₃⁻ (bisulfite). The overall rate in the aqueous phase is:

$$\frac{d[\text{S(VI)}]}{dt} = k_{H_2O_2} \cdot [\text{HSO}_3^-] \cdot [\text{H}_2\text{O}_2(\text{aq})]$$

Reference: Hoffmann and Calvert (1985); Seinfeld and Pandis Chapter 6.

Note: The refactored CAM code uses a slightly different form:
$k = 8.0 \times 10^{4} \cdot \exp(-3650 \cdot (1/T - 1/298)) / (0.1 + [\text{H}^+])$

#### R2: S(IV) + O₃ → S(VI)

The rate depends on both HSO₃⁻ and SO₃²⁻:

$$k_{O_3} = 3.75 \times 10^{5} \cdot \exp(-5530 \cdot \delta_T) \cdot [\text{HSO}_3^-] + 1.59 \times 10^{9} \cdot \exp(-5280 \cdot \delta_T) \cdot [\text{SO}_3^{2-}]$$

where $\delta_T = 1/T - 1/298$.

Reference: Hoffmann and Calvert (1985); Seinfeld and Pandis Chapter 6.

The overall rate is:

$$\frac{d[\text{S(VI)}]}{dt} = k_{O_3} \cdot [\text{O}_3(\text{aq})]$$

#### R3: HO₂(g) → H₂O₂(aq) (Multi-step)

This pseudo-reaction combines:
1. HO₂(g) dissolves: H* = $H_{\text{HO}_2} \cdot (1 + K_{a8}/[\text{H}^+])$
2. Self-reaction and cross-reaction form H₂O₂:

$$k^* = \frac{k_{10a} + k_{10b} \cdot K_{a8}/[\text{H}^+]}{(1 + K_{a8}/[\text{H}^+])^2}$$

$$\frac{d[\text{H}_2\text{O}_2]}{dt} = k^* \cdot [\text{O}_2^{(-1)}]^2$$

Parameters: $k_{10a}$ = 8.6×10⁵ L mol⁻¹ s⁻¹, $k_{10b}$ = 1.0×10⁸ L mol⁻¹ s⁻¹,
$K_{a8}$ = 2.05×10⁻⁵ M. Reference: Schwartz (1984).

### 2.6 Charge Balance (pH Determination)

CAM solves this via bisection on pH ∈ [2, 7]:

$$[\text{H}^+] + [\text{NH}_4^+] = [\text{OH}^-] + [\text{HCO}_3^-] + [\text{NO}_3^-] + [\text{HSO}_3^-] + 2[\text{SO}_3^{2-}] + n_{\text{SO}_4}[\text{SO}_4^{2-}]$$

where $n_{\text{SO}_4}} = 2$ normally or 1 for NH₄HSO₄ aerosol mode.

In MIAM, this becomes a `LinearConstraint`.

---

## 3. Gap Analysis: MIAM Capabilities vs. CAM Needs

### 3.1 What MIAM Can Handle Today

| CAM Feature | MIAM Mechanism | Status |
|-------------|----------------|--------|
| Henry's Law equilibrium partitioning (SO₂, H₂O₂, O₃, HNO₃, NH₃, CO₂) | `HenryLawEquilibriumConstraint` + `HenrysLawConstant` | ✅ Ready |
| Acid dissociation equilibria (E1–E4, E6–E8) | `DissolvedEquilibriumConstraint` + `EquilibriumConstant` | ✅ Ready |
| Charge balance | `LinearConstraint` | ✅ Ready |
| S(IV) + H₂O₂ → S(VI) oxidation | `DissolvedReaction` with custom rate function | ✅ Ready |
| S(IV) + O₃ → S(VI) oxidation | `DissolvedReaction` with custom rate function | ✅ Ready |
| Cloud droplet representation | `SingleMomentMode` | ✅ Ready |
| Temperature-dependent constants | `HenrysLawConstant`, `EquilibriumConstant` | ✅ Ready |
| Mass conservation constraints | `LinearConstraint` | ✅ Ready |

### 3.2 Issues Requiring Attention

#### 3.2.1 Unit Conversion for HLC

MIAM's `HenrysLawConstant` is in **mol m⁻³ Pa⁻¹**. CAM's values are in **M atm⁻¹**.
The conversion factor is:

$$\text{HLC}_{\text{MIAM}} = \text{HLC}_{\text{CAM}} \times \frac{1000}{101325}$$

All parameter values in the example script below include this conversion.

#### 3.2.2 Unit Conversion for Equilibrium Constants

MIAM's `DissolvedEquilibriumConstraint` operates on concentrations in **mol m⁻³**
(the MICM state variables). CAM's equilibrium constants have units involving
mol L⁻¹. Since 1 mol L⁻¹ = 1000 mol m⁻³, the equilibrium constants need to be
adjusted for the concentration scale.

For a reaction like: A(aq) ⇌ B + H⁺, where $K_{eq}$ = [B][H⁺]/[A] in M units:

$$K_{eq,\text{MIAM}} = K_{eq,\text{CAM}} \times 1000^{n_p - n_r}$$

where $n_p$ and $n_r$ are the number of product and reactant species.

However, MIAM's `DissolvedEquilibriumConstraint` normalizes by solvent:

$$G = K_{eq} \cdot \frac{\prod[R_i]}{[S]^{n_r - 1}} - \frac{\prod[P_j]}{[S]^{n_p - 1}} = 0$$

This means the constraint depends on the solvent concentration units as well. The
unit conversion needs to be derived carefully for each equilibrium. **This should be
verified with unit tests for each individual equilibrium.**

#### 3.2.3 NH₃ Protonation Equilibrium

CAM treats NH₃ as a base: NH₃(aq) + H₂O ⇌ NH₄⁺ + OH⁻, using $K_b$ and $K_w$.
In MIAM, this can be expressed directly as a `DissolvedEquilibriumConstraint` with
H₂O as both a reactant and the solvent. The solvent normalization in the constraint
handles the activity of water.

#### 3.2.4 HO₂ Self-Reaction → H₂O₂

The HO₂ self-reaction to form H₂O₂ in the aqueous phase is more complex than a
simple dissolved reaction. It involves:
1. HO₂(g) dissolution (Henry's Law)
2. HO₂(aq) ⇌ O₂⁻ + H⁺ (acid dissociation)
3. HO₂ + HO₂ → H₂O₂ + O₂ (self-reaction)
4. HO₂ + O₂⁻ → H₂O₂ + O₂ (cross-reaction)

**Options**:
- **A (recommended)**: Model as explicit processes in MIAM: Henry's Law constraint
  for HO₂ dissolution, dissociation constraint for HO₂⇌O₂⁻+H⁺, and two
  `DissolvedReaction`s for the self- and cross-reactions.
- **B**: Compute the effective H₂O₂ production rate externally (as CAM does) and
  add it as a source term. This avoids extra species/constraints but is less general.

Option A is preferred because it allows the DAE solver to handle the coupled
equilibria consistently.

#### 3.2.5 Reaction Rate Constants with pH Dependence

The S(IV)+H₂O₂ rate constant depends explicitly on [H⁺]. In MIAM, `DissolvedReaction`
takes a rate constant function `std::function<double(const Conditions&)>`, which does
not have direct access to species concentrations.

**Options**:
- **A**: Make the H⁺-dependent factor part of the rate law by including H⁺ as a
  reactant. Since H⁺ already appears in the system, this would work if the
  `DissolvedReaction` rate law supports it. However, the rate expression
  $k \cdot [\text{H}^+] / (1 + 13[\text{H}^+])$ is not a simple mass-action form.
- **B (recommended)**: Use a `DissolvedReaction` where the reactants are
  HSO₃⁻ and H₂O₂(aq). The pH dependence in the original rate is actually converting
  from total S(IV) to HSO₃⁻ concentration. Since MIAM tracks HSO₃⁻ explicitly
  (via the S(IV) dissociation equilibrium), the rate simplifies to just the
  temperature-dependent Arrhenius factor applied to the actual HSO₃⁻ concentration.
  Similarly, the O₃ reaction has separate terms for HSO₃⁻ and SO₃²⁻.

Option B is cleaner and follows naturally from having explicit dissociation
constraints.

#### 3.2.6 Mass Conservation / Total Species Constraints

CAM works with total mixing ratios (gas + all aqueous forms) and partitions them.
In MIAM's DAE approach, we can either:
- **A**: Initialize all forms and let the solver find equilibrium
- **B (recommended)**: Add `LinearConstraint`s for total sulfur(IV), total nitrogen,
  etc. This ensures mass is conserved and helps the solver.

### 3.3 Summary: No New Process/Constraint Types Needed

All CAM cloud chemistry features can be represented using MIAM's existing 3 process
types and 3 constraint types. The key insight is that by tracking individual ionic
species (HSO₃⁻, SO₃²⁻, etc.) as explicit state variables with equilibrium
constraints, the complex pH-dependent rate expressions in CAM simplify to
straightforward mass-action kinetics.

---

## 4. Species and Phase Definitions

### 4.1 Gas Phase Species

```
SO2_g, H2O2_g, O3_g, HO2_g, NH3_g, HNO3_g, CO2_g
```

### 4.2 Aqueous Phase Species

```
SO2_aq, HSO3m, SO3mm, SO4mm,
H2O2_aq, HO2m_peroxide,
O3_aq,
HO2_aq, O2m,
NH3_aq, NH4p,
HNO3_aq, NO3m,
CO2_aq, HCO3m,
Hp, OHm,
H2O
```

(Naming: `m` = minus/anion, `p` = plus/cation, `mm` = double minus)

---

## 5. Draft MIAM Configuration Script

```cpp
#include <miam/miam.hpp>
#include <miam/processes/constants/henrys_law_constant.hpp>
#include <miam/processes/constants/equilibrium_constant.hpp>
#include <micm/CPU.hpp>

using namespace micm;
using namespace miam;

// ============================================================================
// Unit conversion constants
// ============================================================================
constexpr double M_PER_ATM_TO_MOL_M3_PA = 1000.0 / 101325.0;  // M/atm → mol/(m³·Pa)

// ============================================================================
// 1. Define species
// ============================================================================
// Gas-phase species (need molecular weight for phase transfer)
auto so2_g   = Species{ "SO2",   {{ "molecular weight [kg mol-1]", 0.064 }} };
auto h2o2_g  = Species{ "H2O2",  {{ "molecular weight [kg mol-1]", 0.034 }} };
auto o3_g    = Species{ "O3",    {{ "molecular weight [kg mol-1]", 0.048 }} };
auto ho2_g   = Species{ "HO2",   {{ "molecular weight [kg mol-1]", 0.033 }} };
auto nh3_g   = Species{ "NH3",   {{ "molecular weight [kg mol-1]", 0.017 }} };
auto hno3_g  = Species{ "HNO3",  {{ "molecular weight [kg mol-1]", 0.063 }} };
auto co2_g   = Species{ "CO2",   {{ "molecular weight [kg mol-1]", 0.044 }} };

// Aqueous-phase species (solvent needs density)
auto h2o     = Species{ "H2O",   {{ "molecular weight [kg mol-1]", 0.018 },
                                   { "density [kg m-3]", 1000.0 }} };

// Dissolved neutrals
auto so2_aq  = Species{ "SO2_aq" };
auto h2o2_aq = Species{ "H2O2_aq" };
auto o3_aq   = Species{ "O3_aq" };
auto ho2_aq  = Species{ "HO2_aq" };
auto nh3_aq  = Species{ "NH3_aq" };
auto hno3_aq = Species{ "HNO3_aq" };
auto co2_aq  = Species{ "CO2_aq" };

// Ions
auto hp      = Species{ "Hp" };       // H+
auto ohm     = Species{ "OHm" };      // OH-
auto hso3m   = Species{ "HSO3m" };    // HSO3-
auto so3mm   = Species{ "SO3mm" };    // SO3²⁻
auto so4mm   = Species{ "SO4mm" };    // SO4²⁻
auto ho2m_p  = Species{ "HO2m_p" };   // HO2- (from H2O2 dissociation)
auto o2m     = Species{ "O2m" };       // O2-  (from HO2 dissociation)
auto nh4p    = Species{ "NH4p" };      // NH4+
auto no3m    = Species{ "NO3m" };      // NO3-
auto hco3m   = Species{ "HCO3m" };    // HCO3-

// ============================================================================
// 2. Define phases
// ============================================================================
Phase gas_phase{ "GAS", { so2_g, h2o2_g, o3_g, ho2_g, nh3_g, hno3_g, co2_g } };

Phase aqueous_phase{ "AQUEOUS", {
    h2o,
    so2_aq, hso3m, so3mm, so4mm,
    h2o2_aq, ho2m_p,
    o3_aq,
    ho2_aq, o2m,
    nh3_aq, nh4p,
    hno3_aq, no3m,
    co2_aq, hco3m,
    hp, ohm
} };

// ============================================================================
// 3. Define cloud droplet representation
// ============================================================================
// Typical cloud droplet: geometric mean radius ~10 μm, GSD ~1.4
auto cloud = representation::SingleMomentMode{
    "CLOUD", { aqueous_phase }, 10.0e-6, 1.4 };

// ============================================================================
// 4. Henry's Law equilibrium constraints (gas ⇌ aqueous)
// ============================================================================
// Each constraint replaces the ODE for the condensed-phase species with:
//   HLC * R * T * f_v * [A_g] - [A_aq] = 0

auto hl_so2 = constraint::HenryLawEquilibriumConstraintBuilder()
    .SetGasSpecies(so2_g)
    .SetCondensedSpecies(so2_aq)
    .SetSolvent(h2o)
    .SetCondensedPhase(aqueous_phase)
    .SetHenryLawConstant(process::constant::HenrysLawConstant({
        .HLC_ref_ = 1.23 * M_PER_ATM_TO_MOL_M3_PA,   // SO2: 1.23 M/atm
        .C_ = 3120.0 }))                               // ΔH/R = 3120 K
    .SetMwSolvent(0.018)
    .SetRhoSolvent(1000.0)
    .Build();

auto hl_h2o2 = constraint::HenryLawEquilibriumConstraintBuilder()
    .SetGasSpecies(h2o2_g)
    .SetCondensedSpecies(h2o2_aq)
    .SetSolvent(h2o)
    .SetCondensedPhase(aqueous_phase)
    .SetHenryLawConstant(process::constant::HenrysLawConstant({
        .HLC_ref_ = 7.4e4 * M_PER_ATM_TO_MOL_M3_PA,  // H2O2: 7.4e4 M/atm
        .C_ = 6621.0 }))
    .SetMwSolvent(0.018)
    .SetRhoSolvent(1000.0)
    .Build();

auto hl_o3 = constraint::HenryLawEquilibriumConstraintBuilder()
    .SetGasSpecies(o3_g)
    .SetCondensedSpecies(o3_aq)
    .SetSolvent(h2o)
    .SetCondensedPhase(aqueous_phase)
    .SetHenryLawConstant(process::constant::HenrysLawConstant({
        .HLC_ref_ = 1.15e-2 * M_PER_ATM_TO_MOL_M3_PA, // O3: 1.15e-2 M/atm
        .C_ = 2560.0 }))
    .SetMwSolvent(0.018)
    .SetRhoSolvent(1000.0)
    .Build();

auto hl_ho2 = constraint::HenryLawEquilibriumConstraintBuilder()
    .SetGasSpecies(ho2_g)
    .SetCondensedSpecies(ho2_aq)
    .SetSolvent(h2o)
    .SetCondensedPhase(aqueous_phase)
    .SetHenryLawConstant(process::constant::HenrysLawConstant({
        .HLC_ref_ = 9.0e3 * M_PER_ATM_TO_MOL_M3_PA,  // HO2: 9.0e3 M/atm
        .C_ = 0.0 }))                                  // no T dependence
    .SetMwSolvent(0.018)
    .SetRhoSolvent(1000.0)
    .Build();

auto hl_nh3 = constraint::HenryLawEquilibriumConstraintBuilder()
    .SetGasSpecies(nh3_g)
    .SetCondensedSpecies(nh3_aq)
    .SetSolvent(h2o)
    .SetCondensedPhase(aqueous_phase)
    .SetHenryLawConstant(process::constant::HenrysLawConstant({
        .HLC_ref_ = 58.0 * M_PER_ATM_TO_MOL_M3_PA,   // NH3: 58 M/atm
        .C_ = 4085.0 }))
    .SetMwSolvent(0.018)
    .SetRhoSolvent(1000.0)
    .Build();

auto hl_hno3 = constraint::HenryLawEquilibriumConstraintBuilder()
    .SetGasSpecies(hno3_g)
    .SetCondensedSpecies(hno3_aq)
    .SetSolvent(h2o)
    .SetCondensedPhase(aqueous_phase)
    .SetHenryLawConstant(process::constant::HenrysLawConstant({
        .HLC_ref_ = 2.1e5 * M_PER_ATM_TO_MOL_M3_PA,  // HNO3: 2.1e5 M/atm
        .C_ = 8700.0 }))
    .SetMwSolvent(0.018)
    .SetRhoSolvent(1000.0)
    .Build();

auto hl_co2 = constraint::HenryLawEquilibriumConstraintBuilder()
    .SetGasSpecies(co2_g)
    .SetCondensedSpecies(co2_aq)
    .SetSolvent(h2o)
    .SetCondensedPhase(aqueous_phase)
    .SetHenryLawConstant(process::constant::HenrysLawConstant({
        .HLC_ref_ = 3.1e-2 * M_PER_ATM_TO_MOL_M3_PA, // CO2: 3.1e-2 M/atm
        .C_ = 2423.0 }))
    .SetMwSolvent(0.018)
    .SetRhoSolvent(1000.0)
    .Build();

// ============================================================================
// 5. Acid-base dissociation equilibrium constraints
// ============================================================================
// Each replaces the ODE row for the algebraic species (typically the anion)
// with: K_eq * prod[R] / [S]^(nr-1) - prod[P] / [S]^(np-1) = 0
//
// NOTE: K_eq values from CAM are in mol/L units. Since MIAM state variables
// are in mol/m³, the A parameter for EquilibriumConstant needs adjustment.
// The solvent normalization in the constraint equation handles this partially.
// Exact conversion depends on the constraint's normalization — verify with
// unit tests for each equilibrium.

// E1: H2O ⇌ H+ + OH-
// K_w = [H+][OH-] = 1.0e-14 mol² L⁻²
auto eq_h2o = constraint::DissolvedEquilibriumConstraintBuilder()
    .SetPhase(aqueous_phase)
    .SetReactants({ h2o })
    .SetProducts({ hp, ohm })
    .SetAlgebraicSpecies(ohm)
    .SetSolvent(h2o)
    .SetEquilibriumConstant(process::constant::EquilibriumConstant({
        .A_ = 1.0e-14,  // mol² L⁻² (needs unit conversion verification)
        .C_ = 0.0 }))
    .Build();

// E2: SO2·H2O ⇌ HSO3- + H+
// Ka1 = [HSO3-][H+] / [SO2·H2O] = 1.7e-2 M
auto eq_so2_ka1 = constraint::DissolvedEquilibriumConstraintBuilder()
    .SetPhase(aqueous_phase)
    .SetReactants({ so2_aq })
    .SetProducts({ hso3m, hp })
    .SetAlgebraicSpecies(hso3m)
    .SetSolvent(h2o)
    .SetEquilibriumConstant(process::constant::EquilibriumConstant({
        .A_ = 1.7e-2,   // Ka1 for SO2
        .C_ = 2090.0 }))
    .Build();

// E3: HSO3- ⇌ SO3²⁻ + H+
// Ka2 = [SO3²⁻][H+] / [HSO3-] = 6.0e-8 M
auto eq_so2_ka2 = constraint::DissolvedEquilibriumConstraintBuilder()
    .SetPhase(aqueous_phase)
    .SetReactants({ hso3m })
    .SetProducts({ so3mm, hp })
    .SetAlgebraicSpecies(so3mm)
    .SetSolvent(h2o)
    .SetEquilibriumConstant(process::constant::EquilibriumConstant({
        .A_ = 6.0e-8,
        .C_ = 1120.0 }))
    .Build();

// E4: HNO3(aq) ⇌ NO3- + H+
// Ka = [NO3-][H+] / [HNO3(aq)] = 15.4 M
auto eq_hno3 = constraint::DissolvedEquilibriumConstraintBuilder()
    .SetPhase(aqueous_phase)
    .SetReactants({ hno3_aq })
    .SetProducts({ no3m, hp })
    .SetAlgebraicSpecies(no3m)
    .SetSolvent(h2o)
    .SetEquilibriumConstant(process::constant::EquilibriumConstant({
        .A_ = 15.4,
        .C_ = 0.0 }))
    .Build();

// E5: NH3(aq) + H2O ⇌ NH4+ + OH-
// Kb = [NH4+][OH-] / [NH3(aq)] = 1.7e-5 M
auto eq_nh3 = constraint::DissolvedEquilibriumConstraintBuilder()
    .SetPhase(aqueous_phase)
    .SetReactants({ nh3_aq, h2o })
    .SetProducts({ nh4p, ohm })
    .SetAlgebraicSpecies(nh4p)
    .SetSolvent(h2o)
    .SetEquilibriumConstant(process::constant::EquilibriumConstant({
        .A_ = 1.7e-5,
        .C_ = -4325.0 }))
    .Build();

// E6: CO2(aq) ⇌ HCO3- + H+
// Ka1 = [HCO3-][H+] / [CO2(aq)] = 4.3e-7 M
auto eq_co2 = constraint::DissolvedEquilibriumConstraintBuilder()
    .SetPhase(aqueous_phase)
    .SetReactants({ co2_aq })
    .SetProducts({ hco3m, hp })
    .SetAlgebraicSpecies(hco3m)
    .SetSolvent(h2o)
    .SetEquilibriumConstant(process::constant::EquilibriumConstant({
        .A_ = 4.3e-7,
        .C_ = -913.0 }))
    .Build();

// E7: H2O2(aq) ⇌ HO2- + H+
// Ka = [HO2-][H+] / [H2O2(aq)] = 2.2e-12 M
auto eq_h2o2 = constraint::DissolvedEquilibriumConstraintBuilder()
    .SetPhase(aqueous_phase)
    .SetReactants({ h2o2_aq })
    .SetProducts({ ho2m_p, hp })
    .SetAlgebraicSpecies(ho2m_p)
    .SetSolvent(h2o)
    .SetEquilibriumConstant(process::constant::EquilibriumConstant({
        .A_ = 2.2e-12,
        .C_ = -3730.0 }))
    .Build();

// E8: HO2(aq) ⇌ O2- + H+
// Ka = [O2-][H+] / [HO2(aq)] = 2.05e-5 M
auto eq_ho2 = constraint::DissolvedEquilibriumConstraintBuilder()
    .SetPhase(aqueous_phase)
    .SetReactants({ ho2_aq })
    .SetProducts({ o2m, hp })
    .SetAlgebraicSpecies(o2m)
    .SetSolvent(h2o)
    .SetEquilibriumConstant(process::constant::EquilibriumConstant({
        .A_ = 2.05e-5,
        .C_ = 0.0 }))
    .Build();

// ============================================================================
// 6. Kinetic reactions (S(IV) oxidation)
// ============================================================================

// R1: HSO3- + H2O2(aq) → SO4²⁻ + H2O + H+
// The rate constant from Hoffmann & Calvert (1985) applied to HSO3- directly
// k = 7.45e7 * exp(-4430*(1/T - 1/298))
// (The [H+]/(1+13[H+]) factor in CAM converts total S(IV) to HSO3-;
//  since we track HSO3- explicitly, we use just the Arrhenius part.)
auto rxn_siv_h2o2 = process::DissolvedReactionBuilder()
    .SetPhase(aqueous_phase)
    .SetReactants({ hso3m, h2o2_aq })
    .SetProducts({ so4mm, h2o, hp })
    .SetSolvent(h2o)
    .SetRateConstant([](const Conditions& c) -> double {
        return 7.45e7 * std::exp(-4430.0 * (1.0 / c.temperature_ - 1.0 / 298.0));
    })
    .Build();

// R2a: HSO3- + O3(aq) → SO4²⁻ + H+ (+ O2 implicit)
// k = 3.75e5 * exp(-5530*(1/T - 1/298))
auto rxn_hso3_o3 = process::DissolvedReactionBuilder()
    .SetPhase(aqueous_phase)
    .SetReactants({ hso3m, o3_aq })
    .SetProducts({ so4mm, hp })
    .SetSolvent(h2o)
    .SetRateConstant([](const Conditions& c) -> double {
        return 3.75e5 * std::exp(-5530.0 * (1.0 / c.temperature_ - 1.0 / 298.0));
    })
    .Build();

// R2b: SO3²⁻ + O3(aq) → SO4²⁻ (+ O2 implicit)
// k = 1.59e9 * exp(-5280*(1/T - 1/298))
auto rxn_so3_o3 = process::DissolvedReactionBuilder()
    .SetPhase(aqueous_phase)
    .SetReactants({ so3mm, o3_aq })
    .SetProducts({ so4mm })
    .SetSolvent(h2o)
    .SetRateConstant([](const Conditions& c) -> double {
        return 1.59e9 * std::exp(-5280.0 * (1.0 / c.temperature_ - 1.0 / 298.0));
    })
    .Build();

// R3a: HO2(aq) + HO2(aq) → H2O2(aq) + O2 (implicit O2)
// k = 8.6e5 L mol⁻¹ s⁻¹
auto rxn_ho2_ho2 = process::DissolvedReactionBuilder()
    .SetPhase(aqueous_phase)
    .SetReactants({ ho2_aq, ho2_aq })
    .SetProducts({ h2o2_aq })
    .SetSolvent(h2o)
    .SetRateConstant([](const Conditions&) -> double { return 8.6e5; })
    .Build();

// R3b: HO2(aq) + O2- → H2O2(aq) + O2 (implicit O2)
// k = 1.0e8 L mol⁻¹ s⁻¹
auto rxn_ho2_o2m = process::DissolvedReactionBuilder()
    .SetPhase(aqueous_phase)
    .SetReactants({ ho2_aq, o2m })
    .SetProducts({ h2o2_aq })
    .SetSolvent(h2o)
    .SetRateConstant([](const Conditions&) -> double { return 1.0e8; })
    .Build();

// ============================================================================
// 7. Charge balance (linear constraint)
// ============================================================================
// [H+] + [NH4+] = [OH-] + [HSO3-] + 2[SO3²⁻] + 2[SO4²⁻] + [NO3-]
//                + [HCO3-] + [HO2-] + [O2-]
//
// Rearranged to: [H+] + [NH4+] - [OH-] - [HSO3-] - 2[SO3²⁻] - 2[SO4²⁻]
//              - [NO3-] - [HCO3-] - [HO2-] - [O2-] = 0
//
// The algebraic variable is H+ (its ODE row is replaced by this constraint).
auto charge_balance = constraint::LinearConstraintBuilder()
    .SetAlgebraicSpecies(aqueous_phase, hp)
    .AddTerm(aqueous_phase, hp,     1.0)
    .AddTerm(aqueous_phase, nh4p,   1.0)
    .AddTerm(aqueous_phase, ohm,   -1.0)
    .AddTerm(aqueous_phase, hso3m, -1.0)
    .AddTerm(aqueous_phase, so3mm, -2.0)
    .AddTerm(aqueous_phase, so4mm, -2.0)
    .AddTerm(aqueous_phase, no3m,  -1.0)
    .AddTerm(aqueous_phase, hco3m, -1.0)
    .AddTerm(aqueous_phase, ho2m_p,-1.0)
    .AddTerm(aqueous_phase, o2m,   -1.0)
    .SetConstant(0.0)
    .Build();

// ============================================================================
// 8. Assemble model, system, and solver
// ============================================================================
auto cloud_model = Model{ .name_ = "CLOUD", .representations_ = { cloud } };

cloud_model.AddProcesses({
    rxn_siv_h2o2,
    rxn_hso3_o3,
    rxn_so3_o3,
    rxn_ho2_ho2,
    rxn_ho2_o2m
});

cloud_model.AddConstraints({
    hl_so2, hl_h2o2, hl_o3, hl_ho2, hl_nh3, hl_hno3, hl_co2,
    eq_h2o, eq_so2_ka1, eq_so2_ka2, eq_hno3, eq_nh3, eq_co2, eq_h2o2, eq_ho2,
    charge_balance
});

auto system = System(gas_phase, cloud_model);

auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(
                  RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                  .SetSystem(system)
                  .AddExternalModel(cloud_model)
                  .SetIgnoreUnusedSpecies(true)
                  .Build();

// ============================================================================
// 9. Set initial conditions and solve
// ============================================================================
State state = solver.GetState();

// Environmental conditions (representative mid-troposphere cloud)
state.conditions_[0].temperature_ = 280.0;  // K
state.conditions_[0].pressure_ = 70000.0;   // Pa (~700 hPa)
state.conditions_[0].CalculateIdealAirDensity();

// Gas-phase mixing ratios converted to mol/m³:
// [X] = vmr * air_density (molecules/m³) / Avogadro
double air_density = state.conditions_[0].air_density_;  // mol/m³
state[so2_g]  = 1.0e-9 * air_density;   // ~1 ppb SO2
state[h2o2_g] = 1.0e-9 * air_density;   // ~1 ppb H2O2
state[o3_g]   = 50.0e-9 * air_density;  // ~50 ppb O3
state[ho2_g]  = 10.0e-12 * air_density; // ~10 ppt HO2
state[nh3_g]  = 1.0e-9 * air_density;   // ~1 ppb NH3
state[hno3_g] = 0.5e-9 * air_density;   // ~0.5 ppb HNO3
state[co2_g]  = 420.0e-6 * air_density; // ~420 ppm CO2

// Cloud droplet water content
// Use ~0.3 g/m³ cloud LWC → 300 mol/m³ for liquid water
state[cloud.Species(aqueous_phase, h2o)] = 300.0;

// Initial aqueous concentrations (small seed values; solver will equilibrate)
state[cloud.Species(aqueous_phase, hp)]     = 1.0e-5 * 1000.0;  // pH ~5, in mol/m³
state[cloud.Species(aqueous_phase, ohm)]    = 1.0e-9 * 1000.0;
state[cloud.Species(aqueous_phase, so4mm)]  = 1.0e-6 * 1000.0;  // initial sulfate

// Set default representation parameters (radius, GSD, number conc)
cloud.SetDefaultParameters(state);

// Solve
solver.UpdateStateParameters(state);
double dt = 1800.0;  // 30-minute timestep (typical CAM physics timestep)
auto result = solver.Solve(dt, state);
```

---

## 6. Open Questions and Next Steps

1. **Unit conversion verification**: The equilibrium constant unit conversions
   (mol/L → mol/m³) need careful verification through unit tests for each individual
   equilibrium. The solvent normalization in `DissolvedEquilibriumConstraint` affects
   the required units.

2. **Solver convergence**: The system has 7 algebraic constraints (Henry's Law) +
   8 algebraic constraints (dissociation) + 1 algebraic constraint (charge balance) =
   16 algebraic constraints total, plus 5 kinetic ODE terms. The Rosenbrock solver
   should handle this, but convergence may need tuning (tolerances, initial guesses).

3. **Rate constant forms**: Need to verify that the simplified Arrhenius form for
   S(IV)+H₂O₂ (without the [H+]/(1+13[H+]) factor) is correct when HSO₃⁻ is tracked
   explicitly. The original rate expression from Hoffmann and Calvert may have
   additional nuance.

4. **Cloud droplet number concentration**: `SingleMomentMode` requires specifying
   the geometric mean radius and GSD. The actual number concentration is a state
   parameter. Typical values: N ~ 100-300 cm⁻³, r_g ~ 5-15 μm, σ_g ~ 1.2-1.5.

5. **Mass balance constraints**: Consider adding `LinearConstraint`s for total
   sulfur(IV) conservation, total nitrogen (HNO3+NO3), and total ammonia (NH3+NH4+)
   to help solver convergence and ensure mass conservation.

6. **CAM-specific bugs**: Several bugs documented in the user's notes
   (unit conversion errors, inconsistent Heff usage, pressure assumptions) are
   inherently avoided by MIAM's consistent DAE formulation.

7. **MSA (methanesulfonic acid)**: Present in CAM but not modeled here. Can be
   added later as an additional species if needed.
