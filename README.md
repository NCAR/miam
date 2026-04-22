MIAM
==============

Model-Independent Aerosol Module

[![License](https://img.shields.io/github/license/NCAR/miam.svg)](https://github.com/NCAR/miam/blob/main/LICENSE)

Copyright (C) 2026 University Corporation for Atmospheric Research

> **Note**
> MIAM is currently under development. Expect incomplete and unverified
> functionality across minor revision numbers.

MIAM is a C++20 library for representing and evolving aerosol and cloud
particle populations within atmospheric chemistry models. It is designed to
work alongside [MICM](https://github.com/NCAR/micm) as part of the
[MUSICA](https://github.com/NCAR/musica) (Multi-Scale Infrastructure for
Chemistry and Aerosols) project.

MIAM provides:
- **Representations** — describe particle size distributions
  (single-moment modes, two-moment modes, uniform sections)
- **Processes** — describe condensed-phase chemistry and gas-particle
  mass transfer (dissolved reversible reactions, Henry's Law phase transfer)
- **A common process interface** — add new process types without modifying
  the model core
- **Aerosol property providers** — representation-specific calculations of
  effective radius, number concentration, and phase volume fraction with
  analytical partial derivatives for implicit solvers

## Units and Conventions

All concentrations in MIAM — both gas-phase and condensed-phase — are in units
of **mol m⁻³ of air**. This includes the solvent (e.g. cloud liquid water).

| Parameter | MIAM units | Notes |
|-----------|-----------|-------|
| Gas-phase concentration | mol m⁻³ air | |
| Condensed-phase concentration | mol m⁻³ air | Not mol/L of solution |
| Solvent concentration [S] | mol m⁻³ air | Cloud water: ≈ 0.017 mol m⁻³ for LWC = 0.3 g m⁻³ |
| Rate constant *k* | s⁻¹ | Always s⁻¹ after unit conversion (see below) |
| Equilibrium constant *K*_eq | dimensionless | Always dimensionless after unit conversion |

### Dissolved-phase rate expression

For a dissolved reaction with *n*_r reactants, MIAM computes the rate as:

> *r* = *k* / [S]^(*n*_r − 1) × ∏[R_i]

The solvent-normalization factor [S]^(*n*_r − 1) absorbs the concentration
dimensions, so *k* is always in units of s⁻¹ regardless of reaction order.

### Converting from literature (molar) values

Literature rate constants and equilibrium constants are typically expressed in
molar (mol L⁻¹) units. To convert to MIAM's *k* and *K*_eq:

> *k*_MIAM = *k*_lit × *c*_H₂O^(*n*_r − 1)
>
> *K*_eq,MIAM = *K*_lit / *c*_H₂O^(*n*_p − *n*_r)

where *c*_H₂O = 55.556 mol L⁻¹ is the molar concentration of pure liquid
water (a **fixed** unit-conversion constant, not the cloud water state
variable), and *K*_lit includes all species — including water — in molar
concentration units.

Examples:
- Unimolecular (*n*_r = 1): *k*_MIAM = *k*_lit
- Bimolecular (*n*_r = 2): *k*_MIAM = *k*_lit × 55.556
- Dissociation A → B + C (*n*_r = 1, *n*_p = 2): *K*_eq = *K*_lit / 55.556

# Getting Started

## Installing MIAM locally

MIAM is a header-only library. To build and run the tests, you need a
C++20-capable compiler and CMake ≥ 3.21. MICM and Google Test are fetched
automatically via CMake's `FetchContent`.

```bash
git clone https://github.com/NCAR/miam.git
cd miam
mkdir build && cd build
cmake ..
make -j8
make test
```

If you would later like to uninstall MIAM, you can run `sudo make uninstall`
from the `build/` directory.

## Running a MIAM Docker container

You must have [Docker](https://www.docker.com/get-started) (or
[Podman](https://podman.io/)) installed and running.

```bash
git clone https://github.com/NCAR/miam.git
cd miam
docker build -t miam -f docker/Dockerfile .
docker run -it miam bash
```

Inside the container, run the tests from the `/build/` directory:

```bash
cd /build/
make test
```

# Using the MIAM API

The following example sets up a cloud droplet system and models CO₂
dissolving from the gas phase into aqueous droplets via Henry's Law phase
transfer:

```c++
#include <miam/miam.hpp>
#include <miam/processes/constants/henrys_law_constant.hpp>
#include <micm/CPU.hpp>

#include <iomanip>
#include <iostream>

using namespace micm;
using namespace miam;

int main()
{
  // Define species with physical properties required for mass transfer
  auto co2 = Species{ "CO2",
      { { "molecular weight [kg mol-1]", 0.044 },
        { "density [kg m-3]", 1800.0 } } };
  auto h2o = Species{ "H2O",
      { { "molecular weight [kg mol-1]", 0.018 },
        { "density [kg m-3]", 1000.0 } } };

  // Define phases
  Phase gas_phase{ "GAS", { { co2 } } };
  Phase aqueous_phase{ "AQUEOUS", { { co2 }, { h2o } } };

  // Cloud droplets with a single-moment log-normal distribution
  auto cloud = representation::SingleMomentMode{
    "CLOUD",
    { aqueous_phase },
    5.0e-6,  // geometric mean radius [m]
    1.2      // geometric standard deviation
  };

  // Henry's Law phase transfer: CO2(g) <-> CO2(aq)
  auto co2_transfer = process::HenryLawPhaseTransferBuilder()
    .SetCondensedPhase(aqueous_phase)
    .SetGasSpecies(co2)
    .SetCondensedSpecies(co2)
    .SetSolvent(h2o)
    .SetHenrysLawConstant(process::constant::HenrysLawConstant(
        { .HLC_ref_ = 3.4e-2 }))       // mol m-3 Pa-1 at 298 K
    .SetDiffusionCoefficient(1.5e-5)    // m2 s-1
    .SetAccommodationCoefficient(5.0e-6)
    .Build();

  // Create model and add processes
  auto cloud_model = Model{
    .name_ = "CLOUD",
    .representations_ = { cloud }
  };
  cloud_model.AddProcesses({ co2_transfer });

  // Build solver (MICM Rosenbrock)
  auto system = System(gas_phase, cloud_model);
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(
                    RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                    .SetSystem(system)
                    .AddExternalModel(cloud_model)
                    .SetIgnoreUnusedSpecies(true)
                    .Build();

  State state = solver.GetState();

  // Initial conditions
  state.conditions_[0].temperature_ = 298.15;   // K
  state.conditions_[0].pressure_ = 101325.0;    // Pa
  state.conditions_[0].CalculateIdealAirDensity();

  state[co2] = 1.0e-3;                                  // mol m-3 air
  state[cloud.Species(aqueous_phase, h2o)] = 0.017;      // mol m-3 air (cloud LWC ~ 0.3 g m-3)
  cloud.SetDefaultParameters(state);

  // Integrate
  state.PrintHeader();
  state.PrintState(0);
  for (int i = 1; i <= 10; ++i)
  {
    solver.CalculateRateConstants(state);
    auto result = solver.Solve(0.1, state);  // 100 ms steps
    state.PrintState(i * 100);
  }
}
```

To build and run (assuming the default install location):

```bash
g++ -o cloud_chem cloud_chem.cpp -I/usr/local/miam/include -I/usr/local/micm/include -std=c++20
./cloud_chem
```

Expected output — gas-phase CO₂ dissolves into the cloud droplets, approaching
Henry's Law equilibrium. With realistic cloud liquid water content
(~0.017 mol m⁻³), only a trace fraction of the gas dissolves:

```
  time,                CO2,  CLOUD.AQUEOUS.CO2,  CLOUD.AQUEOUS.H2O
     0,           1.00e-03,           0.00e+00,           1.70e-02
   ...
  1000,           ~1.00e-03,          ~1.0e-08,           1.70e-02
```

See the [MIAM documentation](https://miam.readthedocs.io/) for the full API
reference, science guide, and additional examples including dissolved reversible
reactions and combined gas–aqueous chemistry.

# Citation

MIAM is part of the [MUSICA](https://github.com/NCAR/musica) project and can
be cited by reference to the MUSICA vision paper:

```bibtex
@Article{acom.software.musica-vision,
    author = "Gabriele G. Pfister and Sebastian D. Eastham and Avelino F. Arellano and
              Bernard Aumont and Kelley C. Barsanti and Mary C. Barth and Andrew Conley and
              Nicholas A. Davis and Louisa K. Emmons and Jerome D. Fast and Arlene M. Fiore and
              Benjamin Gaubert and Steve Goldhaber and Claire Granier and Georg A. Grell and
              Marc Guevara and Daven K. Henze and Alma Hodzic and Xiaohong Liu and
              Daniel R. Marsh and John J. Orlando and John M. C. Plane and
              Lorenzo M. Polvani and Karen H. Rosenlof and Allison L. Steiner and
              Daniel J. Jacob and Guy P. Brasseur",
    title = "The Multi-Scale Infrastructure for Chemistry and Aerosols (MUSICA)",
    journal = "Bulletin of the American Meteorological Society",
    year = "2020",
    publisher = "American Meteorological Society",
    address = "Boston MA, USA",
    volume = "101",
    number = "10",
    doi = "10.1175/BAMS-D-19-0331.1",
    pages = "E1743 - E1760",
    url = "https://journals.ametsoc.org/view/journals/bams/101/10/bamsD190331.xml"
}
```

# Community and Contributions

We welcome contributions and feedback from anyone, everything from updating
the content or appearance of the documentation to new and cutting-edge science.

- [Collaboration](https://github.com/NCAR/musica/blob/main/docs/Software%20Development%20Plan.pdf)
  — Anyone interested in scientific collaboration which would add new software
  functionality should read the
  [MUSICA software development plan](https://github.com/NCAR/musica/blob/main/docs/Software%20Development%20Plan.pdf).

# Documentation

Please see the [MIAM documentation](https://miam.readthedocs.io/) for detailed
installation and usage instructions, the science guide, and the full API
reference.

# License

- [Apache 2.0](/LICENSE)

Copyright (C) 2026 University Corporation for Atmospheric Research