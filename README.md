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

The following example sets up a simple cloud droplet system with CO₂
dissolving and undergoing aqueous-phase carbonate chemistry:

```c++
#include <miam/miam.hpp>
#include <micm/CPU.hpp>

#include <iomanip>
#include <iostream>

using namespace micm;
using namespace miam;

int main()
{
  // Define species
  auto co2   = Species{ "CO2" };
  auto h2o   = Species{ "H2O" };
  auto ohm   = Species{ "OH-" };
  auto hp    = Species{ "H+" };
  auto hco3m = Species{ "HCO3-" };
  auto co32m = Species{ "CO32-" };
  auto h2co3 = Species{ "H2CO3" };

  // Define phases
  Phase gas_phase{ "GAS", { { co2, 31.2 } } };
  Phase aqueous_phase{
    "AQUEOUS", { { co2 }, { h2o }, { ohm }, { hp }, { hco3m }, { co32m }, { h2co3 } }
  };

  // Cloud droplets with a single-moment log-normal distribution
  auto cloud_droplets = representation::SingleMomentMode{
    "CLOUD",
    { aqueous_phase },
    1.0e-6,  // geometric mean radius [m]
    1.4      // geometric standard deviation
  };

  // Aqueous equilibrium: H2O <-> OH- + H+
  auto h2o_dissociation = process::DissolvedReversibleReactionBuilder()
    .SetPhase(aqueous_phase)
    .SetReactants({ h2o })
    .SetProducts({ ohm, hp })
    .SetSolvent(h2o)
    .SetEquilibriumConstant(process::constant::EquilibriumConstant(
        { .A_ = 1.14e-2, .C_ = 2300.0, .T0_ = 298.15 }))
    .SetReverseRateConstant(ArrheniusRateConstant(
        { .A_ = 1.4e11, .C_ = 5.1e4 }))
    .Build();

  // Aqueous equilibrium: CO2 + H2O <-> H2CO3
  auto co2_hydration = process::DissolvedReversibleReactionBuilder()
    .SetPhase(aqueous_phase)
    .SetReactants({ co2, h2o })
    .SetProducts({ h2co3 })
    .SetSolvent(h2o)
    .SetEquilibriumConstant(process::constant::EquilibriumConstant(
        { .A_ = 1.70e3, .C_ = 2400.0, .T0_ = 298.15 }))
    .SetReverseRateConstant(ArrheniusRateConstant(
        { .A_ = 1.4e11, .C_ = 5.1e4 }))
    .Build();

  // Create model and add processes
  auto cloud_model = Model{
    .name_ = "CLOUD",
    .representations_ = { cloud_droplets }
  };
  cloud_model.AddProcesses({ h2o_dissociation, co2_hydration });

  // Build solver (MICM Rosenbrock)
  auto system = System(gas_phase, cloud_model);
  auto solver = CpuSolverBuilder<RosenbrockSolverParameters>(
                    RosenbrockSolverParameters::ThreeStageRosenbrockParameters())
                    .SetSystem(system)
                    .AddExternalModelProcesses(cloud_model)
                    .SetIgnoreUnusedSpecies(true)
                    .Build();

  State state = solver.GetState();

  // Initial conditions
  state.conditions_[0].temperature_ = 287.45;  // K
  state.conditions_[0].pressure_ = 101319.9;   // Pa
  state.conditions_[0].CalculateIdealAirDensity();

  state[co2] = 0.2;   // mol m-3
  state[cloud_droplets.Species(aqueous_phase, h2o)] = 55.5;  // mol m-3
  cloud_droplets.SetDefaultParameters(state);

  // Integrate
  state.PrintHeader();
  for (int i = 0; i < 10; ++i)
  {
    solver.CalculateRateConstants(state);
    auto result = solver.Solve(500.0, state);
    state.PrintState(i * 500);
  }
}
```

To build and run (assuming the default install location):

```bash
g++ -o cloud_chem cloud_chem.cpp -I/usr/local/miam/include -I/usr/local/micm/include -std=c++20
./cloud_chem
```

See the [MIAM documentation](https://ncar.github.io/miam/) for the full API
reference, science guide, and additional examples including Henry's Law phase
transfer.

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

Please see the [MIAM documentation](https://ncar.github.io/miam/) for detailed
installation and usage instructions, the science guide, and the full API
reference.

# License

- [Apache 2.0](/LICENSE)

Copyright (C) 2026 University Corporation for Atmospheric Research