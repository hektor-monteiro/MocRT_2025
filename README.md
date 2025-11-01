# MOCASSIN

## Description

MOCASSIN (MOnte CArlo SimulationS of Ionized Nebulae) is a fully 3D or 2D photoionisation and dust radiative transfer code. It employs a Monte Carlo approach to the transfer of radiation through media of arbitrary geometry and density distribution. It was originally developed for the modelling of photoionised regions like HII regions and planetary nebulae and has since expanded and been applied to a variety of astrophysical problems.

## Features

*   3D and 2D photoionisation and dust radiative transfer simulations.
*   Monte Carlo approach for radiation transfer.
*   Supports arbitrary Cartesian grids of variable resolution.
*   Handles complex density fields from SPH calculations.
*   Deals with ionising radiation extending from the Lyman edge to the X-ray.
*   Fully coupled dust and gas microphysics in both radiation transfer and thermal balance.

## System Requirements

*   **Fortran 90 Compiler:** `ifort`, `gfortran`, or `g95`.
*   **MPI Implementation:** An MPI implementation is required for parallel processing.
*   **Operating System:** Tested on various platforms, including Linux and macOS.

## Compilation and Execution

### Compilation

To compile the code, simply run the `make` command in the root directory:

```bash
make
```

This will create the following executables:

*   `mocassin`
*   `mocassinWarm`
*   `mocassinOutput`
*   `mocassinPlot`

### Execution

To run a simulation, you will need an input file. The `examples` directory contains some example input files. To run a simulation, you can use the following command:

```bash
mpirun -np <number_of_processors> mocassin <input_file>
```

## Dust Species

Some notes on dust species:

| Grain Species        | Type               | Size Range (μm) | ρ (g/cm³) | Sublim. Temp (K) | Molecular Formula | Mol. Weight (amu) | ⟨amu/atom⟩ | Work Function (Ry) |
| -------------------- | ------------------ | --------------- | --------- | ---------------- | ----------------- | ----------------- | ---------- | ------------------ |
| Silicates (amorphous)| Mg/Fe silicates    | 0.005–0.25      | 2.5–3.5   | 1200–1500        | Approx. MgFeSiO₄  | ~172              | ~21.5      | ~0.37 (5 eV)       |
| Graphite             | Carbonaceous       | 0.005–0.25      | 2.2       | ~2000            | C                 | 12                | 12.0       | 0.33–0.37 (4.5–5 eV)|
| Olivine              | (Mg,Fe)₂SiO₄       | 0.01–1.0        | 3.3–4.4   | 1400–1500        | Mg₂SiO₄           | 140.7             | 23.5       | ~0.37 (5 eV)       |
| Forsterite           | Mg₂SiO₄ (pure olivine) | 0.1–1.0         | 3.27      | 1400             | Mg₂SiO₄           | 140.7             | 23.5       | ~0.37 (5 eV)       |
| Quartz               | SiO₂               | 0.01–0.5        | 2.65      | 1200–1300        | SiO₂              | 60.08             | 20.0       | ~0.33 (4.5 eV)     |
| Enstatite            | MgSiO₃             | 0.1–1.0         | 3.2       | 1350–1400        | MgSiO₃            | 100.4             | 25.1       | ~0.37 (5 eV)       |

### References & Sources

Here are sources used to compile the physical parameters (need to double check these):

#### Density and Molecular Weights
- **Draine, B. T. (2003)**, *Annual Review of Astronomy and Astrophysics*, 41, 241
- **Henning, T. (2010)**, *ARA&A*, 48, 21
- **Lodders, K. (2003)**, *ApJ*, 591, 1220 (solar abundances)
- **CRC Handbook of Chemistry and Physics** for molecular weights and densities

#### Sublimation Temperatures
- **Gail & Sedlmayr (1999)**, *A&A*, 347, 594
- **Tielens, A. G. G. M. (2005)**, *The Physics and Chemistry of the Interstellar Medium*
- **Duschl et al. (1996)**, *A&A*, 312, 624
- Note: Temperatures vary with ambient pressure and gas density.

#### Work Function
- **Draine & Salpeter (1979)**, *ApJ*, 231, 77 (photoelectric emission from grains)
- **Jenkins (2009)**, *ApJ*, 700, 1299 (dust charging)
- **CRC Handbook** for graphite and silica values
- For silicates, work functions are inferred from electron yield and ionization potentials (~5 eV = 0.37 Ry)

## References

*   **Official MOCASSIN Website:** [https://mocassin.nebulousresearch.org/](https://mocassin.nebulousresearch.org/)

## License

This project is licensed under the GNU Lesser General Public License v3.0. See the [LICENSE.md](LICENSE.md) file for details.
