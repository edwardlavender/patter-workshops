Workshop materials for the
[`patter`](https://github.com/edwardlavender/patter)
[R](https://www.r-project.org/) package
================
Edward Lavender\*

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

<sup>\*</sup>This repository is maintained by Edward Lavender
(<edward.lavender@eawag.ch>).

# Introduction

This repository contains
[`patter`](https://github.com/edwardlavender/patter) workshop materials.
[`patter`](https://github.com/edwardlavender/patter) is an
[R](https://www.r-project.org/) package that fits state-space models to
animal tracking data using sequential Monte Carlo algorithms (particle
filters and smoothers). This wraps the high-performance
[`Julia`](https://julialang.org) package
[`Patter.jl`](https://github.com/edwardlavender/patter.jl).

<figure>
<img src="./fig/particle-maps.png"
alt="patter implements particle filtering and smoothing algorithms. The algorithms reconstruct movements, patterns of space use and residency from animal-tracking data, particularly in acoustic telemetry systems." />
<figcaption aria-hidden="true"><strong><code>patter</code> implements
particle filtering and smoothing algorithms</strong>. The algorithms
reconstruct movements, patterns of space use and residency from
animal-tracking data, particularly in acoustic telemetry
systems.</figcaption>
</figure>

# Dependencies

The project was built in [R](https://www.r-project.org/) (version 4.4.1)
in [RStudio](https://www.rstudio.com/).

# Resources

For `patter` users, these are the most relevant resources:

1.  **`doc/`** contains workshop documentation:

    - **`patter-v.1.0-1/patter-intro.html` is a comprehensive
      introductory presentation for `patter`**
      - Download this file for a comprehensive introduction to `patter`

2.  **`R/`** contains `R` scripts (workshop exercises):

    - **`intro-patter.R` is an introductory exercise for `patter`**
      - Download this file for a practical introduction to the basic
        `patter` workflow

The following directories are relevant for users who clone the project:

1.  **`renv/`** implements local dependency management for `R` packages

2.  **`data/`** contains project management files:

    - `inst/dependencies.rds` is a list of dependencies;
    - `inst/session-info.rds` is a record of information about the `R`
      Session;
    - `inst/tree.rds` is a record of the project directory tree; <br/>

3.  **`dev/`** contains project-management scripts:

    - `01-dev.R` records project set up and development;
    - `02-clone.R` is used to clone the project;
    - `03-utils.R` contains supporting utilities;

4.  **`Julia/`** is the `Julia` project directory

5.  **`fig/`** contains selected figure(s)

# Project clones

To clone the project, follow the instructions in `dev/02-clone.R` to
install packages and directories:

- **Packages.** Work through `dev/02-clone.R` to use
  [`renv`](https://rstudio.github.io/renv/articles/renv.html) to
  regenerate the local project library. Packages can also be manually
  reinstalled via `02-clone.R`.
- **Directories.** Rebuild the project directory tree, via
  `dv::use_template_proj()` and `dv::use_template_tree()`.

Clone the project if you experience issues with package versions on your
own machine and you want to work through the introductory `R` script(s).

# Issues

Please report issues on
[GitHub](https://github.com/edwardlavender/patter/issues).

# Code of conduct

Please note that this project is released with a [Contributor Code of
Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

------------------------------------------------------------------------