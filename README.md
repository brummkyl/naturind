
<!-- README.md is generated from README.Rmd. Please edit that file -->

# naturind

<!-- badges: start -->

[![R-CMD-check](https://github.com/brummkyl/naturind/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/brummkyl/naturind/actions/workflows/R-CMD-check.yaml)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://choosealicense.com/licenses/mit/)
<!-- badges: end -->

The goal of naturind is to support the calculation of performance
indicators that may be used to evaluate the performance of existing
protected areas and expansion priorities identified using spatial
prioritization algorithms.

## Installation

You can install the development version of naturind like so:

``` r
install.packages("remotes")
library(remotes)
remotes::install_github("brummkyl/naturind", update = "never)
library(naturind)
```

## Key functions

Performance indicators (and their associated metrics) can be calculated
using the provided ‘pi.–’ functions. A limited number of helper
functions (‘create.–’) are also available to support the creation of
required input parameters, including transboundary buffers (‘tb.buffers’
parameter –\> ‘pi.transboundary’) and planning unit-by-region lookup
tables (‘–.sub2’ parameters –\> ‘pi.evp’; ‘pi.mta’; ‘pi.rpz’).

## Outputs

After running any of the ‘pi.–’ functions, outputs will automatically be
saved to the current .Rproject, or working directory. These .csv files
are saved using user-defined ‘project.name’ and ‘version’ parameters as
a built-in version control and naming convention. The files are stored
to the “performance_outputs/version” subfolder.
