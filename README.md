# EvoTrace
A collection of R functions to work on Evolutionary Trace (ET) related jobs.

## Install
The lollipop plots in EvoTrace require ggplot v3.4.0. Use
the this command to install that specific version.
```
remotes::install_version("ggplot2", version = "3.4.0")
```

Vignettes requires pymol to build.
```
remotes::install_github("LichtargeLab/EvoTrace", build_vignettes = TRUE)
```

To install without vignettes:
```
remotes::install_github("LichtargeLab/EvoTrace", build_vignettes = FALSE)
```

## Vignettes
```
browseVignettes("EvoTrace")
```
