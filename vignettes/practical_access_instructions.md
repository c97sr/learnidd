---
title: "Instructions for accessing the practicals"
date: "`r Sys.Date()`"
output: html_document
---
# Steps
1. Install Rtools from the website https://cran.r-project.org/bin/windows/Rtools/

2. In RStudio, check that .libPaths() returns two paths: the network path and the shared path
```{r, eval=FALSE}
.libPaths()
```

3. If it does, run
```{r, eval=FALSE}
install.packages("callr", lib = .libPaths()[2])
install.packages("rlang", lib = .libPaths()[2])
devtools::install_github("c97sr/learnidd", build_vignettes = TRUE, force = TRUE, lib = .libPaths()[2])
```
If the vignettes don’t build, repeat `install_github()` but with `build_vignettes = FALSE`.

4. Use R to load the vignettes by running
```{r, eval=FALSE}
browseVignettes("learnidd")
```

5. Click on spatial (HTML)*

* answers are in spatial_ans.
