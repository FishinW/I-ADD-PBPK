# Project setup ---------------------------------------------------------------
# Run scripts from the repository root. This file creates output directories.
rm(list = ls())
curr.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(curr.dir)
setwd("..")
dir.create(file.path("results", "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("results", "figures"), recursive = TRUE, showWarnings = FALSE)
