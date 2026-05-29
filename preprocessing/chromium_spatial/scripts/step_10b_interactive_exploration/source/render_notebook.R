#!/usr/bin/env Rscript

suppressPackageStartupMessages({ library(rmarkdown) })

args <- commandArgs(trailingOnly = TRUE)
notebook_path <- args[1]
output_html <- args[2]
project_root <- args[3]

setwd(project_root)

rmarkdown::render(
  input = notebook_path,
  output_file = basename(output_html),
  output_dir = dirname(output_html),
  params = list(PROJECT_ROOT = project_root),
  envir = new.env()
)
