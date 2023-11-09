#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

use_rstudio = F

list.of.packages <- c("rgl","ggplot2","knitr","rglwidget", "OneR")

len <- length
#
#You should be able to simply reuse the following lines of code as is
#
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
#
if(length(new.packages)) install.packages(new.packages)
#
# By now we have installed the requisite packages. Time to load them .
#
lapply(list.of.packages,function(x){library(x,character.only=TRUE)})

source("windows_override.R")
source("structure_analysis.R")
source("matrix_presets.R")
source("get_red_mats.R")
source("cca_analysis.R")
source("intrinsic_dim_analysis.R")
