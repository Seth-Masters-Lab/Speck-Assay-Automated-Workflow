rm(list=ls())

library(flowCore)
library(flowAI)
library(flowWorkspace)
library(openCyto)
library(ggcyto)
library(gridExtra)
library(drc)
library(writexl)

source('Functions.R')

path <- 'Testing_NLRP3_library_12_sample_1A-H'

fs <- fcsImport(path, T, T)

# Create an empty gating set
gs <- GatingSet(fs)

# QA Step
ggcyto(fs, aes(x = 'FSC.A', y = 'SSC.A')) + geom_hex(bins = 100)

