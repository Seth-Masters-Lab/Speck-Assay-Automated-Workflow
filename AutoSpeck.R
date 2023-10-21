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
ggcyto(gs, subset = 'root', aes(x = 'FSC.A', y = 'SSC.A')) + geom_hex(bins = 100)

# Debris Gate
gate1d(gatingSet = gs, parentPop = 'root', xchannel = 'FSC.A', name = 'debris',
       plot = T, positive = T, range = c(0,1e5))

#Single cell gates
gate1d(gs, 'debris', xchannel = 'SSC.W', range = c(0,5), name = 'single1', positive = F, plot = T)

gate2d(gs, 'single1', xchannel = 'FSC.A', ychannel = 'FSC.H', quantile = 0.98,
       name = 'single2', plot = T)

# ASC Gate

