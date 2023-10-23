library(flowCore)
library(flowAI)
library(flowWorkspace)
library(openCyto)
library(ggcyto)
library(gridExtra)
library(drc)
library(writexl)

rm(list=ls())

source('Functions.R')

path <- 'Testing_NLRP3_library_12_sample_1A-H'

fs <- fcsImport(path, T, T)

# Create an empty gating set
gs <- GatingSet(fs)

# QA Step
ggcyto(gs, subset = 'root', aes(x = 'FSC.A', y = 'SSC.A')) + geom_hex(bins = 100)

# Debris Gate
gate1d(gatingSet = gs, parentPop = 'root', xchannel = 'FSC.A', name = 'debris',
       plot = F, positive = T, range = c(0,1e5), smoothing = 2, peaks = NULL)

#Single cell gates
gate1d(gs, 'debris', xchannel = 'SSC.W', range = c(0,5), name = 'single1', positive = F, plot = T, smoothing = 2, peaks = NULL)

gate2d(gs, 'single1', xchannel = 'FSC.A', ychannel = 'FSC.H', quantile = 0.95,
       name = 'single2', plot = T, kpop = 1)

# ASC Gate
gate2d(gs, 'single2', xchannel = 'FSC.A', ychannel = 'V450.50.A', quantile = 0.95, name = 'asc', plot = F, kpop = 1)

#Speck negative/positive gate
gate1d(gs, 'asc', xchannel = 'V450.50.W', range = c(3.9, 4.05), positive = T,
       name = 'speckNegGate', plot = T, smoothing = 1, peaks = NULL)
gate1d(gs, 'asc', xchannel = 'V450.50.W', range = c(3.9, 4.05), positive = F,
       name = 'speckPosGate', plot = T, smoothing = 1, peaks = NULL)

ggcyto(gs, subset = 'asc', aes(x = 'V450.50.W', y = 'SSC.A')) + geom_hex(bins = 100) + geom_gate()


# Get cell info for speck populations
exportSingleCell('speckPosGate', 'speckNegGate', 'asc', "B530.30.A")


sec50 <- c()
for(i in 1:length(speckName)){
  binData <- stepBin(i, 0.01, speckAll, speckPosRaw, speckNegRaw)
  
  ggplot() +
    geom_line(data = binData, aes(x = bins, y = SpeckPos), col = "red") + 
    geom_line(data = binData, aes(x = bins, y = SpeckNeg), col = "blue")
  
  speck <- data.frame(NLRP3 = binData$bins, speckPositive = binData$SpeckPos/(binData$SpeckPos + binData$SpeckNeg))
  print(ggplot(speck, aes(NLRP3, speckPositive)) + geom_line() + labs(title = speckName[[i]]))
  
}
