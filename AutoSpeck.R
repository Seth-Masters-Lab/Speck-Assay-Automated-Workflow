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

path <- 'Data/20230503_Nlrp3 library 1'


fs <- fcsImport(path, T, T)

# Create an empty gating set
gs <- GatingSet(fs)

# QA Step
ggcyto(gs[1], subset = 'root', aes(x = 'FSC.A', y = 'SSC.A')) + geom_hex(bins = 100)

# Debris Gate
gate1d(gatingSet = gs, parentPop = 'root', xchannel = 'FSC.A', name = 'debris',
       plot = F, positive = T, range = c(0,1e5), smoothing = 2, peaks = NULL)


# ggcyto(gs[1], subset = 'debris', aes(x = 'SSC.W')) + geom_density() +
#     geom_gate('single1')


#Single cell gates
gate1d(gs, 'debris', xchannel = 'SSC.W', range = c(4,4.3), name = 'single1', 
       positive = F, plot = F, smoothing = 2, peaks = NULL)

gate2d(gs, 'single1', xchannel = 'FSC.A', ychannel = 'FSC.H', quantile = 0.95,
       name = 'single2', plot = F, kpop = 1)

# ASC Gate
gate2d(gs, 'single2', xchannel = 'FSC.A', ychannel = 'V450.50.A', 
       quantile = 0.95, name = 'asc', plot = F, kpop = 1)

#Speck negative/positive gate
gate1d(gs, 'asc', xchannel = 'V450.50.W', range = c(3.9, 4.05), positive = T,
       name = 'speckNegGate', plot = F, smoothing = 1, peaks = NULL)
gate1d(gs, 'asc', xchannel = 'V450.50.W', range = c(3.9, 4.05), positive = F,
       name = 'speckPosGate', plot = F, smoothing = 1, peaks = NULL)

ggcyto(gs, subset = 'asc', aes(x = 'V450.50.W', y = 'SSC.A')) + geom_hex(bins = 100) + geom_gate()


# Get cell info for speck populations
exportSingleCell('speckPosGate', 'speckNegGate', 'asc', "B530.30.A")


sec50 <- c()
for(i in 1:length(speckName)){
  
  # Step-binning for data
  binData <- stepBin(i, 0.01, speckAll, speckPosRaw, speckNegRaw)
  
  #Distributions of speck- and speck+ populations
  ggplot() +
    geom_line(data = binData, aes(x = bins, y = SpeckPos), col = "red") + 
    geom_line(data = binData, aes(x = bins, y = SpeckNeg), col = "blue")
  
  # Calculate proportions of speck+ cells
  speck <- data.frame(NLRP3 = binData$bins, speckPositive = binData$SpeckPos/(binData$SpeckPos + binData$SpeckNeg))
  
  print(ggplot(speck, aes(NLRP3, speckPositive)) + geom_line() + labs(title = speckName[[i]]))
  
  
  
  # Curve fitting function - NEEDS MORE WORK
  curve_fit <- drm(
    formula = speck$speckPositive ~ speck$NLRP3,
    data = speck,
    logDose = 10,
    fct = LL.4(names = c("hill", "min_value", "max_value", "ec_50")))
  
  ec50 <- curve_fit$coefficients['ec_50:(Intercept)']
  ec50 <- ec50[[1]]
  
  # Filters out ec50 values which are greater than the range of the curve
  # Only relevant for small noidy curves
  if(10^(max(speck$NLRP3)) < ec50 | ec50 < 1){
    ec50 <- `is.na<-`(ec50)
    cat("Sample", speckName[i], "ec50 out of bounds \n", sep = " ")}
  
  sec50 <- c(sec50, ec50)
  
}

speck50 <- 1 - (sec50/sec50[1])


# Filter out samples with too low of a speck positive population
# Found that ec50 calculation on these samples was wildly inaccurate due to small
# y range of speck proportion

totalSpeck <- c()
for(i in 1:length(speckName)){
  speckpct <- length(speckPosRaw[[i]]) / length(speckAll[[i]]) * 100
  totalSpeck <- c(totalSpeck, speckpct)
  if(speckpct < 5){
    cat("Sample ", speckName[i], " speckpos is too low: ", 
        speckpct, "% \n",
        sep = "")
    sec50[i] = `is.na<-`(sec50[i])
    speck50[i] = `is.na<-`(speck50[i])
  }
}



results <- data.frame(well = speckName, 
                      SpeckTotal = totalSpeck, 
                      ec50 = sec50, 
                      speck50 = speck50)

print(ggplot(results, aes(well, speck50)) + geom_col() + labs(title = path))



# Exports the raw data as an excel sheet along with the graph of relative ec50
ggsave(filename = paste(sub('Data/', x = path, replacement = ""), 'results.png', sep = "_"), device = 'png',
       path = 'Results',
       limitsize = F,
       width = 3840,
       height = 2160,
       units = 'px',
       scale = 2,
       plot = ggplot(results, aes(well, speck50)) + geom_col() + labs(title = path))

write_xlsx(results, path = paste(paste0('Results/',sub('Data/', x = path, replacement = "")),"results.xlsx", sep = "_"))

