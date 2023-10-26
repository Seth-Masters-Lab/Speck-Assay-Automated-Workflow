
rm(list=ls())

source('Functions.R')


path <- 'Data/ASC50 calculation/20230622_Nlrp3 library 12/P1'


fs <- fcsImport(path, T, T)

# Create an empty gating set
gs <- GatingSet(fs)

# QA Step
ggcyto(gs, subset = 'root', aes(x = 'FSC.A', y = 'SSC.A')) + geom_hex(bins = 200)

ggsave(filename = paste0('total', '.png'), device = 'png',
       path = resultDir,
       limitsize = F,
       width = 3840,
       height = 2160,
       units = 'px',
       scale = 2,
       plot = ggcyto(gs, subset = 'root', aes(x = 'FSC.A', y = 'SSC.A')) + geom_hex(bins = 200))

# Debris Gate
gate1d(gatingSet = gs, parentPop = 'root', xchannel = 'FSC.A', name = 'debris',
       plot = F, positive = T, range = c(0.1e5,1e5), smoothing = 2, peaks = NULL,
       save = T)


#Single cell gates
gate2d(gs, 'debris', xchannel = 'SSC.W', ychannel = 'SSC.A', quantile = 0.95,
       name = 'single1', plot = F, kpop = 2, save = T, target = c(4, 3.5))

gate2d(gs, 'single1', xchannel = 'FSC.A', ychannel = 'FSC.H', quantile = 0.95,
       name = 'single2', plot = F, kpop = 1, save = T)


# ASC Gate
gate2d(gs, 'single2', xchannel = 'FSC.A', ychannel = 'V450.50.A', 
       quantile = 0.95, name = 'asc', plot = F, kpop = 1, save = T)

#Speck negative/positive gate
gate1d(gs, 'asc', xchannel = 'V450.50.W', range = c(3.9, 4.02), positive = T,
       name = 'speckNegGate', plot = F, smoothing = 1, peaks = NULL, save = T)
gate1d(gs, 'asc', xchannel = 'V450.50.W', range = c(3.9, 4.02), positive = F,
       name = 'speckPosGate', plot = F, smoothing = 1, peaks = NULL, save = T)


# Get cell info for speck populations
exportSingleCell('speckPosGate', 'speckNegGate', 'asc', "B530.30.A")

# Mean of NLRP3 for each channel to evaluate overall expression level
nlrp3_means <- c()
for(i in 1:length(speckName)){
  sampleMean <- gm_mean(speckAll[[i]])
  nlrp3_means <- c(nlrp3_means, sampleMean)
}

# ggcyto(gs, subset = 'asc', aes(x = "B530.30.A")) + geom_density()

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
  
  ggplot(speck, aes(NLRP3, speckPositive)) + geom_line() + labs(title = speckName[[i]])
  
  
  
  # Curve fitting function
  curve_fit <- drda::drda(
    formula = speck$speckPositive~speck$NLRP3,
    data = speck,
    lower_bound = c(-Inf, -Inf, -Inf, 0),
    upper_bound = c(+Inf, 1, +Inf, +Inf),
    mean_function = 'l4'
  )


  ec50 <- curve_fit$coefficients['phi']
  ec50 <- ec50[[1]]
  
  # Filters out ec50 values which are greater than the range of the curve
  # Only relevant for small noisy curves
  if(max(speck$NLRP3) < ec50){
    ec50 <- `is.na<-`(ec50)
    cat("Sample", speckName[i], "ec50 out of bounds \n", sep = " ")}

  
  # ec50 <- 10^ec50
  sec50 <- c(sec50, ec50)
  
  plot(curve_fit, main = speckName[[i]])
}

speck50 <- 1 - (sec50/sec50[1])


# Filter out samples with too low of a speck positive population
# Found that ec50 calculation on these samples was wildly inaccurate due to small
# y range of speck proportion

totalSpeck <- c()
for(i in 1:length(speckName)){
  speckpct <- length(speckPosRaw[[i]]) / length(speckAll[[i]]) * 100
  totalSpeck <- c(totalSpeck, speckpct)
  if(speckpct < 0){ #Changed to zero to include low speck populations for now
    cat("Sample ", speckName[i], " speckpos is too low: ",
        speckpct, "% \n",
        sep = "")
    sec50[i] = `is.na<-`(sec50[i])
    speck50[i] = `is.na<-`(speck50[i])
  }
}



results <- data.frame(well = speckName, 
                      Speck_Total = totalSpeck,
                      Geometric_Mean_NLRP3 = nlrp3_means,
                      ec50 = sec50, 
                      asc50 = speck50)

# print(ggplot(results, aes(well, speck50)) + geom_col() + labs(title = path))



# Exports the raw data as an excel sheet along with the graph of relative ec50


ggsave(filename = 'graph.png', device = 'png',
       path = resultDir,
       limitsize = F,
       width = 3840,
       height = 2160,
       units = 'px',
       scale = 2,
       plot = ggplot(results, aes(well, speck50)) + geom_col() + labs(title = path))

write_xlsx(results, path = paste0(resultDir,"/raw_data.xlsx"))

print('Done!')
