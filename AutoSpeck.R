
rm(list=ls())

source('Functions.R')

path <- 'Data/20230621_Nlrp3 library 10/P1'

fs <- fcsImportLogicle(path, T, T)

# Create an empty gating set
gs <- GatingSet(fs)

# Pick sample to be used as control for gating
gatingControl <- strtoi(readline("Please enter sample number to be used as control for gating: "))

# QA Step
ggcyto(gs, subset = 'root', aes(x = 'FSC.A', y = 'SSC.A')) + geom_hex(bins = 200)

ggsave(filename = paste0('total', '.png'), device = 'png',
       path = resultDir,
       limitsize = F,
       width = 1920,
       height = 1080,
       units = 'px',
       scale = 4,
       plot = ggcyto(gs, subset = 'root', aes(x = 'FSC.A', y = 'SSC.A')) + geom_hex(bins = 200))

# Debris Gate
gate1dc(gatingSet = gs, parentPop = 'root', xchannel = 'FSC.A', name = 'debris',
       plot = F, positive = T, range = c(0.1e5,1e5), smoothing = 1.5, peaks = NULL,
       save = T, controlSample = 1)

#Single cell gates
gate2dc(gs, 'debris', xchannel = 'SSC.W', ychannel = 'SSC.A', quantile = 0.95,
       name = 'single1', plot = F, kpop = 2, save = T, target = c(4, 3.5), controlSample = gatingControl)

gate2dc(gs, 'single1', xchannel = 'FSC.A', ychannel = 'FSC.H', quantile = 0.95,
       name = 'single2', plot = F, kpop = 1, save = T, controlSample = gatingControl)


# ASC Gate
# gate2dc(gs, 'single2', xchannel = 'FSC.A', ychannel = 'V450.50.A', 
#        quantile = 0.95, name = 'asc', plot = F, kpop = 1, save = T,
#        controlSample = gatingControl)


#Speck negative/positive gate
gate1dc(gs, 'single2', xchannel = 'V450.50.W', range = NULL, positive = T,
       name = 'speckNegGate', plot = F, smoothing = 1.5, peaks = NULL, save = T,
       controlSample = gatingControl)
gate1dc(gs, 'single2', xchannel = 'V450.50.W', range = NULL, positive = F,
       name = 'speckPosGate', plot = F, smoothing = 1.5, peaks = NULL, save = T,
       controlSample = gatingControl)


# Get cell info for speck populations
exportSingleCell('speckPosGate', 'speckNegGate', 'single2', "B530.30.A")

# Mean of NLRP3 for each channel to evaluate overall expression level
nlrp3_means <- c()
for(i in 1:length(speckName)){
  sampleMean <- gm_mean(speckAll[[i]])
  nlrp3_means <- c(nlrp3_means, sampleMean)
}

sec50 <- c()
minimum <- c()
height <- c()
slope <- c()
midpoint <- c()
residual_std_err <- c()
loglik <- c()
aic <- c()
bic <- c()
ec50Plot <- c()

for(i in 1:length(speckName)){
  
  # Step-binning for data
  binData <- stepBin(i, 0.1, speckAll, speckPosRaw, speckNegRaw)
  
  if(!is.null(bin)){
    #Distributions of speck- and speck+ populations
    
    ggplot() +
      geom_point(data = binData, aes(x = bins, y = SpeckPos), col = "red") +
      geom_point(data = binData, aes(x = bins, y = SpeckNeg), col = "blue")
    
    # Calculate proportions of speck+ cells
    speck <- data.frame(NLRP3 = binData$bins, speckPositive = binData$SpeckPos/(binData$SpeckPos + binData$SpeckNeg))
    
    ggplot(speck, aes(NLRP3, speckPositive)) + geom_point() + labs(title = speckName[[i]])
    
    
    
    # Curve fitting function
    
    tryCatch(
      {
        curve_fit <- drda::drda(
          formula = speck$speckPositive~speck$NLRP3,
          data = speck,
          lower_bound = c(0, -Inf, -Inf, 0),
          upper_bound = c(1, 1, +Inf, +Inf),
          mean_function = 'l4'
        )
        
        curve_data <- summary(curve_fit)
        minimum <- c(minimum, curve_data$coefficients[1])
        height <- c(height, curve_data$coefficients[2])
        slope <- c(slope, curve_data$coefficients[3])
        midpoint <- c(midpoint, curve_data$coefficients[4])
        residual_std_err <- c(residual_std_err, curve_data$sigma)
        loglik <- c(loglik, curve_data$loglik)
        aic <- c(aic, curve_data$aic)
        bic <- c(bic, curve_data$bic)
        
        ec50 <- curve_fit$coefficients['phi']
        ec50 <- ec50[[1]]
        
        # Filters out ec50 values which are greater than the range of the curve
        # Only relevant for small noisy curves
        if(max(speck$NLRP3) < ec50){
          ec50 <- `is.na<-`(ec50)
          cat("Sample", speckName[i], "ec50 out of bounds \n", sep = " ")}
        
        trans <- inverseLogicleTransform(logicleTransform())
        ec50 <- trans(ec50)
        sec50 <- c(sec50, ec50)
        
        png(filename = paste0(resultDir, '/ec50_curves/', speckName[[i]],'.png'))
        # plot(curve_fit, main = speckName[[i]])
        plot(curve_fit, main = speckName[[i]], xlim = c(0,4), ylim = c(0,1))
        dev.off()
        
      
        
        
        
      },
      # If fitting fails: 
      error=function(error_message) {
        message("")
        message(speckName[[i]])
        message("Fitting Failed")
        message("ERROR:")
        message(error_message)
        
        ec50 <<- NA
        sec50 <<- c(sec50, ec50)
        minimum <<- c(minimum, NA)
        height <<- c(height, NA)
        slope <<- c(slope, NA)
        midpoint <<- c(midpoint, NA)
        residual_std_err <<- c(residual_std_err, NA)
        loglik <<- c(loglik, NA)
        aic <<- c(aic, NA)
        bic <<- c(bic, NA)
        print(paste0("sample ", speckName[i], " unable to calculate ec50 \n"))
      }
    )
    
    
  } else {
    ec50 <- NA
    sec50 <- c(sec50, ec50)
    minimum <- c(minimum, NA)
    height <- c(height, NA)
    slope <- c(slope, NA)
    midpoint <- c(midpoint, NA)
    residual_std_err <- c(residual_std_err, NA)
    loglik <- c(loglik, NA)
    aic <- c(aic, NA)
    bic <- c(bic, NA)
    print(paste0("sample ", speckName[i], " unable to calculate ec50 \n"))
    }
}


# ASC50 calculation (1 - ratio of sample ec50 vs WT)
speck50 <- 1 - (sec50/sec50[gatingControl])


# Filter out samples with too low of a speck positive population
# Found that ec50 calculation on these samples was inaccurate due to small
# y range of speck proportion

totalSpeck <- c()
for(i in 1:length(speckName)){
  speckpct <- length(speckPosRaw[[i]]) / length(speckAll[[i]]) * 100
  totalSpeck <- c(totalSpeck, speckpct)
  if(speckpct < 0 | is.na(speckpct)){ #Changed to zero to include low speck populations for now
    cat("Sample ", speckName[i], " speckpos is too low: ",
        speckpct, "% \n",
        sep = "")
    sec50[i] = `is.na<-`(sec50[i])
    speck50[i] = `is.na<-`(speck50[i])
  }
}



results <- data.frame(well = speckName, 
                      asc50 = speck50,
                      ec50 = sec50,
                      Speck_Total = totalSpeck,
                      Geometric_Mean_NLRP3 = nlrp3_means,
                      minimum = minimum,
                      height = height,
                      slope = slope,
                      midpoint = midpoint,
                      residual_std_err = residual_std_err,
                      loglik = loglik,
                      aic = aic,
                      bic = bic)

print(ggplot(results, aes(well, speck50)) + geom_col() + labs(title = path))

# Exports the raw data as an excel sheet along with the graph of relative ec50
ggsave(filename = 'graph.png', device = 'png',
       path = resultDir,
       limitsize = F,
       width = 1920,
       height = 1080,
       units = 'px',
       scale = 4,
       plot = ggplot(results, aes(well, speck50)) + geom_col() + labs(title = path))

write_xlsx(results, path = paste0(resultDir,"/raw_data.xlsx"))
print(paste0(path, ' Done!'))
