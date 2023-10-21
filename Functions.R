



################################################################################


## Step Gating ##

# This function takes single cell fluorescence values for speck+ and speck-
# cell populations, applying a single bin that incrementally steps through the
# data. This allows for the generation of a smoothed frequency curve that can
# be used for modelling an EC50

stepBin <- function(index, stepLen, speckAll, speckPosRaw, speckNegRaw){
  
  plotRange <- c(min(speckAll[[index]]), max(speckAll[[index]]))
  binStart <- plotRange[1]
  binSize <- (plotRange[2] - plotRange[1]) / 5
  bin <<- c()
  speckPosCounts <<- c()
  speckNegCounts <<- c()
  binEnd <- 0
  
  while(binEnd < plotRange[2]){
    binEnd <- binStart + binSize
    posBinCount <- sum(speckPosRaw[[index]] >= binStart & speckPosRaw[[index]] <= binEnd)
    negBinCount <- sum(speckNegRaw[[index]] >= binStart & speckNegRaw[[index]] <= binEnd)
    bin <<- c(bin, binEnd)
    speckPosCounts <<- c(speckPosCounts, posBinCount)
    speckNegCounts <<- c(speckNegCounts, negBinCount)
    binStart <- binStart + stepLen}
}

################################################################################


