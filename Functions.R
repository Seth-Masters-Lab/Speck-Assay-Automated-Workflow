library(flowCore)
library(flowAI)
library(flowWorkspace)
library(openCyto)
library(ggcyto)
library(gridExtra)
library(drc)
library(writexl)



################################################################################



## Import FCS files and perform cleanup if desired ##

fcsImport <- function(path, clean, logTrans){
  
  # Load Data
  myfiles <- list.files(path = path, pattern = ".fcs", ignore.case = T)
  fs <- read.flowSet(myfiles, path = path, alter.names = T)
  
  # Assign well ID to samples
  pData(fs)$well <- gsub(".*_.*_(.*)_.*.fcs","\\1",sampleNames(fs))
  
  # Data Cleaning
  if(clean == T){
    fs <- flow_auto_qc(fs, mini_report = F, html_report = F, fcs_QC = F, folder_results = F)
    }
  
  # Log transformation
  if(logTrans == T){
    print(fs[[1]]@parameters@data[1])
    logRange <- readline("Please input the channels to be converted to log scale, (eg. 5:10): ")
    lower <- gsub(x = logRange, pattern = ":.*", replacement = "")
    lower <- strtoi(lower)
    upper <- gsub(x = logRange, pattern = "*.:", replacement = "")
    upper <- strtoi(upper)
    
    trans <- estimateLogicle(fs[[1]], colnames(fs[,lower:upper]))
    fs <- transform(fs, trans)
  }
  return(fs)
}



################################################################################


## 2D Gate

gate2d <- function(gatingSet, parentPop, xchannel, ychannel, quantile, name, plot, kpop){
  setData <- gs_pop_get_data(gatingSet, parentPop)
  gate <- fsApply(setData, function(fr) openCyto::gate_flowclust_2d(
    fr,
    xChannel = xchannel,
    yChannel = ychannel,
    K = kpop,
    quantile = quantile))
  
  gs_pop_add(gatingSet, gate, parent = parentPop, name = name)
  recompute(gatingSet)
  if(plot == T){
    print(ggcyto(gatingSet, subset = parentPop, aes(x = {{xchannel}}, y = {{ychannel}})) + geom_hex(bins = 100) +
      geom_gate(name))
  }
}


## 1D Gate

gate1d <- function(gatingSet, parentPop, xchannel, range, name, plot, positive, smoothing, peaks){
  setData <- gs_pop_get_data(gatingSet)
  gate <- fsApply(setData, function(fr) openCyto::gate_mindensity(
    fr,
    channel = xchannel,
    gate_range = range,
    positive = positive,
    adjust = smoothing,
    peaks = peaks
  ))
  gs_pop_add(gatingSet, gate, parent = parentPop, name = name)
  recompute(gatingSet)
  if(plot == T){
    print(ggcyto(gatingSet, subset = parentPop, aes(x = {{xchannel}})) + geom_density() +
            geom_gate(name))
  }
}


################################################################################

## Export single cell data for marker of interest (eg. NLRP3) ##

exportSingleCell <- function(speckPosGate, speckNegGate, ascGate, facsChannel){
  
  # Get cell info for speck populations
  speckPosData <- gs_pop_get_data(gs, speckPosGate)
  speckNegData <- gs_pop_get_data(gs, speckNegGate)
  totalCellData <- gs_pop_get_data(gs, ascGate)
  
  # Export single cell data based on 
  speckName <<- c()
  speckPosRaw <<- c()
  speckNegRaw <<- c()
  speckAll <<- c()
  
  for(i in 1:length(speckPosData)){
    # Export well names
    temp <- pData(speckPosData[i])[,2]
    speckName <<- c(speckName, temp)
    
    # Export Speck positive population
    temp <- exprs(speckPosData[[i]])
    temp <- temp[,facsChannel]
    temp <- list(temp)
    speckPosRaw <<- c(speckPosRaw, temp)
    
    # Export speck negative population
    temp <- exprs(speckNegData[[i]])
    temp <- temp[,facsChannel]
    temp <- list(temp)
    speckNegRaw <<- c(speckNegRaw, temp)
    
    # Export all asc cells
    temp <- exprs(totalCellData[[i]])
    temp <- temp[,facsChannel]
    temp <- list(temp)
    speckAll <<- c(speckAll, temp)
  }
  cat("Exported: \n", "speckName \n", "speckPosRaw \n", "speckNegRaw \n", "speckAll", sep = "")
}



################################################################################

## Step Gating ##

# This function takes single cell fluorescence values for speck+ and speck-
# cell populations, applying a single bin that incrementally steps through the
# data. This allows for the generation of a smoothed frequency curve that can
# be used for modelling an EC50

stepBin <- function(index, stepLen, speckAll, speckPosRaw, speckNegRaw){
  
  plotRange <- c(min(speckAll[[index]]), max(speckAll[[index]]))
  binStart <- plotRange[1]
  binSize <- (plotRange[2] - plotRange[1]) / 4
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
  
  return(data.frame(bins = bin, SpeckPos = speckPosCounts, SpeckNeg = speckNegCounts))
}

################################################################################


