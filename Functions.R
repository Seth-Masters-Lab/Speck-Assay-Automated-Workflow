library(flowCore)
library(flowAI)
library(flowWorkspace)
library(openCyto)
library(ggcyto)
library(gridExtra)
library(drc)
library(writexl)
library(drda)


################################################################################

## List Channel Info

channelInfo <- function(path){
  myfiles <- list.files(path = path, pattern = ".fcs", ignore.case = T)
  fs <- read.flowSet(myfiles[1], path = path, alter.names = T)
  print(fs[[1]]@parameters@data[1])
}


## Import FCS files and perform cleanup if desired ##

fcsImportLogicle <- function(path, clean, logTrans){
  
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
    
    trans <- logicleTransform()
    transformat <- transformList(colnames(fs[,lower:upper]), trans)
    fs <- transform(fs, transformat)
    
  }
  resultDir <<- sub("Data", 'Results', path)
  dir.create(resultDir, recursive = T)
  dir.create(paste0(resultDir, '/ec50_curves/'))
  return(fs)
}

fcsImportLogicleLib <- function(path, clean, logTrans){
  
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
    trans <- logicleTransform()
    transformat <- transformList(colnames(fs[,lower:upper]), trans)
    fs <- transform(fs, transformat)
  }
  resultDir <<- sub("Data", 'Results', path)
  dir.create(resultDir, recursive = T)
  dir.create(paste0(resultDir, '/ec50_curves/'))
  return(fs)
}


fcsImportLog <- function(path, clean, logTrans){
  
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
    
    trans <- logTransform(logbase = 10, r = 1, d = 1)
    t2 <- transformList(colnames(fs[,lower:upper]), trans)
    
    fs <- suppressWarnings(transform(fs, t2))
  }
  resultDir <<- sub("Data", 'Results', path)
  dir.create(sub("Data", 'Results', path), recursive = T)
  return(fs)
}


################################################################################


## 2D Gate

gate2d <- function(gatingSet, parentPop, xchannel, ychannel, quantile, name, plot, kpop, save, target){
  if(missing(target)){
    target = NULL
  }
  setData <- gs_pop_get_data(gatingSet, parentPop)
  
  gate <- fsApply(setData, function(fr) openCyto::gate_flowclust_2d(
    fr,
    xChannel = xchannel,
    yChannel = ychannel,
    K = kpop,
    target = target,
    quantile = quantile))
  gs_pop_add(gatingSet, gate, parent = parentPop, name = name)
  recompute(gatingSet)
  if(plot == T){
    print(ggcyto(gatingSet, subset = parentPop, aes(x = {{xchannel}}, y = {{ychannel}})) + geom_hex(bins = 100) +
      geom_gate(name))
  }
  if(save == T){
    ggsave(filename = paste0(name, '.png'), device = 'png',
           path = resultDir,
           limitsize = F,
           width = 1920,
           height = 1080,
           units = 'px',
           scale = 4,
           plot = ggcyto(gatingSet, subset = parentPop, aes(x = {{xchannel}}, y = {{ychannel}})) + geom_hex(bins = 200) +
             geom_gate(name))
  }
}


## 2D Gate control

gate2dc <- function(gatingSet, parentPop, xchannel, ychannel, quantile, name, plot, kpop, save, target, controlSample){
  if(missing(target)){
    target = NULL
  }
  fr <- gh_pop_get_data(gs[[controlSample]], parentPop, returnType = 'flowFrame')
  gate <- openCyto::gate_flowclust_2d(fr,
                                     xChannel = xchannel,
                                     yChannel = ychannel,
                                     K = kpop,
                                     target = target,
                                     quantile = quantile)
  gs_pop_add(gatingSet, gate, parent = parentPop, name = name)
  recompute(gatingSet)
  if(plot == T){
    print(ggcyto(gatingSet, subset = parentPop, aes(x = {{xchannel}}, y = {{ychannel}})) + geom_hex(bins = 100) +
            geom_gate(name))
  }
  if(save == T){
    ggsave(filename = paste0(name, '.png'), device = 'png',
           path = resultDir,
           limitsize = F,
           width = 1920,
           height = 1080,
           units = 'px',
           scale = 4,
           plot = ggcyto(gatingSet, subset = parentPop, aes(x = {{xchannel}}, y = {{ychannel}})) + geom_hex(bins = 200) +
             geom_gate(name))
  }
}


## 1D Gate

gate1d <- function(gatingSet, parentPop, xchannel, range, name, plot, positive, smoothing, peaks, save){
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
  if(save == T){
    ggsave(filename = paste0(name, '.png'), device = 'png',
           path = resultDir,
           limitsize = F,
           width = 1920,
           height = 1080,
           units = 'px',
           scale = 4,
           plot = ggcyto(gatingSet, subset = parentPop, aes(x = {{xchannel}})) + geom_density() +
             geom_gate(name))
  }
}

## 1D Gate control

gate1dc <- function(gatingSet, parentPop, xchannel, range, name, plot, positive, smoothing, peaks, save, controlSample){
  fr <- gh_pop_get_data(gs[[controlSample]], parentPop, returnType = 'flowFrame')
  gate <- openCyto:::.mindensity(fr, 
                                 channels = xchannel,
                                 gate_range = range,
                                 positive = positive,
                                 adjust = smoothing,
                                 peaks = peaks)
  gs_pop_add(gatingSet, gate, parent = parentPop, name = name)
  recompute(gatingSet)
  if(plot == T){
    print(ggcyto(gatingSet, subset = parentPop, aes(x = {{xchannel}})) + geom_density() +
            geom_gate(name))
  }
  if(save == T){
    ggsave(filename = paste0(name, '.png'), device = 'png',
           path = resultDir,
           limitsize = F,
           width = 1920,
           height = 1080,
           units = 'px',
           scale = 4,
           plot = ggcyto(gatingSet, subset = parentPop, aes(x = {{xchannel}})) + geom_density() +
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
    temp <- temp[is.finite(temp)]
    temp <- list(temp)
    speckPosRaw <<- c(speckPosRaw, temp)
    
    # Export speck negative population
    temp <- exprs(speckNegData[[i]])
    temp <- temp[,facsChannel]
    temp <- temp[is.finite(temp)]
    temp <- list(temp)
    speckNegRaw <<- c(speckNegRaw, temp)
    
    # Export all asc cells
    temp <- exprs(totalCellData[[i]])
    temp <- temp[,facsChannel]
    temp <- temp[is.finite(temp)]
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
  
  plotRange <- c(0, 5)
  binStart <- plotRange[1]
  binSize <- (plotRange[2] - plotRange[1]) / 3
  bin <<- c()
  speckPosCounts <<- c()
  speckNegCounts <<- c()
  binEnd <- 0
  
  while(binEnd < plotRange[2]){
    binEnd <- binStart + binSize
    posBinCount <- sum(speckPosRaw[[index]] >= binStart & speckPosRaw[[index]] <= binEnd)
    negBinCount <- sum(speckNegRaw[[index]] >= binStart & speckNegRaw[[index]] <= binEnd)
    # if(posBinCount+negBinCount < 100){
    #   break
    # }
    bin <<- c(bin, binEnd)
    speckPosCounts <<- c(speckPosCounts, posBinCount)
    speckNegCounts <<- c(speckNegCounts, negBinCount)
    binStart <- binStart + stepLen}
  
  return(data.frame(bins = bin, SpeckPos = speckPosCounts, SpeckNeg = speckNegCounts))
}

################################################################################

## Geometric Mean ##
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
