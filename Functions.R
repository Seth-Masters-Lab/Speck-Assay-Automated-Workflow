################################################################################



## Import FSC files and perform cleanup if desired

fcsImport <- function(path, clean, logTrans){
  
  # Load Data
  myfiles <- list.files(path = path, pattern = ".fcs", ignore.case = T)
  fs <<- read.flowSet(myfiles, path = path, alter.names = T)
  
  # Assign well ID to samples
  pData(fs)$well <- gsub(".*_.*_(.*)_.*.fcs","\\1",sampleNames(fs))
  
  # Data Cleaning
  if(clean == T){
    fs <<- flow_auto_qc(fs, mini_report = F, html_report = F, fcs_QC = F, folder_results = F)}
  
  # Log transformation
  if(logTrans == T){
    print(fs[[1]]@parameters@data[1])
    logRange <- readline("Please input the channels to be converted to log scale, (eg. 5:10): ")
    lower <- gsub(x = logRange, pattern = ":.*", replacement = "")
    upper <- gsub(x = logRange, pattern = "*.:", replacement = "")
    
    trans <- estimateLogicle(fs[[1]], colnames(fs[,strtoi(lower):strtoi(upper)]))
    fs <<- transform(fs, trans)}
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


