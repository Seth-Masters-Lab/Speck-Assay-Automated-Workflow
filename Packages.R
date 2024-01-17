if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

install.packages("shiny")

BiocManager::install("flowAI")
BiocManager::install("flowCore")
BiocManager::install("flowWorkspace")
BiocManager::install("openCyto")
BiocManager::install("ggcyto")
BiocManager::install("gridExtra")
