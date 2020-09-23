suppressWarnings(suppressMessages(library(miceadds)))
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(breakpointR)))
suppressWarnings(suppressMessages(library(GenomicRanges)))
#setwd("/")
#setwd("/Users/zeidh/Desktop/BreakpointerGithub/")
args = commandArgs(trailingOnly=TRUE)
n=0
#wrapper function for filtering and plotting functions
master <- function(inputDirectory, outputDirectory,plottingDirectory){
	#file=list.files(inputDirectory,full.names = T)[1]
	for (file in list.files(inputDirectory,full.names = T)){
		n=n+1
		perc = round((n/length(list.files(inputDirectory,full.names = T)))*100,2)
		cat(paste0("Filtering ",tail(strsplit(as.character(file),"/")[[1]],1)," ... ",perc,"% done \n"))
		blacklistBreakpoint(as.character(file),outputDirectory)
		
	}
	cat("\n Done filtering out breakpoints ... plotting out chromsome ideograms \n\n")
	#filteredOutputDir <- list.files(outputDirectory)
	plotBreakpoints(outputDirectory,plottingDirectory)
}

removeCentromereSCEs <- function(data,bedfile){
	breaks <- data$breaks
	confint <- data$confint
	centromeres <- read.table(bedfile,header=F) #%>% select(-c(V4))
	centromeres <-centromeres %>% dplyr::rename("seqnames"=V1,"start"=V2,"end"=V3)
	
	centroGRange <- GRanges(centromeres)
	
	breakGRanges <-breaks[-queryHits(findOverlaps(breaks, centroGRange, type="any")),]
	confintGRanges <-confint[-queryHits(findOverlaps(confint, centroGRange, type="any")),]
	
	data$breaks <- breakGRanges
	data$confint <- confintGRanges
	return(data)
}

blacklistBreakpoint <- function(file,outputDirectory){
	#setwd("/")
	data <- load.Rdata2(file)
	
	breakpoints <- data$breaks
	
	data <- removeCentromereSCEs(data, "Input/Centromeres/centromeres2.txt")
	#write.table(as.character(length(data$breaks)),"/Users/zeidh/Downloads/test2.txt",quote = F,row.names = F,col.names = F,append = T)
	class.breakpoint <- "BreakPoint"
	class(data) <- class.breakpoint
	outputFileName = paste0("bl_",data$ID)
	outputDir <- paste0(outputDirectory,outputFileName,".RData")
	save(data, file=outputDir)
}

#plot breakpoints
plotBreakpoints <- function(fileDir,plotDir){
	datapath <- file.path(fileDir)
	plotspath <- file.path(plotDir)
	files2plot <- list.files(datapath, pattern = ".RData", full.names = TRUE)
	breakpointR::plotBreakpoints(files2plot=files2plot, file=file.path(plotspath, 'breaksPlot***BL_all.pdf')) -> beQuiet
}

master(args[1],args[2],args[3])

#args=c("Input/RData_good/","Input/RData_blacklisted/","Output/BPR_breaksPlots/")
#inputDirectory =  args[1]
#list.files(inputDirectory)
#bedfile="/Users/zeidh/Desktop/Coding/for_ERIBA/blacklist/centromeres2.txt"
#outputDirectory = args[2]
#outputDir=outputDirectory 
#plottingDirectory=args[3]
#plotDir=args[3]
#file=inputDirectory[20]


