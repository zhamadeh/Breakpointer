#This script is for collecting Breakpoints from multiple libraries
#Also for creating seperate files for different cell lines (ie. blm, recq5 ..)
setwd("/")
args = commandArgs(trailingOnly=TRUE)

suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(miceadds)))

#This is the master function
master <- function(inputDir,outputDir,metricsDir){
	bind <- masterCollectWidths(list.files(inputDir,full.names = T))
	cat("\nCollecting breakpoints from",args[1], "...\n\n")
	cat("Annotating breakpoinst with metadata from",args[3],"...\n\n")
	bind <- annotatingBreakpoints(bind,metricsDir)
	bind <- removingInversions(bind)
	cat("Saving breakpoints to",args[2],"...\n\n")
	write.table(bind,paste0(outputDir,args[4]),append = F,row.names = F,col.names = T,quote=F,sep="\t")
	return(bind)
}

masterCollectWidths <- function(inputDirectory){
	bind=data.frame()
	#file=inputDirectory[1]
	for (file in inputDirectory){
		data <- load.Rdata2(file)
		breakpoints <- as.data.frame(data$breaks)
		#cat(paste0(ncol(breakpoints),data$ID))
		if (nrow(breakpoints)!= 0 ){
			breakpoints$library = data$ID
			bind <- rbind(breakpoints,bind)
			
		}
		
	}
	return(bind)
}

annotatingBreakpoints <- function(bind,metricsDir){
	bind$gene <- "gene"
	bind$Library <- bind$library
	suppressWarnings(bind <- bind %>% separate(library, c("a","b","c","d","e","f"), "[_-]+"))
	

	
	
	for (row in 1:nrow(bind)){
		for (letter in c("a","b","c","d","e","f")){
			#cat(bind[1,letter])
			
			if (bind[row,letter]=="WT" | bind[row,letter]=="wt"){
				bind[row,"gene"]="WT"
			}
			else if (bind[row,letter]=="blm" | bind[row,letter]=="BLM" ) {
				bind[row,"gene"]="BLM"
			}
			
			else if (bind[row,letter]=="RECQL5" | bind[row,letter]=="recql5" | bind[row,letter]=="RECQ5" | bind[row,letter]=="recq5" ) {
				if (bind[row,"gene"]=="BLM"){
					bind[row,"gene"]="BLM/RECQL5"
				}
				else{
					bind[row,"gene"]="RECQL5"
				}
			}
		}
		
	}
	bind <- select(bind,-c(a,b,c,d,e,f, deltaW, genoT, strand))
	bind$Library <- str_split_fixed(bind$Library,"[.]",2)[,1]
	metrics <- collectMetrics(metricsDir)
	bind <- merge(bind,metrics,by.x="Library")
	bind <- bind[,c(2,3,4,5,6,7,8,9,1)]
	return(bind)
}

collectMetrics <- function(metricsDir){
	files <- list.files(metricsDir, full.names=T)
	metrics <- data.frame()
	if (length(files)!=0){
		#file=files[2]
		for (file in files){
			met <- read.table(file,header=T,fill=T) 
			if ("Reads_aligned_postfiltering" %in% colnames(met)){
				met <- met %>% select(Library, Reads_aligned_postfiltering,Background,Reads_per_Mb)
				names(met)[names(met) == "Reads_aligned_postfiltering"] <- "Reads"
			}
			else if ("Postfiltering_reads_aligned" %in% colnames(met)){
				met <- met %>% select(Library, Postfiltering_reads_aligned , Background,Reads_per_Mb)
				names(met)[names(met) == "Postfiltering_reads_aligned"] <- "Reads"
			}
			metrics <- rbind(met,metrics)
			 
		}
	}
	return(metrics)
}

removingInversions <- function(bind){
	countedSCEs <- bind#[-c(3,22:24,27,33,35,118,123:126,178,182,191,224,226,251,279,287,293,302:308,320,321,334,342,348,362,395,399,400,402,405,412,443,448),]
	rows <- c()
	#removing chr15 and chr16 inversions at 60-90Mb and 21-23Mb
	#removes 258 inversion calls (129) from 1062 leaving behind 804 SCEs
	for (row in 1:nrow(countedSCEs)){
		if (countedSCEs[row,1]=="chr15"){
			if (countedSCEs[row,2] <= 60000000 & countedSCEs[row,3] >= 90000000 ){
				rows <- append(row,rows)
			}
			else if (countedSCEs[row,2] >= 60000000 & countedSCEs[row,2] <= 90000000 & countedSCEs[row,3] <= 90000000 & countedSCEs[row,3] >= 60000000){
				rows <- append(row,rows)
			}
			else if (countedSCEs[row,2] >= 60000000 & countedSCEs[row,2] <= 90000000 & countedSCEs[row,3] >= 90000000){
				rows <- append(row,rows)
			}
			else if (countedSCEs[row,2] <= 60000000 & countedSCEs[row,3] <= 90000000 & countedSCEs[row,3] >= 60000000){
				rows <- append(row,rows)
			}
		}
		else if (countedSCEs[row,1]=="chr16"){
			if (countedSCEs[row,2] <= 21000000 & countedSCEs[row,3] >= 23000000 ){
				rows <- append(row,rows)
			}
			else if (countedSCEs[row,2] >= 21000000 & countedSCEs[row,2] <= 23000000 & countedSCEs[row,3] <= 23000000 & countedSCEs[row,3] >= 21000000){
				rows <- append(row,rows)
			}
			else if (countedSCEs[row,2] >= 21000000 & countedSCEs[row,2] <= 23000000 & countedSCEs[row,3] >= 23000000){
				rows <- append(row,rows)
			}
			else if (countedSCEs[row,2] <= 21000000 & countedSCEs[row,3] <= 23000000 & countedSCEs[row,3] >= 21000000){
				rows <- append(row,rows)
			}
		}
		else if (countedSCEs[row,1]=="chr8"){
			if (countedSCEs[row,2] >= 70569000 & countedSCEs[row,2] <= 70800000 & countedSCEs[row,3] <= 70800000 & countedSCEs[row,3] >= 70569000){
				rows <- append(row,rows)
			}
			else if (countedSCEs[row,2] >= 70569000 & countedSCEs[row,2] <= 70800000 & countedSCEs[row,3] >= 70800000){
				rows <- append(row,rows)
			}
			else if (countedSCEs[row,2] <= 70569000 & countedSCEs[row,3] <= 70800000 & countedSCEs[row,3] >= 70569000){
				rows <- append(row,rows)
			}
		}
	}
	countedSCEs <- (countedSCEs[-c(rows),])
	return(countedSCEs)
}

bind <- master(args[1],args[2],args[3])


#write.table(bind,"/Users/zeidh/Desktop/unfilt.bed",sep="\t",quote=F,row.names = F,col.names = F)
#write.table(bind2,"/Users/zeidh/Desktop/filt.bed",sep="\t",quote=F,row.names = F,col.names = F)
#args2=c("/Users/zeidh/Desktop/Data/RDATA/ALL_BL/","/Users/zeidh/Desktop/BreakpointerGithub/Breakpoints/","/Users/zeidh/Desktop/Data/METRICS")
#inputDir="/Users/zeidh/Desktop/Data/RDATA/ALL_BL/"

#inputDir=args[1]
#inputDirectory=list.files(inputDir,full.names = T)
#outputDir=args[2]
#metricsDir=args[3]
#args=c("/Users/zeidh/Desktop/BreakpointerGithub/Input/RData_blacklisted/","/Users/zeidh/Desktop/BreakpointerGithub/Output/Breakpoints/","/Users/zeidh/Desktop/BreakpointerGithub/Input/Metrics/")


#inputDirectory=list.files(inputDir,full.names = T)
#(file=inputDirectory[2])
#utputDir="/Users/zeidh/Desktop/BreakpointerGithub/Breakpoints/"
#metricsDir="/Users/zeidh/Desktop/Data/METRICS"


								   