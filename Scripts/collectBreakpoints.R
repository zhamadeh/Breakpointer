#This script is for collecting Breakpoints from multiple libraries
#Also for creating seperate files for different cell lines (ie. blm, recq5 ..)
setwd("/")
args = commandArgs(trailingOnly=TRUE)

suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(miceadds)))

#args=c("/Users/zeidh/Desktop/Data/RDATA/ALL_BL/","/Users/zeidh/Desktop/BreakpointerGithub/Breakpoints/","/Users/zeidh/Desktop/Data/METRICS")
#args=c("/Users/zeidh/Desktop/Data/RDATA/Aug28-2020/blacklistedForPlotting","/Users/zeidh/Desktop/BreakpointerGithub/Breakpoints/","/Users/zeidh/Desktop/Data/METRICS")
#This is the master function
master <- function(inputDir,outputDir,metricsDir){
	bind <- masterCollectWidths(list.files(inputDir,full.names = T))
	cat("\nCollecting breakpoints from",args[1], "...\n\n")
	cat("Annotating breakpoinst with metadata from",args[3],"...\n\n")
	bind <- annotatingBreakpoints(bind,metricsDir)
	cat("Saving breakpoints to",args[2],"...\n\n")
	write.table(bind,paste0(outputDir,"genomicBreakpoints_BL_feb25_mar6_aug18_aug28.txt"),append = F,row.names = F,col.names = T,quote=F,sep="\t")
	return(bind)
}

masterCollectWidths <- function(inputDirectory){
	bind=data.frame()
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
	
	feb25 <- read.table(files[3],header=T,fill=T) %>% select(Library, Reads_aligned_postfiltering,Background,Reads_per_Mb) 
	#feb25 <- rename(feb25, c("Reads" = "Reads_aligned_postfiltering"  , "Background" =  "Background" , "ReadsPerMb" = "Reads_per_Mb" ))
	names(feb25)[names(feb25) == "Reads_aligned_postfiltering"] <- "Reads"
	names(feb25)[names(feb25) == "Reads_per_Mb"] <- "ReadsPerMb"
	
	mar6<- read.table(files[4],header=T,fill=T)%>% select(Library, Postfiltering_reads_aligned,Background,Reads_per_Mb)
	#mar6<- rename(mar6, c("Reads"="Postfiltering_reads_aligned" , "Background" =  "Background" , "ReadsPerMb" = "Reads_per_Mb" ))
	names( mar6)[names( mar6) == "Postfiltering_reads_aligned"] <- "Reads"
	names( mar6)[names( mar6) == "Reads_per_Mb"] <- "ReadsPerMb"
	
	aug2020<- read.table(files[1],header=T,fill=T)%>% select(Library, Postfiltering_reads_aligned,Background,Reads_per_Mb)
	#aug2020<- rename(aug2020, c("Reads"="Postfiltering_reads_aligned" , "Background" =  "Background" , "ReadsPerMb" = "Reads_per_Mb" ))
	names(aug2020)[names(aug2020) == "Postfiltering_reads_aligned"] <- "Reads"
	names(aug2020)[names(aug2020) == "Reads_per_Mb"] <- "ReadsPerMb"
	
	aug28<- read.table(files[2],header=T,fill=T)%>% select(Library, Postfiltering_reads_aligned,Background,Reads_per_Mb)
	#au28<- rename(au28, c("Reads"="Postfiltering_reads_aligned" , "Background" =  "Background" , "ReadsPerMb" = "Reads_per_Mb" ))
	names(aug28)[names(aug28) == "Postfiltering_reads_aligned"] <- "Reads"
	names(aug28)[names(aug28) == "Reads_per_Mb"] <- "ReadsPerMb"
	
	nrow(feb25) + nrow(mar6) + nrow(aug2020) + nrow(aug28)
	
	metrics <- rbind(aug28,feb25,mar6,aug2020)
	str(metrics$Library)
	
	return(metrics)
}

bind <- master(args[1],args[2],args[3])


#write.table(bind,"/Users/zeidh/Desktop/unfilt.bed",sep="\t",quote=F,row.names = F,col.names = F)
#write.table(bind2,"/Users/zeidh/Desktop/filt.bed",sep="\t",quote=F,row.names = F,col.names = F)
#args2=c("/Users/zeidh/Desktop/Data/RDATA/ALL_BL/","/Users/zeidh/Desktop/BreakpointerGithub/Breakpoints/","/Users/zeidh/Desktop/Data/METRICS")
#inputDir="/Users/zeidh/Desktop/Data/RDATA/ALL_BL/"

#inputDir="/Users/zeidh/Desktop/Data/RDATA/Aug28-2020/blacklistedForPlotting"
#inputDirectory=list.files(inputDir,full.names = T)
#outputDir="/Users/zeidh/Desktop/BreakpointerGithub/Breakpoints/"
#metricsDir="/Users/zeidh/Desktop/Data/METRICS"

#inputDirectory=list.files(inputDir,full.names = T)
#(file=inputDirectory[2])
#utputDir="/Users/zeidh/Desktop/BreakpointerGithub/Breakpoints/"
#metricsDir="/Users/zeidh/Desktop/Data/METRICS"


								   