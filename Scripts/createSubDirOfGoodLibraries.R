#setwd("/Users/zeidh/Desktop/BreakpointerGithub/")
args = commandArgs(trailingOnly=TRUE)
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(miceadds)))
#args= c("Input/Metrics","Input/RData_unfiltered/","Input/RData_good/")
cat("\n")
cat("Filtering out poor quality libraries ....\n\n")

collectMetrics <- function(metricsDir){
	files <- list.files(metricsDir, full.names=T)
	metrics <- data.frame()
	if (length(files)!=0){
		file=files[1]
		for (file in files){
			met <- read.table(file,header=T,fill=T) 
			if ("Reads_aligned_postfiltering" %in% colnames(met)){
				met <- met %>% select(Library, Reads_aligned_postfiltering,Background,Reads_per_Mb,Quality)
				names(met)[names(met) == "Reads_aligned_postfiltering"] <- "Reads"
			}
			else if ("Postfiltering_reads_aligned" %in% colnames(met)){
				met <- met %>% select(Library, Postfiltering_reads_aligned , Background,Reads_per_Mb,Quality)
				names(met)[names(met) == "Postfiltering_reads_aligned"] <- "Reads"
			}
			metrics <- rbind(met,metrics)
			
		}
	}
	return(metrics)
}
libraries <- collectMetrics(args[1])



good <- filter(libraries, libraries$Quality=="g")
good$Library <- paste0(good$Library,".trimmed.mdup.bam")
goodFound <- intersect(paste0(good[,1],".RData"),list.files(args[2]))

n=0

for (row in 1:nrow(good)){
	#row=1
	n=n+1
	if (paste0(good[row,1],".RData") %in% list.files(args[2])){
		cat(paste0("Adding ",good[row,1]," to good libraries directory ... ",round((n/length(goodFound)*100),2),"%\n"))
		filename <- paste0(args[2],good[row,1],".RData")
		file <- load.Rdata2(filename)
		fileDir = paste0(args[3],good[row,1],".RData")
		save(file=fileDir, file)
	}
}

write.table(intersect(paste0(good[,1],".RData"),list.files(args[2])),"Output/Good_libraries/August28-2020.txt",sep="\t",quote=F,row.names = F,col.names = F)

#summary quality stats
q <- libraries
q$library=q$Library

suppressWarnings(q <- q %>% separate(library, c("a","b","c","d","e","f"), "[_-]+"))
q$gene <- "gene"

for (row in 1:nrow(q)){

	for (letter in c("a","b","c","d","e","f")){
		if (is.na(q[row,letter])!=T){
			if (q[row,letter]=="WT" | q[row,letter]=="wt"){
				q[row,"gene"]="WT"
			}
			else if (q[row,letter]=="blm" | q[row,letter]=="BLM" ) {
				q[row,"gene"]="BLM"
			}
			
			else if (q[row,letter]=="RECQL5" | q[row,letter]=="recql5" | q[row,letter]=="RECQ5" | q[row,letter]=="recq5" ) {
				if (q[row,"gene"]=="BLM"){
					q[row,"gene"]="BLM/RECQL5"
				}
				else{
					q[row,"gene"]="RECQL5"
				}
			}
		}
	}
	
}

q <- select(q,-c(a,b,c,d,e,f))
q <- q[paste0(q[,1],".trimmed.mdup.bam.RData") %in% list.files(args[2]),]
#intersect(paste0(q[,1],".trimmed.mdup.bam.RData"),list.files(args[2]))
suppressMessages(t <- group_by(q,Quality,gene) %>% dplyr::summarize(n()))
file <- tail(strsplit(args[1],"/")[[1]],1)

perc <- (sum(filter(as.data.frame(t),Quality=="g")$`n()`)/nrow(q))*100

cat(paste0("\n There were ",nrow(good)," good quality libraries in the metrics file but only ", length(intersect(paste0(good[,1],".RData"),list.files(args[2]))), " were found in the Input folder\n\n"))
cat(paste0("Those ",length(intersect(paste0(good[,1],".RData"),list.files(args[2])))," libraries have been added to ",args[3],"\n\n"))

print(as.data.frame(t))


