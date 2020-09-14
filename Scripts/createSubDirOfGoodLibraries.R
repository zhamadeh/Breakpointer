#setwd("/Users/zeidh/Desktop/BreakpointerGithub/")
args = commandArgs(trailingOnly=TRUE)
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(miceadds)))
args= c("Input/Metrics/Aug28.txt","Input/RData_unfiltered/","Input/RData_good/")
cat("\n")
cat("Filtering out poor quality libraries ....\n\n")

libraries <- read.table(args[1],header=T)
good <- filter(libraries, libraries$Quality=="g")
good$Library <- paste0(good$Library,".trimmed.mdup.bam")

n=0
for (row in 1:nrow(good)){
	n=n+1
	cat(paste0("Adding ",good[row,1]," to good libraries directory ... ",round((n/nrow(good)*100),2),"%\n"))
	filename <- paste0(args[2],good[row,1],".RData")
	file <- load.Rdata2(filename)
	fileDir = paste0(args[3],good[row,1],".RData")
	save(file=fileDir, file)
}

write.table(good$Library,"/Users/zeidh/Desktop/BreakpointerGithub/Output/Good_libraries/August28-2020.txt",sep="\t",quote=F,row.names = F,col.names = F)

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
suppressMessages(t <- group_by(q,Quality,gene) %>% dplyr::summarize(n()))
file <- tail(strsplit(args[1],"/")[[1]],1)
perc <- (sum(filter(as.data.frame(t),Quality=="g")$`n()`)/nrow(q))*100

cat(paste0("\n There were ",nrow(good)," good quality libraries have been filtered from ",file,". There are ",round(perc, digits = 2),"% good libraries and they were moved to ",args[3],"\n\n"))
bl <- c()
for (row in 1:nrow(t)){
	if (t[row,1] != "p" & t[row,1] != "g" ){
		bl <- append(row,bl)
	}
}
t <- t[-c(bl),]
print(as.data.frame(t))
