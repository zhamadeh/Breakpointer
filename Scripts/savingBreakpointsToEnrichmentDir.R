args = commandArgs(trailingOnly=TRUE)
args=c("Output/Breakpoints/Aug28_breakpoints.txt")

countedSCEs <- read.table(args[1],header=T)

bl <- c()
r=countedSCEs
n=0
for (row in 1:nrow(countedSCEs)){
	if(as.numeric(as.character(countedSCEs[row,3])) - as.numeric(as.character(countedSCEs[row,2])) <= 10000){
		#print("small")
	}
	else if (as.numeric(as.character(countedSCEs[row,3])) - as.numeric(as.character(countedSCEs[row,2])) >= 10000){
		#print("big")
		bl <- append(bl,row)
		n=n+1
	}
}

cat("\n\n There were ",n," breakpoints with an interval greater than 10kb, leaving behind ", (nrow(countedSCEs)-n)," out of",nrow(countedSCEs),"\n\n")

q <- as.data.frame(r[-c(bl),])
countedSCEs=q

cat(nrow(dplyr::filter(countedSCEs,countedSCEs$gene=="BLM"))," ... BLM \n")
cat(nrow(dplyr::filter(countedSCEs,countedSCEs$gene=="BLM/RECQL5"))," ... BLM/RECQL5 \n")
cat(nrow(dplyr::filter(countedSCEs,countedSCEs$gene=="RECQL5"))," ... RECQL5 \n")
cat(nrow(dplyr::filter(countedSCEs,countedSCEs$gene=="WT"))," ... WT \n")

write.table(dplyr::filter(countedSCEs,countedSCEs$gene=="BLM"),"Output/BED_enrichment_analysis/blm.bed",quote=F,col.names = F,row.names = F,sep="\t")
write.table(dplyr::filter(countedSCEs,countedSCEs$gene=="RECQL5"),"Output/BED_enrichment_analysis/recql5.bed",quote=F,col.names = F,row.names = F,sep="\t")
write.table(dplyr::filter(countedSCEs,countedSCEs$gene=="BLM/RECQL5"),"Output/BED_enrichment_analysis/blm-recql5.bed",quote=F,col.names = F,row.names = F,sep="\t")
write.table(dplyr::filter(countedSCEs,countedSCEs$gene=="WT"),"Output/BED_enrichment_analysis/wt.bed",quote=F,col.names = F,row.names = F,sep="\t")
