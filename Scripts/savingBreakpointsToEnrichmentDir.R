args = commandArgs(trailingOnly=TRUE)


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
nrow(countedSCEs)
n

if (length(bl)!=0){
	q <- r[-c(bl),]
}

countedSCEs=q

write.table(filter(countedSCEs,gene=="BLM"),"Output/BED_enrichment_analysis/blm.bed",quote=F,col.names = F,row.names = F,sep="\t")
write.table(filter(countedSCEs,gene=="RECQL5"),"Output/BED_enrichment_analysis/recql5.bed",quote=F,col.names = F,row.names = F,sep="\t")
write.table(filter(countedSCEs,gene=="BLM/RECQL5"),"Output/BED_enrichment_analysis/blm-recql5.bed",quote=F,col.names = F,row.names = F,sep="\t")
write.table(filter(countedSCEs,gene=="WT"),"Output/BED_enrichment_analysis/wt.bed",quote=F,col.names = F,row.names = F,sep="\t")
