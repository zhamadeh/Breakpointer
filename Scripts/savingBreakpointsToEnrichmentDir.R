

write.table(filter(countedSCEs,gene=="BLM"),"/Users/zeidh/Desktop/EnricherGithub/BED-Input/blm.bed",quote=F,col.names = F,row.names = F,sep="\t")
write.table(filter(countedSCEs,gene=="RECQL5"),"/Users/zeidh/Desktop/EnricherGithub/BED-Input/recql5.bed",quote=F,col.names = F,row.names = F,sep="\t")
write.table(filter(countedSCEs,gene=="BLM/RECQL5"),"/Users/zeidh/Desktop/EnricherGithub/BED-Input/blm-recql5.bed",quote=F,col.names = F,row.names = F,sep="\t")
write.table(filter(countedSCEs,gene=="WT"),"/Users/zeidh/Desktop/EnricherGithub/BED-Input/wt.bed",quote=F,col.names = F,row.names = F,sep="\t")
