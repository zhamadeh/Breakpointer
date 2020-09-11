#setwd("/")
library(plyr)
library(tidyverse)
args = commandArgs(trailingOnly=TRUE)
args=c("Output/Breakpoints/Aug28_breakpoints.txt")

countedSCEs <- read.table(args[1],header=T)
countedSCEs$gene <- as.factor(countedSCEs$gene)
countedSCEs$Library <- as.factor(countedSCEs$Library)




#sces per chroomosome
all <- as.data.frame(countedSCEs %>% group_by(seqnames) %>% dplyr::summarize(ALL=n()))
b <- as.data.frame(countedSCEs %>% filter(gene=="BLM")  %>% group_by(seqnames) %>% dplyr::summarize(BLM=n()))
br <- as.data.frame(countedSCEs %>% filter(gene=="BLM/RECQL5")  %>% group_by(seqnames)%>% dplyr::summarize("BLM/RECQL5"= n()))
r <- as.data.frame(countedSCEs %>% filter(gene=="RECQL5")  %>% group_by(seqnames)%>% dplyr::summarize(RECQL5=n()))
w <- as.data.frame(countedSCEs %>% filter(gene=="WT")  %>% group_by(seqnames)%>% dplyr::summarize(WT=n()))
byChr <- merge(b,br,by="seqnames")
byChr <- merge(byChr,r,by="seqnames",all=T)
byChr <- merge(byChr,w,by="seqnames",all=T)
byChr <- merge(byChr,all,by="seqnames",all=T)
byChr[is.na(byChr)] <- 0
write.table(byChr,"Output/Tables/perChrom.txt",quote=F,row.names = F,col.names = T,sep="\t")


numOfLibsPerGene <- data.frame(gene=character(),n=numeric())
b <- (as.data.frame(countedSCEs %>% filter(gene=="BLM")  %>% droplevels()))# group_by(seqnames)%>% summarize(n()))
numOfLibsPerGene= add_row(numOfLibsPerGene, gene="BLM",n=as.numeric(length(levels(b$Library))))
br <- (as.data.frame(countedSCEs %>% filter(gene=="BLM/RECQL5")  %>% droplevels()))# group_by(seqnames)%>% summarize(n()))
numOfLibsPerGene= add_row(numOfLibsPerGene, gene="BLM.RECQL5",n=as.numeric(length(levels(br$Library))))
r <- (as.data.frame(countedSCEs %>% filter(gene=="RECQL5")  %>% droplevels()))# group_by(seqnames)%>% summarize(n()))
numOfLibsPerGene= add_row(numOfLibsPerGene, gene="RECQL5",n=as.numeric(length(levels(r$Library))))
w <- (as.data.frame(countedSCEs %>% filter(gene=="WT")  %>% droplevels()))# group_by(seqnames)%>% summarize(n()))
numOfLibsPerGene= add_row(numOfLibsPerGene, gene="WT",n=as.numeric(length(levels(w$Library))))
numOfLibsPerGene= add_row(numOfLibsPerGene, gene="ALL",n=as.numeric(length(levels(countedSCEs$Library))))
write.table(numOfLibsPerGene,"Output/Tables/numOfLibsPerGene.txt",quote=F,row.names = F,col.names = F,sep="\t")


numOfLibsPerGene=read.table("Output/Tables/numOfLibsPerGene.txt",header=F)
lengths <- read.table("Output/Tables/chrLengths.txt",header=T)
byChr<- read.table("Output/Tables/perChrom.txt",header=T)

for (i in c("BLM","BLM.RECQL5","RECQL5","WT","ALL")){
	#print(byChr[,i])
	for (j in 1:nrow(numOfLibsPerGene)){
		if (numOfLibsPerGene[j,1]==i){
			#print(numOfLibsPerGene[j,2])
			byChr[,i] = byChr[,i]/ numOfLibsPerGene[j,2]
		}
	}
}

scePerChrPerGeneVsLength <- merge(lengths,byChr,by.x="Chromosome",by.y="seqnames")
tidy <- gather(scePerChrPerGeneVsLength, gene,sce_per_chr,BLM:ALL)
all <- filter(tidy,gene=="ALL")
tidy <- filter(tidy,gene!="ALL")

ggplot(tidy) + geom_point(aes(Length,sce_per_chr,group=gene,color=gene))+
	geom_smooth(aes(Length,sce_per_chr,group=gene,color=gene),se=F,method="lm")+
	theme_classic(base_size = 18) +
	ylab("SCEs/chr/lib")+
	xlab("Chromosome Length") +
	ggsave("Output/Plots/scePerChrPerGeneVsLength.png")

ggplot(tidy) + geom_point(aes(Length,sce_per_chr,group=gene,color=gene))+
	geom_smooth(all,mapping=aes(Length,sce_per_chr),color="black",method="lm",size=2)+
	theme_classic(base_size = 18) +
	ylab("SCEs/chr/lib")+
	xlab("Chromosome Length") +
	ggsave("Output/Plots/scePerChrPerGeneVsLength_ALL.png")




test<-as.data.frame(countedSCEs %>%
	group_by(Library) %>%
	dplyr::summarize(n()))
test$gene <- "gene"
test$library <- test$Library
test <- test %>% separate(Library, c("a","b","c","d","e","f"), "[_-]+")

for (row in 1:nrow(test)){
	for (letter in c("a","b","c","d","e","f")){
		#print(test[1,letter])
		
		if (test[row,letter]=="WT" | test[row,letter]=="wt"){
			test[row,"gene"]="WT"
		}
		else if (test[row,letter]=="blm" | test[row,letter]=="BLM" ) {
			test[row,"gene"]="BLM"
		}
		
		else if (test[row,letter]=="RECQL5" | test[row,letter]=="recql5" | test[row,letter]=="RECQ5" | test[row,letter]=="recq5" ) {
			if (test[row,"gene"]=="BLM"){
				test[row,"gene"]="BLM/RECQL5"
			}
			else{
				test[row,"gene"]="RECQL5"
			}
		}
	}
	
}

test <- select(test,c("n()","gene"))
test<-dplyr::rename(test,c("sces"="n()"))

ggplot(test) + geom_jitter(aes(gene,sces, color=gene))+ geom_boxplot(aes(gene,sces),width=0.1,coef = 5) +
	theme_classic()+
	theme(text=element_text(size=15)) +
	ggsave("Output/Plots/SCEperGene.png")

#my_comparisons <- list( c("WT", "RECQL5"), c("WT", "BLM/RECQL5"), c("WT", "BLM") )
#ggboxplot(test, x = "gene", y = "sces",
#		  color = "black",  add = "jitter",width=0.25, add.params = list(color = "gene"),
#		  xlab="Gene",ylab="SCEs/library") +
#	stat_compare_means(comparisons = my_comparisons,label.y = c(31, 34, 37)) +
#	stat_compare_means(label = "p.signif", method = "t.test",
#					   ref.group = "WT") 






ggplot(countedSCEs) + geom_smooth(aes(Reads_per_Mb, width,color=gene),se=F) +
	scale_y_log10() +
	theme_classic() +
	theme(text=element_text(size=15))+
	geom_hline(yintercept=10000, linetype="dashed", color = "red") +
	ggsave("Output/Plots/resolutionVsDepth.png")



sce_summary=data.frame(gene=character(),SCE=numeric(),mean_resolution=numeric(),median_resolution=numeric())
b<- filter(countedSCEs, gene=="BLM")
sce_summary<- add_row(sce_summary,gene="BLM",SCE=nrow(b),mean_resolution=mean(b$width),median_resolution=median(b$width))
r<- filter(countedSCEs, gene=="RECQL5")
sce_summary<- add_row(sce_summary,gene="RECQL5",SCE=nrow(r),mean_resolution=mean(r$width),median_resolution=median(r$width))
br<- filter(countedSCEs, gene=="BLM/RECQL5")
sce_summary<- add_row(sce_summary,gene="BLM/RECQL5",SCE=nrow(br),mean_resolution=mean(br$width),median_resolution=median(br$width))
w<- filter(countedSCEs, gene=="WT")
sce_summary<- add_row(sce_summary,gene="WT",SCE=nrow(w),mean_resolution=mean(w$width),median_resolution=median(w$width))
write.table(sce_summary,"Output/Tables/SCE_summary.txt",quote=F,row.names = F,col.names = T,sep="\t")




#plot here
ggplot(countedSCEs) + stat_ecdf(aes(width,color=gene)) +
	scale_x_log10() +
	theme_classic() +
	ylab("SCEs Mapped (%)") +
	xlab("Resolution") +
	annotation_logticks(sides = "b") +
	theme(text = element_text(size=15))+
	geom_density(aes(width),size=1.1)+
	geom_vline(xintercept=median(countedSCEs$width), linetype="dashed", color = "red") +
	geom_text(aes(x=5000, label=paste0("Median\n",median(countedSCEs$width)," bp"), y=0.8))  +
	ggsave("Output/Plots/breakpointResolution.png")


