countedSCEs$gene <- as.factor(countedSCEs$gene)
countedSCEs$Library <- as.factor(countedSCEs$Library)
group_by(countedSCEs,gene) 

library(plyr)

#sces per chroomosome
write.table(as.data.frame(group_by(countedSCEs,seqnames) %>% summarize(n())),"Users/zeidh/Desktop/perChrom.txt",quote=F,row.names = F,col.names = F,sep="\t")
b <- as.data.frame(countedSCEs %>% filter(gene=="BLM")  %>% group_by(seqnames) %>% dplyr::summarize(n()))
br <- as.data.frame(countedSCEs %>% filter(gene=="BLM/RECQL5")  %>% group_by(seqnames)%>% dplyr::summarize(n()))
r <- as.data.frame(countedSCEs %>% filter(gene=="RECQL5")  %>% group_by(seqnames)%>% dplyr::summarize(n()))
w <- as.data.frame(countedSCEs %>% filter(gene=="WT")  %>% group_by(seqnames)%>% dplyr::summarize(n()))

b <- (as.data.frame(countedSCEs %>% filter(gene=="BLM")  %>% droplevels()))# group_by(seqnames)%>% summarize(n()))
length(levels(b$Library))
br <- (as.data.frame(countedSCEs %>% filter(gene=="BLM/RECQL5")  %>% droplevels()))# group_by(seqnames)%>% summarize(n()))
length(levels(br$Library))
r <- (as.data.frame(countedSCEs %>% filter(gene=="RECQL5")  %>% droplevels()))# group_by(seqnames)%>% summarize(n()))
length(levels(r$Library))
w <- (as.data.frame(countedSCEs %>% filter(gene=="WT")  %>% droplevels()))# group_by(seqnames)%>% summarize(n()))
length(levels(w$Library))


lengths <- read.table("Users/zeidh/Desktop/chromLen.txt",header=T)
tidy <- gather(lengths, gene,sce_per_chr,BLM:WT)

ggplot(tidy) + geom_point(aes(Chromosome.1,sce_per_chr,group=gene,color=gene))+
	geom_smooth(aes(Chromosome.1,sce_per_chr,group=gene,color=gene),se=F,method="lm")+
	theme_classic(base_size = 18) +
	ylab("SCEs/chr/lib")+
	xlab("Chromosome Length")


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
test<-rename(test,c("n()"="sces"))


ggplot(test) + geom_jitter(aes(gene,sces, color=gene))+ geom_boxplot(aes(gene,sces),width=0.1,coef = 5) 

ggplot(countedSCEs) + geom_smooth(aes(ReadsPerMb, width,color=gene),se=F) +
	scale_y_log10() +
	theme_classic() +
	theme(text=element_text(size=15)) +
	geom_vline(xintercept=150, linetype="dashed", color = "red") 

my_comparisons <- list( c("WT", "RECQL5"), c("WT", "BLM/RECQL5"), c("WT", "BLM") )
ggboxplot(test, x = "gene", y = "sces",
		  color = "black",  add = "jitter",width=0.25, add.params = list(color = "gene"),
		  xlab="Gene",ylab="SCEs/library") +
	stat_compare_means(comparisons = my_comparisons,label.y = c(31, 34, 37)) +
	stat_compare_means(label = "p.signif", method = "t.test",
					   ref.group = "WT") 


b<- filter(countedSCEs, gene=="BLM")
mean(b$width)
median(b$width)
b<- filter(countedSCEs, gene=="RECQL5")
mean(b$width)
median(b$width)
b<- filter(countedSCEs, gene=="BLM/RECQL5")
mean(b$width)
median(b$width)
b<- filter(countedSCEs, gene=="WT")
mean(countedSCEs$width)
median(countedSCEs$width)


library(tidyverse)
library(plyr)
nonWidthSizeCutOff = countedSCEs

#plot here
ggplot(countedSCEs) + stat_ecdf(aes(width,color=gene)) +
	scale_x_log10() +
	theme_classic() +
	ylab("SCEs Mapped (%)") +
	xlab("Resolution") +
	annotation_logticks(sides = "b") +
	theme(text = element_text(size=15))+
	geom_density(aes(width),adjust=adj,size=1.1)+
	geom_vline(xintercept=median(countedSCEs$width), linetype="dashed", color = "red") +
	geom_text(aes(x=5000, label=paste0("Median\n","22.6Kb"), y=0.8))
mean<-mean(countedSCEs$width)
sort(countedSCEs$width)
ggsave("Plots/truePosCumlative.png")
ggplot(q) + stat_ecdf(aes(width,color=gene)) +
	scale_x_log10() +
	theme_classic() +
	ylab("SCEs Mapped (%)") +
	xlab("Resolution") +
	annotation_logticks(sides = "b") +
	theme(text = element_text(size=20))+
	geom_density(aes(width,color=gene),adjust=adj, alpha=0.7)

smooth_ecd = function(adj = 1) {
	
	# Fake data
	dat = data.frame(x=countedSCEs$width)
	
	# Extend range of density estimate beyond data
	e =  0.3 * diff(range(dat$x))
	
	# Kernel density estimate of fake data
	#adj=1
	dens = density(dat$x, adjust=adj, from=min(dat$x)-e, to=max(dat$x) +e)
	dens = data.frame(x=dens$x, y=dens$y)
	
	# Plot kernel density (blue), ecdf (red) and smoothed ecdf (black)
	ggplot(dat, aes(x)) + 
		geom_density(adjust=adj, colour="blue", alpha=0.7) +
		geom_line(data=dens, aes(x=x, y=cumsum(y)/sum(y)), size=0.7, colour='grey30') +
		stat_ecdf(colour="red", size=0.6, alpha=0.6) +
		theme_classic() +
		labs(title=paste0("adj=",adj))+
		scale_x_log10()
}
smooth_ecd(adj=1)
smooth_ecd(adj=0.3)
smooth_ecd(adj=0.1)




bl <- c()
r=countedSCEs
n=0
for (row in 1:nrow(countedSCEs)){
	if(as.numeric(as.character(countedSCEs[row,3])) - as.numeric(as.character(countedSCEs[row,2])) <= 1000000){
		print("small")
	}
	else if (as.numeric(as.character(countedSCEs[row,3])) - as.numeric(as.character(countedSCEs[row,2])) >= 1000000){
		print("big")
		bl <- append(bl,row)
		n=n+1
	}
}
#nrow(countedSCEs)
n
mean(countedSCEs$width)
mean(q$width)
if (length(bl)!=0){
	q <- r[-c(bl),]
	
}

countedSCEs=q


