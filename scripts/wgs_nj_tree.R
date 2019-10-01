library(vcfR);library(ape);library(progress);library(data.table)
library(plyr);library(ggplot2);library(ggrepel);library(zoo);library(cowplot)
setwd("~/selasphorus_evolution/")
vcf <- read.vcfR("~/Dropbox/selasphorus/called_snps/selasphorus_het_depth_mq30_biallel_nomd.recode.vcf.gz")

#concatenated NJ tree
dna <- vcfR2DNAbin(vcf,unphased_as_NA = F,consensus = T,extract.haps = F)
dist <- dist.dna(dna,model="K80")
nj <- bionj(dist)
nj <- root(nj,c("Scal1","Scal2","Scal3","Scal4","Scal5","Scal6","Scal7"))
plot(nj,cex=0.5)

#get local nj trees in nonoverlapping windows
trees <- list()
pb <- progress_bar$new(total = length(unique(vcf@fix[,1])))
for(contig in unique(vcf@fix[,1])){
  a <- vcf[vcf@fix[,1]==contig]
  start <- 1;step <- 2e4
  if(max(as.numeric(a@fix[,2]))>step){ #for contigs longer than step
    for(i in seq(step,max(as.numeric(a@fix[,2])),step)){
      b <- a[as.numeric(a@fix[,2])>start & as.numeric(a@fix[,2])<i]
      nsnps <- nrow(b@gt)
      c <- vcfR2DNAbin(b,unphased_as_NA = F,consensus = T,extract.haps = F,verbose=F)
      dist <- dist.dna(c,model="K80")
      if(length(c)>=20 & !is.infinite(max(dist)) & !is.na(max(dist))){ #require >50 SNPs
        nj <- nj(dist)
        write.tree(nj,paste0("analysis/twisst/20kb_trees.tre"),append=T)
        write(x=paste0(contig,"\t",start,"\t",i),file="analysis/twisst/20kb_tree_regions.txt",append=T)
        start <- start+step
      }
    }
  } else { #for single-window contigs
    start <- 1;i <- step
    b <- a[as.numeric(a@fix[,2])>start & as.numeric(a@fix[,2])<i]
    nsnps <- nrow(b@gt)
    c <- vcfR2DNAbin(b,unphased_as_NA = F,consensus = T,extract.haps = F,verbose=F)
    dist <- dist.dna(c,model="K80")
    if(length(c)>=20 & !is.infinite(max(dist)) & !is.na(max(dist))){
      nj <- nj(dist)
      write.tree(nj,paste0("analysis/twisst/20kb_trees.tre"),append=T)
      write(x=paste0(contig,"\t",start,"\t",i),file="analysis/twisst/20kb_tree_regions.txt",append=T)
      start <- start+step
    }
  }
  pb$tick()
}

#run twisst 
system(paste0("cd ~/selasphorus_evolution/analysis/twisst/;\
              /anaconda3/bin/python2 twisst.py -t 50kb_trees.tre -w 50kb.weights.csv -g rufus -g sasin -g sedentarius -g calliope --method complete --groupsFile pop_groups.tsv"))


################ plotting #################
tw <- fread("analysis/twisst/50kb.weights.csv")
regions <- fread("analysis/twisst/50kb_tree_regions.txt")
tw <- cbind(tw,regions)

chralign <- fread("data/wgs/CannaMUGM01_to_Tgut2_mummer_summary.txt")

tw <- merge(tw,chralign,by.x="V1",by.y="qName") %>% data.frame()
tw$topo1prop <- tw$topo1/2800
tw$topo2prop <- tw$topo2/2800
tw$topo3prop <- tw$topo3/2800

#chromosome order for pretty plots
chr_order <- c("1","1A","1B","2","3","4","4A",as.character(5:28),"Z","M","NA")
tw$chr <- factor(gsub("chr","",tw$refName),levels=chr_order)
tw <- arrange(tw,chr,refStart,V2)
tw$x <- 1:nrow(tw)
tw$rolltopo1 <- rollmean(tw$topo1prop,40,fill = NA)
tw$rolltopo2 <- rollmean(tw$topo2prop,40,fill = NA)
tw$rolltopo3 <- rollmean(tw$topo3prop,40,fill = NA)
mtw <- tw[,c("V1","V2","refName","refStart","topo1prop","topo2prop","topo3prop","chr","x","rolltopo1","rolltopo2","rolltopo3")]
mtw <- melt(mtw,id.vars=c("V1","V2","refName","refStart","chr","x"))
mtw$topology <- NA
mtw$topology[mtw$variable %in% c("rolltopo1","topo1prop")] <- " rufus    sasin   sedentarius"
mtw$topology[mtw$variable %in% c("rolltopo2","topo2prop")] <- " rufus sedentarius   sasin"
mtw$topology[mtw$variable %in% c("rolltopo3","topo3prop")] <- " rufus     sasin  sedentarius"
mtw$topology <- factor(mtw$topology,levels=c(" rufus    sasin   sedentarius"," rufus sedentarius   sasin"," rufus     sasin  sedentarius"))

chr_labels <- ddply(mtw,.(chr),summarize,mid=median(x),start=min(x),stop=max(x))
chr_labels$chr <- as.character(chr_labels$chr)
chr_labels$chr[chr_labels$chr %in% as.character(21:28)] <- "21-28"
chr_labels$mid[chr_labels$chr=="21-28"] <- median(chr_labels$mid[chr_labels$chr=="21-28"],na.rm=T)
chr_labels$start[chr_labels$chr=="21-28"] <- min(chr_labels$start[chr_labels$chr=="21-28"],na.rm=T)
chr_labels$stop[chr_labels$chr=="21-28"] <- max(chr_labels$stop[chr_labels$chr=="21-28"],na.rm=T)
chr_labels$start[chr_labels$chr=="1B"] <- 0
chr_labels$stop[chr_labels$chr=="1B"] <- 0
chr_labels$value <- -0.1 #y position in the plot
chr_labels <- subset(chr_labels,!is.na(chr) & !duplicated(chr))

t1 <- function(){
        treestr <- "((rufus,sasin_sasin),sasin_sedentarius);"
        tree <- read.tree(text=treestr)
        plot.phylo(tree,show.tip.label = F,direction="upwards",edge.color="orangered")
}
t2 <- function(){
  treestr <- "((rufus,sasin_sasin),sasin_sedentarius);"
  tree <- read.tree(text=treestr)
  plot.phylo(tree,show.tip.label = F,direction="upwards",edge.color="gold")
}
t3 <- function(){
  treestr <- "(rufus,(sasin_sasin,sasin_sedentarius));"
  tree <- read.tree(text=treestr)
  plot.phylo(tree,show.tip.label = F,direction="upwards",edge.color="steelblue2")
}

png("fig/nj_tree_scan_50kb.png",width=6.5,height=2,res=600,units = "in")
a <- ggplot(data=mtw,aes(x=x,y=value,color=topology))+
  theme_classic()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=8),
        axis.text.y=element_text(size=6),
        legend.title=element_text(size=8),
        legend.text=element_text(size=6),
        legend.position = "bottom",
        legend.box.margin = margin(-20,0,0,0),
        legend.key.size = unit(0,"mm"),
        legend.background = element_blank())+
  ylab("proportion of subtrees")+
  ylim(-0.4,1)+
  scale_color_manual(values = c("orangered","gold","steelblue2"))+
  #facet_wrap(~topology,ncol=1)+
  geom_step(data=subset(mtw,variable %in% c("topo1prop","topo2prop","topo3prop")),lwd=0.1,alpha=0.5)+
  geom_step(data=subset(mtw,variable %in% c("rolltopo1","rolltopo2","rolltopo3")),lwd=0.4)+
  #geom_smooth(se=F,method="loess",span=0.05,lwd=0.5)+
  geom_segment(data=chr_labels,aes(x=start+10,xend=stop-10,y=value,yend=value,col=NA),col="black")+
  geom_text_repel(data=chr_labels,aes(label=chr,x=mid,y=value,col=NA),force=2,ylim=c(-.45,-0.1),
                  col="black",size=2,angle=0,direction="y",box.padding = 0.15,
                  segment.size=0.2)
ggdraw()+
  draw_plot(a,0,.05,1,0.95)+
  draw_plot(t1,0.32,0,0.14,0.15)+
  draw_plot(t2,0.49,0,0.14,0.15)+
  draw_plot(t3,0.65,0,0.14,0.15)

dev.off()


#get stats on top topology by window
toptree <- c();maxsupport <- c()
for(i in 1:nrow(tw)){
  row <- tw[i,]
  maxtopo=max(row$topo1,row$topo2,row$topo3)
  toptree[i] <- c("topo1","topo2","topo3")[row[,c("topo1","topo2","topo3")]==maxtopo]
  maxsupport[i] <- maxtopo
}
summary(factor(toptree[tw$chr=="Z"]))/nrow(tw[tw$chr=="Z",])
summary(factor(toptree))/nrow(tw)

###consensus tree
# trees <- read.tree("analysis/twisst/50kb_trees.tre")
# consensus_tree <- root(consensus(trees,p=0.2),outgroup=grep("cal",trees[[1]]$tip.label,value=T))
# plot(consensus_tree,cex=0.5)
