#selasphorus Fst scan visualization
#get vcftools windowed fst output
setwd("~/Dropbox/selasphorus/")
library(data.table);library(plyr);library(ggplot2);library(pbapply)
library(magrittr);library(foreach);library(ggrepel);library(cowplot)
cr <- fread("window_stats/ruf_cal_angsd_fstwindow.txt")[,c(2,3,4,5)]
cr$comp <- "rufus x\ncalliope"
rs <- fread("window_stats/ruf_sas_angsd_fstwindow.txt")[,c(2,3,4,5)]
rs$comp <- "rufus x\nsasin"
ew <- fread("window_stats/east_west_angsd_fstwindow.txt")[,c(2,3,4,5)]
ew$comp <- "east x west\nrufus"
ss <- fread("window_stats/sas_sed_angsd_fstwindow.txt")[,c(2,3,4,5)]
ss$comp <- "sasin x\nsedentarius"
rsd <- fread("window_stats/ruf_sed_angsd_fstwindow.txt")[,c(2,3,4,5)]
rsd$comp <- "rufus x\nsedentarius"
fst <- rbind(cr,rs,ew,ss,rsd)
colnames(fst) <- c("contig","window_start","snps","Fst","comp")
fst$comp <- factor(fst$comp,levels=c("rufus x\ncalliope","rufus x\nsedentarius","rufus x\nsasin",
                                     "sasin x\nsedentarius","east x west\nrufus"))
fst <- subset(fst,! contig %in% c("MUGM01000982.1","MUGM01000053.1")) #drop 2 contigs with alignment issues - seem to be paralogs in there


gff <- read.table("~/Dropbox/augustus/Canna_MUGM01_fwd.gff",comment.char = "#",skip = 19,sep="\t",stringsAsFactors = F)
colnames(gff) <- c("contig","whytho","feature","start","end","score","strand","frame","attribute")

#summarize mummer output
#read in and clean up formatting for coords files
files <- list.files("~/Dropbox/mummer_alignments/mummer_out/Canna_MUGM01/coords/",full.names = T)
files <- files[grepl("\\.coords",files)]
i <- 1
for(f in files){
  if(i==1){
    dat <- fread(f,sep=" ",data.table=F)
    i=i+1
  } else {
    tmp <- fread(f,sep=" ",data.table=F)
    dat <- rbind(tmp,dat)
    i=i+1
  }
}
colnames(dat) <- c("refStart","refStop","sep1","qStart","qStop","sep2","refLength","qLength","sep3",
                   "p.identity","sep4","names")
dat$refName <- strsplit(dat$names,"\t") %>% sapply(function(e) unlist(e)[1])
dat$qName <- strsplit(dat$names,"\t") %>% sapply(function(e) unlist(e)[2])
dat <- arrange(dat,refName,refStart)
# contig_lengths <- data.frame(qName=unique(dat$qName),
#                              contig_length=pbsapply(unique(dat$qName),
#                                                     function(e) max(subset(dat,qName==e)$qStop)))
#dat <- merge(dat,contig_lengths,by="qName")
sum <- ddply(dat,
             .(qName,refName),
             summarize,
             totalMatch=sum(qLength),
             refStart=min(refStart))
sum <- arrange(sum,refName,refStart,totalMatch)
sum <- subset(sum,totalMatch>10000)                                                 
sum <- ddply(sum,.(qName),function(e){      #get contigs with hits to only one chromosome > 1000bp
  #a <- subset(e,refName==e$refName[e$totalMatch==max(totalMatch)])
  if(nrow(e)==1){
    e
  }
})

#chromosome order for pretty plots
chr_order <- c("1","1A","1B","2","3","4","4A",as.character(5:28),"Z","M","NA")

#merge mummer info with angsd windowed Fst's
fst <- merge(fst,sum,by.x="contig",by.y="qName",all.x=T,all.y=F)
fst$chr <- gsub("chr","",fst$refName)
fst$chr[!fst$chr %in% chr_order] <- "NA"
fst$chr <- factor(fst$chr,levels=chr_order)
fst <- ddply(fst,.(comp),function(e) {
  e <- arrange(e,chr,refStart)
  e$row <- 1:nrow(e)
  return(e)})
fst <- subset(fst,chr != "M")

#build second dataset for zebra finch chromosome plotting
chr_labels <- ddply(fst,.(chr),summarize,mid=median(row),start=min(row),stop=max(row))
chr_labels$chr <- as.character(chr_labels$chr)
chr_labels$chr[chr_labels$chr %in% as.character(21:28)] <- "21-28"
chr_labels$mid[chr_labels$chr=="21-28"] <- median(chr_labels$mid[chr_labels$chr=="21-28"],na.rm=T)
chr_labels$start[chr_labels$chr=="21-28"] <- min(chr_labels$start[chr_labels$chr=="21-28"],na.rm=T)
chr_labels$stop[chr_labels$chr=="21-28"] <- max(chr_labels$stop[chr_labels$chr=="21-28"],na.rm=T)
chr_labels$start[chr_labels$chr=="1B"] <- 0
chr_labels$stop[chr_labels$chr=="1B"] <- 0
chr_labels$Fst <- -0.1 #y position in the plot
chr_labels <- subset(chr_labels,!is.na(chr) & !duplicated(chr))
chr_labels$comp <- "east x west\nrufus"
chr_labels$comp <- factor(chr_labels$comp,levels=c("rufus x\ncalliope","rufus x\nsedentarius","rufus x\nsasin",
                                                   "sasin x\nsedentarius","east x west\nrufus"))

#line segments for labeling chromosomes
chr_segments <- ddply(fst,.(chr),summarize,mid=median(row),start=min(row),stop=max(row))
chr_segments$chr <- as.character(chr_segments$chr)
chr_segments$chr[chr_segments$chr %in% as.character(21:28)] <- "21-28"
chr_segments$mid[chr_segments$chr=="21-28"] <- median(chr_segments$mid[chr_segments$chr=="21-28"],na.rm=T)
chr_segments$start[chr_segments$chr=="21-28"] <- min(chr_segments$start[chr_segments$chr=="21-28"],na.rm=T)
chr_segments$stop[chr_segments$chr=="21-28"] <- max(chr_segments$stop[chr_segments$chr=="21-28"],na.rm=T)
chr_segments$start[chr_segments$chr=="1B"] <- 0
chr_segments$stop[chr_segments$chr=="1B"] <- 0
chr_segments$Fst <- -0.1
chr_segments <- subset(chr_segments,!is.na(chr) & !duplicated(chr))
chr_segments <- rbind(chr_segments,chr_segments,chr_segments,chr_segments,chr_segments)
chr_segments$comp <- c(rep("east x west\nrufus",nrow(chr_segments)/5),
                       rep("rufus x\nsasin",nrow(chr_segments)/5),
                       rep("rufus x\ncalliope",nrow(chr_segments)/5),
                       rep("sasin x\nsedentarius",nrow(chr_segments)/5),
                       rep("rufus x\nsedentarius",nrow(chr_segments)/5))
chr_segments$comp <- factor(chr_segments$comp,levels=c("rufus x\ncalliope","rufus x\nsedentarius","rufus x\nsasin",
                                                       "sasin x\nsedentarius","east x west\nrufus"))

#set negative Fst's to 0
fst$Fst[fst$Fst<0] <- 0

#rolling means
fst$rollmean <- unlist(dlply(fst,.(comp),function(e) zoo::rollmean(e$Fst,20,fill=c(NA,NA,NA))))

#get outlier windows
outliers <- ddply(fst,.(comp),function(e) subset(e,Fst>quantile(e$Fst,0.995,na.rm=T)))

#blank layer for ylims including labels
ylims <- data.frame(min=c(0,0,0,0,-0.6),max=c(1,.75,.75,.5,.25),comp=c("rufus x\ncalliope","rufus x\nsasin","rufus x\nsedentarius",
                                                           "sasin x\nsedentarius","east x west\nrufus"),x=c(1,1,1,1,1))

#plot multi fst comps
png(file="~/Dropbox/selasphorus/ms/fig/angsd_fst_fig.png",res=600,width=6.5,height=4,units="in")
#pdf(file="~/Dropbox/selasphorus/ms/fig/angsd_fst_fig.pdf",width=6.5,height=4,useDingbats = F)
ggplot(data=fst,aes(x=row,y=Fst,col=contig))+theme_classic()+facet_grid(comp~.,scales="free")+
  theme(text=element_text(size=8),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x=element_blank(),
        axis.title.x=element_blank(),
        strip.text=element_text(size=7),
        panel.grid.major.y=element_line(color="grey",size = .25),
        legend.position="none")+
  scale_y_continuous(breaks=c(0,.25,.5,.75,1))+
  scale_color_manual(values=rep(c("grey55","grey80"),length(unique(fst$contig))/2+1))+
  geom_point(size=0.4,shape=21,stroke=0.2)+
  geom_line(aes(y=rollmean),lwd=0.15,col="black")+
  geom_segment(data=chr_segments,aes(x=start+50,xend=stop-50,y=Fst,yend=Fst,col=NA),col="black")+
  geom_point(data=outliers,size=0.4,shape=21,stroke=0.4,col="red")+
  geom_text_repel(data=chr_labels,aes(label=chr,x=mid,y=Fst,col=NA),force=2,ylim=c(-.45,-0.1),
                  col="black",size=2,angle=0,direction="y",box.padding = 0.15,
                  segment.size=0.2)+
  geom_hline(aes(yintercept=0))+
  geom_blank(data=ylims,aes(x=x,y=min,comp=comp,col=NA))+
  geom_blank(data=ylims,aes(x=x,y=max,comp=comp,col=NA))
dev.off()

#z vs autosome fst comps
fst$z <- fst$chr=="Z"
zs <- ddply(fst,.(comp),function(e){
  mean_fst_aut <- mean(e$Fst[e$z==FALSE])
  mean_fst_z <- mean(e$Fst[e$z==TRUE])
  ratio <- mean_fst_z/mean_fst_aut
  outliers <- subset(e,Fst>quantile(e$Fst,0.99))
  nz <- nrow(subset(outliers,chr=="Z"))
  naut <- nrow(subset(outliers,chr!="Z"))
  z_outlier_proportion <- nz/(nz+naut)
  data.frame(mean_fst_z,mean_fst_aut,ratio,z_outlier_proportion)
})
ggplot(data=fst,aes(x=Fst,fill=z))+facet_wrap(~comp,scales="free")+
  theme(legend.position = c(0.8,0.25))+
  #xlim(0,0.5)+
  geom_density(alpha=0.5)



#pairwise correlation in Fst
#get n standard deviations from mean values
fst <- arrange(fst,comp)
zfst <- c()
for(i in unique(fst$comp)){
  a <- subset(fst,comp==i)
  b <- (a$Fst-mean(a$Fst,na.rm=T))/sd(a$Fst,na.rm=T)
  zfst <- append(zfst,b)
}
fst$zfst <- zfst

shared <- merge(rc,rs,by=c("contig","window_start"),all.x=F,all.y=F)[,c("contig","window_start")] %>% 
  merge(rsd,by=c("contig","window_start"),all.x=F,all.y=F) %>% .[,c("contig","window_start")] %>% 
  merge(ss,by=c("contig","window_start"),all.x=F,all.y=F) %>% .[,c("contig","window_start")] %>% 
  merge(ew,by=c("contig","window_start"),all.x=F,all.y=F) %>% .[,c("contig","window_start")]

zfst <- subset(fst,paste(fst$contig,fst$window_start) %in% paste(shared$contig,shared$window_start))

rc <- subset(zfst,comp=="rufus x\ncalliope")
rs <- subset(zfst,comp=="rufus x\nsasin")
rsd <- subset(zfst,comp=="rufus x\nsedentarius")
ss <- subset(zfst,comp=="sasin x\nsedentarius")
ew <- subset(zfst,comp=="east x west\nrufus")

zfst <- data.frame(shared,rc=rc$Fst,rs=rs$Fst,rsd=rsd$Fst,ss=ss$Fst,ew=ew$Fst)
plot_data <- data.frame(x=c(zfst$rc,zfst$rc,zfst$rc,zfst$rc,zfst$rs,zfst$rs,zfst$rs,zfst$rsd,zfst$rsd,zfst$ss),
                        y=c(zfst$rs,zfst$rsd,zfst$ss,zfst$ew,zfst$rsd,zfst$ss,zfst$ew,zfst$ss,zfst$ew,zfst$ew),
                        xsp=c(rep("rufus x\n calliope",nrow(rc)*4),
                              rep("rufus x\n sasin",nrow(rc)*3),
                              rep("rufus x\n sedentarius",nrow(rc)*2),
                              rep("sasin x\n sedentarius",nrow(rc))),
                        ysp=c(rep("rufus x\n sasin",nrow(rc)),
                              rep("rufus x\n sedentarius",nrow(rc)),
                              rep("sasin x\n sedentarius",nrow(rc)),
                              rep("east x\n west rufus",nrow(rc)),
                              rep("rufus x\n sedentarius",nrow(rc)),
                              rep("sasin x\n sedentarius",nrow(rc)),
                              rep("east x\n west rufus",nrow(rc)),
                              rep("sasin x\n sedentarius",nrow(rc)),
                              rep("east x\n west rufus",nrow(rc)),
                              rep("east x\n west rufus",nrow(rc))))
plot_data$ysp <- factor(plot_data$ysp,levels=c("rufus x\n sasin","rufus x\n sedentarius","sasin x\n sedentarius","east x\n west rufus"))

p <- ggplot(data=plot_data,aes(x=x,y=y))+facet_grid(ysp~xsp,switch="both",scales="free")+
  theme(strip.background = element_blank(),
        legend.position=c(0.83,0.76),
        legend.text = element_text(size=7),
        legend.title = element_text(size=8),
        strip.text=element_text(size=8),
        strip.placement = 'outside',
        axis.line=element_blank(),
        axis.title=element_blank(),
        panel.background = element_rect(color="black",size = 1),
        axis.text=element_text(size=7),
        axis.text.x=element_text(size=7,angle=45,hjust=1))+
  xlab("z(Fst)")+ylab("z(Fst)")+
  scale_fill_distiller(palette = "RdYlBu")+
  stat_bin_hex(aes(fill=log(..count..)))+
  geom_smooth(method="lm",col="black",lwd=0.75)
p <- p+guides(fill=guide_colorbar(barheight = unit(23,"mm"),barwidth = unit(5,"mm")))
pdf("ms/fig/angsd_fst_matrix.pdf",width=6.5,height=6)
print(p)
dev.off()




