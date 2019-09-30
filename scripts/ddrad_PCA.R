library(data.table);library(ggplot2);library(ggridges)
library(vcfR);library(pegas);library(plyr);library(magrittr);library(raster);library(viridis)
library(RColorBrewer);library(foreach);library(doMC);library(hierfstat);library(cowplot)
registerDoMC(cores=8)
setwd("~/selasphorus_evolution/")

#genetic data
vcf <- read.vcfR("data/ddrad/ddRAD_rufus+sasin+calliope_populations.snps.vcf")
dna <- vcfR2DNAbin(vcf,unphased_as_NA = F,consensus = T,extract.haps = F)
#str <- fread("~/Dropbox/ruhu/rad/stacks/outfiles/ruhu_alhu/batch_1.structure.tsv")
rufous <- dna[-grep("sasin|cal",rownames(dna)),]
ru <- DNAbin2genind(dna)
ru2 <- ru[grepl("OR|WA|AK|MT|CA|NM|UT|TX|AZ",rownames(ru@tab))]
#ru2 <- ru2[!grepl("AK5|AK1|AK3|MT1",rownames(ru2@tab))]
ru3 <- ru[grepl("OR|WA|AK|MT|CA|NM|UT|TX|AZ|sasin",rownames(ru@tab))]
#ru3 <- ru3[!grepl("AK5|AK1|AK3",rownames(ru3@tab))]
a <- fread("data/ddrad/ruf_sas_mac3.str")

#specimen data
dat <- fread("data/samples/ddrad_sample_data.csv",data.table=F)
dat$loc_clust <- dat[,c("Longitude","Latitude")] %>% SpatialPoints() %>% 
  sp::spDists() %>% as.dist() %>% hclust() %>% cutree(h=2)
dat$tmpclust <- gsub("[[:digit:]]","",dat$SampleID)
dat$tmpclust <- mapvalues(factor(dat$tmpclust),
                          from=c("AK","AZ","CA","LOU","MT","NM","OR","TX","UT","WA","S.sasin"),
                          to=c("AK+MT","E. Migrants","W. Migrants","LA","AK+MT","E. Migrants","OR","E. Migrants","E. Migrants","WA","S. sasin"))
dat$tmpclust <- factor(dat$tmpclust,levels=c("AK+MT","WA","OR","W. Migrants","E. Migrants","LA","S. sasin"))
dat <- subset(dat,dat$SampleID %in% rownames(ru3@tab))

########## multivariate ordination & clustering ##########
#replace missing data with mean allele frequency
rumat2 <- scaleGen(ru2,NA.method="mean",center=T,scale=F)
rumat2 <- rumat2[,!apply(rumat2,2,function(e) any(is.na(e)))]

rumat3 <- scaleGen(ru3,NA.method="mean",center=T,scale=F)
rumat3 <- rumat3[,!apply(rumat3,2,function(e) any(is.na(e)))]

#pca over all samples
pc_ruf_sas <- prcomp(rumat2,center=F,scale=F)$x %>% data.frame() %>% .[1:5]
pc_ruf_sas$SampleID <- rownames(pc_ruf_sas)
pc_ruf <- prcomp(rumat3,center=F,scale=F)$x %>% data.frame() %>% .[1:5]
pc_ruf$SampleID <- rownames(pc_ruf)
dat <- merge(dat,pc_ruf_sas,by="SampleID",all.x=T)
dat <- merge(dat,pc_ruf,by="SampleID",all.x=T)

#k-means clustering
dat$clust <- kmeans(rumat3,4)$cluster
ggplot(data=dat,aes(x=PC1.y,y=PC2.y,col=factor(clust),label=SampleID))+theme_minimal()+
  geom_hline(aes(yintercept=0),lwd=0.25)+geom_vline(aes(xintercept=0),lwd=0.25)+
  geom_text()

############# mantel tests for IBD ##############
# breeding <- dna[grepl("OR|WA|MT|AK",labels(dna)),]
# gen_dist <- dist.dna(breeding,as.matrix=T,pairwise.deletion = T)
# dat2 <- fread("~/Dropbox/ruhu/specimen_data.csv",data.table=F)
# dat2 <- dat2[match(labels(breeding),dat2$SampleID),] 
# geo_dist <- spDists(as.matrix(dat2[,c("Long","Lat")]),longlat = T)
# m <- data.frame(gen=c(gen_dist),geo=c(geo_dist))
# ggplot(data=m,aes(x=geo,y=gen))+
#   geom_point(shape=1)+geom_smooth(method="lm")
# mantel.test(gen_dist,geo_dist)
# 
# migrants <- dna[grepl("CA|AZ|UT|TX|NM",labels(dna)),]
# gen_dist <- dist.dna(migrants,as.matrix=T,pairwise.deletion = T)
# dat2 <- fread("~/Dropbox/ruhu/specimen_data.csv",data.table=F)
# dat2 <- dat2[match(labels(migrants),dat2$SampleID),] 
# geo_dist <- spDists(as.matrix(dat2[,c("Long","Lat")]),longlat = T)
# m <- data.frame(gen=c(gen_dist),geo=c(geo_dist))
# ggplot(data=m,aes(x=geo,y=gen))+
#   geom_point(shape=1)+geom_smooth(method="lm")
# mantel.test(gen_dist,geo_dist) 
# 
# breeding <- dna[grepl("OR|WA|MT|AK|sasin",labels(dna)),]
# #breeding <- breeding[!labels(breeding) %in% c("S.sasin11","S.sasin6","S.sasin5","S.sasin8","S.sasin2"),]
# gen_dist <- dist.dna(breeding,as.matrix=T,pairwise.deletion = T)
# dat2 <- fread("~/Dropbox/ruhu/specimen_data.csv",data.table=F)
# dat2 <- dat2[match(labels(breeding),dat2$SampleID),] 
# geo_dist <- spDists(as.matrix(dat2[,c("Long","Lat")]),longlat = T)
# m <- data.frame(gen=c(gen_dist),geo=c(geo_dist))
# ggplot(data=m,aes(x=geo,y=gen))+
#   geom_point(shape=1)+geom_smooth(method="lm")
# mantel.test(gen_dist,geo_dist) 
# 
# #correlations bw lat/long and PC1 
# #breeding range
# dat2 <- subset(dat,Species=="Selasphorus rufus")
# pclm <- lm(PC1.y~poly(Latitude,2),data=dat2[grepl("OR|WA|MT|AK|sasin",dat$SampleID),])
# summary(pclm) #strong correlation bw PC1 and latitude
# ggplot(data=dat[grepl("OR|WA|MT|AK|sasin",dat2$SampleID),],aes(x=PC1.y,y=Latitude))+
#   theme_minimal()+
#   xlab("Genotype PC1")+ylab("Latitude")+
#   geom_point(shape=1)+
#   geom_text(aes(label=SampleID))+
#   geom_smooth(method="lm",formula=y~poly(x,2),col="black",fill="grey85",lwd=0.75)
# 
# #migrants
# dat2 <- dat#subset(dat,!SampleID %in% c("CA15","CA12","CA11"))
# pclm <- lm(PC1.y~poly(Longitude,1),data=dat2[grepl("CA|NM|UT|TX",dat2$SampleID),])
# summary(pclm)
# ggplot(data=dat2[grepl("CA|NM|UT|TX",dat2$SampleID),],aes(x=PC1.y,y=Longitude))+
#   theme_minimal()+
#   xlab("Genotype PC1")+ylab("Longitude")+
#   geom_point(shape=1)+
#   geom_text(aes(label=SampleID))+
#   geom_smooth(method="lm",formula=y~poly(x,2),col="black",fill="grey85",lwd=0.75)

