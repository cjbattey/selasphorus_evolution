#d-test summaries
library(pbapply);library(data.table);library(ggridges)
setwd("~/selasphorus_evolution/analysis/abba_baba/")
d <- fread("rsc_D_out.txt",data.table = F)
d$H1[grepl("sasin2|sasin5|sasin6|sasin11",d$H1)] <- "sed"
d$H2[grepl("sasin2|sasin5|sasin6|sasin11",d$H2)] <- "sed"
d$H3[grepl("sasin2|sasin5|sasin6|sasin11",d$H3)] <- "sed"
d$H1[grepl("sasin4|sasin7|sasin9|sasin10",d$H1)] <- "sas"
d$H2[grepl("sasin4|sasin7|sasin9|sasin10",d$H2)] <- "sas"
d$H3[grepl("sasin4|sasin7|sasin9|sasin10",d$H3)] <- "sas"
d$H1[grepl("Scal",d$H1)] <- "cal"
d$H2[grepl("Scal",d$H2)] <- "cal"
d$H3[grepl("Scal",d$H3)] <- "cal"
# d$H1[grepl("CA|NM|WA|OR|UT|AK",d$H1)] <- "ruf"
# d$H2[grepl("CA|NM|WA|OR|UT|AK",d$H2)] <- "ruf"
# d$H3[grepl("CA|NM|WA|OR|UT|AK",d$H3)] <- "ruf"

trees <- pbapply(d[,1:3],1,function(e){
  loc_tree <- gsub(".bam","",e) %>% gsub("[[:digit:]]","",.)
  sp_tree <- gsub(".bam","",e) %>% gsub("[[:digit:]]","",.) %>% 
              gsub("AK|NM|WA|OR|CA|UT|AZ","ruf",.)
  c(paste0("((",loc_tree[[1]],",",loc_tree[[2]],"),",loc_tree[[3]],")"),
    paste0("((",sp_tree[[1]],",",sp_tree[[2]],"),",sp_tree[[3]],")"))
}) %>% t() %>% as.data.frame()

d <- cbind(d,trees)
d <- subset(d,V2 %in% c("((ruf,ruf),sas)",
                        "((ruf,ruf),sed)",
                        "((ruf,ruf),cal)",
                        "((sas,sas),ruf)",
                        "((sas,sas),cal)",
                        "((sas,sas),sed)",
                        "((sas,sed),ruf)",
                        "((sas,sed),cal)",
                        "((sed,sed),ruf)",
                        "((sed,sas),ruf)",
                        "((sed,sed),cal)",
                        "((sed,sas),cal)",
                        "((cal,cal),sas)",
                        "((cal,cal),sed)",
                        "((cal,cal),ruf)"))

d$V2[d$V2=="((sed,sas),ruf)"] <- "((sas,sed),ruf)"
d$V2[d$V2=="((sed,sas),cal)"] <- "((sas,sed),cal)"


dtable <- ddply(d,.(V2),summarize,
      #nABBA=round(median(nABBA)),
      #nBABA=round(median(nBABA)),
      median_D=format(signif(median(abs(Dstat)),3),scientific=T),
      #D_CI=paste(format(signif(quantile(abs(Dstat),0.025),3),scientific=T),
       #          format(signif(quantile(abs(Dstat),0.975),3),scientific=T),sep=" - "),
      median_Z=signif(median(abs(Z)),3),
      #Z_CI=paste(signif(quantile(abs(Z),0.025),3),signif(quantile(abs(Z),0.975),3),sep=" - "),
      prop_sig=signif(length(Z[abs(Z)>4.3])/length(Z),3)
)

pdf("~/selasphorus_evolution/fig/fig4_dstat.pdf",width=3.5,height=3)
ggplot(data=d,aes(x=abs(Z),y=V2,fill=..x..))+
  theme(axis.title=element_text(size=8),
        axis.text=element_text(size=8))+
  ylab("Topology")+xlab("Z")+
  scale_fill_distiller(palette = "RdYlBu",guide=F)+
  geom_vline(aes(xintercept=4.3),col="red",lwd=0.5)+
  geom_density_ridges_gradient(scale=1.5,lwd=0.35)+
  geom_point(data=dtable,shape=8,size=.75,aes(x=-2,y=as.integer(factor(V2))+.6,col=median_Z>=4.31))+
  scale_color_manual(values=rep(c("white","black"),nrow(dtable)),guide=F)
dev.off()

d2 <- subset(d,V2 %in% c("((ruf,ruf),sas)","((ruf,ruf),sed)"))







