library(data.table);library(ggplot2);library(magrittr);library(plyr)
setwd("~/Dropbox/selasphorus/moments/")
theme_set(theme_classic()+theme(axis.text=element_text(size=7,angle=45,hjust=1),
                                axis.title=element_blank(),
                                strip.text=element_text(size=8),
                                strip.background = element_blank()))

#plot model parameters assuming mu=2.3e-9 substitutions/bp/year and L=9.5e8
rp <- fread("realparams_boots.txt")
names(rp) <- c("nref","ncal","nruf","nsas","nrufsas","Trs","Trsc","mcrs","mrs","mcr","mcs","ll")
#rp <- subset(rp,ll>quantile(rp$ll,0.75))
mrp <- melt(rp)
best <- subset(rp,ll==max(rp$ll))
mbest <- melt(best)
#mbest <- subset(mbest,value>0)
#mrp <- subset(mrp,value>0)
mrp <- subset(mrp,variable != "ll")
medians <- ddply(mrp,.(variable),summarize,median=median(value),low=quantile(value,0.025),high=quantile(value,0.975))
pdf("~/selasphorus_evolution/fig/moments_param_distributions.pdf",width=6,height=4,useDingbats = F)
ggplot(data=mrp,aes(x=value))+
  facet_wrap(~variable,scales="free")+
  #scale_x_log10()+
  #geom_density(fill="grey",lwd=0.1)+
  geom_histogram(fill="grey",lwd=0.1,bins = 20)+
  geom_vline(data=medians,aes(xintercept=median),col="red",lwd=0.5)+
  geom_vline(data=medians,aes(xintercept=low),col="orange",lwd=0.2)+
  geom_vline(data=medians,aes(xintercept=high),col="red",lwd=0.2)
dev.off()

#check for parameter values near upper/lower limits
mp <- fread("modelparams_boots.txt")
names(mp) <- c("ncal","nruf","nsas","nrufsas","Trs","Trsc","mcrs","mrs","mcr","mcs","ll","theta")
#mp <- subset(mp,ll>quantile(mp$ll,0.5))
mmp <- melt(mp)
mbest2 <- subset(mp,ll==max(mp$ll)) %>% melt()
ggplot(data=mmp,aes(x=value))+
  facet_wrap(~variable,scales="free")+
  geom_histogram(fill="grey",lwd=0.1)+
  geom_vline(data=mbest2,aes(xintercept=value),col="red",lwd=0.3)

nrow(rp)
print(mbest)
