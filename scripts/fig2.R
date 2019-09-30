#fig 1
setwd("~/selasphorus_evolution/")
source("scripts/ddrad_PCA.R")
pal <- "RdYlBu"

#load birdlife international range map shapefiles
ruhu_range <- shapefile("data/range_maps/Selasphorus_rufus_22688296.shp")
cahu_range <- shapefile("data/range_maps/Selasphorus_calliope_22688232.shp")
alhu_range <- shapefile("data/range_maps/Selasphorus_sasin_22688299.shp")

#create data frame with separate breeding/nonbreeding range polygons
range_list <- list(ruhu_range,cahu_range,alhu_range)
ranges <- ldply(range_list,function(e){
  brd <- e[e@data$SEASONAL %in% c(1,2),] %>% fortify()
  brd$season <- "breeding"
  brd$species <- e@data$SCINAME[1]
  wnt <- e[e@data$SEASONAL %in% c(1,3),] %>% fortify()
  wnt$season <- "nonbreeding"
  wnt$species <- e@data$SCINAME[1]
  rbind(brd,wnt)
})

#load and summarize multivariate data
dat2 <- subset(dat,Species %in% c("Selasphorus rufus","Selasphorus sasin sasin","Selasphorus sasin sedentarius") & !is.na(PC1.y))
sum_breeding <- ddply(subset(dat2,grepl("AK|MT|WA|OR|sasin",dat2$SampleID)),
                      .(loc_clust,Species),
                      summarize,
                      long=mean(Longitude),lat=mean(Latitude),
                      Samples=length(Longitude),
                      PC1=mean(PC1.y,na.rm=T))
sum_all_rufous <- ddply(subset(dat2,grepl("AK|MT|WA|OR|CA|NM|UT|AZ|TX|LOU",dat2$SampleID) & !is.na(PC1.y)),
                        .(loc_clust,Species),
                        summarize,
                        long=mean(Longitude),lat=mean(Latitude),
                        Samples=length(Longitude),
                        PC1=mean(PC1.y,na.rm=T))
#brd <- shapefile("~/Dropbox/BirdlifeIntnl_Range_Maps/Selasphorus_rufus_22688296.shp")
brd <- subset(ranges,species %in% c("Selasphorus rufus") & season=="breeding")
wnt <- subset(ranges,species %in% c("Selasphorus rufus") & season=="nonbreeding")
map <- map_data("world")
states <- map_data("state")

mapplot_sasin <- ggplot()+coord_map()+theme_bw()+
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border=element_blank(),
        plot.margin = unit(rep(1,4),"mm"),
        axis.text=element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = c(0.11,0.43),
        legend.margin = unit(3,"mm"),
        #legend.box = "horizontal",
        #legend.box.just = "bottom",
        legend.background = element_blank(),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8))+
  xlim(-149,-110)+ylim(30,62)+
  scale_color_distiller(palette = pal,direction=1,name="Genotype\nPC1")+
  scale_size_continuous(breaks=c(1,5,10),range = c(1,8))+
  geom_polygon(data=subset(ranges,species %in% c("Selasphorus sasin","Selasphorus rufus") & season=="breeding"),
               aes(x=long,y=lat,group=group),fill="grey85")+
  #geom_path(data=states,aes(x=long,y=lat,group=group),lwd=0.05,col="black")+
  geom_path(data=map,aes(x=long,y=lat,group=group),lwd=0.25)+
  geom_point(data=sum_breeding,aes(y=lat,x=long,col=PC1,size=Samples))+
  geom_point(data=sum_breeding,aes(y=lat,x=long,size=Samples),shape=1,stroke=.4,col="black")+
  annotate(geom="segment",x=-138,xend=-128,y=42,yend=42,colour="black")+
  annotate(geom="text",x=-133,y=43.5,label="S. rufus",size=2.75)+
  annotate(geom="text",x=-133,y=40.5,label="S. sasin",size=2.75)+
  annotate(geom="segment",x=-133,xend=-133,y=45,yend=48,arrow=arrow(length = unit(2,"mm")))+
  annotate(geom="segment",x=-133,xend=-133,y=38.5,yend=35,arrow=arrow(length = unit(2,"mm")))
mapplot_sasin <- mapplot_sasin+guides(color=guide_colorbar(order=2,
                                                           barheight = unit(24,"mm"),
                                                           barwidth = unit(4,"mm")),
                                      size=guide_legend(order=1))

pcplot_facets_sasin <- ggplot(data=subset(dat2,grepl("AK|MT|WA|OR|sasin",dat2$SampleID)),aes(x=PC1.y,y=PC2.y))+
  theme(text=element_text(size=8),
        plot.margin = unit(c(1,3,1,1),"mm"),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size=8),
        axis.text=element_text(size=7),
        axis.title = element_text(size=8),
        axis.ticks = element_line(size=0.25),
        legend.position = "none")+
  xlab("PC1")+ylab("PC2")+
  scale_color_distiller(palette = pal,direction=1)+
  facet_wrap(~Species,ncol=1)+
  geom_vline(aes(xintercept=0),lwd=0.25)+geom_hline(aes(yintercept=0),lwd=0.25)+
  geom_point(data=subset(dat2,grepl("AK|MT|WA|OR|sasin",dat2$SampleID))[,colnames(dat) != "Species"],col="grey60",shape=21,stroke=.2)+
  geom_point(aes(col=PC1.y))

mapplot_rufus <- ggplot()+coord_map()+theme_bw()+
  theme(panel.grid = element_blank(),
        panel.border=element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        plot.margin = unit(rep(1,4),"mm"),
        axis.text=element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = c(0.18,0.26),
        legend.margin = unit(1,"mm"),
        legend.box = "horizontal",
        legend.box.just = "bottom",
        legend.background = element_blank(),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8))+
  xlim(-149,-100)+ylim(29,62)+
  scale_color_distiller(palette = pal,direction=1,name="Genotype\nPC1")+
  scale_size_continuous(breaks=c(1,5,10,15),range = c(1,8))+
  geom_polygon(data=brd,aes(x=long,y=lat,group=group),fill="grey85")+
  #geom_polygon(data=wnt,aes(x=long,y=lat,group=group),fill="grey85")+
  #geom_path(data=states,aes(x=long,y=lat,group=group),lwd=0.05,col="black")+
  geom_path(data=map,aes(x=long,y=lat,group=group),lwd=0.25)+
  geom_point(data=sum_all_rufous,aes(y=lat,x=long,col=PC1,size=Samples))+
  geom_point(data=sum_all_rufous,aes(y=lat,x=long,size=Samples),shape=1,stroke=.4,col="black")+
  annotate(geom="segment",x=-121,xend=-118,y=41,yend=37,arrow=arrow(length=unit(2,"mm")))+
  annotate(geom="segment",x=-115,xend=-108,y=45,yend=38,arrow=arrow(length=unit(2,"mm")))
mapplot_rufus <- mapplot_rufus+guides(color=guide_colorbar(order=1,
                                                             barwidth = unit(4,"mm"),
                                                             barheight = unit(24,"mm")),
                                        size=guide_legend(order=2))

pcplot_facets_rufus <- ggplot(data=dat2[!is.na(dat2$tmpclust),],aes(x=PC1.y,y=PC2.y))+
  theme(panel.grid = element_blank(),
        plot.margin = unit(c(1,3,1,1),"mm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size=8),
        #axis.text = element_text(size=7),
        axis.text = element_blank(),
        axis.title = element_text(size=8),
        axis.ticks = element_blank(),
        legend.position = "none")+
  xlab("PC1")+ylab("PC2")+
  scale_color_distiller(palette = pal,direction=1)+
  xlim(-9,18)+
  facet_wrap(~tmpclust)+
  geom_vline(aes(xintercept=0),lwd=0.25)+geom_hline(aes(yintercept=0),lwd=0.25)+
  geom_point(data=dat2[colnames(dat2) != "tmpclust"],col="grey60",shape=21,stroke=.2)+
  geom_point(aes(col=PC1.y))

dat2 <- subset(dat,Species=="Selasphorus rufus" & !grepl("LOU",dat$SampleID))
dat2$breeding <- grepl("OR|WA|MT|AK",dat2$SampleID)
pclm <- lm(PC1.y~poly(Latitude,2),data=dat2[dat2$breeding==T,])
summary(pclm) #strong correlation bw PC1 and latitude
lat_plot <- ggplot(data=dat2[dat2$breeding==T,],aes(x=PC1.y,y=Latitude))+ggtitle("Breeding")+
  theme(text=element_text(size=8),
        plot.title=element_text(size=10,face = "plain"),
        axis.text=element_text(size=7),
        axis.title.x=element_text(hjust = 2.1),
        axis.ticks = element_line(size=0.25),
        plot.background = element_blank(),
        panel.background = element_blank())+
  xlab("Genotype PC1")+ylab("Latitude")+
  geom_point(shape=21,stroke=.2)+
  #geom_text(aes(label=SampleID))+
  geom_smooth(method="lm",col="black",fill="grey85",lwd=0.5,se=F)

dat3 <- subset(dat2,!is.na(PC1.y) & breeding==F & !SampleID %in% c("CA15","CA12"))
pclm <- lm(Longitude~poly(PC1.y,2),data=dat3)
summary(pclm)
long_plot <- ggplot(data=dat3,aes(x=PC1.y,y=Longitude))+ggtitle("Fall Migration")+
  theme(text=element_text(size=8),
        plot.title=element_text(size=10,face = "plain"),
        axis.text=element_text(size=7),
        axis.title.y=element_text(margin=unit(c(0,0,0,0),"mm")),
        axis.ticks = element_line(size=0.25),
        plot.background = element_blank(),
        panel.background = element_blank())+
  xlab(" ")+ylab("Longitude")+
  geom_point(data=dat3,shape=21,stroke=.2)+
  #geom_text(aes(label=SampleID))+
  geom_smooth(data=dat3,method="lm",col="black",fill="grey85",lwd=0.5,se = F)

source("scripts/str2df.R")
str_mig <- str2df("analysis/structure/out/ruf_sas_mig_mac3.str_k2_run01_f")
str_mig <- merge(str_mig,dat,by.x="id",by.y="SampleID",all.x=T,all.y=F)
str_mig$id <- factor(str_mig$id,levels=unique(arrange(str_mig,Longitude)$id))
strplot_mig <- ggplot(data=str_mig,aes(y=value,x=id,fill=variable))+
  theme_classic()+theme(axis.text=element_text(angle=45,hjust=1),
                        plot.background = element_blank(),
                        panel.background = element_blank(),
                        text=element_text(size=8),
                        axis.text.y=element_blank(),
                        axis.ticks.y=element_blank(),
                        axis.title=element_blank(),
                        legend.position = "none")+
  scale_fill_brewer(palette=pal)+
  geom_bar(stat="identity")

str_brd <- str2df("analysis/structure/out/ruf_sas_brd_mac3.str_k2_run01_f")
str_brd <- subset(str_brd,id!="S.sp1")
str_brd <- merge(str_brd,dat,by.x="id",by.y="SampleID",all.x=T,all.y=F)
str_brd$id <- factor(str_brd$id,levels=unique(arrange(str_brd,Latitude)$id))
strplot_brd_vert <- ggplot(data=str_brd,aes(y=value,x=id,fill=variable))+
  coord_flip()+
  theme_classic()+theme(axis.text=element_text(),
                        text=element_text(size=8),
                        axis.text.x=element_blank(),
                        axis.ticks.x=element_blank(),
                        axis.title=element_blank(),
                        legend.position = "none")+
  scale_fill_brewer(palette=pal)+
  scale_x_discrete(position="top")+
  geom_bar(stat="identity")
strplot_brd <- ggplot(data=str_brd,aes(y=value,x=id,fill=variable))+
                  theme(axis.line = element_blank(),
                        text=element_text(size=8),
                        axis.text=element_blank(),
                        axis.ticks=element_blank(),
                        axis.title=element_blank(),
                        panel.border = element_blank(),
                        legend.position = "none")+
  coord_cartesian(clip="off")+
  scale_fill_manual(values = c(brewer.pal(3,"RdYlBu")[3],brewer.pal(3,"RdYlBu")[1]))+
  geom_bar(stat="identity")+
  annotate(geom="text",x=2.5,y=-0.24,label="S. sasin\nsedentarius",size=2.5)+
  annotate(geom="segment",x=0.5,xend=5,y=-0.03,yend=-0.03)+
  annotate(geom="text",x=9.5,y=-0.24,label="S. sasin\nsasin",size=2.5)+
  annotate(geom="segment",x=5.5,xend=11,y=-0.03,yend=-0.03)+
  annotate(geom="text",x=20,y=-0.24,label="S. rufus\nOR",size=2.5)+
  annotate(geom="segment",x=11.5,xend=30,y=-0.03,yend=-0.03)+
  annotate(geom="text",x=31.5,y=-0.24,label="S. rufus\nWA",size=2.5)+
  annotate(geom="segment",x=30.5,xend=35,y=-0.03,yend=-0.03)+
  annotate(geom="text",x=38.5,y=-0.24,label="S. rufus\nMT",size=2.5)+
  annotate(geom="segment",x=35.5,xend=40,y=-0.03,yend=-0.03)+
  annotate(geom="text",x=45.5,y=-0.24,label="S. rufus\nWA",size=2.5)+
  annotate(geom="segment",x=40.5,xend=49,y=-0.03,yend=-0.03)+
  annotate(geom="text",x=55.5,y=-0.24,label="S. rufus\nAK",size=2.5)+
  annotate(geom="segment",x=49.5,xend=61,y=-0.03,yend=-0.03)


pdf("fig/fig_2.pdf",width=6.5,height=6.5,useDingbats = F)
ggdraw()+
  draw_plot(mapplot_sasin,width=0.4,x=0,y=0.25)+
  draw_plot(mapplot_rufus,width=0.47,x=0,y=-0.24)+
  draw_plot(pcplot_facets_sasin,x=0.63,y=0.50,width=0.36,height=0.5)+
  draw_plot(lat_plot,x=0.45,y=0,width=0.27,height=0.27)+
  draw_plot(long_plot,x=0.72,y=0,width=0.27,height=0.27)+
  draw_plot(strplot_brd,x=0.4,y=0.265,width=0.6,height=0.2)+
  draw_image("data/Selasphorus_Illustrations/Rufus_Hummingbird_Final_RGB.png",
             scale=0.3,x=0.02,y=0.35)+
  draw_image("data/Selasphorus_Illustrations/Allen's_Hummingbird_Final_RGB.png",
             scale=0.3,x=0.02,y=0.11)+
  draw_plot_label(c("A","B","C","D","E"),x=c(0.005,0.4,0.63,0.005,0.4),y=c(1,1,1,.5,.5),fontface="plain")
dev.off()


