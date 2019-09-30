library(raster);library(ggplot2);library(magrittr);library(plyr)
setwd("~/Dropbox/range_map_demo/")

#load birdlife international range map shapefiles
ruhu_range <- shapefile("Selasphorus_rufus_22688296.shp")
cahu_range <- shapefile("Selasphorus_calliope_22688232.shp")
brdthu_range <- shapefile("Selasphorus_platycercus_22688293.shp")
alhu_range <- shapefile("Selasphorus_sasin_22688299.shp")

#create data frame with separate breeding/nonbreeding range polygons
range_list <- list(ruhu_range,cahu_range,brdthu_range,alhu_range)
ranges <- ldply(range_list,function(e){
  brd <- e[e@data$SEASONAL %in% c(1,2),] %>% fortify()
  brd$season <- "breeding"
  brd$species <- e@data$SCINAME[1]
  wnt <- e[e@data$SEASONAL %in% c(1,3),] %>% fortify()
  wnt$season <- "nonbreeding"
  wnt$species <- e@data$SCINAME[1]
  rbind(brd,wnt)
})

#get country/state line basemaps
world <- map_data("world")
state <- map_data("state")

#plot
ggplot(data=ranges,aes(x=long,y=lat,group=group))+coord_map()+
  facet_wrap(~species)+
  xlim(-150,-80)+ylim(5,65)+xlab("")+ylab("")+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.background = element_blank(),
        panel.grid=element_blank(),
        legend.position = "right")+
  geom_path(data=state,lwd=0.25,col="grey")+
  geom_path(data=world,lwd=0.25,col="black")+
  geom_polygon(aes(fill=season))

