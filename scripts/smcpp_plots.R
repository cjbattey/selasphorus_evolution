#plot smcpp results
library(magrittr);library(plyr);library(pbapply);library(ggplot2);library(RColorBrewer)
setwd("~/Dropbox/selasphorus/smcpp/")

smc <- read.csv("plotdata.csv",stringsAsFactors = F)
smc$label[smc$label=="sasin"] <- "sasin\nsasin"
smc$label[smc$label=="sedentarius"] <- "sasin\nsedentarius"
smc$label <- factor(smc$label,levels=c("rufus","sasin\nsasin","sasin\nsedentarius","calliope"))
pdf("~/Dropbox/selasphorus/ms/fig/smcpp_fig.pdf",width=3.25,height=4)
ggplot(data=smc,aes(x=x,y=y))+facet_grid(label~.)+
  theme_classic()+
  theme(axis.title=element_text(size=8),
        strip.text=element_text(size=8),
        axis.text=element_text(size=8),
        panel.grid.major=element_line(color="grey",size=0.25),
        strip.background = element_blank())+
  xlab("Generations")+ylab("Ne")+
  ylim(0,1.25e5)+xlim(2e3,1e5)+
  # scale_color_manual(values=c(rep(brewer.pal(4,"RdYlBu")[2],11),
  #                             rep(brewer.pal(4,"RdYlBu")[3],11),
  #                             rep(brewer.pal(4,"RdYlBu")[4],11),
  #                             rep("violet",11)),guide=F)+
  geom_rect(aes(xmin=18000,xmax=26500,ymin=0,ymax=1.25e5),fill="grey",alpha=0.5)+
  geom_path(lwd=.35,alpha=0.8)
dev.off()


# model <- readLines("models/rufus/model.final.json")
# rho <- model[grep("rho",model)] %>% strsplit(": |,") %>% unlist() %>% .[2] %>% as.numeric()
# theta <- model[grep("theta",model)] %>% strsplit(": |,") %>% unlist() %>% .[2] %>% as.numeric()
# alpha <- model[grep("alpha",model)] %>% strsplit(": |,") %>% unlist() %>% .[2] %>% as.numeric()
# N0 <- model[grep("N0",model)] %>% strsplit(": |,") %>% unlist() %>% .[2] %>% as.numeric()
# knots <- model[(grep("knots",model)+1):(grep("knots",model)+10)] %>% 
#           sapply(function(e) gsub(" |,","",e)) %>%
#             as.numeric()
# hidden_states <- model[(grep("hidden_states",model)+2):(grep("model",model)-3)] %>% 
#                   sapply(function(e) gsub(" |,","",e)) %>%
#                     as.numeric()
# y <- model[(grep("y",model)+1):(grep("rho",model)-3)] %>% 
#       sapply(function(e) gsub(" |,","",e)) %>%
#         as.numeric()
# 
# df <- data.frame(knots,y)
# df$times <- df$knots*(2*N0)*2
# df$endtimes <- c(df$times[2:length(df$times)],df$times[length(df$times)])
# df$pops <- df$y*N0
# 
# ggplot(data=df,aes(x=times,y=pops,xend=endtimes,yend=pops))+
#   scale_x_log10()+
#   ylim(-1e4,1.5e5)+
#   geom_segment()