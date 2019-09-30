#convert vcftools 012 output to treemix
library(pbapply);library(magrittr)
input <- "~/Dropbox/selasphorus/called_snps/selasphorus_het_depth_mq30_biallel_nomd.012"
samples <- read.table("~/Dropbox/selasphorus/called_snps/selasphorus_het_depth_mq30_biallel_nomd.012.indv",header=F)
pops <- list(rufus=c(1,2,17,18,37,38,39,40),#c(1:18,34:40),
             calliope=c(19:25),
             sasin_sasin=c(26,29,33,34),
             sasin_sedentarius=c(27,28,30,31))
data <- fread(input,data.table=F)[,-1]
rownames(data) <- samples[,1]

tmix <- pbapply(data,2,function(e){
  rufus <- e[pops$rufus] %>% sum()
  calliope <- e[pops$calliope] %>% sum()
  sasin_sasin <- e[pops$sasin_sasin] %>% sum()
  sasin_sedentarius <- e[pops$sasin_sedentarius] %>% sum()
  out <- c(paste(80-rufus,rufus,sep=","),
           paste(80-calliope,calliope,sep=","),
           paste(80-sasin_sasin,sasin_sasin,sep=","),
           paste(80-sasin_sedentarius,sasin_sedentarius,sep=","))
})
tmix <- data.frame(t(tmix))
colnames(tmix) <- c("rufus","calliope","sasin","sedentarius")

write.table(tmix,"~/Dropbox/selasphorus/treemix/data/selasphorus_tmix_t10k.tmix",
            row.names = F,col.names = T,quote = F,sep=" ")

source("/Applications/treemix-1.13/src/plotting_funcs.R")
plot_tree("~/Dropbox/selasphorus/treemix/1M")

