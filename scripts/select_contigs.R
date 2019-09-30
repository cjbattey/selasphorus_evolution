library(plyr);library(magrittr);library(data.table);library(ggplot2)
#select contigs by length for use in smc++
gff <- read.table("~/Dropbox/augustus/Canna_MUGM01_fwd.gff",comment.char = "#",skip = 19,sep="\t",stringsAsFactors = F)

#get contigs over 1e6 bases
ddply(gff,.(V1),summarize,length=max(V5)) %>%
  subset(length>1e6) %>% 
  .[,1] %>%
  write.table("~/Dropbox/selasphorus/smcpp/anhu_mugm_contigs.txt",
              row.names = F,col.names = F,quote = F)


contigs <- ddply(gff,.(V1),summarize,length=max(V5)) %>%
  subset(length>1e6) %>% 
  .[,1]

#rufus bootstraps
files <- list.files("~/Dropbox/selasphorus/smcpp/data/rufus_sasin/")
rufboot <- files[grepl("rufus_M",files)] %>% .[!grepl("sasin",.)]
for(i in 1:10){
  boot <- sample(contigs,length(contigs),replace=T)
  set <- grep(paste(boot,collapse = "|"),rufboot,value=T)
  set <- paste0("/home/ubuntu/Dropbox/selasphorus/smcpp/data/rufus_sasin/",set)
  print(set)
  write.table(set,paste0("~/Dropbox/selasphorus/smcpp/boot/ruf_boot",i,".txt"),row.names=F,col.names=F,quote=F)
}

#sasin sedentarius bootstraps
files <- list.files("~/Dropbox/selasphorus/smcpp/data/sas_sed/")
sedboot <- files[grepl("sedentarius_M",files)]
for(i in 1:10){
  boot <- sample(contigs,length(contigs),replace=T)
  set <- grep(paste(boot,collapse = "|"),sedboot,value=T)
  set <- paste0("/home/ubuntu/Dropbox/selasphorus/smcpp/data/sas_sed/",set)
  print(set)
  write.table(set,paste0("~/Dropbox/selasphorus/smcpp/boot/sed_boot",i,".txt"),row.names=F,col.names=F,quote=F)
}

#sasin sasin bootstraps
sasboot <- files[grepl("sasin_M",files)]
for(i in 1:10){
  boot <- sample(contigs,length(contigs),replace=T)
  set <- grep(paste(boot,collapse = "|"),sasboot,value=T)
  set <- paste0("/home/ubuntu/Dropbox/selasphorus/smcpp/data/sas_sed/",set)
  print(set)
  write.table(set,paste0("~/Dropbox/selasphorus/smcpp/boot/sas_boot",i,".txt"),row.names=F,col.names=F,quote=F)
}

#calliope bootstraps
calboot <- list.files("~/Dropbox/selasphorus/smcpp/data/calliope/")
for(i in 1:10){
  boot <- sample(contigs,length(contigs),replace=T)
  set <- grep(paste(boot,collapse = "|"),calboot,value=T)
  set <- paste0("/home/ubuntu/Dropbox/selasphorus/smcpp/data/calliope/",set)
  print(set)
  write.table(set,paste0("~/Dropbox/selasphorus/smcpp/boot/cal_boot",i,".txt"),row.names=F,col.names=F,quote=F)
}
