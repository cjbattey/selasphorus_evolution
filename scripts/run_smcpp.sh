#!/bin/bash
cd ~/Dropbox/selasphorus/smcpp/

####get list of longest contigs from gff in R
#library(magrittr);library(plyr)
#gff <- read.table("~/Dropbox/augustus/Canna_MUGM01_fwd.gff",comment.char = "#",skip = 19,sep="\t",stringsAsFactors = F)
#ddply(gff,.(V1),summarize,length=max(V5)) %>% arrange(desc(length)) %>% .[1:10,] %>% .[,1] %>% write.table("~/Dropbox/selasphorus/smcpp/anhu_mugm_contigs.txt",col.names=F,row.names = F,sep="\t",quote = F)

parallel smc++ vcf2smc -d {2} {2} \
-c 20000 \
../called_snps/selasphorus_mq30_nomissing_biallelic.recode.vcf.gz \
data/rufus_sasin/rufus\_{1}\_{2}.smc.gz \
{1} \
rufus:AK1,AK3,WA22,WA31,WA35,WA36,CA10,CA15,CA19,CA22,UT1,UT2,UT3,CA33,CA5,CA6,NM1,NM2,NM3,NM4,NM5,NM6,NM7 \
::: `cat anhu_mugm_contigs.txt` ::: WA22 WA31 WA35 WA36 AK3 CA33 NM1 NM2

parallel smc++ vcf2smc -d {2} {2} \
--missing-cutoff 20000 \
../called_snps/selasphorus_mq30_nomissing_biallelic.recode.vcf.gz \
data/rufus_sasin/sasin\_{1}\_{2}.smc.gz \
{1} \
sasin:Ssasin2,Ssasin4,Ssasin5,Ssasin6,Ssasin7,Ssasin9,Ssasin10,Ssasin11 \
::: `cat anhu_mugm_contigs.txt` ::: Ssasin4 Ssasin5 Ssasin6 Ssasin7 Ssasin9 Ssasin10 Ssasin11

parallel smc++ vcf2smc -d {2} {2} \
--missing-cutoff 20000 \
../called_snps/selasphorus_mq30_nomissing_biallelic.recode.vcf.gz \
data/sas_sed/sedentarius\_{1}\_{2}.smc.gz \
{1} \
sedentarius:Ssasin2,Ssasin5,Ssasin6,Ssasin11 \
::: `cat anhu_mugm_contigs.txt` ::: Ssasin5 Ssasin6 Ssasin2 Ssasin11

parallel smc++ vcf2smc -d {2} {2} \
-c 20000 \
../called_snps/selasphorus_mq30_nomissing_biallelic.recode.vcf.gz \
data/sas_sed/sasin\_{1}\_{2}.smc.gz \
{1} \
sasin:Ssasin4,Ssasin7,Ssasin9,Ssasin10 \
::: `cat anhu_mugm_contigs.txt` ::: Ssasin4 Ssasin7 Ssasin9 Ssasin10

parallel smc++ vcf2smc -d {2} {2} \
-c 20000 \
../called_snps/selasphorus_mq30_nomissing_biallelic.recode.vcf.gz \
data/sas_sed/sas_sed\_{1}\_{2}.smc.gz \
{1} \
sasin:Ssasin4,Ssasin7,Ssasin9,Ssasin10 sedentarius:Ssasin2,Ssasin5,Ssasin6,Ssasin11 \
::: `cat anhu_mugm_contigs.txt` ::: Ssasin4 Ssasin7 Ssasin5 Ssasin11

parallel smc++ vcf2smc -d {2} {2} \
-c 20000 \
../called_snps/selasphorus_mq30_nomissing_biallelic.recode.vcf.gz \
data/sas_sed/sed_sas\_{1}\_{2}.smc.gz \
{1} \
sedentarius:Ssasin2,Ssasin5,Ssasin6,Ssasin11 sasin:Ssasin4,Ssasin7,Ssasin9,Ssasin10 \
::: `cat anhu_mugm_contigs.txt` ::: Ssasin4 Ssasin7 Ssasin5 Ssasin11

parallel smc++ vcf2smc -d {2} {2} \
-c 20000 \
../called_snps/selasphorus_mq30_nomissing_biallelic.recode.vcf.gz \
data/calliope/calliope\_{1}\_{2}.smc.gz \
{1} \
calliope:Scal1,Scal2,Scal3,Scal4,Scal5,Scal6,Scal7 \
::: `cat anhu_mugm_contigs.txt` ::: Scal1 Scal2 Scal3 Scal4 Scal5 Scal6 Scal7

parallel smc++ vcf2smc -d {2} {2} \
-c 20000 \
../called_snps/selasphorus_mq30_nomissing_biallelic.recode.vcf.gz \
data/rufus_sasin/rufus_sasin\_{1}\_{2}.smc.gz \
{1} \
rufus:AK1,AK3,CA10,CA15,CA19,CA22,UT1,UT2,UT3,WA22,WA31,WA35,WA36,CA33,CA5,CA6,NM1,NM2,NM3,NM4,NM5,NM6,NM7,OR10,OR2 \
sasin:Ssasin2,Ssasin4,Ssasin5,Ssasin6,Ssasin7,Ssasin9,Ssasin10,Ssasin11 \
::: `cat anhu_mugm_contigs.txt` ::: OR2 WA22 WA35 AK3 Ssasin5 Ssasin9 Ssasin11

parallel smc++ vcf2smc -d {2} {2} \
-c 20000 \
../called_snps/selasphorus_mq30_nomissing_biallelic.recode.vcf.gz \
data/rufus_sasin/sasin_rufus\_{1}\_{2}.smc.gz \
{1} \
sasin:Ssasin2,Ssasin4,Ssasin5,Ssasin6,Ssasin7,Ssasin9,Ssasin10,Ssasin11 \
rufus:AK1,AK3,CA10,CA15,CA19,CA22,UT1,UT2,UT3,WA22,WA31,WA35,WA36,CA33,CA5,CA6,NM1,NM2,NM3,NM4,NM5,NM6,NM7,OR10,OR2 \
::: `cat anhu_mugm_contigs.txt` ::: OR2 WA22 WA35 AK3 Ssasin5 Ssasin9 Ssasin11

smc++ estimate -o models/rufus/ -c 30000 --cores 8 --unfold \
 --thinning 400 4.6e-9 data/rufus_sasin/rufus_M*

smc++ estimate -o models/sasin/ -c 30000 --cores 8 --unfold \
 --thinning 400 4.6e-9 data/rufus_sasin/sasin_M*

smc++ estimate -o models/calliope/ -c 30000 --cores 8 --unfold \
 --thinning 400 4.6e-9 data/calliope/calliope_M*

smc++ estimate -o models/sas/ -c 30000 --cores 8 --unfold \
 --thinning 400 4.6e-9 data/sas_sed/sasin_M*

smc++ estimate -o models/sed/ -c 30000 --cores 8 --unfold \
 --thinning 400 4.6e-9 data/sas_sed/sedentarius_M*

smc++ split -o models/ruf_sas/ -c 30000 --cores 8 --unfold --thinning 400 \
models/rufus/model.final.json models/sasin/model.final.json data/rufus_sasin/*

smc++ split -o models/sas_sed/ -c 30000 --cores 8 --unfold --thinning 400 \
models/sas/model.final.json models/sed/model.final.json data/sas_sed/*

