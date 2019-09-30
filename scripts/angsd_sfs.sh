#!/bin/bash
cd ~/Dropbox/selasphorus/

##estimate SFS per species with angsd
angsd -b bamlists/rufus_angsd_sfs.txt -anc anhu_alignment/Canna_MUGM01_haploid.fa \
-dosaf 1 -remove_bads -unique_only -minMapQ 20 -minQ 20 -only_proper_pairs 1 -trim 0 \
-gl 1 -P 60 -out ~/selasphorus/data/wgs/sfs/ruf

angsd -b bamlists/sasin_sasin.txt -anc anhu_alignment/Canna_MUGM01_haploid.fa \
-dosaf 1 -remove_bads -unique_only -minMapQ 20 -minQ 20 -only_proper_pairs 1 -trim 0 \
-gl 1 -P 60 -out ~/selasphorus/data/wgs/sfs/sas

angsd -b bamlists/sasin_sedentarius.txt -anc anhu_alignment/Canna_MUGM01_haploid.fa \
-dosaf 1 -remove_bads -unique_only -minMapQ 20 -minQ 20 -only_proper_pairs 1 -trim 0 \
-gl 1 -P 60 -out ~/selasphorus/data/wgs/sfs/sed

angsd -b bamlists/calliope.txt -anc anhu_alignment/Canna_MUGM01_haploid.fa \
-dosaf 1 -remove_bads -unique_only -minMapQ 20 -minQ 20 -only_proper_pairs 1 -trim 0 \
-gl 1 -P 60 -out ~/selasphorus/data/wgs/sfs/cal

angsd -b bamlists/sasin.txt -anc anhu_alignment/Canna_MUGM01_haploid.fa \
-dosaf 1 -remove_bads -unique_only -minMapQ 20 -minQ 20 -only_proper_pairs 1 -trim 0 \
-gl 1 -P 60 -out ~/selasphorus/data/wgs/sfs/sasin

#note this requires ~300GB RAM (!)
realSFS -P 60 ~/selasphorus/data/wgs/sfs/cal.saf.idx \
~/selasphorus/data/wgs/sfs/ruf.saf.idx \
~/selasphorus/data/wgs/sfs/sas.saf.idx \
~/selasphorus/data/wgs/sfs/sed.saf.idx > ~/selasphorus/data/wgs/sfs/cal_ruf_sas_sed.sfs