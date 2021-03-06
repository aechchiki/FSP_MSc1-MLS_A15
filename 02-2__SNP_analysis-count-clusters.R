#####
#
#  Author: Amina Echchiki [mailto:amina.echchiki@unil.ch]
#  Affiliation: MSc student, MLS bioinformatics, University of Lausanne
#  Date: 2015-12-10
#  Name: SNP_analysis-count-clusters.R
#  Purpose: Analyze SNP counts extracted by Nu6_SNP-analysis.sh and Rhiir2_SNP-analysis.sh, and clusterize isolates using SNP count as distance
#  Comments: Requires files of mono/poly SNP count (run SNP_count-data_download-script.sh)
#  Testing platforms: desktop Linux 3.16.0-50-generic (Debian 7.6); server Linux 2.6.32-504.el6.x86_64
#  R version: R version 3.0.2 (2013-09-25) -- "Frisbee Sailing"
#
#####


#
# step 0. retrieve SNP count data
#

# follow instruction in SNP_count-data_download-script.sh


#
# step 1. analyze Nu6 SNP data 
#


## prepare working directory 
wdir="/home/dovah/Documents/MSc/MSc-A15/FSP/09-12_NEW-SCRIPTS/0R_analysis/" # change string to own dwd directory
setwd(wdir)
library(lattice) #graphical package  
library(cluster) #clustering package 

## define genome length
Nu6_len=115079203; Nu6_len_kb=Nu6_len/1000 #Nu6 length
Rhiir2_len=136807476; Rhiir2_len_kb=Rhiir2_len/1000 #Rhiir2 length


## import files 
monoNu6=read.table("isolateNu6_monoRS-sort", h=T) #read Nu6 monoSNPs file
polyNu6=read.table("isolateNu6_polyRS-sort", h=T) #read Nu6 polySNPs file 


## SNP count processing
#mono
monoNu6_SNPcount=monoNu6[2:5] #retrieve monoSNP count
isonames=monoNu6[1] #retrieve isolates order
#poly
polyNu6_SNPcount=polyNu6[2:5] #retrieve polySNP count
#global
globalNu6_SNPcount=monoNu6_SNPcount+polyNu6_SNPcount #retrieve globalSNP count 
globalNu6=data.frame(isonames, globalNu6_SNPcount) #add isolates names
write.table(globalNu6, "globalNu6") #write to file: import for final report 
#isolates_list
rownames(isonames)=isonames$isolate


## SNP count results
#mono
monoNu6_genome=rowSums(monoNu6_SNPcount)
mean(monoNu6_genome); sd(monoNu6_genome)
min(monoNu6_genome); isonames$isolate[which.min(monoNu6_genome)] #minimal monoSNP count
min(monoNu6_genome[-which.min(monoNu6_genome)]); isonames$isolate[which.min(monoNu6_genome[-which.min(monoNu6_genome)])] #minimal monoSNP count except DAOM
max(monoNu6_genome); isonames$isolate[which.max(monoNu6_genome)] #maximal monoSNP count
#poly
polyNu6_genome=rowSums(polyNu6_SNPcount)
mean(polyNu6_genome); sd(polyNu6_genome)
min(polyNu6_genome); isonames$isolate[which.min(polyNu6_genome)] #minimal polySNP count
min(polyNu6_genome[-which.min(polyNu6_genome)]); isonames$isolate[which.min(polyNu6_genome[-which.min(polyNu6_genome)])] #minimal polySNP count except DAOM
max(polyNu6_genome); isonames$isolate[which.max(polyNu6_genome)] #maximal polySNP count
#global
globalNu6_genome=rowSums(globalNu6_SNPcount)
mean(globalNu6_genome); sd(globalNu6_genome)
mean(globalNu6_genome)/Nu6_len_kb; sd(globalNu6_genome/Nu6_len_kb)
min(globalNu6_genome); isonames$isolate[which.min(globalNu6_genome)] #minimal globalSNP count
#report
#min(globalNu6_genome)/Nu6_len_kb
min(globalNu6_genome[-which.min(globalNu6_genome)]); isonames$isolate[which.min(globalNu6_genome[-which.min(globalNu6_genome)])] #minimal globalSNP count except DAOM
#report
#min(globalNu6_genome[-which.min(globalNu6_genome)])/Nu6_len_kb
max(globalNu6_genome); isonames$isolate[which.max(globalNu6_genome)] #maximal globalSNP count
#report
#max(globalNu6_genome)/Nu6_len_kb
# SNP frequency
globalNu6_SNPfreq=globalNu6_SNPcount/globalNu6_genome
globalNu6_SNPfreq_noDAOM=globalNu6_SNPfreq[-18,]
max(as.matrix(globalNu6_SNPfreq))
min(as.matrix(globalNu6_SNPfreq))
min(as.matrix(globalNu6_SNPfreq_noDAOM))
max(as.matrix(globalNu6_SNPfreq_noDAOM))
# 
mean(monoNu6_genome/globalNu6_genome); sd(monoNu6_genome/globalNu6_genome)
mean(monoNu6_genome)/Nu6_len_kb; sd(monoNu6_genome/Nu6_len_kb)
mean(polyNu6_genome)/Nu6_len_kb; sd(polyNu6_genome/Nu6_len_kb)
# correlation 
poly_N_dens=polyNu6_genome/Nu6_len_kb; mono_N_dens=monoNu6_genome/Nu6_len_kb
dfN=data.frame(mono_N_dens, poly_N_dens); rownames(dfN)=isonames$isolate
cor.test(mono_N_dens, poly_N_dens)
plot(poly_N_dens~mono_N_dens , main="SNP density: monoallelic vs polyallelic SNPs (Nu6 reference)", xlab="Monoallelic SNPs density", ylab="Polyallelic SNPs density")
mtext("cor=0.53, p.value=0.002")
with(dfN, text(poly_N_dens~mono_N_dens, labels=rownames(dfN), pos=4, cex=0.8))


## SNP count graphics
par(mfrow=c(1,3))
boxplot(monoNu6_SNPcount, ylim=c(0, 4000), main="monoSNP count against N6 reference"); boxplot(polyNu6_SNPcount, ylim=c(0, 4000),  main="polySNP count against N6 reference"); boxplot(globalNu6_SNPcount, ylim=c(0, 4000),  main="SNP count against N6 reference")
par(mfrow=c(1,1))
data_monoNu6=prop.table(as.matrix(monoNu6_SNPcount), 1) 
barchart(data_monoNu6, scales=list(y=list(labels=monoNu6$isolate)), auto.key=list(space='right'), ylab="Isolates", main="Distribution of monoSNPs in Nu6", axes=FALSE)
data_polyNu6=prop.table(as.matrix(polyNu6_SNPcount), 1) 
barchart(data_polyNu6, scales=list(y=list(labels=polyNu6$isolate)), auto.key=list(space='right'), ylab="Isolates", main="Distribution of polySNPs in Nu6", axes=FALSE)
data_globalNu6=prop.table(as.matrix(globalNu6_SNPcount), 1) 
barchart(data_globalNu6, scales=list(y=list(labels=globalNu6$isolate)), auto.key=list(space='right'), ylab="Isolates", main="Distribution of SNPs in Nu6", axes=FALSE)


## SNP density processing
# mono
monoNu6_SNPdensity_region=monoNu6_SNPcount/Nu6_len_kb #Nu6 SNP density per genome region
monoNu6_SNPdensity_genome=rowSums(monoNu6_SNPdensity_region) #Nu6 SNP density genomewide
monoNu6_SNPdensity_genometable=data.frame(isonames,monoNu6_SNPdensity_genome) #add isolates names
write.table(monoNu6_SNPdensity_genometable, "monoNu6_density") #write to file: import for final report 
#poly
polyNu6_SNPdensity_region=polyNu6_SNPcount/Nu6_len_kb #Nu6 SNP density per genome region
polyNu6_SNPdensity_genome=rowSums(polyNu6_SNPdensity_region) #Nu6 SNP density genomewide
polyNu6_SNPdensity_genometable=data.frame(isonames,polyNu6_SNPdensity_genome) #add isolates names
write.table(polyNu6_SNPdensity_genometable, "polyNu6_density") #write to file: import for final report 
#global
globalNu6_SNPdensity_region=monoNu6_SNPdensity_region+polyNu6_SNPdensity_region
globalNu6_SNPdensity_genome=monoNu6_SNPdensity_genome+polyNu6_SNPdensity_genome
globalNu6_SNPdensity_genometable=data.frame(isonames, globalNu6_SNPdensity_genome)
write.table(globalNu6_SNPdensity_genometable, "globalNu6_density")


## SNP density results 
#mono
mean(monoNu6_SNPdensity_genome); sd(monoNu6_SNPdensity_genome)
min(monoNu6_SNPdensity_genome); isonames$isolate[which.min(monoNu6_SNPdensity_genome)]
min(monoNu6_SNPdensity_genome[-which.min(monoNu6_SNPdensity_genome)]); isonames$isolate[which.min(monoNu6_SNPdensity_genome[-which.min(monoNu6_SNPdensity_genome)])]
max(monoNu6_SNPdensity_genome); isonames$isolate[which.max(monoNu6_SNPdensity_genome)]
#poly
mean(polyNu6_SNPdensity_genome); sd(polyNu6_SNPdensity_genome)
min(polyNu6_SNPdensity_genome); isonames$isolate[which.min(polyNu6_SNPdensity_genome)]
min(polyNu6_SNPdensity_genome[-which.min(polyNu6_SNPdensity_genome)]); isonames$isolate[which.min(polyNu6_SNPdensity_genome[-which.min(polyNu6_SNPdensity_genome)])]
max(polyNu6_SNPdensity_genome); isonames$isolate[which.max(polyNu6_SNPdensity_genome)]#global
#global
mean(globalNu6_SNPdensity_genome); sd(globalNu6_SNPdensity_genome)
min(globalNu6_SNPdensity_genome); isonames$isolate[which.min(globalNu6_SNPdensity_genome)]
min(globalNu6_SNPdensity_genome[-which.min(globalNu6_SNPdensity_genome)]); isonames$isolate[which.min(globalNu6_SNPdensity_genome[-which.min(globalNu6_SNPdensity_genome)])]
max(globalNu6_SNPdensity_genome); isonames$isolate[which.max(globalNu6_SNPdensity_genome)]


## SNP clusters
#k-means on global
fit_globalNu6 <- kmeans(globalNu6_SNPcount, 6)
rownames(globalNu6_SNPcount)=isonames$isolate
clusplot(globalNu6_SNPcount, fit_globalNu6$cluster,  color=TRUE,  labels=3, lines=0, main="Isolate clustering using Nu6 SNP count")


#
# step 2. analyze Rhiir2 SNP data 
#


## import files 
monoRhiir2=read.table("isolateRhiir2_monoRS-sort", h=T) #read Rhiir2 monoSNPs file
polyRhiir2=read.table("isolateRhiir2_polyRS-sort", h=T) #read Rhiir2 polySNPs file 


## SNP count processing
#mono
monoRhiir2_SNPcount=monoRhiir2[2:5] #retrieve monoSNP count
#poly
polyRhiir2_SNPcount=polyRhiir2[2:5] #retrieve polySNP count
#global
globalRhiir2_SNPcount=monoRhiir2_SNPcount+polyRhiir2_SNPcount #retrieve globalSNP count 
globalRhiir2=data.frame(isonames, globalRhiir2_SNPcount) #add isolates names
write.table(globalRhiir2, "globalRhiir2") #write to file: import for final report 
#isolates_list
isonames=monoRhiir2[1]; isonames #retrieve isolates order
rownames(isonames)=isonames$isolate


## SNP count results
#mono
monoRhiir2_genome=rowSums(monoRhiir2_SNPcount)
mean(monoRhiir2_genome); sd(monoRhiir2_genome)
min(monoRhiir2_genome); isonames$isolate[which.min(monoRhiir2_genome)] #minimal monoSNP count
min(monoRhiir2_genome[-which.min(monoRhiir2_genome)]); isonames$isolate[which.min(monoRhiir2_genome[-which.min(monoRhiir2_genome)])] #minimal monoSNP count except DAOM
max(monoRhiir2_genome); isonames$isolate[which.max(monoRhiir2_genome)] #maximal monoSNP count
#poly
polyRhiir2_genome=rowSums(polyRhiir2_SNPcount)
mean(polyRhiir2_genome); sd(polyRhiir2_genome)
min(polyRhiir2_genome); isonames$isolate[which.min(polyRhiir2_genome)] #minimal polySNP count
min(polyRhiir2_genome[-which.min(polyRhiir2_genome)]); isonames$isolate[which.min(polyRhiir2_genome[-which.min(polyRhiir2_genome)])] #minimal polySNP count except DAOM
max(polyRhiir2_genome); isonames$isolate[which.max(polyRhiir2_genome)] #maximal polySNP count
#global
globalRhiir2_genome=rowSums(globalRhiir2_SNPcount)
mean(globalRhiir2_genome); sd(globalRhiir2_genome)
#report
#mean(globalRhiir2_genome)/Rhiir2_len_kb; sd(globalRhiir2_genome/Rhiir2_len_kb)
min(globalRhiir2_genome); isonames$isolate[which.min(globalRhiir2_genome)] #minimal globalSNP count
min(globalRhiir2_genome)/Rhiir2_len_kb
min(globalRhiir2_genome[-which.min(globalRhiir2_genome)]); isonames$isolate[which.min(globalRhiir2_genome[-which.min(globalRhiir2_genome)])] #minimal globalSNP count except DAOM
min(globalRhiir2_genome[-which.min(globalRhiir2_genome)])/Rhiir2_len_kb
max(globalRhiir2_genome); isonames$isolate[which.max(globalRhiir2_genome)] #maximal globalSNP count
max(globalRhiir2_genome)/Rhiir2_len_kb
# SNP frequency 
mean(globalRhiir2_genome); sd(globalRhiir2_genome)
globalRhiir2_SNPfreq=globalRhiir2_SNPcount/globalRhiir2_genome
globalRhiir2_SNPfreq_noDAOM=globalRhiir2_SNPfreq[-18,]
max(as.matrix(globalRhiir2_SNPfreq))
min(as.matrix(globalRhiir2_SNPfreq))
min(as.matrix(globalRhiir2_SNPfreq_noDAOM))
max(as.matrix(globalRhiir2_SNPfreq_noDAOM))
# % of monoallelic SNPs
mean(monoRhiir2_genome/globalRhiir2_genome); sd(monoRhiir2_genome/globalRhiir2_genome)
mean(monoRhiir2_genome)/Rhiir2_len_kb; sd(monoRhiir2_genome/Rhiir2_len_kb)
mean(polyRhiir2_genome)/Rhiir2_len_kb; sd(polyRhiir2_genome/Rhiir2_len_kb)
#correlation 
poly_R_dens=polyRhiir2_genome/Rhiir2_len_kb; mono_R_dens=monoRhiir2_genome/Rhiir2_len_kb
dfR=data.frame(mono_R_dens, poly_R_dens); rownames(dfR)=isonames$isolate
cor.test(mono_R_dens, poly_R_dens)
plot(poly_R_dens~mono_R_dens , main="SNP density: monoallelic vs polyallelic SNPs (Rhiir2 reference)", xlab="Monoallelic SNPs density", ylab="Polyallelic SNPs density")
mtext("cor=0.74, p.value<0.001")
with(dfR, text(poly_R_dens~mono_R_dens, labels=rownames(dfR), pos=4, cex=0.8))

## SNP count graphics
par(mfrow=c(1,3))
boxplot(monoRhiir2_SNPcount, ylim=c(0, 4000), main="monoSNP count against Rhiir2 reference"); boxplot(polyRhiir2_SNPcount, ylim=c(0, 4000),  main="polySNP count against Rhiir2 reference"); boxplot(globalRhiir2_SNPcount, ylim=c(0, 4000),  main="SNP count against Rhiir2 reference")
par(mfrow=c(1,1))
data_monoRhiir2=prop.table(as.matrix(monoRhiir2_SNPcount), 1) 
barchart(data_monoRhiir2, scales=list(y=list(labels=monoRhiir2$isolate)), auto.key=list(space='right'), ylab="Isolates", main="Distribution of monoSNPs in Rhiir2", axes=FALSE)
data_polyRhiir2=prop.table(as.matrix(polyRhiir2_SNPcount), 1) 
barchart(data_polyRhiir2, scales=list(y=list(labels=polyRhiir2$isolate)), auto.key=list(space='right'), ylab="Isolates", main="Distribution of polySNPs in Rhiir2", axes=FALSE)
data_globalRhiir2=prop.table(as.matrix(globalRhiir2_SNPcount), 1) 
barchart(data_globalRhiir2, scales=list(y=list(labels=globalRhiir2$isolate)), auto.key=list(space='right'), ylab="Isolates", main="Distribution of SNPs in Rhiir2", axes=FALSE)


## SNP density processing
# mono
monoRhiir2_SNPdensity_region=monoRhiir2_SNPcount/Rhiir2_len_kb #Rhiir2 SNP density per genome region
monoRhiir2_SNPdensity_genome=rowSums(monoRhiir2_SNPdensity_region) #Rhiir2 SNP density genomewide
monoRhiir2_SNPdensity_genometable=data.frame(isonames,monoRhiir2_SNPdensity_genome) #add isolates names
write.table(monoRhiir2_SNPdensity_genometable, "monoRhiir2_density") #write to file: import for final report 
#poly
polyRhiir2_SNPdensity_region=polyRhiir2_SNPcount/Rhiir2_len_kb #Rhiir2 SNP density per genome region
polyRhiir2_SNPdensity_genome=rowSums(polyRhiir2_SNPdensity_region) #Rhiir2 SNP density genomewide
polyRhiir2_SNPdensity_genometable=data.frame(isonames,polyRhiir2_SNPdensity_genome) #add isolates names
write.table(polyRhiir2_SNPdensity_genometable, "polyRhiir2_density") #write to file: import for final report 
#global
globalRhiir2_SNPdensity_region=monoRhiir2_SNPdensity_region+polyRhiir2_SNPdensity_region
globalRhiir2_SNPdensity_genome=monoRhiir2_SNPdensity_genome+polyRhiir2_SNPdensity_genome
globalRhiir2_SNPdensity_genometable=data.frame(isonames, globalRhiir2_SNPdensity_genome)
write.table(globalRhiir2_SNPdensity_genometable, "globalRhiir2_density")


## SNP density results 
#mono
mean(monoRhiir2_SNPdensity_genome); sd(monoRhiir2_SNPdensity_genome)
min(monoRhiir2_SNPdensity_genome); isonames$isolate[which.min(monoRhiir2_SNPdensity_genome)]
min(monoRhiir2_SNPdensity_genome[-which.min(monoRhiir2_SNPdensity_genome)]); isonames$isolate[which.min(monoRhiir2_SNPdensity_genome[-which.min(monoRhiir2_SNPdensity_genome)])]
max(monoRhiir2_SNPdensity_genome); isonames$isolate[which.max(monoRhiir2_SNPdensity_genome)]
#poly
mean(polyRhiir2_SNPdensity_genome); sd(polyRhiir2_SNPdensity_genome)
min(polyRhiir2_SNPdensity_genome); isonames$isolate[which.min(polyRhiir2_SNPdensity_genome)]
min(polyRhiir2_SNPdensity_genome[-which.min(polyRhiir2_SNPdensity_genome)]); isonames$isolate[which.min(polyRhiir2_SNPdensity_genome[-which.min(polyRhiir2_SNPdensity_genome)])]
max(polyRhiir2_SNPdensity_genome); isonames$isolate[which.max(polyRhiir2_SNPdensity_genome)]#global
#global
mean(globalRhiir2_SNPdensity_genome); sd(globalRhiir2_SNPdensity_genome)
min(globalRhiir2_SNPdensity_genome); isonames$isolate[which.min(globalRhiir2_SNPdensity_genome)]
min(globalRhiir2_SNPdensity_genome[-which.min(globalRhiir2_SNPdensity_genome)]); isonames$isolate[which.min(globalRhiir2_SNPdensity_genome[-which.min(globalRhiir2_SNPdensity_genome)])]
max(globalRhiir2_SNPdensity_genome); isonames$isolate[which.max(globalRhiir2_SNPdensity_genome)]


## SNP clusters
#k-means on global
fit_globalRhiir2 <- kmeans(globalRhiir2_SNPcount, 6)
rownames(globalRhiir2_SNPcount)=isonames$isolate
clusplot(globalRhiir2_SNPcount, fit_globalRhiir2$cluster, color=TRUE,  labels=3, lines=0, main="Isolate clustering using Rhiir2 SNP count")
