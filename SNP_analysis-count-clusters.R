#####
#
#  Author: Amina Echchiki [mailto:amina.echchiki@unil.ch]
#  Affiliation: MSc student, MLS bioinformatics, University of Lausanne
#  Date: 2015-11-20 
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


## import files 
monoNu6=read.table("isolateNu6_monoRS-sort", h=T) #read Nu6 monoSNPs file
polyNu6=read.table("isolateNu6_polyRS-sort", h=T) #read Nu6 polySNPs file 


## SNP count processing
#mono
monoNu6_SNPcount=monoNu6[2:5] #retrieve monoSNP count
#poly
polyNu6_SNPcount=polyNu6[2:5] #retrieve polySNP count
#global
globalNu6_SNPcount=monoNu6_SNPcount+polyNu6_SNPcount #retrieve globalSNP count 
globalNu6=data.frame(isonames, globalNu6_SNPcount) #add isolates names
write.table(globalNu6, "globalNu6") #write to file: import for final report 
#isolates_list
isonames=monoNu6[1]; isonames #retrieve isolates order
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
min(globalNu6_genome); isonames$isolate[which.min(globalNu6_genome)] #minimal globalSNP count
min(globalNu6_genome[-which.min(globalNu6_genome)]); isonames$isolate[which.min(globalNu6_genome[-which.min(globalNu6_genome)])] #minimal globalSNP count except DAOM
max(globalNu6_genome); isonames$isolate[which.max(globalNu6_genome)] #maximal globalSNP count


## SNP count graphics
par(mfrow=c(1,3))
boxplot(monoNu6_SNPcount, ylim=c(0, 4000), main="monoSNP count against N6 reference"); boxplot(polyNu6_SNPcount, ylim=c(0, 4000),  main="polySNP count against N6 reference"); boxplot(globalNu6_SNPcount, ylim=c(0, 4000),  main="SNP count against N6 reference")
data_monoNu6=prop.table(as.matrix(monoNu6_SNPcount), 1) 
par(mfrow=c(1,1))
barchart(data_monoNu6, scales=list(y=list(labels=monoNu6$isolate)), auto.key=list(space='right'), ylab="Isolates", main="Distribution of monoSNPs in Nu6", axes=FALSE)
data_polyNu6=prop.table(as.matrix(polyNu6_SNPcount), 1) 
barchart(data_polyNu6, scales=list(y=list(labels=polyNu6$isolate)), auto.key=list(space='right'), ylab="Isolates", main="Distribution of polySNPs in Nu6", axes=FALSE)
data_globalNu6=prop.table(as.matrix(globalNu6_SNPcount), 1) 
barchart(data_globalNu6, scales=list(y=list(labels=globalNu6$isolate)), auto.key=list(space='right'), ylab="Isolates", main="Distribution of SNPs in Nu6", axes=FALSE)


## SNP density processing
Nu6_len=115079203; Nu6_len_kb=Nu6_len/1000 #Nu6 length
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
#mono
distance_monoNu6=dist(monoNu6_SNPcount, method="euclidean")
cluster_monoNu6=hclust(distance_monoNu6, method="average")
par(mfrow=c(1,1)); plot(cluster_monoNu6, hang=1, label=monoNu6$isolates , main="Isolate clustering using Nu6 SNP count")
#poly
distance_polyNu6=dist(polyNu6_SNPcount, method="euclidean")
cluster_polyNu6=hclust(distance_polyNu6, method="average")
par(mfrow=c(1,1)); plot(cluster_polyNu6, hang=1, label=polyNu6$isolates , main="Isolate clustering using Nu6 SNP count")
#global
distance_globalNu6=dist(globalNu6_SNPcount, method="euclidean")
cluster_globalNu6=hclust(distance_globalNu6, method="average")
par(mfrow=c(1,1)); plot(cluster_globalNu6, hang=1, label=globalNu6$isolates , main="Isolate clustering using Nu6 SNP count")
#k-means on global
fit_globalNu6 <- kmeans(globalNu6_SNPcount, 6)
rownames(globalNu6_SNPcount)=isonames$isolate
clusplot(globalNu6_SNPcount, fit_globalNu6$cluster,  color=TRUE,  labels=3, lines=0, main="Isolate clustering using Nu6 SNP count")


#
# step 2. analyze Rhiir2 SNP data 
#


## prepare working directory 
wdir="/home/dovah/Documents/MSc/MSc-A15/FSP/09-12_NEW-SCRIPTS/0R_analysis/" # change string to own dwd directory
setwd(wdir)
library(lattice) #graphical package  
library(cluster) #clustering package 


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
min(globalRhiir2_genome); isonames$isolate[which.min(globalRhiir2_genome)] #minimal globalSNP count
min(globalRhiir2_genome[-which.min(globalRhiir2_genome)]); isonames$isolate[which.min(globalRhiir2_genome[-which.min(globalRhiir2_genome)])] #minimal globalSNP count except DAOM
max(globalRhiir2_genome); isonames$isolate[which.max(globalRhiir2_genome)] #maximal globalSNP count

barchart(data_monoNu6, scales=list(y=list(labels=monoNu6$isolate)), auto.key=list(space='right'), ylab="Isolates", main="Distribution of monoSNPs in Nu6", axes=FALSE)


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
Rhiir2_len=115079203; Rhiir2_len_kb=Rhiir2_len/1000 #Rhiir2 length
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
#mono
distance_monoRhiir2=dist(monoRhiir2_SNPcount, method="euclidean")
cluster_monoRhiir2=hclust(distance_monoRhiir2, method="average")
par(mfrow=c(1,1)); plot(cluster_monoRhiir2, hang=1, label=monoRhiir2$isolates , main="Isolate clustering using Rhiir2 SNP count")
#poly
distance_polyRhiir2=dist(polyRhiir2_SNPcount, method="euclidean")
cluster_polyRhiir2=hclust(distance_polyRhiir2, method="average")
par(mfrow=c(1,1)); plot(cluster_polyRhiir2, hang=1, label=polyRhiir2$isolates , main="Isolate clustering using Rhiir2 SNP count")
#global
distance_globalRhiir2=dist(globalRhiir2_SNPcount, method="euclidean")
cluster_globalRhiir2=hclust(distance_globalRhiir2, method="average")
par(mfrow=c(1,1)); plot(cluster_globalRhiir2, hang=1, label=globalRhiir2$isolates , main="Isolate clustering using Rhiir2 SNP count")
#k-means on global
fit_globalRhiir2 <- kmeans(globalRhiir2_SNPcount, 5)
rownames(globalRhiir2_SNPcount)=isonames$isolate
clusplot(globalRhiir2_SNPcount, fit_globalRhiir2$cluster,  color=TRUE,  labels=3, lines=0, main="Isolate clustering using Rhiir2 SNP count")
