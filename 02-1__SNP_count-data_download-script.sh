#!/bin/bash


#####
#
#  Author: Amina Echchiki [mailto:amina.echchiki@unil.ch]
#  Affiliation: MSc student, MLS bioinformatics, University of Lausanne
#  Date: 2015-12-10
#  Name: SNP_count-data_download-script.sh
#  Purpose: Retrieve SNP count data from remote, prior to R script launching
#  Comments: /
#  Testing platforms: desktop Linux 3.16.0-50-generic (Debian 7.6); server Linux 2.6.32-504.el6.x86_64
#  Shell version: GNU bash, version 4.3.11(1)-release (x86_64-pc-linux-gnu)
#
#####


#get formatted SNP count data per genomic region and data for analysis

mkdir ./SNP_count-data; mv SNP_count-data
wget https://copy.com/vnXBpjpowyx0k0Yv
mv vnXBpjpowyx0k0Yv vnXBpjpowyx0k0Yv.zip
unzip vnXBpjpowyx0k0Yv.zip 
ls

## ls should read: 
# isolateNu6_monoRS-sort
# isolateNu6_polyRS-sort
# isolateRhiir2_monoRS-sort
# isolateRhiir2_polyRS-sort
