#!/bin/bash

#####
#
#  Author: Amina Echchiki [mailto:amina.echchiki@unil.ch]
#  Affiliation: MSc student, MLS bioinformatics, University of Lausanne
#  Date: 2015-11-24 
#  Name: Nu6mono_SNP-density.sh
#  Purpose: Analyse SNPs in R.irregularis (ref Nu6 and Rhiir2 genome) to calculate SNP distribution over the 20 most significant scaffolds
#  Comments: /
#  Testing platforms: desktop Linux 3.16.0-50-generic (Debian 7.6); server Linux 2.6.32-504.el6.x86_64
#  Shell version: GNU bash, version 4.3.11(1)-release (x86_64-pc-linux-gnu)
#
#####


#
# step 1.0: get monoNu6 files, prepare wdir
#

# prepare working directory
mkdir ./20-scaffold_SNP-density
mkdir ./Nu6mono_SNP-density; cd ./Nu6mono_SNP-density
# prepare table file 
echo -e "Snb\tREPCOD\tREPNCOD\tNREPCOD\tNREPNCOD" > Nu6mono

# download data 
wget https://copy.com/amb9NDjJWXFhtWHm
mv amb9NDjJWXFhtWHm amb9NDjJWXFhtWHm.zip
unzip amb9NDjJWXFhtWHm.zip
 
# download script 
wget https://copy.com/nWBcXYj4gUtwC9Wj
mv nWBcXYj4gUtwC9Wj nWBcXYj4gUtwC9Wj.zip
unzip nWBcXYj4gUtwC9Wj.zip
# prepare execution scripts 
cp SNP_density_Nu6mono.sh SNP_density_Nu6mono_edit.sh 
chmod +x SNP_density_Nu6mono.sh; chmod +x SNP_density_Nu6mono_edit.sh


#
# step 1.1: retrieve monoNu6 SNP count
#

echo "ANALYZING Nu6 monoSNP..."; echo "---"
./SNP_density_Nu6mono.sh; echo "ANALYSING SCAFFOLD 1..."
sed -i 's/1/2/g' SNP_density_Nu6mono_edit.sh; echo "ANALYSING SCAFFOLD 2..."
./SNP_density_Nu6mono_edit.sh
sed -i 's/2/3/g' SNP_density_Nu6mono_edit.sh ; echo "ANALYSING SCAFFOLD 3..."
./SNP_density_Nu6mono_edit.sh 
sed -i 's/Scaffold3/Scaffold4/g' SNP_density_Nu6mono_edit.sh ; sed -i 's/S3/S4/g' SNP_density_Nu6mono_edit.sh ; echo "ANALYSING SCAFFOLD 4..."
./SNP_density_Nu6mono_edit.sh 
sed -i 's/4/5/g' SNP_density_Nu6mono_edit.sh; echo "ANALYSING SCAFFOLD 5..."
./SNP_density_Nu6mono_edit.sh
sed -i 's/5/6/g' SNP_density_Nu6mono_edit.sh; echo "ANALYSING SCAFFOLD 6..."
./SNP_density_Nu6mono_edit.sh
sed -i 's/Scaffold6/Scaffold7/g' SNP_density_Nu6mono_edit.sh ; sed -i 's/S6/S7/g' SNP_density_Nu6mono_edit.sh; echo "ANALYSING SCAFFOLD 7..."
./SNP_density_Nu6mono_edit.sh
sed -i 's/7/8/g' SNP_density_Nu6mono_edit.sh; echo "ANALYSING SCAFFOLD 8..."
./SNP_density_Nu6mono_edit.sh
sed -i 's/8/9/g' SNP_density_Nu6mono_edit.sh; echo "ANALYSING SCAFFOLD 9..."
./SNP_density_Nu6mono_edit.sh
sed -i 's/9/10/g' SNP_density_Nu6mono_edit.sh; echo "ANALYSING SCAFFOLD 10..."
./SNP_density_Nu6mono_edit.sh
sed -i 's/10/11/g' SNP_density_Nu6mono_edit.sh; echo "ANALYSING SCAFFOLD 11..."
./SNP_density_Nu6mono_edit.sh
sed -i 's/11/12/g' SNP_density_Nu6mono_edit.sh; echo "ANALYSING SCAFFOLD 12..."
./SNP_density_Nu6mono_edit.sh
sed -i 's/12/13/g' SNP_density_Nu6mono_edit.sh; echo "ANALYSING SCAFFOLD 13..."
./SNP_density_Nu6mono_edit.sh
sed -i 's/13/14/g' SNP_density_Nu6mono_edit.sh; echo "ANALYSING SCAFFOLD 14..."
./SNP_density_Nu6mono_edit.sh
sed -i 's/14/15/g' SNP_density_Nu6mono_edit.sh; echo "ANALYSING SCAFFOLD 15..."
./SNP_density_Nu6mono_edit.sh
sed -i 's/15/16/g' SNP_density_Nu6mono_edit.sh; echo "ANALYSING SCAFFOLD 16..."
./SNP_density_Nu6mono_edit.sh
sed -i 's/16/17/g' SNP_density_Nu6mono_edit.sh; echo "ANALYSING SCAFFOLD 17..."
./SNP_density_Nu6mono_edit.sh
sed -i 's/17/18/g' SNP_density_Nu6mono_edit.sh; echo "ANALYSING SCAFFOLD 18..."
./SNP_density_Nu6mono_edit.sh
sed -i 's/18/19/g' SNP_density_Nu6mono_edit.sh; echo "ANALYSING SCAFFOLD 19..."
./SNP_density_Nu6mono_edit.sh
sed -i 's/19/20/g' SNP_density_Nu6mono_edit.sh; echo "ANALYSING SCAFFOLD 20..."
./SNP_density_Nu6mono_edit.sh
echo "...Nu6 monoSNP ANALYZED!"; echo "---"

# copy count file for analysis
cp ./Nu6mono ../20-scaffold_SNP-density
cd ../


#
# step 2.0: get polyNu6 files, prepare wdir
#

# prepare working directory

mkdir ./Nu6poly_SNP-density; cd ./Nu6poly_SNP-density
# prepare table file 
echo -e "Snb\tREPCOD\tREPNCOD\tNREPCOD\tNREPNCOD" > Nu6poly

# download data 
wget https://copy.com/OLoKMUbHq6ol0EcT
mv OLoKMUbHq6ol0EcT OLoKMUbHq6ol0EcT.zip
unzip OLoKMUbHq6ol0EcT.zip
 
# download script 
wget https://copy.com/M4EtH2uqwLJmfT6l
mv M4EtH2uqwLJmfT6l M4EtH2uqwLJmfT6l.zip
unzip M4EtH2uqwLJmfT6l.zip

# prepare execution scripts 
cp SNP_density_Nu6poly.sh SNP_density_Nu6poly_edit.sh 
chmod +x SNP_density_Nu6poly.sh; chmod +x SNP_density_Nu6poly_edit.sh

#
# step 2.1: retrieve polyNu6 SNP count 
#

echo "ANALYZING Nu6 polySNP..."; echo "---"
./SNP_density_Nu6poly.sh; echo "ANALYSING SCAFFOLD 1..."
sed -i 's/1/2/g' SNP_density_Nu6poly_edit.sh; echo "ANALYSING SCAFFOLD 2..."
./SNP_density_Nu6poly_edit.sh
sed -i 's/2/3/g' SNP_density_Nu6poly_edit.sh; echo "ANALYSING SCAFFOLD 3..."
./SNP_density_Nu6poly_edit.sh 
sed -i 's/Scaffold3/Scaffold4/g' SNP_density_Nu6poly_edit.sh ; sed -i 's/S3/S4/g' SNP_density_Nu6poly_edit.sh; echo "ANALYSING SCAFFOLD 4..."
./SNP_density_Nu6poly_edit.sh 
sed -i 's/4/5/g' SNP_density_Nu6poly_edit.sh; echo "ANALYSING SCAFFOLD 5..."
./SNP_density_Nu6poly_edit.sh
sed -i 's/5/6/g' SNP_density_Nu6poly_edit.sh; echo "ANALYSING SCAFFOLD 6..."
./SNP_density_Nu6poly_edit.sh
sed -i 's/Scaffold6/Scaffold7/g' SNP_density_Nu6poly_edit.sh ; sed -i 's/S6/S7/g' SNP_density_Nu6poly_edit.sh ; echo "ANALYSING SCAFFOLD 7..."
./SNP_density_Nu6poly_edit.sh
sed -i 's/7/8/g' SNP_density_Nu6poly_edit.sh; echo "ANALYSING SCAFFOLD 8..."
./SNP_density_Nu6poly_edit.sh
sed -i 's/8/9/g' SNP_density_Nu6poly_edit.sh; echo "ANALYSING SCAFFOLD 9..."
./SNP_density_Nu6poly_edit.sh
sed -i 's/9/10/g' SNP_density_Nu6poly_edit.sh; echo "ANALYSING SCAFFOLD 10..."
./SNP_density_Nu6poly_edit.sh
sed -i 's/10/11/g' SNP_density_Nu6poly_edit.sh; echo "ANALYSING SCAFFOLD 11..."
./SNP_density_Nu6poly_edit.sh
sed -i 's/11/12/g' SNP_density_Nu6poly_edit.sh; echo "ANALYSING SCAFFOLD 12..."
./SNP_density_Nu6poly_edit.sh
sed -i 's/12/13/g' SNP_density_Nu6poly_edit.sh; echo "ANALYSING SCAFFOLD 13..."
./SNP_density_Nu6poly_edit.sh
sed -i 's/13/14/g' SNP_density_Nu6poly_edit.sh; echo "ANALYSING SCAFFOLD 14..."
./SNP_density_Nu6poly_edit.sh
sed -i 's/14/15/g' SNP_density_Nu6poly_edit.sh; echo "ANALYSING SCAFFOLD 15..."
./SNP_density_Nu6poly_edit.sh
sed -i 's/15/16/g' SNP_density_Nu6poly_edit.sh; echo "ANALYSING SCAFFOLD 16..."
./SNP_density_Nu6poly_edit.sh
sed -i 's/16/17/g' SNP_density_Nu6poly_edit.sh; echo "ANALYSING SCAFFOLD 17..."
./SNP_density_Nu6poly_edit.sh
sed -i 's/17/18/g' SNP_density_Nu6poly_edit.sh; echo "ANALYSING SCAFFOLD 18..."
./SNP_density_Nu6poly_edit.sh
sed -i 's/18/19/g' SNP_density_Nu6poly_edit.sh; echo "ANALYSING SCAFFOLD 19..."
./SNP_density_Nu6poly_edit.sh
sed -i 's/19/20/g' SNP_density_Nu6poly_edit.sh; echo "ANALYSING SCAFFOLD 20..."
./SNP_density_Nu6poly_edit.sh
echo "...Nu6 polySNP ANALYZED!"; echo "---"

# copy count file for analysis
cp ./Nu6poly ../20-scaffold_SNP-density
cd ../


#
# step 3.0: get monoRhiir2 files, prepare wdir
#

# prepare working directory
mkdir ./Rhiir2mono_SNP-density; cd ./Rhiir2mono_SNP-density
# prepare table file 
echo -e "Snb\tREPCOD\tREPNCOD\tNREPCOD\tNREPNCOD" > Rhiir2mono

# download data 
wget https://copy.com/ozazXhF8PbSxyBVm
mv ozazXhF8PbSxyBVm ozazXhF8PbSxyBVm.zip
unzip ozazXhF8PbSxyBVm.zip

# download script 
wget https://copy.com/rGKZFZvXlB8lfrpa
mv rGKZFZvXlB8lfrpa rGKZFZvXlB8lfrpa.zip
unzip rGKZFZvXlB8lfrpa.zip

# prepare execution scripts 
cp SNP_density_Rhiir2mono.sh SNP_density_Rhiir2mono_edit.sh 
chmod +x SNP_density_Rhiir2mono.sh; chmod +x SNP_density_Rhiir2mono_edit.sh


#
# step 3.1: retrieve SNP count
#

echo "ANALYZING Rhiir2 monoSNP...."; echo "---"
./SNP_density_Rhiir2mono.sh; echo "ANALYSING SCAFFOLD 1..."
sed -i 's/1/2/g' SNP_density_Rhiir2mono_edit.sh; echo "ANALYSING SCAFFOLD 2..."
./SNP_density_Rhiir2mono_edit.sh
sed -i 's/scaffold2/scaffold3/g' SNP_density_Rhiir2mono_edit.sh; sed -i 's/S2/S3/g' SNP_density_Rhiir2mono_edit.sh; echo "ANALYSING SCAFFOLD 3..."
./SNP_density_Rhiir2mono_edit.sh 
sed -i 's/scaffold3/scaffold4/g' SNP_density_Rhiir2mono_edit.sh ; sed -i 's/S3/S4/g' SNP_density_Rhiir2mono_edit.sh ; echo "ANALYSING SCAFFOLD 4..."
./SNP_density_Rhiir2mono_edit.sh 
sed -i 's/4/5/g' SNP_density_Rhiir2mono_edit.sh; echo "ANALYSING SCAFFOLD 5..."
./SNP_density_Rhiir2mono_edit.sh
sed -i 's/5/6/g' SNP_density_Rhiir2mono_edit.sh; echo "ANALYSING SCAFFOLD 6..."
./SNP_density_Rhiir2mono_edit.sh
sed -i 's/6/7/g' SNP_density_Rhiir2mono_edit.sh; echo "ANALYSING SCAFFOLD 7..."
./SNP_density_Rhiir2mono_edit.sh
sed -i 's/7/8/g' SNP_density_Rhiir2mono_edit.sh; echo "ANALYSING SCAFFOLD 8..."
./SNP_density_Rhiir2mono_edit.sh
sed -i 's/8/9/g' SNP_density_Rhiir2mono_edit.sh; echo "ANALYSING SCAFFOLD 9..."
./SNP_density_Rhiir2mono_edit.sh
sed -i 's/9/10/g' SNP_density_Rhiir2mono_edit.sh; echo "ANALYSING SCAFFOLD 10..."
./SNP_density_Rhiir2mono_edit.sh
sed -i 's/10/11/g' SNP_density_Rhiir2mono_edit.sh; echo "ANALYSING SCAFFOLD 11..."
./SNP_density_Rhiir2mono_edit.sh
sed -i 's/11/12/g' SNP_density_Rhiir2mono_edit.sh; echo "ANALYSING SCAFFOLD 12..."
./SNP_density_Rhiir2mono_edit.sh
sed -i 's/12/13/g' SNP_density_Rhiir2mono_edit.sh; echo "ANALYSING SCAFFOLD 13..."
./SNP_density_Rhiir2mono_edit.sh
sed -i 's/13/14/g' SNP_density_Rhiir2mono_edit.sh; echo "ANALYSING SCAFFOLD 14..."
./SNP_density_Rhiir2mono_edit.sh
sed -i 's/14/15/g' SNP_density_Rhiir2mono_edit.sh; echo "ANALYSING SCAFFOLD 15..."
./SNP_density_Rhiir2mono_edit.sh
sed -i 's/15/16/g' SNP_density_Rhiir2mono_edit.sh; echo "ANALYSING SCAFFOLD 16..."
./SNP_density_Rhiir2mono_edit.sh
sed -i 's/16/17/g' SNP_density_Rhiir2mono_edit.sh; echo "ANALYSING SCAFFOLD 17..."
./SNP_density_Rhiir2mono_edit.sh
sed -i 's/17/18/g' SNP_density_Rhiir2mono_edit.sh; echo "ANALYSING SCAFFOLD 18..."
./SNP_density_Rhiir2mono_edit.sh
sed -i 's/18/19/g' SNP_density_Rhiir2mono_edit.sh; echo "ANALYSING SCAFFOLD 19..."
./SNP_density_Rhiir2mono_edit.sh
sed -i 's/19/20/g' SNP_density_Rhiir2mono_edit.sh; echo "ANALYSING SCAFFOLD 20..."
./SNP_density_Rhiir2mono_edit.sh
echo "...Rhiir2 monoSNP ANALYZED!"; echo "---"

# copy count file for analysis
cp ./Rhiir2mono ../20-scaffold_SNP-density
cd ../


#
# step 4.0: get polyRhiir2 files, prepare wdir
#

# prepare working directory
mkdir ./Rhiir2poly_SNP-density; cd ./Rhiir2poly_SNP-density
# prepare table file 
echo -e "Snb\tREPCOD\tREPNCOD\tNREPCOD\tNREPNCOD" > Rhiir2poly

# download data 
wget https://copy.com/pLe9dgAOL2tn7AJX
mv pLe9dgAOL2tn7AJX pLe9dgAOL2tn7AJX.zip
unzip pLe9dgAOL2tn7AJX.zip

# download script 
wget https://copy.com/LLdiH87fFjuGtTNC
mv LLdiH87fFjuGtTNC LLdiH87fFjuGtTNC.zip
unzip LLdiH87fFjuGtTNC.zip

# prepare execution scripts 
cp SNP_density_Rhiir2poly.sh SNP_density_Rhiir2poly_edit.sh 
chmod +x SNP_density_Rhiir2poly.sh; chmod +x SNP_density_Rhiir2poly_edit.sh


#
# step 4.1: retrieve polyRhiir2 SNP count
#

echo "ANALYZING Rhiir2 polySNP...."; echo "---" 
./SNP_density_Rhiir2poly.sh; echo "ANALYSING SCAFFOLD 1..."
sed -i 's/1/2/g' SNP_density_Rhiir2poly_edit.sh;echo "ANALYSING SCAFFOLD 2..."
./SNP_density_Rhiir2poly_edit.sh
sed -i 's/scaffold2/scaffold3/g' SNP_density_Rhiir2poly_edit.sh; sed -i 's/S2/S3/g' SNP_density_Rhiir2poly_edit.sh; echo "ANALYSING SCAFFOLD 3..."
./SNP_density_Rhiir2poly_edit.sh 
sed -i 's/scaffold3/scaffold4/g' SNP_density_Rhiir2poly_edit.sh ; sed -i 's/S3/S4/g' SNP_density_Rhiir2poly_edit.sh ; echo "ANALYSING SCAFFOLD 4..."
./SNP_density_Rhiir2poly_edit.sh 
sed -i 's/4/5/g' SNP_density_Rhiir2poly_edit.sh; echo "ANALYSING SCAFFOLD 5..."
./SNP_density_Rhiir2poly_edit.sh
sed -i 's/5/6/g' SNP_density_Rhiir2poly_edit.sh; echo "ANALYSING SCAFFOLD 6..."
./SNP_density_Rhiir2poly_edit.sh
sed -i 's/6/7/g' SNP_density_Rhiir2poly_edit.sh; echo "ANALYSING SCAFFOLD 7..."
./SNP_density_Rhiir2poly_edit.sh
sed -i 's/7/8/g' SNP_density_Rhiir2poly_edit.sh; echo "ANALYSING SCAFFOLD 8..."
./SNP_density_Rhiir2poly_edit.sh
sed -i 's/8/9/g' SNP_density_Rhiir2poly_edit.sh; echo "ANALYSING SCAFFOLD 9..."
./SNP_density_Rhiir2poly_edit.sh
sed -i 's/9/10/g' SNP_density_Rhiir2poly_edit.sh; echo "ANALYSING SCAFFOLD 10..."
./SNP_density_Rhiir2poly_edit.sh
sed -i 's/10/11/g' SNP_density_Rhiir2poly_edit.sh; echo "ANALYSING SCAFFOLD 11..."
./SNP_density_Rhiir2poly_edit.sh
sed -i 's/11/12/g' SNP_density_Rhiir2poly_edit.sh; echo "ANALYSING SCAFFOLD 12..."
./SNP_density_Rhiir2poly_edit.sh
sed -i 's/12/13/g' SNP_density_Rhiir2poly_edit.sh; echo "ANALYSING SCAFFOLD 13..."
./SNP_density_Rhiir2poly_edit.sh
sed -i 's/13/14/g' SNP_density_Rhiir2poly_edit.sh; echo "ANALYSING SCAFFOLD 14..."
./SNP_density_Rhiir2poly_edit.sh
sed -i 's/14/15/g' SNP_density_Rhiir2poly_edit.sh; echo "ANALYSING SCAFFOLD 15..."
./SNP_density_Rhiir2poly_edit.sh
sed -i 's/15/16/g' SNP_density_Rhiir2poly_edit.sh; echo "ANALYSING SCAFFOLD 16..."
./SNP_density_Rhiir2poly_edit.sh
sed -i 's/16/17/g' SNP_density_Rhiir2poly_edit.sh; echo "ANALYSING SCAFFOLD 17..."
./SNP_density_Rhiir2poly_edit.sh
sed -i 's/17/18/g' SNP_density_Rhiir2poly_edit.sh; echo "ANALYSING SCAFFOLD 18..."
./SNP_density_Rhiir2poly_edit.sh
sed -i 's/18/19/g' SNP_density_Rhiir2poly_edit.sh; echo "ANALYSING SCAFFOLD 19..."
./SNP_density_Rhiir2poly_edit.sh
sed -i 's/19/20/g' SNP_density_Rhiir2poly_edit.sh; echo "ANALYSING SCAFFOLD 20..."
./SNP_density_Rhiir2poly_edit.sh

# copy count file for analysis
cp ./Rhiir2poly ../20-scaffold_SNP-density
cd ../


