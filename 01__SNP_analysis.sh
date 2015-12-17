#!/bin/bash


#####
#
#  Author: Amina Echchiki [mailto:amina.echchiki@unil.ch]
#  Affiliation: MSc student, MLS bioinformatics, University of Lausanne
#  Date: 2015-12-10
#  Name: SNP-analysis.sh
#  Purpose: Analyse mono- and poly-allelic SNPs in R.irregularis isolates using Nu6 and Rhiir2 as reference genome
#  Comments: requires *.subVCF files of each isolate; requires pwd free-space of aprox. 5GB
#  Testing platforms: desktop Linux 3.16.0-50-generic (Debian 7.6); server Linux 2.6.32-504.el6.x86_64
#  Shell version: GNU bash, version 4.3.11(1)-release (x86_64-pc-linux-gnu)
#
#####


#
# step 1.0: get Nu6 .subVCF files
#      


mkdir ./Nu6_analysis
cd ./Nu6_analysis
# cp /home/dovah/Documents/MSc/MSc-A15/FSP/subVCF/Nu6_subVCF-files.zip .
# unzip Nu6_subVCF-files.zip

wget https://copy.com/CwpIrl9rA4gzC9vx
mv CwpIrl9rA4gzC9vx CwpIrl9rA4gzC9vx.zip
unzip CwpIrl9rA4gzC9vx.zip

#
# step 1.1: prepare N6 subVCF for SNP analysis 
#

echo "ANALYZING Nu6 DATA...."; echo "..."
echo "PREPARING subVCF FOR ANALYSIS..."
# define coverage threshold 
cov=10
# filter subVCF files for coverage ($13) value above threshold > *.cov
ls *.subVCF | while read FN; do awk -v cov="$cov" -F "\t" 'BEGIN{OFS="\t"}{if ($13 > cov) print$0; }' $FN > ${FN/subVCF/cov}; done

# define wanted polymorphisms ($10): S=substitutions, R=same as reference
# filter cov files with R and or S for mutation type ($10) > *.RS
ls *.cov | while read FN; do awk -F "\t" 'BEGIN{OFS="\t"}{if ($10 ~ /[RS]/) print$0; }'  $FN > ${FN/cov/RS}; done
# define unwanted polymorphisms ($10): I=insertions; D=deletions; MNP=multi-nucleotide-polymorphism
# clean *.RS files from unwanted polymorpisms > *.cleanRS
ls *.RS | while read FN; do awk -F "\t" 'BEGIN{OFS="\t"}{if ($10 !~ /[IDM]/) print$0; }'  $FN > ${FN/RS/cleanRS}; done
echo "...DONE!"


#
# step 1.2: define mono and polyallelic SNPs
#


echo "DEFINING MONO/POLYALLELIC SNPs..."
# define polyallelic SNPs: will have more than 1 mutation event reported in $10, field separator within $10: comma (',')
# look for polyallelic SNPs in *.clanRS > *.polyRS 
ls *.cleanRS | while read FN; do awk -F "\t" 'BEGIN{OFS="\t"}{if ($10 ~ /,/) print$0; }'  $FN > ${FN/cleanRS/polyRS}; done
# define monoallelic SNPs: will have exactly 1 mutation event reported in $10, no field separator (',') needed within $10, or no mutation events reported in $10
# look for monoallelic SNPs in *.cleanRS > *.monoRS 
ls *.cleanRS | while read FN; do awk -F "\t" 'BEGIN{OFS="\t"}{if ($10 !~ /,/) print$0; }'  $FN > ${FN/cleanRS/monoRS}; done
echo "...DONE!"


#
# step 1.3: analyse monoallelic SNPs
#


echo "ANALYSING MONOALLELIC SNPs...."
# select interesting info from *.catmonoRS: $1_isolate; $2_scaffoldnb; $4_SNPposition; $7_rep(>=1)-nonrep(0); $8_cod(1)-noncod(0); $9_ref-nt; $10_mutation-type; $11_mutated-nt > *.selectedmonoRS
ls *.monoRS | while read FN; do awk -F "\t" 'BEGIN{OFS="\t"}{print $1, $2, $4, $7, $8, $9, $10, $11 }'  $FN > ${FN/monoRS/selectedmonoRS}; done

# cut 1st column to compare replicates per isolate
ls *.selectedmonoRS | while read FN; do cut -f2- $FN > ${FN/selectedmonoRS/cutmonoRS}; done

# split digit from scaffold string
ls *.cutmonoRS | while read FN; do sed 's/Scaffold/&\t/' $FN > ${FN/cutmonoRS/splitmonoRS}; done

# keep only the SNPs which are conserved among  all the replicates within an isolate 
# uniq: to keep the SNPs once, if it's recorded more than once; sort is a requirement for uniq to be executed; an additional column will be printed in *.allmonoRS $1, depending on how many times the position was found
sort A1*.splitmonoRS |uniq -c > A1.allmonoRS; awk '{if ($1 == "3") print $0;}' A1.allmonoRS > A1.conservedmonoRS
sort A2*.splitmonoRS |uniq -c > A2.allmonoRS; awk '{if ($1 == "3") print $0;}' A2.allmonoRS > A2.conservedmonoRS
sort A3*.splitmonoRS |uniq -c > A3.allmonoRS; awk '{if ($1 == "3") print $0;}' A3.allmonoRS > A3.conservedmonoRS
sort A4*.splitmonoRS |uniq -c > A4.allmonoRS; awk '{if ($1 == "3") print $0;}' A4.allmonoRS > A4.conservedmonoRS
sort A5*.splitmonoRS |uniq -c > A5.allmonoRS; awk '{if ($1 == "3") print $0;}' A5.allmonoRS > A5.conservedmonoRS
sort B10*.splitmonoRS |uniq -c > B10.allmonoRS; awk '{if ($1 == "3") print $0;}' B10.allmonoRS > B10.conservedmonoRS
sort B11*.splitmonoRS |uniq -c > B11.allmonoRS; awk '{if ($1 == "3") print $0;}' B11.allmonoRS > B11.conservedmonoRS
sort B12*.splitmonoRS |uniq -c > B12.allmonoRS; awk '{if ($1 == "3") print $0;}' B12.allmonoRS > B12.conservedmonoRS
sort B14*.splitmonoRS |uniq -c > B14.allmonoRS; awk '{if ($1 == "2") print $0;}' B14.allmonoRS > B14.conservedmonoRS 
sort B15*.splitmonoRS |uniq -c > B15.allmonoRS; awk '{if ($1 == "3") print $0;}' B15.allmonoRS > B15.conservedmonoRS
sort B17*.splitmonoRS |uniq -c > B17.allmonoRS; awk '{if ($1 == "3") print $0;}' B17.allmonoRS > B17.conservedmonoRS
sort B1-*.splitmonoRS |uniq -c > B1.allmonoRS; awk '{if ($1 == "3") print $0;}' B1.allmonoRS > B1.conservedmonoRS
sort B2*.splitmonoRS |uniq -c > B2.allmonoRS; awk '{if ($1 == "3") print $0;}' B2.allmonoRS > B2.conservedmonoRS
sort B3*.splitmonoRS |uniq -c > B3.allmonoRS; awk '{if ($1 == "3") print $0;}' B3.allmonoRS > B3.conservedmonoRS
sort B4*.splitmonoRS |uniq -c > B4.allmonoRS; awk '{if ($1 == "3") print $0;}' B4.allmonoRS > B4.conservedmonoRS
sort B7*.splitmonoRS |uniq -c > B7.allmonoRS; awk '{if ($1 == "3") print $0;}' B7.allmonoRS > B7.conservedmonoRS
sort B8*.splitmonoRS |uniq -c > B8.allmonoRS; awk '{if ($1 == "3") print $0;}' B8.allmonoRS > B8.conservedmonoRS
sort C1*.splitmonoRS |uniq -c > C1.allmonoRS; awk '{if ($1 == "3") print $0;}' C1.allmonoRS > C1.conservedmonoRS
sort C2*.splitmonoRS |uniq -c > C2.allmonoRS; awk '{if ($1 == "3") print $0;}' C2.allmonoRS > C2.conservedmonoRS
sort C3*.splitmonoRS |uniq -c > C3.allmonoRS; awk '{if ($1 == "3") print $0;}' C3.allmonoRS > C3.conservedmonoRS
sort C4*.splitmonoRS |uniq -c > C4.allmonoRS; awk '{if ($1 == "3") print $0;}' C4.allmonoRS > C4.conservedmonoRS
sort C5*.splitmonoRS |uniq -c > C5.allmonoRS; awk '{if ($1 == "3") print $0;}' C5.allmonoRS > C5.conservedmonoRS
sort CAN*.splitmonoRS |uniq -c > CAN.allmonoRS; awk '{if ($1 == "3") print $0;}' CAN.allmonoRS > CAN.conservedmonoRS
sort D1*.splitmonoRS |uniq -c > D1.allmonoRS; awk '{if ($1 == "3") print $0;}' D1.allmonoRS > D1.conservedmonoRS
sort D2*.splitmonoRS |uniq -c > D2.allmonoRS; awk '{if ($1 == "3") print $0;}' D2.allmonoRS > D2.conservedmonoRS
sort D3*.splitmonoRS |uniq -c > D3.allmonoRS; awk '{if ($1 == "3") print $0;}' D3.allmonoRS > D3.conservedmonoRS
sort D4*.splitmonoRS |uniq -c > D4.allmonoRS; awk '{if ($1 == "3") print $0;}' D4.allmonoRS > D4.conservedmonoRS
sort E1*.splitmonoRS |uniq -c > E1.allmonoRS; awk '{if ($1 == "3") print $0;}' E1.allmonoRS > E1.conservedmonoRS
sort F1*.splitmonoRS |uniq -c > F1.allmonoRS; awk '{if ($1 == "3") print $0;}' F1.allmonoRS > F1.conservedmonoRS
sort G1*.splitmonoRS |uniq -c > G1.allmonoRS; awk '{if ($1 == "3") print $0;}' G1.allmonoRS > G1.conservedmonoRS

# sort by scaffold number ($2), position ($3) 
ls *.conservedmonoRS | while read FN; do sort -t$'\t' -k2,2n -k3,3n  $FN > ${FN/conservedmonoRS/sortedconservedmonoRS}; done 
# paste together "scaffold" string-number (convert first tab into no-character)
ls *.sortedconservedmonoRS | while read FN; do sed 's/\t//' $FN > ${FN/sortedconservedmonoRS/consmonoRS}; done 

# extract SNPs within coding (1) regions ($4)
ls *.consmonoRS | while read FN; do awk -F "\t" 'BEGIN{OFS="\t"}{if ($4==1) print$0; }'  $FN > ${FN/consmonoRS/conscodmonoRS}; done
# define repeat threshold 
rep=0
# extract SNPs within coding repeated (>1) regions ($3)  
ls *.conscodmonoRS | while read FN; do awk -v rep="$rep" -F "\t" 'BEGIN{OFS="\t"}{if ($3 > rep) print$0; }'  $FN > ${FN/conscodmonoRS/consrepcodmonoRS}; done
# extract SNPs within coding non-repeated (0) regions ($3) 
ls *.conscodmonoRS | while read FN; do awk -F "\t" 'BEGIN{OFS="\t"}{if ($3==0) print$0; }'  $FN > ${FN/conscodmonoRS/consnrepcodmonoRS}; done

# extract SNPs within non-coding (1) regions ($4)
ls *.consmonoRS | while read FN; do awk -F "\t" 'BEGIN{OFS="\t"}{if ($4==0) print$0; }'  $FN > ${FN/consmonoRS/consncodmonoRS}; done
# extract SNPs within non-coding repeated (>1) regions ($3)  
ls *.consncodmonoRS | while read FN; do awk -v rep="$rep" -F "\t" 'BEGIN{OFS="\t"}{if ($3 > rep) print$0; }'  $FN > ${FN/consncodmonoRS/consrepncodmonoRS}; done
# extract SNPs within non-coding non-repeated (0) regions ($3) 
ls *.consncodmonoRS | while read FN; do awk -F "\t" 'BEGIN{OFS="\t"}{if ($3==0) print$0; }'  $FN > ${FN/consncodmonoRS/consnrepncodmonoRS}; done
echo "...DONE!"

echo "EXTRACTING COUNT FILES..."
# extract count infos
ls *.consrepcodmonoRS| while read FN; do wc -l $FN | awk '{print $1}'  >> repcodmono; done
ls *.consrepncodmonoRS| while read FN; do wc -l $FN  | awk '{print $1}' >> repncodmono; done
ls *.consnrepcodmonoRS| while read FN; do wc -l $FN  | awk '{print $1}' >> nrepcodmono; done
ls *.consnrepncodmonoRS| while read FN; do wc -l $FN  | awk '{print $1}' >> nrepncodmono; done

# extract isolates names
ls *.consrepcodmonoRS > isolate_order
sed 's/\./\t/g' isolate_order > isolate_name
awk '{print $1}' isolate_name > isolates

# write count file: for R analysis
echo -e "isolate\tREPCOD\tREPNCOD\tNREPCOD\tNREPNCOD" > isolateNu6_monoRS
paste isolates repcodmono repncodmono nrepcodmono nrepncodmono >> isolateNu6_monoRS

# reformat for R analysis 
sed -i 's/./&\t/' isolateNu6_monoRS # split isolates to be sorted
sort -k1,1 -k2,2n  isolateNu6_monoRS > isolateNu6_monoRS-sort #sort file by 1st col alpha and 2nd col numerically
sed -i '1h;1d;$!H;$!d;G' isolateNu6_monoRS-sort #put last line top
sed -i 's/\t//' isolateNu6_monoRS-sort # reset isolates names

echo "...DONE!"


#
# step 1.4: analyse polyallelic SNPs
#


echo "ANALYSING POLYALLELIC SNPs..."
# select interesting info from *.catmonoRS: $1_isolate; $2_scaffoldnb; $4_SNPposition; $7_rep(>=1)-nonrep(0); $8_cod(1)-noncod(0); $9_ref-nt; $10_mutation-type; $11_mutated-nt > *.selectedmonoRS
ls *.polyRS | while read FN; do awk -F "\t" 'BEGIN{OFS="\t"}{print $1, $2, $4, $7, $8, $9, $10, $11 }'  $FN > ${FN/polyRS/selectedpolyRS}; done

# cut 1st column to compare replicates per isolate
ls *.selectedpolyRS | while read FN; do cut -f2- $FN > ${FN/selectedpolyRS/cutpolyRS}; done

# split digit from scaffold string
ls *.cutpolyRS | while read FN; do sed 's/Scaffold/&\t/' $FN > ${FN/cutpolyRS/splitpolyRS}; done

# keep only the SNPs which are conserved among  all the replicates within an isolate 
# uniq: to keep the SNPs once, if it's recorded more than once; sort is a requirement for uniq to be executed; an additional column will be printed in *.allmonoRS $1, depending on how many times the position was found
sort A1*.splitpolyRS |uniq -c > A1.allpolyRS; awk '{if ($1 == "3") print $0;}' A1.allpolyRS > A1.conservedpolyRS
sort A2*.splitpolyRS |uniq -c > A2.allpolyRS; awk '{if ($1 == "3") print $0;}' A2.allpolyRS > A2.conservedpolyRS
sort A3*.splitpolyRS |uniq -c > A3.allpolyRS; awk '{if ($1 == "3") print $0;}' A3.allpolyRS > A3.conservedpolyRS
sort A4*.splitpolyRS |uniq -c > A4.allpolyRS; awk '{if ($1 == "3") print $0;}' A4.allpolyRS > A4.conservedpolyRS
sort A5*.splitpolyRS |uniq -c > A5.allpolyRS; awk '{if ($1 == "3") print $0;}' A5.allpolyRS > A5.conservedpolyRS
sort B10*.splitpolyRS |uniq -c > B10.allpolyRS; awk '{if ($1 == "3") print $0;}' B10.allpolyRS > B10.conservedpolyRS
sort B11*.splitpolyRS |uniq -c > B11.allpolyRS; awk '{if ($1 == "3") print $0;}' B11.allpolyRS > B11.conservedpolyRS
sort B12*.splitpolyRS |uniq -c > B12.allpolyRS; awk '{if ($1 == "3") print $0;}' B12.allpolyRS > B12.conservedpolyRS
sort B14*.splitpolyRS |uniq -c > B14.allpolyRS; awk '{if ($1 == "2") print $0;}' B14.allpolyRS > B14.conservedpolyRS
sort B15*.splitpolyRS |uniq -c > B15.allpolyRS; awk '{if ($1 == "3") print $0;}' B15.allpolyRS > B15.conservedpolyRS
sort B17*.splitpolyRS |uniq -c > B17.allpolyRS; awk '{if ($1 == "3") print $0;}' B17.allpolyRS > B17.conservedpolyRS
sort B1-*.splitpolyRS |uniq -c > B1.allpolyRS; awk '{if ($1 == "3") print $0;}' B1.allpolyRS > B1.conservedpolyRS
sort B2*.splitpolyRS |uniq -c > B2.allpolyRS; awk '{if ($1 == "3") print $0;}' B2.allpolyRS > B2.conservedpolyRS
sort B3*.splitpolyRS |uniq -c > B3.allpolyRS; awk '{if ($1 == "3") print $0;}' B3.allpolyRS > B3.conservedpolyRS
sort B4*.splitpolyRS |uniq -c > B4.allpolyRS; awk '{if ($1 == "3") print $0;}' B4.allpolyRS > B4.conservedpolyRS
sort B7*.splitpolyRS |uniq -c > B7.allpolyRS; awk '{if ($1 == "3") print $0;}' B7.allpolyRS > B7.conservedpolyRS
sort B8*.splitpolyRS |uniq -c > B8.allpolyRS; awk '{if ($1 == "3") print $0;}' B8.allpolyRS > B8.conservedpolyRS
sort C1*.splitpolyRS |uniq -c > C1.allpolyRS; awk '{if ($1 == "3") print $0;}' C1.allpolyRS > C1.conservedpolyRS
sort C2*.splitpolyRS |uniq -c > C2.allpolyRS; awk '{if ($1 == "3") print $0;}' C2.allpolyRS > C2.conservedpolyRS
sort C3*.splitpolyRS |uniq -c > C3.allpolyRS; awk '{if ($1 == "3") print $0;}' C3.allpolyRS > C3.conservedpolyRS
sort C4*.splitpolyRS |uniq -c > C4.allpolyRS; awk '{if ($1 == "3") print $0;}' C4.allpolyRS > C4.conservedpolyRS
sort C5*.splitpolyRS |uniq -c > C5.allpolyRS; awk '{if ($1 == "3") print $0;}' C5.allpolyRS > C5.conservedpolyRS
sort CAN*.splitpolyRS |uniq -c > CAN.allpolyRS; awk '{if ($1 == "3") print $0;}' CAN.allpolyRS > CAN.conservedpolyRS
sort D1*.splitpolyRS |uniq -c > D1.allpolyRS; awk '{if ($1 == "3") print $0;}' D1.allpolyRS > D1.conservedpolyRS
sort D2*.splitpolyRS |uniq -c > D2.allpolyRS; awk '{if ($1 == "3") print $0;}' D2.allpolyRS > D2.conservedpolyRS
sort D3*.splitpolyRS |uniq -c > D3.allpolyRS; awk '{if ($1 == "3") print $0;}' D3.allpolyRS > D3.conservedpolyRS
sort D4*.splitpolyRS |uniq -c > D4.allpolyRS; awk '{if ($1 == "3") print $0;}' D4.allpolyRS > D4.conservedpolyRS
sort E1*.splitpolyRS |uniq -c > E1.allpolyRS; awk '{if ($1 == "3") print $0;}' E1.allpolyRS > E1.conservedpolyRS
sort F1*.splitpolyRS |uniq -c > F1.allpolyRS; awk '{if ($1 == "3") print $0;}' F1.allpolyRS > F1.conservedpolyRS
sort G1*.splitpolyRS |uniq -c > G1.allpolyRS; awk '{if ($1 == "3") print $0;}' G1.allpolyRS > G1.conservedpolyRS

# sort by scaffold number ($2), position ($3) 
ls *.conservedpolyRS | while read FN; do sort -t$'\t' -k2,2n -k3,3n  $FN > ${FN/conservedpolyRS/sortedconservedpolyRS}; done 
# paste together "scaffold" string-number (convert first tab into no-character)
ls *.sortedconservedpolyRS | while read FN; do sed 's/\t//' $FN > ${FN/sortedconservedpolyRS/conspolyRS}; done 


# extract SNPs within coding (1) regions ($4)
ls *.conspolyRS | while read FN; do awk -F "\t" 'BEGIN{OFS="\t"}{if ($4==1) print$0; }'  $FN > ${FN/conspolyRS/conscodpolyRS}; done
# define repeat threshold 
rep=0
# extract SNPs within coding repeated (>1) regions ($3) 
ls *.conscodpolyRS | while read FN; do awk  -v rep="$rep" -F "\t" 'BEGIN{OFS="\t"}{if ($3 > rep) print$0; }'  $FN > ${FN/conscodpolyRS/consrepcodpolyRS}; done
# extract SNPs within coding non-repeated (0) regions ($3) 
ls *.conscodpolyRS | while read FN; do awk -F "\t" 'BEGIN{OFS="\t"}{if ($3==0) print$0; }'  $FN > ${FN/conscodpolyRS/consnrepcodpolyRS}; done

# extract SNPs within non-coding (1) regions ($4)
ls *.conspolyRS | while read FN; do awk -F "\t" 'BEGIN{OFS="\t"}{if ($4==0) print$0; }'  $FN > ${FN/conspolyRS/consncodpolyRS}; done
# extract SNPs within non-coding repeated (>1) regions ($3)  
ls *.consncodpolyRS | while read FN; do awk -v rep="$rep" -F "\t" 'BEGIN{OFS="\t"}{if ($3 > rep) print$0; }'  $FN > ${FN/consncodpolyRS/consrepncodpolyRS}; done
# extract SNPs within non-coding non-repeated (0) regions ($3) 
ls *.consncodpolyRS | while read FN; do awk -F "\t" 'BEGIN{OFS="\t"}{if ($3==0) print$0; }'  $FN > ${FN/consncodpolyRS/consnrepncodpolyRS}; done

echo "...DONE!"

echo "EXTRACTING COUNT FILES..."
# extract count infos
ls *.consrepcodpolyRS| while read FN; do wc -l $FN | awk '{print $1}' >> repcodpoly; done
ls *.consrepncodpolyRS| while read FN; do wc -l $FN | awk '{print $1}' >> repncodpoly; done
ls *.consnrepcodpolyRS| while read FN; do wc -l $FN | awk '{print $1}' >> nrepcodpoly; done
ls *.consnrepncodpolyRS| while read FN; do wc -l $FN | awk '{print $1}' >> nrepncodpoly; done

# write count file: for R analysis
# isolates file same as extracted in mono
echo -e "isolate\tREPCOD\tREPNCOD\tNREPCOD\tNREPNCOD" > isolateNu6_polyRS
paste isolates repcodpoly repncodpoly nrepcodpoly nrepncodpoly >> isolateNu6_polyRS

# reformat for R analysis 
sed -i 's/./&\t/' isolateNu6_polyRS # split isolates to be sorted
sort -k1,1 -k2,2n  isolateNu6_polyRS > isolateNu6_polyRS-sort #sort file by 1st col alpha and 2nd col numerically
sed -i '1h;1d;$!H;$!d;G' isolateNu6_polyRS-sort #put last line top
sed -i 's/\t//' isolateNu6_polyRS-sort  # reset isolates names

echo "...DONE!"
echo "Nu6 DATA ANALYSIS COMPLETED!"


###


#
# step 2.0: get Rhiir2 .subVCF files
#      


cd ../
mkdir ./Rhiir2_analysis
cd ./Rhiir2_analysis
# cp /home/dovah/Documents/MSc/MSc-A15/FSP/subVCF/Rhiir2_subVCF-files.zip . 
# unzip Rhiir2_subVCF-files.zip

wget https://copy.com/Pc4b9J3AgpPC7E7S
mv Pc4b9J3AgpPC7E7S Pc4b9J3AgpPC7E7S.zip
unzip Pc4b9J3AgpPC7E7S.zip


#
# step 2.1: prepare subVCF for SNP analysis 
#

echo "ANALYSING Rhiir2 DATA..."
echo "PREPARING subVCF FOR ANALYSIS..."
# define coverage threshold 
cov=10
# filter subVCF files for coverage ($13) value above threshold > *.cov
ls *.subVCF | while read FN; do awk -v cov="$cov" -F "\t" 'BEGIN{OFS="\t"}{if ($13 > cov) print$0; }' $FN > ${FN/subVCF/cov}; done

# define wanted polymorphisms ($10): S=substitutions, R=same as reference
# filter cov files with R and or S for mutation type ($10) > *.RS
ls *.cov | while read FN; do awk -F "\t" 'BEGIN{OFS="\t"}{if ($10 ~ /[RS]/) print$0; }'  $FN > ${FN/cov/RS}; done

# define unwanted polymorphisms ($10): I=insertions; D=deletions; MNP=multi-nucleotide-polymorphism
# clean *.RS files from unwanted polymorpisms > *.cleanRS
ls *.RS | while read FN; do awk -F "\t" 'BEGIN{OFS="\t"}{if ($10 !~ /[IDM]/) print$0; }'  $FN > ${FN/RS/cleanRS}; done


#
# step 2.2: define mono and polyallelic SNPs
#


echo "DEFINING MONO/POLYALLELIC SNPs..."
# define polyallelic SNPs: will have more than 1 mutation event reported in $10, field separator within $10: comma (',')
# look for polyallelic SNPs in *.clanRS > *.polyRS 
ls *.cleanRS | while read FN; do awk -F "\t" 'BEGIN{OFS="\t"}{if ($10 ~ /,/) print$0; }'  $FN > ${FN/cleanRS/polyRS}; done

# define monoallelic SNPs: will have exactly 1 mutation event reported in $10, no field separator (',') needed within $10, or no mutation events reported in $10
# look for monoallelic SNPs in *.cleanRS > *.monoRS 
ls *.cleanRS | while read FN; do awk -F "\t" 'BEGIN{OFS="\t"}{if ($10 !~ /,/) print$0; }'  $FN > ${FN/cleanRS/monoRS}; done
echo "...DONE!"


#
# step 2.3: analyse monoallelic SNPs
#


echo "ANALYSING MONOALLELIC SNPs...."
# select interesting info from *.catmonoRS: $1_isolate; $2_scaffoldnb; $4_SNPposition; $7_rep(>=1)-nonrep(0); $8_cod(1)-noncod(0); $9_ref-nt; $10_mutation-type; $11_mutated-nt > *.selectedmonoRS
ls *.monoRS | while read FN; do awk -F "\t" 'BEGIN{OFS="\t"}{print $1, $2, $4, $7, $8, $9, $10, $11 }'  $FN > ${FN/monoRS/selectedmonoRS}; done

# cut 1st column to compare replicates per isolate
ls *.selectedmonoRS | while read FN; do cut -f2- $FN > ${FN/selectedmonoRS/cutmonoRS}; done

# split digit from scaffold string
ls *.cutmonoRS | while read FN; do sed 's/scaffold/&\t/' $FN > ${FN/cutmonoRS/spltmonoRS}; done
ls *.spltmonoRS | while read FN; do sed 's/_//g' $FN > ${FN/spltmonoRS/splitmonoRS}; done

# keep only the SNPs which are conserved among  all the replicates within an isolate 
# uniq: to keep the SNPs once, if it's recorded more than once; sort is a requirement for uniq to be executed; an additional column will be printed in *.allmonoRS $1, depending on how many times the position was found
sort A1*.splitmonoRS |uniq -c > A1.allmonoRS; awk '{if ($1 == "3") print $0;}' A1.allmonoRS > A1.conservedmonoRS
sort A2*.splitmonoRS |uniq -c > A2.allmonoRS; awk '{if ($1 == "3") print $0;}' A2.allmonoRS > A2.conservedmonoRS
sort A3*.splitmonoRS |uniq -c > A3.allmonoRS; awk '{if ($1 == "3") print $0;}' A3.allmonoRS > A3.conservedmonoRS
sort A4*.splitmonoRS |uniq -c > A4.allmonoRS; awk '{if ($1 == "3") print $0;}' A4.allmonoRS > A4.conservedmonoRS
sort A5*.splitmonoRS |uniq -c > A5.allmonoRS; awk '{if ($1 == "3") print $0;}' A5.allmonoRS > A5.conservedmonoRS
sort B10*.splitmonoRS |uniq -c > B10.allmonoRS; awk '{if ($1 == "3") print $0;}' B10.allmonoRS > B10.conservedmonoRS
sort B11*.splitmonoRS |uniq -c > B11.allmonoRS; awk '{if ($1 == "3") print $0;}' B11.allmonoRS > B11.conservedmonoRS
sort B12*.splitmonoRS |uniq -c > B12.allmonoRS; awk '{if ($1 == "3") print $0;}' B12.allmonoRS > B12.conservedmonoRS
sort B14*.splitmonoRS |uniq -c > B14.allmonoRS; awk '{if ($1 == "3") print $0;}' B14.allmonoRS > B14.conservedmonoRS 
sort B15*.splitmonoRS |uniq -c > B15.allmonoRS; awk '{if ($1 == "3") print $0;}' B15.allmonoRS > B15.conservedmonoRS
sort B17*.splitmonoRS |uniq -c > B17.allmonoRS; awk '{if ($1 == "3") print $0;}' B17.allmonoRS > B17.conservedmonoRS
sort B1-*.splitmonoRS |uniq -c > B1.allmonoRS; awk '{if ($1 == "3") print $0;}' B1.allmonoRS > B1.conservedmonoRS
sort B2*.splitmonoRS |uniq -c > B2.allmonoRS; awk '{if ($1 == "3") print $0;}' B2.allmonoRS > B2.conservedmonoRS
sort B3*.splitmonoRS |uniq -c > B3.allmonoRS; awk '{if ($1 == "3") print $0;}' B3.allmonoRS > B3.conservedmonoRS
sort B4*.splitmonoRS |uniq -c > B4.allmonoRS; awk '{if ($1 == "3") print $0;}' B4.allmonoRS > B4.conservedmonoRS
sort B7*.splitmonoRS |uniq -c > B7.allmonoRS; awk '{if ($1 == "3") print $0;}' B7.allmonoRS > B7.conservedmonoRS
sort B8*.splitmonoRS |uniq -c > B8.allmonoRS; awk '{if ($1 == "3") print $0;}' B8.allmonoRS > B8.conservedmonoRS
sort C1*.splitmonoRS |uniq -c > C1.allmonoRS; awk '{if ($1 == "3") print $0;}' C1.allmonoRS > C1.conservedmonoRS
sort C2*.splitmonoRS |uniq -c > C2.allmonoRS; awk '{if ($1 == "3") print $0;}' C2.allmonoRS > C2.conservedmonoRS
sort C3*.splitmonoRS |uniq -c > C3.allmonoRS; awk '{if ($1 == "3") print $0;}' C3.allmonoRS > C3.conservedmonoRS
sort C4*.splitmonoRS |uniq -c > C4.allmonoRS; awk '{if ($1 == "3") print $0;}' C4.allmonoRS > C4.conservedmonoRS
sort C5*.splitmonoRS |uniq -c > C5.allmonoRS; awk '{if ($1 == "3") print $0;}' C5.allmonoRS > C5.conservedmonoRS
sort CAN*.splitmonoRS |uniq -c > CAN.allmonoRS; awk '{if ($1 == "3") print $0;}' CAN.allmonoRS > CAN.conservedmonoRS
sort D1*.splitmonoRS |uniq -c > D1.allmonoRS; awk '{if ($1 == "3") print $0;}' D1.allmonoRS > D1.conservedmonoRS
sort D2*.splitmonoRS |uniq -c > D2.allmonoRS; awk '{if ($1 == "3") print $0;}' D2.allmonoRS > D2.conservedmonoRS
sort D3*.splitmonoRS |uniq -c > D3.allmonoRS; awk '{if ($1 == "3") print $0;}' D3.allmonoRS > D3.conservedmonoRS
sort D4*.splitmonoRS |uniq -c > D4.allmonoRS; awk '{if ($1 == "3") print $0;}' D4.allmonoRS > D4.conservedmonoRS
sort E1*.splitmonoRS |uniq -c > E1.allmonoRS; awk '{if ($1 == "3") print $0;}' E1.allmonoRS > E1.conservedmonoRS
sort F1*.splitmonoRS |uniq -c > F1.allmonoRS; awk '{if ($1 == "3") print $0;}' F1.allmonoRS > F1.conservedmonoRS
sort G1*.splitmonoRS |uniq -c > G1.allmonoRS; awk '{if ($1 == "3") print $0;}' G1.allmonoRS > G1.conservedmonoRS

# sort by scaffold number ($2), position ($3) 
ls *.conservedmonoRS | while read FN; do sort -t$'\t' -k2,2n -k3,3n  $FN > ${FN/conservedmonoRS/sortedconservedmonoRS}; done 
# paste together "scaffold" string-number (convert first tab into no-character)
ls *.sortedconservedmonoRS | while read FN; do sed 's/\t//' $FN > ${FN/sortedconservedmonoRS/consmonoRS}; done 

# extract SNPs within coding (1) regions ($4)
ls *.consmonoRS | while read FN; do awk -F "\t" 'BEGIN{OFS="\t"}{if ($4==1) print$0; }'  $FN > ${FN/consmonoRS/conscodmonoRS}; done
# define repeat threshold 
rep=0
# extract SNPs within coding repeated (>1) regions ($3)  
ls *.conscodmonoRS | while read FN; do awk -v rep="$rep" -F "\t" 'BEGIN{OFS="\t"}{if ($3 > rep) print$0; }'  $FN > ${FN/conscodmonoRS/consrepcodmonoRS}; done
# extract SNPs within coding non-repeated (0) regions ($3) 
ls *.conscodmonoRS | while read FN; do awk -F "\t" 'BEGIN{OFS="\t"}{if ($3==0) print$0; }'  $FN > ${FN/conscodmonoRS/consnrepcodmonoRS}; done

# extract SNPs within non-coding (1) regions ($4)
ls *.consmonoRS | while read FN; do awk -F "\t" 'BEGIN{OFS="\t"}{if ($4==0) print$0; }'  $FN > ${FN/consmonoRS/consncodmonoRS}; done
# extract SNPs within non-coding repeated (>1) regions ($3)  
ls *.consncodmonoRS | while read FN; do awk -v rep="$rep" -F "\t" 'BEGIN{OFS="\t"}{if ($3 > rep) print$0; }'  $FN > ${FN/consncodmonoRS/consrepncodmonoRS}; done
# extract SNPs within non-coding non-repeated (0) regions ($3) 
ls *.consncodmonoRS | while read FN; do awk -F "\t" 'BEGIN{OFS="\t"}{if ($3==0) print$0; }'  $FN > ${FN/consncodmonoRS/consnrepncodmonoRS}; done
echo "...DONE!"

echo "EXTRACTING COUNT FILES..."
# extract count infos
ls *.consrepcodmonoRS| while read FN; do wc -l $FN | awk '{print $1}'  >> repcodmono; done
ls *.consrepncodmonoRS| while read FN; do wc -l $FN  | awk '{print $1}' >> repncodmono; done
ls *.consnrepcodmonoRS| while read FN; do wc -l $FN  | awk '{print $1}' >> nrepcodmono; done
ls *.consnrepncodmonoRS| while read FN; do wc -l $FN  | awk '{print $1}' >> nrepncodmono; done

# extract isolates names
ls *.consrepcodmonoRS > isolate_order
sed 's/\./\t/g' isolate_order > isolate_name
awk '{print $1}' isolate_name > isolates

# writecount file: for R analysis
echo -e "isolate\tREPCOD\tREPNCOD\tNREPCOD\tNREPNCOD" > isolateRhiir2_monoRS
paste isolates repcodmono repncodmono nrepcodmono nrepncodmono >> isolateRhiir2_monoRS

# reformat for R analysis 
sed -i 's/./&\t/' isolateRhiir2_monoRS #split alpha-num isolate name 
sort -k1,1 -k2,2n  isolateRhiir2_monoRS > isolateRhiir2_monoRS-sort #sort file by 1st col alpha and 2nd col numerically
sed -i '1h;1d;$!H;$!d;G' isolateRhiir2_monoRS-sort #put last line top
sed -i 's/\t//' isolateRhiir2_monoRS-sort #restore isolate name 
echo "...DONE!"


#
# step 2.4: analyse polyallelic SNPs
#


echo "ANALYSING POLYALLELIC SNPs..."
# select interesting info from *.catmonoRS: $1_isolate; $2_scaffoldnb; $4_SNPposition; $7_rep(>=1)-nonrep(0); $8_cod(1)-noncod(0); $9_ref-nt; $10_mutation-type; $11_mutated-nt > *.selectedmonoRS
ls *.polyRS | while read FN; do awk -F "\t" 'BEGIN{OFS="\t"}{print $1, $2, $4, $7, $8, $9, $10, $11 }'  $FN > ${FN/polyRS/selectedpolyRS}; done

# cut 1st column to compare replicates per isolate
ls *.selectedpolyRS | while read FN; do cut -f2- $FN > ${FN/selectedpolyRS/cutpolyRS}; done

# split digit from scaffold string
ls *.cutpolyRS | while read FN; do sed 's/scaffold/&\t/' $FN > ${FN/cutpolyRS/spltpolyRS}; done
ls *.spltpolyRS | while read FN; do sed 's/_//g' $FN > ${FN/spltpolyRS/splitpolyRS}; done

# keep only the SNPs which are conserved among  all the replicates within an isolate 
# uniq: to keep the SNPs once, if it's recorded more than once; sort is a requirement for uniq to be executed; an additional column will be printed in *.allmonoRS $1, depending on how many times the position was found
sort A1*.splitpolyRS |uniq -c > A1.allpolyRS; awk '{if ($1 == "3") print $0;}' A1.allpolyRS > A1.conservedpolyRS
sort A2*.splitpolyRS |uniq -c > A2.allpolyRS; awk '{if ($1 == "3") print $0;}' A2.allpolyRS > A2.conservedpolyRS
sort A3*.splitpolyRS |uniq -c > A3.allpolyRS; awk '{if ($1 == "3") print $0;}' A3.allpolyRS > A3.conservedpolyRS
sort A4*.splitpolyRS |uniq -c > A4.allpolyRS; awk '{if ($1 == "3") print $0;}' A4.allpolyRS > A4.conservedpolyRS
sort A5*.splitpolyRS |uniq -c > A5.allpolyRS; awk '{if ($1 == "3") print $0;}' A5.allpolyRS > A5.conservedpolyRS
sort B10*.splitpolyRS |uniq -c > B10.allpolyRS; awk '{if ($1 == "3") print $0;}' B10.allpolyRS > B10.conservedpolyRS
sort B11*.splitpolyRS |uniq -c > B11.allpolyRS; awk '{if ($1 == "3") print $0;}' B11.allpolyRS > B11.conservedpolyRS
sort B12*.splitpolyRS |uniq -c > B12.allpolyRS; awk '{if ($1 == "3") print $0;}' B12.allpolyRS > B12.conservedpolyRS
sort B14*.splitpolyRS |uniq -c > B14.allpolyRS; awk '{if ($1 == "3") print $0;}' B14.allpolyRS > B14.conservedpolyRS 
sort B15*.splitpolyRS |uniq -c > B15.allpolyRS; awk '{if ($1 == "3") print $0;}' B15.allpolyRS > B15.conservedpolyRS
sort B17*.splitpolyRS |uniq -c > B17.allpolyRS; awk '{if ($1 == "3") print $0;}' B17.allpolyRS > B17.conservedpolyRS
sort B1-*.splitpolyRS |uniq -c > B1.allpolyRS; awk '{if ($1 == "3") print $0;}' B1.allpolyRS > B1.conservedpolyRS
sort B2*.splitpolyRS |uniq -c > B2.allpolyRS; awk '{if ($1 == "3") print $0;}' B2.allpolyRS > B2.conservedpolyRS
sort B3*.splitpolyRS |uniq -c > B3.allpolyRS; awk '{if ($1 == "3") print $0;}' B3.allpolyRS > B3.conservedpolyRS
sort B4*.splitpolyRS |uniq -c > B4.allpolyRS; awk '{if ($1 == "3") print $0;}' B4.allpolyRS > B4.conservedpolyRS
sort B7*.splitpolyRS |uniq -c > B7.allpolyRS; awk '{if ($1 == "3") print $0;}' B7.allpolyRS > B7.conservedpolyRS
sort B8*.splitpolyRS |uniq -c > B8.allpolyRS; awk '{if ($1 == "3") print $0;}' B8.allpolyRS > B8.conservedpolyRS
sort C1*.splitpolyRS |uniq -c > C1.allpolyRS; awk '{if ($1 == "3") print $0;}' C1.allpolyRS > C1.conservedpolyRS
sort C2*.splitpolyRS |uniq -c > C2.allpolyRS; awk '{if ($1 == "3") print $0;}' C2.allpolyRS > C2.conservedpolyRS
sort C3*.splitpolyRS |uniq -c > C3.allpolyRS; awk '{if ($1 == "3") print $0;}' C3.allpolyRS > C3.conservedpolyRS
sort C4*.splitpolyRS |uniq -c > C4.allpolyRS; awk '{if ($1 == "3") print $0;}' C4.allpolyRS > C4.conservedpolyRS
sort C5*.splitpolyRS |uniq -c > C5.allpolyRS; awk '{if ($1 == "3") print $0;}' C5.allpolyRS > C5.conservedpolyRS
sort CAN*.splitpolyRS |uniq -c > CAN.allpolyRS; awk '{if ($1 == "3") print $0;}' CAN.allpolyRS > CAN.conservedpolyRS
sort D1*.splitpolyRS |uniq -c > D1.allpolyRS; awk '{if ($1 == "3") print $0;}' D1.allpolyRS > D1.conservedpolyRS
sort D2*.splitpolyRS |uniq -c > D2.allpolyRS; awk '{if ($1 == "3") print $0;}' D2.allpolyRS > D2.conservedpolyRS
sort D3*.splitpolyRS |uniq -c > D3.allpolyRS; awk '{if ($1 == "3") print $0;}' D3.allpolyRS > D3.conservedpolyRS
sort D4*.splitpolyRS |uniq -c > D4.allpolyRS; awk '{if ($1 == "3") print $0;}' D4.allpolyRS > D4.conservedpolyRS
sort E1*.splitpolyRS |uniq -c > E1.allpolyRS; awk '{if ($1 == "3") print $0;}' E1.allpolyRS > E1.conservedpolyRS
sort F1*.splitpolyRS |uniq -c > F1.allpolyRS; awk '{if ($1 == "3") print $0;}' F1.allpolyRS > F1.conservedpolyRS
sort G1*.splitpolyRS |uniq -c > G1.allpolyRS; awk '{if ($1 == "3") print $0;}' G1.allpolyRS > G1.conservedpolyRS

# sort by scaffold number ($2), position ($3) 
ls *.conservedpolyRS | while read FN; do sort -t$'\t' -k2,2n -k3,3n  $FN > ${FN/conservedpolyRS/sortedconservedpolyRS}; done 
# paste together "scaffold" string-number (convert first tab into no-character)
ls *.sortedconservedpolyRS | while read FN; do sed 's/\t//' $FN > ${FN/sortedconservedpolyRS/conspolyRS}; done 

# extract SNPs within coding (1) regions ($4)
ls *.conspolyRS | while read FN; do awk -F "\t" 'BEGIN{OFS="\t"}{if ($4==1) print$0; }'  $FN > ${FN/conspolyRS/conscodpolyRS}; done
# define repeat threshold 
rep=0
# extract SNPs within coding repeated (>1) regions ($3) 
ls *.conscodpolyRS | while read FN; do awk  -v rep="$rep" -F "\t" 'BEGIN{OFS="\t"}{if ($3 > rep) print$0; }'  $FN > ${FN/conscodpolyRS/consrepcodpolyRS}; done
# extract SNPs within coding non-repeated (0) regions ($3) 
ls *.conscodpolyRS | while read FN; do awk -F "\t" 'BEGIN{OFS="\t"}{if ($3==0) print$0; }'  $FN > ${FN/conscodpolyRS/consnrepcodpolyRS}; done

# extract SNPs within non-coding (1) regions ($4)
ls *.conspolyRS | while read FN; do awk -F "\t" 'BEGIN{OFS="\t"}{if ($4==0) print$0; }'  $FN > ${FN/conspolyRS/consncodpolyRS}; done
# extract SNPs within non-coding repeated (>1) regions ($3)  
ls *.consncodpolyRS | while read FN; do awk -v rep="$rep" -F "\t" 'BEGIN{OFS="\t"}{if ($3 > rep) print$0; }'  $FN > ${FN/consncodpolyRS/consrepncodpolyRS}; done
# extract SNPs within non-coding non-repeated (0) regions ($3) 
ls *.consncodpolyRS | while read FN; do awk -F "\t" 'BEGIN{OFS="\t"}{if ($3==0) print$0; }'  $FN > ${FN/consncodpolyRS/consnrepncodpolyRS}; done

echo "...DONE!"
echo "EXTRACTING COUNT FILES..."

# extract count infos
ls *.consrepcodpolyRS| while read FN; do wc -l $FN | awk '{print $1}' >> repcodpoly; done
ls *.consrepncodpolyRS| while read FN; do wc -l $FN | awk '{print $1}' >> repncodpoly; done
ls *.consnrepcodpolyRS| while read FN; do wc -l $FN | awk '{print $1}' >> nrepcodpoly; done
ls *.consnrepncodpolyRS| while read FN; do wc -l $FN | awk '{print $1}' >> nrepncodpoly; done

# write count file: for R analysis
# isolates file same as extracted in mono
echo -e "isolate\tREPCOD\tREPNCOD\tNREPCOD\tNREPNCOD" > isolateRhiir2_polyRS
paste isolates repcodpoly repncodpoly nrepcodpoly nrepncodpoly >> isolateRhiir2_polyRS

#reformat file for R analysis 
sed -i 's/./&\t/' isolateRhiir2_polyRS #split alpha-num isolate name 
sort -k1,1 -k2,2n  isolateRhiir2_polyRS > isolateRhiir2_polyRS-sort #sort file by 1st col alpha and 2nd col numerically
sed -i '1h;1d;$!H;$!d;G' isolateRhiir2_polyRS-sort #put last line top
sed -i 's/\t//' isolateRhiir2_polyRS-sort #restore original isolate name 

echo "...DONE!"
echo "Rhiir2 DATA ANALYSIS COMPLETED!"
