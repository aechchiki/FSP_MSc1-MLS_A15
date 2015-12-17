#!/bin/bash


#####
#
#  Author: Amina Echchiki [mailto:amina.echchiki@unil.ch]
#  Affiliation: MSc student, MLS bioinformatics, University of Lausanne
#  Date: 2015-12-10
#  Name: 00__bam_to_subVCF.sh
#  Purpose: Retrieve *.subVCF files from *.bam inputs using N6 and Rhiir2 as reference genome
#  Comments: requires *.bam files of input isolates
#  Testing platforms: server Linux 2.6.32-504.el6.x86_64
#  Shell version: GNU bash, version 4.3.11(1)-release (x86_64-pc-linux-gnu)
#
#####


#
# N6 genome
#


# from *.bam files, generate *.vcf using FreeBayes
#
a=0; for i in $(ls *.bam); do echo $i;a=$((a + 1)); bsub -q dee-hugemem -L /bin/bash -J freebayes$a -u amina.echchiki@unil.ch -N " export PATH=/software/bin:$PATH; module add UHTS/Analysis/freebayes/0.9.9.2; freebayes -f /scratch/cluster/monthly/twyss/Amina/Genome_files/Glomus_Nu6_genome.fa -p 10 -J -K -F 0.1 -0 -u -v $(echo $i | cut -d'.' -f1)'.vcf' -b $i"; done


# from *.bam files, generate *.coverage using SamTools
#
b=0; for i in $(ls *.bam); do echo $i; a=$(echo $i | cut -d'.' -f1 | cut -d'_' -f3); echo $a; b=$((b + 1)); bsub -q dee-hugemem -L /bin/bash -J samtools$b -u amina.echchiki@unil.ch -N " export PATH=/software/bin:$PATH; module use /software/module/;module add UHTS/Analysis/samtools/0.1.19; samtools depth $i > $a'.coverage' "; done

# from *.vcf files, generate high-quality *Q30.vcf using vcffilter
#
for i in $(ls *RG.vcf); do echo $i; bsub -q dee-hugemem -L /bin/bash -J samtools$b -u amina.echchiki@unil.ch -N " /scratch/cluster/monthly/twyss/Amina/vcffilter -f 'QUAL > 30' $i > $(echo $i | cut -d'.' -f1)'Q30.vcf'" ; done

# merge variants in *combinedvcf using FM perl script n.1
#
bsub -q dee-hugemem -L /bin/bash -J samtools$b -u amina.echchiki@unil.ch -M 6000000 -N "perl /scratch/cluster/monthly/twyss/Amina/KIT_to_make_subVCF/PP25_00_1-9_Haplotyper_VCF_aggregator_record_all_positions.pl /scratch/cluster/monthly/twyss/Amina/KIT_to_make_subVCF/Predicted-RADseqFragments_V4-2015-Nu6__CluRepCod_TR3.txt"


# merge coverage and variants in *.subVCF using FM perl script n.2
#
for i in $(ls *.coverage); do echo $i; a=$(echo $i | cut -d'.' -f1); b=$(ls *_$(printf %q "${a}")_*.vcf); echo $b; bsub -q dee-hugemem -L /bin/bash -J subVCF_$a -u amina.echchiki@unil.ch -N -R "rusage[mem=6000]" -M 6000000 "perl /scratch/cluster/monthly/twyss/Amina/KIT_to_make_subVCF/PP25_00_2-9_Generate_lists_of_variants__subVCF_maker.pl All_positions_with_variants20151015PRV3.combinedVCF $i $b "; done


#
# Rhiir2 genome
#


# from *.bam files, generate *.vcf using FreeBayes
#
a=0; for i in $(ls *.bam); do echo $i;a=$((a + 1)); bsub -q dee-hugemem -L /bin/bash -J freebayes$a -u amina.echchiki@unil.ch -N " export PATH=/software/bin:$PATH; module add UHTS/Analysis/freebayes/0.9.9.2; freebayes -f /scratch/cluster/monthly/twyss/Amina/Genome_files/Rhiir2_DAOM197198_AssemblyScaffolds.fasta -p 10 -J -K -F 0.1 -0 -u -v $(echo $i | cut -d'.' -f1)'.vcf' -b $i"; done

# from *.bam files, generate *.coverage using SamTools
#
b=0; for i in $(ls *.bam); do echo $i; a=$(echo $i | cut -d'.' -f1 | cut -d'_' -f3); echo $a; b=$((b + 1)); bsub -q dee-hugemem -L /bin/bash -J samtools$b -u amina.echchiki@unil.ch -N " export PATH=/software/bin:$PATH; module use /software/module/;module add UHTS/Analysis/samtools/0.1.19; samtools depth $i > $a'.coverage' "; done

# from *.vcf files, generate high-quality *Q30.vcf using vcffilter
#
for i in $(ls *RG.vcf); do echo $i; bsub -q dee-hugemem -L /bin/bash -J samtools$b -u amina.echchiki@unil.ch -N " /scratch/cluster/monthly/twyss/Amina/vcffilter -f 'QUAL > 30' $i > $(echo $i | cut -d'.' -f1)'Q30.vcf'" ; done

# merge variants in *combinedvcf using FM perl script n.1
#
bsub -q dee-hugemem -L /bin/bash -J samtools$b -u amina.echchiki@unil.ch -M 6000000 -N "perl /scratch/cluster/monthly/twyss/Amina/KIT_to_make_subVCF/PP25_00_1-9_Haplotyper_VCF_aggregator_record_all_positions.pl Predicted-RADseqFragments_V3__CluRepCod_TR3_Rhii2.txt"

# merge coverage and variants in *.subVCF using FM perl script n.2
#
for i in $(ls *.coverage); do echo $i; a=$(echo $i | cut -d'.' -f1); b=$(ls *_$(printf %q "${a}")_*.vcf); echo $b; bsub -q dee-hugemem -L /bin/bash -J subVCF_$a -u amina.echchiki@unil.ch -N -R "rusage[mem=6000]" -M 6000000 "perl /scratch/cluster/monthly/twyss/Amina/KIT_to_make_subVCF/PP25_00_2-9_Generate_lists_of_variants__subVCF_maker.pl All_positions_with_variants*.combinedVCF $i $b "; done