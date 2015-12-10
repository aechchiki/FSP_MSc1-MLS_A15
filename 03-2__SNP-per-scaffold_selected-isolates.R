setwd("/home/dovah/Documents/MSc/MSc-A15/FSP/09-12_NEW-SCRIPTS/20-scaffold_SNP-density/")
x=c(1:20)

## isolate d1
# genome N6
d1mono_nu6=read.table("d1Nu6mono", h=T); d1poly_nu6=read.table("d1Nu6poly", h=T)
d1tot_nu6=d1mono_nu6[2:5]+d1poly_nu6[2:5]; d1tot_genome_nu6=rowSums(d1tot_nu6)
barplot(d1tot_genome_nu6, xlab="Scaffold number", ylab="SNP number", names.arg = x, main="SNP per scaffold: isolate D1 (N6)")
# genome rhiir2
d1mono_rh2=read.table("d1Rhiir2mono", h=T); d1poly_rh2=read.table("d1Rhiir2poly", h=T)
d1tot_rh2=d1mono_rh2[2:5]+d1poly_rh2[2:5]; d1tot_genome_rh2=rowSums(d1tot_rh2)
barplot(d1tot_genome_rh2, xlab="Scaffold number", ylab="SNP number", names.arg = x, main="SNP per scaffold: isolate D1 (Rhiir2)")

## isolate daom
# genome N6
canmono_nu6=read.table("canNu6mono", h=T); canpoly_nu6=read.table("canNu6poly", h=T)
cantot_nu6=canmono_nu6[2:5]+canpoly_nu6[2:5]; cantot_genome_nu6=rowSums(cantot_nu6)
barplot(cantot_genome_nu6, xlab="Scaffold number", ylab="SNP number", names.arg = x, main="SNP per scaffold: isolate DAOM (N6)")
# genome rhiir2
canmono_rh2=read.table("canRhiir2mono", h=T); canpoly_rh2=read.table("canRhiir2poly", h=T)
cantot_rh2=canmono_rh2[2:5]+canpoly_rh2[2:5]; cantot_genome_rh2=rowSums(cantot_rh2)
barplot(cantot_genome_rh2, xlab="Scaffold number", ylab="SNP number", names.arg = x, main="SNP per scaffold: isolate DAOM (Rhiir2)")

## isolate C3
# genome N6
c3mono_nu6=read.table("c3Nu6mono", h=T); c3poly_nu6=read.table("c3Nu6poly", h=T)
c3tot_nu6=c3mono_nu6[2:5]+c3poly_nu6[2:5]; c3tot_genome_nu6=rowSums(c3tot_nu6)
barplot(c3tot_genome_nu6, xlab="Scaffold number", ylab="SNP number", names.arg = x, main="SNP per scaffold: isolate C3 (N6)")
# genome rhiir2
c3mono_rh2=read.table("c3Rhiir2mono", h=T); c3poly_rh2=read.table("c3Rhiir2poly", h=T)
c3tot_rh2=c3mono_rh2[2:5]+c3poly_rh2[2:5]; c3tot_genome_rh2=rowSums(c3tot_rh2)
barplot(c3tot_genome_rh2, xlab="Scaffold number", ylab="SNP number", names.arg = x, main="SNP per scaffold: isolate C3 (Rhiir2)")
## isolate c2
# genome N6
c2mono_nu6=read.table("c2Nu6mono", h=T); c2poly_nu6=read.table("c2Nu6poly", h=T)
c2tot_nu6=c2mono_nu6[2:5]+c2poly_nu6[2:5]; c2tot_genome_nu6=rowSums(c2tot_nu6)
barplot(c2tot_genome_nu6, xlab="Scaffold number", ylab="SNP number", names.arg = x, main="SNP per scaffold: isolate C2 (N6)")
# genome rhiir2
c2mono_rh2=read.table("c2Rhiir2mono", h=T); c2poly_rh2=read.table("c2Rhiir2poly", h=T)
c2tot_rh2=c2mono_rh2[2:5]+c2poly_rh2[2:5]; c2tot_genome_rh2=rowSums(c2tot_rh2)
barplot(c2tot_genome_rh2, xlab="Scaffold number", ylab="SNP number", names.arg = x, main="SNP per scaffold: isolate C2 (Rhiir2)")

## isolate b14
# genome N6
b14mono_nu6=read.table("b14Nu6mono", h=T); b14poly_nu6=read.table("b14Nu6poly", h=T)
b14tot_nu6=b14mono_nu6[2:5]+b14poly_nu6[2:5]; b14tot_genome_nu6=rowSums(b14tot_nu6)
barplot(b14tot_genome_nu6, xlab="Scaffold number", ylab="SNP number", names.arg = x, main="SNP per scaffold: isolate B14 (N6)")
# genome rhiir2
b14mono_rh2=read.table("b14Rhiir2mono", h=T); b14poly_rh2=read.table("b14Rhiir2poly", h=T)
b14tot_rh2=b14mono_rh2[2:5]+b14poly_rh2[2:5]; b14tot_genome_rh2=rowSums(b14tot_rh2)
barplot(b14tot_genome_rh2, xlab="Scaffold number", ylab="SNP number", names.arg = x, main="SNP per scaffold: isolate B14 (Rhiir2)")

## isolate b11
# genome N6
b11mono_nu6=read.table("b11Nu6mono", h=T); b11poly_nu6=read.table("b11Nu6poly", h=T)
b11tot_nu6=b11mono_nu6[2:5]+b11poly_nu6[2:5]; b11tot_genome_nu6=rowSums(b11tot_nu6)
barplot(b11tot_genome_nu6, xlab="Scaffold number", ylab="SNP number", names.arg = x, main="SNP per scaffold: isolate B11 (N6)")
# genome rhiir2
b11mono_rh2=read.table("b11Rhiir2mono", h=T); b11poly_rh2=read.table("b11Rhiir2poly", h=T)
b11tot_rh2=b11mono_rh2[2:5]+b11poly_rh2[2:5]; b11tot_genome_rh2=rowSums(b11tot_rh2)
barplot(b11tot_genome_rh2, xlab="Scaffold number", ylab="SNP number", names.arg = x, main="SNP per scaffold: isolate B11 (Rhiir2)")

## isolate b3
# genome N6
b3mono_nu6=read.table("b3Nu6mono", h=T); b3poly_nu6=read.table("b3Nu6poly", h=T)
b3tot_nu6=b3mono_nu6[2:5]+b3poly_nu6[2:5]; b3tot_genome_nu6=rowSums(b3tot_nu6)
barplot(b3tot_genome_nu6, xlab="Scaffold number", ylab="SNP number", names.arg = x, main="SNP per scaffold: isolate B3 (N6)")
# genome rhiir2
b3mono_rh2=read.table("b3Rhiir2mono", h=T); b3poly_rh2=read.table("b3Rhiir2poly", h=T)
b3tot_rh2=b3mono_rh2[2:5]+b3poly_rh2[2:5]; b3tot_genome_rh2=rowSums(b3tot_rh2)
barplot(b3tot_genome_rh2, xlab="Scaffold number", ylab="SNP number", names.arg = x, main="SNP per scaffold: isolate B3 (Rhiir2)")

## isolate a2
# genome N6
a2mono_nu6=read.table("a2Nu6mono", h=T); a2poly_nu6=read.table("a2Nu6poly", h=T)
a2tot_nu6=a2mono_nu6[2:5]+a2poly_nu6[2:5]; a2tot_genome_nu6=rowSums(a2tot_nu6)
barplot(a2tot_genome_nu6, xlab="Scaffold number", ylab="SNP number", names.arg = x, main="SNP per scaffold: isolate A2 (N6)")
# genome rhiir2
a2mono_rh2=read.table("a2Rhiir2mono", h=T); a2poly_rh2=read.table("a2Rhiir2poly", h=T)
a2tot_rh2=a2mono_rh2[2:5]+a2poly_rh2[2:5]; a2tot_genome_rh2=rowSums(a2tot_rh2)
barplot(a2tot_genome_rh2, xlab="Scaffold number", ylab="SNP number", names.arg = x, main="SNP per scaffold: isolate A2 (Rhiir2)")

## isolate a1
# genome N6
a1mono_nu6=read.table("a1Nu6mono", h=T); a1poly_nu6=read.table("a1Nu6poly", h=T)
a1tot_nu6=a1mono_nu6[2:5]+a1poly_nu6[2:5]; a1tot_genome_nu6=rowSums(a1tot_nu6)
barplot(a1tot_genome_nu6, xlab="Scaffold number", ylab="SNP number", names.arg = x, main="SNP per scaffold: isolate A1 (N6)")
# genome rhiir2
a1mono_rh2=read.table("a1Rhiir2mono", h=T); a1poly_rh2=read.table("a1Rhiir2poly", h=T)
a1tot_rh2=a1mono_rh2[2:5]+a1poly_rh2[2:5]; a1tot_genome_rh2=rowSums(a1tot_rh2)
barplot(a1tot_genome_rh2, xlab="Scaffold number", ylab="SNP number", names.arg = x, main="SNP per scaffold: isolate A1 (Rhiir2)")


###
# genome N6
averagemono_nu6=read.table("Nu6mono", h=T); averagepoly_nu6=read.table("Nu6poly", h=T)
averagetot_nu6=averagemono_nu6[2:5]+averagepoly_nu6[2:5]; averagetot_genome_nu6=rowSums(averagetot_nu6)
barplot(averagetot_genome_nu6, xlab="Scaffold number", ylab="SNP number", names.arg = x, main="SNP per scaffold: average (N6)")
# genome rhiir2
averagemono_rh2=read.table("Rhiir2mono", h=T); averagepoly_rh2=read.table("Rhiir2poly", h=T)
averagetot_rh2=averagemono_rh2[2:5]+averagepoly_rh2[2:5]; averagetot_genome_rh2=rowSums(averagetot_rh2)
barplot(averagetot_genome_rh2, xlab="Scaffold number", ylab="SNP number", names.arg = x, main="SNP per scaffold: average (Rhiir2)")





