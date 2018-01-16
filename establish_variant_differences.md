Work performed in directory
```bash
/home/groups/harrisonlab/project_files/venturia
```

After calling SNPs and stuctural variants, it may be useful to contrast variant population
frequencies among different groups of individuals (populations) to e.g. further zoom in on
candidate effector genes likely involved in resistance in different cultivars.

## V. inequalis SNPs
First, create a cut-down VCF file containing only individuals of interest 
(Worcester and Bramley populations) and then filter to remove SNPs with too many missing genotypes.

```bash
source /home/sobczm/bin/marias_profile

vcflib=/home/sobczm/bin/vcflib/bin

$vcflib/vcfremovesamples SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3.vcf 083 096 097 098 101 106 119 >Ash_farm_172_pacbio_contigs_unmasked_3_bw.vcf 
$vcflib/vcfremovesamples SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3.vcf 049 172 173 182 190 196 197 202 >Ash_farm_172_pacbio_contigs_unmasked_3_bc.vcf
$vcflib/vcfremovesamples SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3.vcf 007 024 025 030 044 199 >Ash_farm_172_pacbio_contigs_unmasked_3_cw.vcf

mv Ash_farm_172_pacbio_contigs_unmasked_3_bw.vcf SNP_calling/
mv Ash_farm_172_pacbio_contigs_unmasked_3_bc.vcf SNP_calling/
mv Ash_farm_172_pacbio_contigs_unmasked_3_cw.vcf SNP_calling/

vcftools=/home/sobczm/bin/vcftools/bin

$vcftools/vcftools --vcf SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_bw.vcf  --max-missing 0.95 --recode --out SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_bw_filtered
$vcftools/vcftools --vcf SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_bc.vcf  --max-missing 0.95 --recode --out SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_bc_filtered
$vcftools/vcftools --vcf SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_cw.vcf  --max-missing 0.95 --recode --out SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_cw_filtered
```

## Create custom SnpEff genome database

```bash
SnpEff=/home/sobczm/bin/snpEff
nano $SnpEff/snpEff.config
```

Add the following lines to the section with databases:
```bash
#---
# EMR Databases
#----
# Fus2 genome
Fus2v1.0.genome : Fus2
# Bc16 genome
Bc16v1.0.genome: BC-16
# P414 genome
P414v1.0.genome: 414
# 172_pacbio genome
172_pacbiov1.0.genome: 172_pacbio
```

# Collect input files

```bash
Reference=$(ls repeat_masked/v.inaequalis/172_pacbio/filtered_contigs_repmask/172_pacbio_contigs_unmasked.fa)
Gff=$(ls gene_pred/codingquary/v.inaequalis/172_pacbio/final/final_genes_appended.gff3)
SnpEff=/home/sobczm/bin/snpEff
mkdir $SnpEff/data/172_pacbiov1.0
cp $Reference $SnpEff/data/172_pacbiov1.0/sequences.fa
cp $Gff $SnpEff/data/172_pacbiov1.0/genes.gff
```

# Build database using GFF3 annotation
```bash
java -jar $SnpEff/snpEff.jar build -gff3 -v 172_pacbiov1.0
```

## Annotate VCF files

```bash
for a in $(ls SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_bw_filtered.recode.vcf); do
    echo $a
    filename=$(basename "$a")
    SnpEff=/home/sobczm/bin/snpEff
    java -Xmx4g -jar $SnpEff/snpEff.jar -v -ud 0 172_pacbiov1.0 $a > ${filename%.vcf}_annotated.vcf
    mv snpEff_genes.txt SNP_calling/snpEff_genes_${filename%.vcf}.txt
    mv snpEff_summary.html SNP_calling/snpEff_summary_${filename%.vcf}.html
    mv *_filtered* SNP_calling/.
done
```

<!--
```bash
vcf=SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_bw_filtered.recode.vcf
genome_name=172_pacbio

snpeff=/home/sobczm/bin/snpEff
scripts=/home/sobczm/bin/popgen

#Create subsamples of SNPs containing those in a given category

#genic (includes 5', 3' UTRs)
java -jar $snpeff/SnpSift.jar filter "(ANN[*].EFFECT has 'missense_variant') || (ANN[*].EFFECT has 'nonsense_variant') || (ANN[*].EFFECT has 'synonymous_variant') || (ANN[*].EFFECT has 'intron_variant') || (ANN[*].EFFECT has '5_prime_UTR_variant') || (ANN[*].EFFECT has '3_prime_UTR_variant')" ${vcf%.vcf}_annotated.vcf > ${vcf%.vcf}_gene.vcf
#coding
java -jar $snpeff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant') || (ANN[0].EFFECT has 'synonymous_variant')" ${vcf%.vcf}_annotated.vcf > ${vcf%.vcf}_coding.vcf
#non-synonymous
java -jar $snpeff/SnpSift.jar filter "(ANN[0].EFFECT has 'missense_variant') || (ANN[0].EFFECT has 'nonsense_variant')" ${vcf%.vcf}_annotated.vcf > ${vcf%.vcf}_nonsyn.vcf
#synonymous
java -jar $snpeff/SnpSift.jar filter "(ANN[0].EFFECT has 'synonymous_variant')" ${vcf%.vcf}_annotated.vcf > ${vcf%.vcf}_syn.vcf
#Four-fold degenrate sites (output file suffix: 4fd)
python $scripts/summary_stats/parse_snpeff_synonymous.py ${vcf%.vcf}_syn.vcf
```
-->

Groups of isolates from different cultivars described, 8 isolates for Worcester (pop1 below) and 6 isolates for Bramley (pop2 below); two Bramley isolates lost due to poor sequencing (036 and 057)

Important: check script options below, in order to use the correct ones in
a given analysis (e.g. "ply" argument for ploidy).

```bash
scripts=/home/sobczm/bin/popgen/summary_stats
python $scripts/vcf_find_difference_pop.py --vcf SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_bw_filtered.recode.vcf --out SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_bw_filtered_fixed.vcf --ply 1 --pop1 202,,182,,173,,190,,172,,197,,196,,049 --pop2 024,,030,,007,,025,,044,,199 --thr 0.95
```

# V. inequalis structural variants
```bash
scripts=/home/sobczm/bin/popgen/summary_stats
input=/home/sobczm/popgen/snp/sv_calling
python $scripts/vcf_find_difference_pop.py --vcf $input/vinequalis/vinequalis_struc_variants.vcf --out SNP_calling/Ash_farm_struc_variants_fixed.vcf --ply 1 --pop1 202,,182,,173,,190,,172,,197,,196,,049 --pop2 024,,030,,007,,025,,044,,199 --thr 0.95
