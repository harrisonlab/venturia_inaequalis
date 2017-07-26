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

mv Ash_farm_172_pacbio_contigs_unmasked_3_bw.vcf SNP_calling/

vcftools=/home/sobczm/bin/vcftools/bin
$vcftools/vcftools --vcf SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_bw.vcf  --max-missing 0.95 --recode --out SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_bw_filtered
```

Create custom SnpEff genome database
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

Collect input files

```bash
Reference=$(ls repeat_masked/v.inaequalis/172_pacbio/filtered_contigs_repmask/172_pacbio_contigs_unmasked.fa)
Gff=$(ls gene_pred/codingquary/v.inaequalis/172_pacbio/final/final_genes_appended.gff3)
SnpEff=/home/sobczm/bin/snpEff
mkdir $SnpEff/data/172_pacbiov1.0
cp $Reference $SnpEff/data/172_pacbiov1.0/sequences.fa
cp $Gff $SnpEff/data/172_pacbiov1.0/genes.gff
```

#Build database using GFF3 annotation
java -jar $SnpEff/snpEff.jar build -gff3 -v 172_pacbiov1.0


<!--
Groups of isolates from different cultivars described, 8 isolates for Worcester (pop1 below) and 6 isolates for Bramley (pop2 below); two Bramley isolates lost due to poor sequencing (036 and 057)

Important: check script options below, in order to use the correct ones in
a given analysis (e.g. "ply" argument for ploidy).

```bash
scripts=/home/sobczm/bin/popgen/summary_stats
python $scripts/vcf_find_difference_pop.py --vcf Ash_farm_172_pacbio_contigs_unmasked_bw_filtered.recode.vcf --out Ash_farm_172_pacbio_contigs_unmasked_bw_filtered_fixed.vcf --ply 1 --pop1 202,,182,,173,,190,,172,,197,,196,,049 --pop2 024,,030,,007,,025,,044,,199 --thr 0.95
```

##V. inequalis structural variants
```bash
input=/home/sobczm/popgen/snp/sv_calling
python $scripts/vcf_find_difference_pop.py --vcf $input/vinequalis/vinequalis_struc_variants.vcf --out $input/vinequalis/vinequalis_struc_variants_fixed.vcf --ply 1 --pop1 202,,182,,173,,190,,172,,197,,196 --pop2 057,,024,,030,,007,,025,,044 --thr 0.95
```