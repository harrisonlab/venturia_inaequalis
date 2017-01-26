#$ -S /bin/bash
#$ -cwd
#$ -pe smp 6
#$ -l h_vmem=4G
#$ -l h=blacklace01.blacklace|blacklace02.blacklace|blacklace04.blacklace|blacklace05.blacklace|blacklace06.blacklace|blacklace07.blacklace|blacklace08.blacklace|blacklace09.blacklace|blacklace10.blacklace|blacklace12.blacklace

# Testing parallelisation of GATk HaplotypeCaller - may crash. (It did not! Resulted in 2x speedup)
# NOTE: this is a haploid organism. For diploid organism, change "ploidy" argument to 2.
# Changes required in the script:
# VARIABLES
# Reference - the genome reference used in read mapping.
# INSIDE THE GATK command:
# To specify which BAM mapping files (output from pre_SNP_calling_cleanup.sh, filename ending with "_rg" -> that is, with 
# read group added) are to be used in SNP calling, use the -I argument with full path to each file following after that. 
# Each new BAM file has to be specified after a separate -I

input=/home/groups/harrisonlab/project_files/venturia/alignment/bowtie/v.inaequalis
reference=/home/groups/harrisonlab/project_files/venturia/repeat_masked/v.inaequalis/172_pacbio/filtered_contigs_repmask/172_pacbio_contigs_unmasked.fa

filename=$(basename "$reference")
output="${filename%.*}_temp.vcf"
output2="${filename%.*}.vcf"

gatk=/home/sobczm/bin/GenomeAnalysisTK-3.6

java -jar $gatk/GenomeAnalysisTK.jar \
     -R $reference \
     -T HaplotypeCaller \
     -ploidy 1 \
     -nct 6 \
     --allow_potentially_misencoded_quality_scores \
     -I $input/007/vs_172_PacBio/007_172_pacbio_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/024/vs_172_PacBio/024_172_pacbio_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/025/vs_172_PacBio/025_172_pacbio_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/030/vs_172_PacBio/030_172_pacbio_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/036/vs_172_PacBio/036_172_pacbio_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/044/vs_172_PacBio/044_172_pacbio_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/049/vs_172_PacBio/049_172_pacbio_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/057/vs_172_PacBio/057_172_pacbio_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/083/vs_172_PacBio/083_172_pacbio_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/096/vs_172_PacBio/096_172_pacbio_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/097/vs_172_PacBio/097_172_pacbio_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/098/vs_172_PacBio/098_172_pacbio_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/101/vs_172_PacBio/101_172_pacbio_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/106/vs_172_PacBio/106_172_pacbio_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/118/vs_172_PacBio/118_172_pacbio_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/119/vs_172_PacBio/119_172_pacbio_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/172/vs_172_PacBio/172_172_pacbio_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/173/vs_172_PacBio/173_172_pacbio_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/182/vs_172_PacBio/182_172_pacbio_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/190/vs_172_PacBio/190_172_pacbio_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/196/vs_172_PacBio/196_172_pacbio_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/197/vs_172_PacBio/197_172_pacbio_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/199/vs_172_PacBio/199_172_pacbio_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/202/vs_172_PacBio/202_172_pacbio_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -I $input/saturn/vs_172_PacBio/saturn_172_pacbio_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam \
     -o $output

#Break down complex SNPs into primitive ones with VariantsToAllelicPrimitives
#This tool will take an MNP (e.g. ACCCA -> TCCCG) and break it up into separate records for each component part (A-T and A->G).
#This tool modifies only bi-allelic variants.

java -jar $gatk/GenomeAnalysisTK.jar \
   -T VariantsToAllelicPrimitives \
   -R $input/$reference \
   -V $output \
   -o $output2 \


#####################################
# Notes on GATK parallelisation
#####################################
# http://gatkforums.broadinstitute.org/gatk/discussion/1975/how-can-i-use-parallelism-to-make-gatk-tools-run-faster
