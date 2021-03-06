# SNP comparison of V. inaequalis isolates from different cultivars in same mixed orchard

Work performed in directory:
```
/home/groups/harrisonlab/project_files/venturia
```

## Alignment of trimmed MiSeq reads of isolates to Assembled PacBio genome

```bash
Assembly=$(ls repeat_masked/v.inaequalis/172_pacbio/filtered_contigs_repmask/172_pacbio_contigs_unmasked.fa)
	for IsolateDir in $(ls -d qc_dna/paired/v.inaequalis/*); do
		Jobs=$(qstat | grep 'sub_bo' | grep 'qw' | wc -l)
		while [ $Jobs -gt 1 ]; do
		sleep 10
		printf "."
		Jobs=$(qstat | grep 'sub_bo' | grep 'qw' | wc -l)
		done
	Organism=$(echo $IsolateDir | rev | cut -d '/' -f2 | rev)
	Isolate=$(echo $IsolateDir | rev | cut -d '/' -f1 | rev)
	echo "$Organism - $Strain"
	Read_F=$(ls $IsolateDir/F/*.fq.gz)
	Read_R=$(ls $IsolateDir/R/*.fq.gz)
	echo $Read_F
	echo $Read_R
	OutDir=alignment/bowtie/$Organism/$Isolate/vs_172_PacBio
	ProgDir=/home/passet/git_repos/tools/seq_tools/genome_alignment/bowtie/
	qsub $ProgDir/sub_bowtie.sh $Assembly $Read_F $Read_R $OutDir
	done
```
## Pre SNP calling cleanup

Setting up correct formatting for SNP analysis

```
input=alignment/bowtie/
scripts=home/passet/git_repos/scripts/popgen/snp
```

### Rename input mapping files in each folder by prefixing with the strain ID

```bash
	cd $input/*/007/*
	for filename in 172_pacbio_contigs_unmasked.fa_aligned.sam
	do
		mv "$filename" "007_$filename"
	done

	cd $input/*/024/*
	for filename in 172_pacbio_contigs_unmasked.fa_aligned.sam
	do
		mv "$filename" "024_$filename"
	done

	cd $input/*/025/*
	for filename in 172_pacbio_contigs_unmasked.fa_aligned.sam
	do
		mv "$filename" "025_$filename"
	done

	cd ../../../../..
	cd $input/*/030/*
	for filename in 172_pacbio_contigs_unmasked.fa_aligned.sam
	do
		mv "$filename" "030_$filename"
	done

	cd ../../../../..
	cd $input/*/036/*
	for filename in 172_pacbio_contigs_unmasked.fa_aligned.sam
	do
		mv "$filename" "036_$filename"
	done

	cd ../../../../..
	cd $input/*/044/*
	for filename in 172_pacbio_contigs_unmasked.fa_aligned.sam
	do
		mv "$filename" "044_$filename"
	done

	cd ../../../../..
	cd $input/*/049/*
	for filename in 172_pacbio_contigs_unmasked.fa_aligned.sam
	do
		mv "$filename" "049_$filename"
	done

	cd ../../../../..
	cd $input/*/057/*
	for filename in 172_pacbio_contigs_unmasked.fa_aligned.sam
	do
		mv "$filename" "057_$filename"
	done

	cd ../../../../..
	cd $input/*/083/*
	for filename in 172_pacbio_contigs_unmasked.fa_aligned.sam
	do
		mv "$filename" "083_$filename"
	done

	cd ../../../../..
	cd $input/*/096/*
	for filename in 172_pacbio_contigs_unmasked.fa_aligned.sam
	do
		mv "$filename" "096_$filename"
	done

	cd ../../../../..
	cd $input/*/097/*
	for filename in 172_pacbio_contigs_unmasked.fa_aligned.sam
	do
		mv "$filename" "097_$filename"
	done

	cd ../../../../..
	cd $input/*/098/*
	for filename in 172_pacbio_contigs_unmasked.fa_aligned.sam
	do
		mv "$filename" "098_$filename"
	done

	cd ../../../../..
	cd $input/*/101/*
	for filename in 172_pacbio_contigs_unmasked.fa_aligned.sam
	do
		mv "$filename" "101_$filename"
	done

	cd ../../../../..
	cd $input/*/106/*
	for filename in 172_pacbio_contigs_unmasked.fa_aligned.sam
	do
		mv "$filename" "106_$filename"
	done

	cd ../../../../..
	cd $input/*/118/*
	for filename in 172_pacbio_contigs_unmasked.fa_aligned.sam
	do
		mv "$filename" "118_$filename"
	done

	cd ../../../../..
	cd $input/*/119/*
	for filename in 172_pacbio_contigs_unmasked.fa_aligned.sam
	do
		mv "$filename" "119_$filename"
	done

	cd ../../../../..
	cd $input/*/172/*
	for filename in 172_pacbio_contigs_unmasked.fa_aligned.sam
	do
		mv "$filename" "172_$filename"
	done

	cd ../../../../..
	cd $input/*/173/*
	for filename in 172_pacbio_contigs_unmasked.fa_aligned.sam
	do
		mv "$filename" "173_$filename"
	done

	cd ../../../../..
	cd $input/*/182/*
	for filename in 172_pacbio_contigs_unmasked.fa_aligned.sam
	do
		mv "$filename" "182_$filename"
	done

	cd ../../../../..
	cd $input/*/190/*
	for filename in 172_pacbio_contigs_unmasked.fa_aligned.sam
	do
		mv "$filename" "190_$filename"
	done

	cd ../../../../..
	cd $input/*/196/*
	for filename in 172_pacbio_contigs_unmasked.fa_aligned.sam
	do
		mv "$filename" "196_$filename"
	done

	cd ../../../../..
	cd $input/*/197/*
	for filename in 172_pacbio_contigs_unmasked.fa_aligned.sam
	do
		mv "$filename" "197_$filename"
	done

	cd ../../../../..
	cd $input/*/199/*
	for filename in 172_pacbio_contigs_unmasked.fa_aligned.sam
	do
		mv "$filename" "199_$filename"
	done

	cd ../../../../..
	cd $input/*/202/*
	for filename in 172_pacbio_contigs_unmasked.fa_aligned.sam
	do
		mv "$filename" "202_$filename"
	done

	cd ../../../../..
	cd $input/*/saturn/*
	for filename in 172_pacbio_contigs_unmasked.fa_aligned.sam
	do
		mv "$filename" "saturn_$filename"
	done
```


### Calculate read coverage over each contig

This wasnt required for SNP calling, but important for reporting alignment statistics
```bash
for Bam in $(ls alignment/bowtie/v.inaequalis/*/vs_172_PacBio/*aligned_sorted.bam); do
	Strain=$(echo $Bam | rev | cut -f3 -d '/' | rev)
	Organism=$(echo $Bam | rev | cut -f4 -d '/' | rev)
	Subject=$(echo $Bam | rev | cut -f2 -d '/' | rev)
	echo "$Organism - $Strain"
	# OutDir=$(dirname $Bam)
	OutDir=/data/scratch/armita/venturia/$(dirname $Bam)
	mkdir -p $OutDir
	samtools depth -aa $Bam > $OutDir/${Organism}_${Strain}_${Subject}_depth.tsv
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/genome_alignment/coverage_analysis
	$ProgDir/cov_by_window.py --cov $OutDir/${Organism}_${Strain}_${Subject}_depth.tsv > $OutDir/${Organism}_${Strain}_${Subject}_depth_10kb.tsv
	sed -i "s/$/\t$Strain/g" $OutDir/${Organism}_${Strain}_${Subject}_depth_10kb.tsv
done
OutDir=/data/scratch/armita/venturia/analysis/genome_alignment/bowtie/grouped
mkdir -p $OutDir
cat /data/scratch/armita/venturia/analysis/genome_alignment/bowtie/*/*/*_${Subject}_unmasked/*_*_${Subject}_depth_10kb.tsv > analysis/genome_alignment/bowtie/grouped/vs_12008_grouped_depth.tsv


for Cov in $(ls /data/scratch/armita/venturia/alignment/bowtie/v.inaequalis/*/vs_172_PacBio/v.inaequalis_*_vs_172_PacBio_depth.tsv); do
	echo $(basename $Cov)
	echo "Total Depth"
	cat $Cov | cut -f3 | sort | uniq -c | sort -nr | head -n10
	echo "Depth per 10 kb"
	cat $Cov | cut -f3 | sort | uniq -c | sort -nr | head -n5
	echo ""
done > /data/scratch/armita/venturia/analysis/genome_alignment/bowtie/grouped/summary_depth.txt

mkdir /data/scratch/armita/venturia/alignment/bowtie/v.inaequalis/coverage
cp -s /data/scratch/armita/venturia/alignment/bowtie/v.inaequalis/*/vs_172_PacBio/v.inaequalis_*_vs_172_PacBio_depth.tsv /data/scratch/armita/venturia/alignment/bowtie/v.inaequalis/coverage

R
# df0 <- data.frame()
files <- list.files(path="/data/scratch/armita/venturia/alignment/bowtie/v.inaequalis/coverage", pattern="*_vs_172_PacBio_depth.tsv", full.names=TRUE, recursive=FALSE)
lapply(files, function(x) {
		out <- gsub("_PacBio_depth.tsv", "_PacBio_median.txt", x)
    t <- read.table(x, header=FALSE, sep="\t") # load file
		m1 <- median(t$V3)
		t2 <- subset(t, V3!=0)
		m2 <- median(t2$V3)
		# df1 <- df0
		df2 <- data.frame(x, m1, m2)
		write.table(df2, out, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
		# df0 <- rbind(df1, df2)
})
# write.table(df0, "/data/scratch/armita/venturia/alignment/bowtie/v.inaequalis/coverage/summarised_cov.tsv", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
```


### Remove multimapping reads, discordant reads. PCR and optical duplicates, and add read group and sample name to each mapped read (preferably, the shortest ID possible)

Checked with one isolate first

```bash
	for Alignment in $(ls alignment/bowtie/v.inaequalis/*/*/*.fa_aligned.sam | grep -w "007"); do
		Isolate=$(echo $Alignment | rev | cut -d '/' -f3 | rev)
		echo $Alignment
		echo $Isolate
		ProgDir=/home/passet/git_repos/scripts/popgen/snp
		qsub $ProgDir/sub_pre_snp_calling.sh $Alignment $Isolate
	done
```

Ran cleanup of rest of alignments
```bash
	for Alignment in $(ls alignment/bowtie/v.inaequalis/*/*/*.fa_aligned.sam | grep -v "007"); do
			Jobs=$(qstat | grep 'sub_pr' | grep 'qw' | wc -l)
			while [ $Jobs -gt 1 ]; do
			sleep 10
			printf "."
			Jobs=$(qstat | grep 'sub_pr' | grep 'qw' | wc -l)
			done
		Isolate=$(echo $Alignment | rev | cut -d '/' -f3 | rev)
		echo $Alignment
		echo $Isolate
		ProgDir=/home/passet/git_repos/scripts/popgen/snp
		qsub $ProgDir/sub_pre_snp_calling.sh $Alignment $Isolate
	done
```
Copy cleanup outputs into alignment folders
```bash
	for Isolate in 007 024 025 030 036 044 049 057 083 096 097 098 101 106 118 119 172 173 182 190 196 197 199 202 saturn
	do
		Bam="$Isolate"_172_pacbio_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam
		Bai="$Isolate"_172_pacbio_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup_rg.bam.bai
		Txt="$Isolate"_172_pacbio_contigs_unmasked.fa_aligned_nomulti_proper_sorted_nodup.txt
		Directory=alignment/bowtie/*/$Isolate/vs_172_PacBio
		mv $Bam $Directory
		mv $Bai $Directory
		mv $Txt $Directory
	done
```

## SNP calling
Run a SNP calling script developed by Maria

To change in each analysis:

```bash
input=/home/groups/harrisonlab/project_files/venturia/alignment/bowtie
reference=/home/groups/harrisonlab/project_files/venturia/repeat_masked/v.inaequalis/172_pacbio/filtered_contigs_repmask/172_pacbio_contigs_unmasked.fa


filename=$(basename "$reference")
output="${filename%.*}.dict"
```

### Prepare genome reference indexes required by GATK
```bash
java -jar /home/sobczm/bin/picard-tools-2.5.0/picard.jar CreateSequenceDictionary R=$reference O=$input/$output
samtools faidx $reference
```
.dict file saved in wrong place
```bash
cp alignment/bowtie/172_pacbio_contigs_unmasked.dict repeat_masked/v.inaequalis/172_pacbio/filtered_contigs_repmask/
```

Copy index file to same folder as BAM alignments

```bash
for Isolate in 007 024 025 030 036 044 049 057 083 096 097 098 101 106 118 119 172 173 182 190 196 197 199 202 saturn
	do
    Index=repeat_masked/v.inaequalis/172_pacbio/filtered_contigs_repmask/172_pacbio_contigs_unmasked.fa.fai
    Directory=alignment/bowtie/*/$Isolate/vs_172_PacBio/
    cp $Index $Directory
done
```

Move to the directory where the output of SNP calling should be placed

```bash
mkdir -p /home/groups/harrisonlab/project_files/venturia/SNP_calling
cd /home/groups/harrisonlab/project_files/venturia/SNP_calling
```

### Start SNP calling with GATK
The submission script required needs to be custom-prepared for each analysis, depending on what samples are being analysed.
See inside the submission script below.

```bash
scripts=/home/passet/git_repos/scripts/venturia_inaequalis
qsub $scripts/sub_SNP_calling_multithreaded.sh
```

## Determining genentic structure


### Removal of Isolates for analysis
Want to be able to run analysis without:
Isolate 036 due to poor quality of sequencing (due to initial low library concentration on to MiSeq)
Saturn isolate as not from the same orchard (and therefore not currently of interest)

```bash
source /home/sobczm/bin/marias_profile
vcftools=/home/sobczm/bin/vcftools/bin
vcflib=/home/sobczm/bin/vcflib/bin

$vcflib/vcfremovesamples SNP_calling/172_pacbio_contigs_unmasked.vcf 036 saturn >SNP_calling/Ash_farm_172_pacbio_contigs_unmasked.vcf
```
Isolate 118 also has poor quality sequencing and therefore also removed

```bash
source /home/sobczm/bin/marias_profile
vcftools=/home/sobczm/bin/vcftools/bin
vcflib=/home/sobczm/bin/vcflib/bin

$vcflib/vcfremovesamples SNP_calling/Ash_farm_172_pacbio_contigs_unmasked.vcf 118 >SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_2.vcf
```
Isolate 057 also noticed as being poorly sequenced after removal of isolate 118 so also removed isolate 057 to form a 21 isolate group

```bash
source /home/sobczm/bin/marias_profile
vcftools=/home/sobczm/bin/vcftools/bin
vcflib=/home/sobczm/bin/vcflib/bin

$vcflib/vcfremovesamples SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_2.vcf 057 >SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3.vcf
```


### Only retain biallelic high-quality SNPS with no missing data for genetic analyses.
```bash
	for vcf in $(ls SNP_calling/*_contigs_unmasked.vcf)
	do
		echo $vcf
		script=/home/passet/git_repos/scripts/popgen/snp
		qsub $script/sub_vcf_parser.sh $vcf
	done
```
Repeated on revised Ash Farm only
```bash
	for vcf in $(ls SNP_calling/Ash_farm_172_pacbio_contigs_unmasked.vcf)
	do
		echo $vcf
		script=/home/passet/git_repos/scripts/popgen/snp
		qsub $script/sub_vcf_parser.sh $vcf
	done
```

```bash
	for vcf in $(ls SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_2.vcf)
	do
		echo $vcf
		script=/home/passet/git_repos/scripts/popgen/snp
		qsub $script/sub_vcf_parser.sh $vcf
	done
```

Repeated for 21 isolate group
```bash
	for vcf in $(ls SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3.vcf)
	do
		echo $vcf
		script=/home/passet/git_repos/scripts/popgen/snp
		qsub $script/sub_vcf_parser.sh $vcf
	done
```


### In some organisms, may want to thin (subsample) SNPs in high linkage diseqilibrium down to
1 SNP  per e.g. 10 kbp just for the population structure analyses. Files had to be renamed after running
```bash
input=SNP_calling
vcftools=/home/sobczm/bin/vcftools/bin
$vcftools/vcftools --vcf $input/172_pacbio_contigs_unmasked.vcf --thin 10000 --recode --out ${input/vcf%.vcf}172_pacbio_contigs_unmasked_thinned
mv SNP_calling172_pacbio_contigs_unmasked_thinned.log SNP_calling/172_pacbio_contigs_unmasked_thinned.log
mv SNP_calling172_pacbio_contigs_unmasked_thinned.recode.vcf SNP_calling/172_pacbio_contigs_unmasked_thinned.recode.vcf

input=SNP_calling
vcftools=/home/sobczm/bin/vcftools/bin
$vcftools/vcftools --vcf $input/172_pacbio_contigs_unmasked_filtered.vcf --thin 10000 --recode --out ${input/vcf%.vcf}172_pacbio_contigs_unmasked_filtered_thinned
mv SNP_calling172_pacbio_contigs_unmasked_filtered_thinned.log SNP_calling/172_pacbio_contigs_unmasked_filtered_thinned.log
mv SNP_calling172_pacbio_contigs_unmasked_filtered_thinned.recode.vcf SNP_calling/172_pacbio_contigs_unmasked_filtered_thinned.recode.vcf

input=SNP_calling
vcftools=/home/sobczm/bin/vcftools/bin
$vcftools/vcftools --vcf $input/Ash_farm_172_pacbio_contigs_unmasked.vcf --thin 10000 --recode --out ${input/vcf%.vcf}Ash_farm_172_pacbio_contigs_unmasked_thinned
mv SNP_callingAsh_farm_172_pacbio_contigs_unmasked_thinned.log SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_thinned.log
mv SNP_callingAsh_farm_172_pacbio_contigs_unmasked_thinned.recode.vcf SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_thinned.recode.vcf

input=SNP_calling
vcftools=/home/sobczm/bin/vcftools/bin
$vcftools/vcftools --vcf $input/Ash_farm_172_pacbio_contigs_unmasked_filtered.vcf --thin 10000 --recode --out ${input/vcf%.vcf}Ash_farm_172_pacbio_contigs_unmasked_filtered_thinned
mv SNP_callingAsh_farm_172_pacbio_contigs_unmasked_filtered_thinned.log SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_filtered_thinned.log
mv SNP_callingAsh_farm_172_pacbio_contigs_unmasked_filtered_thinned.recode.vcf SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_filtered_thinned.recode.vcf
```

Re-ran but with less severe thinning to 1 SNP per 1000 bp (as opposed to 10000 above)

```bash
input=SNP_calling
vcftools=/home/sobczm/bin/vcftools/bin
$vcftools/vcftools --vcf $input/172_pacbio_contigs_unmasked.vcf --thin 1000 --recode --out ${input/vcf%.vcf}172_pacbio_contigs_unmasked_thinned_1000
mv SNP_calling172_pacbio_contigs_unmasked_thinned_1000.log SNP_calling/172_pacbio_contigs_unmasked_thinned_1000.log
mv SNP_calling172_pacbio_contigs_unmasked_thinned_1000.recode.vcf SNP_calling/172_pacbio_contigs_unmasked_thinned_1000.recode.vcf

input=SNP_calling
vcftools=/home/sobczm/bin/vcftools/bin
$vcftools/vcftools --vcf $input/172_pacbio_contigs_unmasked_filtered.vcf --thin 1000 --recode --out ${input/vcf%.vcf}172_pacbio_contigs_unmasked_filtered_thinned_1000
mv SNP_calling172_pacbio_contigs_unmasked_filtered_thinned_1000.log SNP_calling/172_pacbio_contigs_unmasked_filtered_thinned_1000.log
mv SNP_calling172_pacbio_contigs_unmasked_filtered_thinned_1000.recode.vcf SNP_calling/172_pacbio_contigs_unmasked_filtered_thinned_1000.recode.vcf

input=SNP_calling
vcftools=/home/sobczm/bin/vcftools/bin
$vcftools/vcftools --vcf $input/Ash_farm_172_pacbio_contigs_unmasked_2.vcf --thin 1000 --recode --out ${input/vcf%.vcf}Ash_farm_172_pacbio_contigs_unmasked_2_thinned_1000
mv SNP_callingAsh_farm_172_pacbio_contigs_unmasked_2_thinned_1000.log SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_2_thinned_1000.log
mv SNP_callingAsh_farm_172_pacbio_contigs_unmasked_2_thinned_1000.recode.vcf SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_thinned_2_1000.recode.vcf

input=SNP_calling
vcftools=/home/sobczm/bin/vcftools/bin
$vcftools/vcftools --vcf $input/Ash_farm_172_pacbio_contigs_unmasked_2_filtered.vcf --thin 1000 --recode --out ${input/vcf%.vcf}Ash_farm_172_pacbio_contigs_unmasked_2_filtered_thinned_1000
mv SNP_callingAsh_farm_172_pacbio_contigs_unmasked_2_filtered_thinned_1000.log SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_2_filtered_thinned_1000.log
mv SNP_callingAsh_farm_172_pacbio_contigs_unmasked_2_filtered_thinned_1000.recode.vcf SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_2_filtered_thinned_1000.recode.vcf
```

```
This kept ~26000 of ~914000 variants
```

<!--
```bash
qlogin
cd /projects/oldhome/groups/harrisonlab/project_files/venturia
vcftools=/home/sobczm/bin/vcftools/bin
$vcftools/vcftools --vcf SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3.vcf --max-missing 1 --remove-indels --recode --out SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_no-missing
```

```
After filtering, kept 526661 out of a possible 914651 Sites
Run Time = 64.00 seconds
```
-->

Ran thinning to 1 SNP per 1000 bp on 21 isolate group

```bash
input=SNP_calling
vcftools=/home/sobczm/bin/vcftools/bin
$vcftools/vcftools --vcf $input/Ash_farm_172_pacbio_contigs_unmasked_3.vcf --thin 1000 --recode --out ${input/vcf%.vcf}Ash_farm_172_pacbio_contigs_unmasked_3_thinned_1000
mv SNP_callingAsh_farm_172_pacbio_contigs_unmasked_3_thinned_1000.log SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_thinned_1000.log
mv SNP_callingAsh_farm_172_pacbio_contigs_unmasked_3_thinned_1000.recode.vcf SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_thinned_3_1000.recode.vcf

input=SNP_calling
vcftools=/home/sobczm/bin/vcftools/bin
$vcftools/vcftools --vcf Ash_farm_172_pacbio_contigs_unmasked_3_filtered.vcf --thin 1000 --recode --out ${input/vcf%.vcf}Ash_farm_172_pacbio_contigs_unmasked_3_filtered_thinned_1000
mv SNP_callingAsh_farm_172_pacbio_contigs_unmasked_3_filtered_thinned_1000.log SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_filtered_thinned_1000.log
mv SNP_callingAsh_farm_172_pacbio_contigs_unmasked_3_filtered_thinned_1000.recode.vcf SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_filtered_thinned_1000.recode.vcf
```


### General VCF stats (remember that vcftools needs to have the PERL library exported)
```bash
export PERL5LIB=/home/sobczm/bin/vcftools/share/perl/5.14.2

source /home/sobczm/bin/marias_profile
perl /home/sobczm/bin/vcftools/bin/vcf-stats \
SNP_calling/172_pacbio_contigs_unmasked.vcf >SNP_calling/172_pacbio_contigs_unmasked.stat

perl /home/sobczm/bin/vcftools/bin/vcf-stats \
SNP_calling/172_pacbio_contigs_unmasked_filtered.vcf >SNP_calling/172_pacbio_contigs_unmasked_filtered.stat

perl /home/sobczm/bin/vcftools/bin/vcf-stats \
SNP_calling/Ash_farm_172_pacbio_contigs_unmasked.vcf >SNP_calling/Ash_farm_172_pacbio_contigs_unmasked.stat

perl /home/sobczm/bin/vcftools/bin/vcf-stats \
SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_filtered.vcf >SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_filtered.stat

perl /home/sobczm/bin/vcftools/bin/vcf-stats \
SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_2.vcf >SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_2.stat

perl /home/sobczm/bin/vcftools/bin/vcf-stats \
SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_2_filtered.vcf >SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_2_filtered.stat

perl /home/sobczm/bin/vcftools/bin/vcf-stats \
SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3.vcf >SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3.stat

perl /home/sobczm/bin/vcftools/bin/vcf-stats \
SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_filtered.vcf >SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_filtered.stat
```


### Calculate the index for percentage of shared SNP alleles between the individuals
```bash
	for vcf in $(ls SNP_calling/*_contigs_unmasked_filtered.vcf)
	do
		scripts=/home/passet/git_repos/scripts/popgen/snp
		$scripts/similarity_percentage.py $vcf
	done
```
For the 22 isolates
```bash
	for vcf in $(ls SNP_calling/*_contigs_unmasked_2_filtered.vcf)
	do
		scripts=/home/passet/git_repos/scripts/popgen/snp
		$scripts/similarity_percentage.py $vcf
	done
```

For the 21 isolate group
```bash
	for vcf in $(ls SNP_calling/*_contigs_unmasked_3_filtered.vcf)
	do
		scripts=/home/passet/git_repos/scripts/popgen/snp
		$scripts/similarity_percentage.py $vcf
	done
```


Using R version 3.2.2 installed locally:
```bash
export PATH=/home/armita/prog/R/R-3.2.2/bin:${PATH}
```
And libraries stored in
```bash
export R_LIBS=/home/sobczm/R/x86_64-pc-linux-gnu-library/3.2:$R_LIBS
```

Visualise the output as heatmap and clustering dendrogram
```bash
scripts=/home/passet/git_repos/scripts/popgen/snp
Rscript --vanilla $scripts/distance_matrix.R SNP_calling/172_pacbio_contigs_unmasked_filtered_distance.log

Rscript --vanilla $scripts/distance_matrix.R SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_filtered_distance.log

Rscript --vanilla $scripts/distance_matrix.R SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_2_filtered_distance.log

Rscript --vanilla $scripts/distance_matrix.R SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_filtered_distance.log
```
<!--
Final plot for publication:
```bash
ProgDir=/projects/oldhome/armita/git_repos/emr_repos/scripts/venturia_inaequalis
Rscript --vanilla $ProgDir/plot_SNP_distance_matrix.R SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_filtered_distance.log
``` -->

The final distance plot used for publication was made after filtering the vcf
file to remove any missing data, including where it was present for a single
isolate. These commands are documented below.

### Carry out PCA and plot the results
```bash
scripts=/home/passet/git_repos/scripts/popgen/snp

Rscript --vanilla $scripts/pca.R SNP_calling/172_pacbio_contigs_unmasked_filtered.vcf

Rscript --vanilla $scripts/pca.R SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_filtered.vcf

Rscript --vanilla $scripts/pca.R SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_2_filtered.vcf

Rscript --vanilla $scripts/pca.R SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_filtered.vcf

```

Create R plot but with sample cultivar identifier and same scale on both axis

```bash
scripts=/home/passet/git_repos/scripts/popgen/snp

Rscript --vanilla $scripts/pca_axis_scale.R SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_filtered_label.vcf

```

### Calculate an NJ tree based on all the SNPs. Outputs a basic diplay of the tree, plus a Newick file to be used for displaying the tree in FigTree and beautifying it.

```bash
scripts=/home/passet/git_repos/scripts/popgen/snp
$scripts/nj_tree.sh SNP_calling/172_pacbio_contigs_unmasked_filtered.vcf

$scripts/nj_tree.sh SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_filtered.vcf
```
No tree produced above as no missing data allowed


Need to filter SNPs to retain those with no missing data in any individual
<!-- ```bash
scripts=/home/sobczm/bin/popgen/snp
qsub $scripts/sub_vcf_parser.sh SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_2_filtered.vcf 40 30 10 30 1 Y
```

And for final 21 isolates: -->
<!-- ```bash
scripts=/home/armita/git_repos/emr_repos/scripts/popgen/snp
qsub $scripts/sub_vcf_parser.sh SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_filtered.vcf 40 30 10 30 1 Y
``` -->

This was performed for the final set of 21 isolates.

```bash
	vcftools=/home/sobczm/bin/vcftools/bin
	$vcftools/vcftools --vcf SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_filtered.vcf --max-missing 1 --remove-indels --recode --out SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_filtered_no-missing
```

The final similarity matrix was made from this vcf file:

```bash
for vcf in $(ls SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3*_no-missing.recode.vcf); do
scripts=/home/passet/git_repos/scripts/popgen/snp
$scripts/similarity_percentage.py $vcf
done
```

Final plot for publication:

```bash
ProgDir=/projects/oldhome/armita/git_repos/emr_repos/scripts/venturia_inaequalis
# Rscript --vanilla $ProgDir/plot_SNP_distance_matrix.R SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_filtered_no-missing.recode_distance.log
Rscript --vanilla $ProgDir/plot_SNP_distance_matrix.R SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_no-missing.recode_distance.log
```



```bash
scripts=/home/passet/git_repos/scripts/popgen/snp
# $scripts/nj_tree.sh SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_filtered.vcf 1
# $scripts/nj_tree.sh SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_filtered_no-missing.recode.vcf 1
$scripts/nj_tree.sh SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_no-missing.recode.vcf 1
```

Run the following from Rstudio on my local computer:

```r
setwd("/Users/armita/Downloads/scab")
#===============================================================================
#       Load libraries
#===============================================================================

library(ape)
library(ggplot2)

library(ggtree)
library(phangorn)
library(treeio)

t <- read.tree("/Users/armita/Downloads/scab/Ash_farm_172_pacbio_contigs_unmasked_3_filtered_no-missing.recode_nj.nwk")
layout(matrix(1:4, 2, 2))
ape::plot.phylo(t, type = "unrooted")
ape::plot.phylo(t, type = "unrooted", use.edge.length = FALSE)
ape::plot.phylo(t, type = "fan")
ape::plot.phylo(t, type = "fan", use.edge.length = FALSE)
layout(matrix(1))

p1 <- ggtree(t, layout="unrooted") + geom_tiplab()
p2 <- ggtree(t, layout="unrooted", branch.length="none") + geom_tiplab()
p3 <- ggtree(t, layout="circular") + geom_tiplab2()
p4 <- ggtree(t, layout="circular", branch.length="none") + geom_tiplab2()
multiplot(p1, p2, p3, p4, ncol = 2)

ggtree(t, layout="circular") + geom_tiplab2(aes(angle=angle), color='blue')

# show node labels
p3 + geom_text(aes(label=node))

# Format nodes by values
nodes <- data.frame(p3$data)
nodes$label <- as.numeric(nodes$label)
as.numeric(nodes$label)
nodes$support[nodes$isTip] <- '≥ 80%'
nodes$support[(!nodes$isTip) & (nodes$label > 80)] <- '≥ 80%'
nodes$support[(!nodes$isTip) & (nodes$label < 80)] <- '< 80%'
nodes$support[(!nodes$isTip) & (nodes$label == '')] <- '≥ 80%'
p5 <- p3 + aes(linetype=nodes$support)
nodes$label[nodes$label > 80] <- ''
p5 <- p5 + geom_nodelab(data=nodes, size=3, hjust=-0.05)
hilight(node=21, fill="steelblue", alpha=.6)

p3 + geom_hilight(node=30, fill="steelblue", alpha=.6)

nodes$support <- factor(nodes$support, levels = c("≥ 80%", '< 80%'))

support <- nodes$support

cls <- list(cox=c("106", "083", "101", "098", "097", "119", "096"),
            bramley=c("030", "024", "007", "025", "044", "199"),
            worcester=c("172", "182", "202", "173", "049", "197", "190", "196")
            )

t <- groupOTU(t, cls)




library(RColorBrewer)
cols <- brewer.pal(3, "Set2")
p6 <- ggtree(t, layout="circular", aes(color=group)) + geom_tiplab2(offset=+0.01) + scale_color_manual(values=brewer.pal(3, "Dark2")) + theme(legend.position="right") + aes(linetype=support)

# Format nodes by values
nodes2 <- data.frame(p6$data)
nodes2$label <- as.numeric(nodes2$label)
as.numeric(nodes2$label)
nodes2$support[nodes2$isTip] <- '≥ 80%'
nodes2$support[(!nodes2$isTip) & (nodes2$label > 80)] <- '≥ 80%'
nodes2$support[(!nodes2$isTip) & (nodes2$label < 80)] <- '< 80%'
nodes2$support[(!nodes2$isTip) & (nodes2$label == '')] <- '≥ 80%'
nodes2$label[nodes2$label > 80] <- ''

nodes2$support <- factor(nodes2$support, levels = c("≥ 80%", '< 80%'))

support <- nodes2$support

p6 <- p6 + geom_nodelab(data=nodes2, size=2.5, nudge_x=+0.005)

ggsave("final_tree_support.pdf", p6)
ggsave("final_tree_support.tiff", p6)

```

<!--
```bash
scripts=/home/sobczm/bin/popgen/snp
qsub $scripts/sub_vcf_parser.sh SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_filtered.vcf 40 30 10 30 1 Y
```

Prepare tree
```bash
scripts=/home/passet/git_repos/scripts/popgen/snp
$scripts/nj_tree.sh SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_2_filtered.vcf 1

scripts=/home/passet/git_repos/scripts/popgen/snp
$scripts/nj_tree.sh SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_filtered.vcf 1
```
-->

### AMOVA analysis
```bash
scripts=/home/passet/git_repos/scripts/venturia_inaequalis
Rscript --vanilla $scripts/Ash_farm_amova_dapc.R SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_filtered.vcf
```
Move and rename AMOVA file to SNP_calling
```bash
mv amova.txt SNP_calling/Ash_farm_amova.txt
```
<!--
Rscript --vanilla $scripts/amova_dapc.R SNP_calling/172_pacbio_contigs_unmasked_filtered.vcf
```
-->

Ran structure analysis on Ash Farm only samples using the structure_analysis.sh script - ran individual parts of the script rather than submitting running whole script in one

Re-ran structure after removing isolate 118 using structure_analysis_2.sh

Re-ran structure with 21 isolate group using structure_analysis_3.sh; Changed iterations to 1M (from 10M) and burin to 100k (from 1M) - Errors in script so wouldn't work. Submitted jobs again, as below, after Matia created initial input file:

```bash
input=/home/sobczm/popgen/other/passey/Structure_4
scripts=/home/passet/git_repos/scripts/popgen/snp

#Run replicate STRUCTURE runs, with K from 1 to 5
qsub $scripts/execute_structure.sh $input/Ash_farm_172_pacbio_contigs_unmasked_3_filtered_thinned_1000.recode.struc 1 1 5
qsub $scripts/execute_structure.sh $input/Ash_farm_172_pacbio_contigs_unmasked_3_filtered_thinned_1000.recode.struc 1 2 5
qsub $scripts/execute_structure.sh $input/Ash_farm_172_pacbio_contigs_unmasked_3_filtered_thinned_1000.recode.struc 1 3 5
qsub $scripts/execute_structure.sh $input/Ash_farm_172_pacbio_contigs_unmasked_3_filtered_thinned_1000.recode.struc 1 4 5
qsub $scripts/execute_structure.sh $input/Ash_farm_172_pacbio_contigs_unmasked_3_filtered_thinned_1000.recode.struc 1 5 5

#Tidy working directory
mkdir SNP_calling/Structure_3
mv structure* SNP_calling/Structure_3
mv execute_structure.sh.* SNP_calling/Structure_3

#Analyze STRUCTURE output
# Generate a folder containing all STRUCTURE output files for all K analyzed
cd SNP_calling/Structure_3
mkdir structureHarvester_3
for d in $PWD/*
do
mv $d/*_f $PWD/structureHarvester_3
done
cd ../..


# structureHarvester - summarise the results
input=/home/groups/harrisonlab/project_files/venturia/SNP_calling
harvester=/home/sobczm/bin/structureHarvester/structureHarvester.py
$harvester --dir=$input/Structure_3/structureHarvester_3 --out=$input/Structure_3/structureHarvester_3 --evanno --clumpp
# CLUMPP - permute the results
cd SNP_calling/Structure_3/structureHarvester_3
clumpp=/home/sobczm/bin/CLUMPP_Linux64.1.1.2
cp $clumpp/paramfile_ind ./
mv paramfile_ind paramfile

#Options fed to CLUMPP
#-i: indfile from StructureHarvester output
#-p: popfile from StructureHarvester output
#-o: output Q matrix for distruct input
#-k: K value (number of clusters tested)

###!!! Options to be changed in each analysis manually
#c: number of individuals (change according to STRUCTURE mainparam file)
#r: number of replicate runs
#s: minimum number of population clusters (K) tested
#f: maximum number of population clusters (K) tested
c=21
r=5
s=1
f=5
for i in $(seq $s $f) #input range of K values tested
do
$clumpp/CLUMPP -i K$i.indfile -p K$i.popfile -o K$i.indivq -k $i -c $c -r $r
done

cp $clumpp/paramfile_pop ./
mv paramfile_pop paramfile

for i in $(seq $s $f) #input range of K values tested
do
$clumpp/CLUMPP -i K$i.indfile -p K$i.popfile -o K$i.popq -k $i -c $c -r $r
done

#Key options in the paramfile
# DISTRUCT to visualise the results
###!!!! Options to be changed in each analysis manually
#-M number of populations assigned in the Structure input file
#-N number of individuals
m=21
n=21
#-K K value
#-p input file (population q's)
#-i input file (individual q's)
#-a input file (labels atop figure)
#-b input file (labels below figure)
#-o output file
distruct=/home/sobczm/bin/distruct1.1
cp $distruct/drawparams ./
for i in $(seq $s $f) #input range of K values tested
do
$distruct/distructLinux1.1 -i K$i.indivq -p K$i.popq -a $names -o K$i.ps -k $i -M $m -N $n -K $i
done

#Output is a number of PostScript files showing the average proportion of each
#individual's genome belonging to a given cluster and allowing some ancestry inference based on
#the most likely true number of population clusters as summarised by StructureHarvester.

```


Re-ran with 10M and 1M iterations and burnin respec. as genetics journals might request
```bash
input=/home/sobczm/popgen/other/passey/Structure_4
scripts=/home/sobczm/bin/popgen/snp

#Run replicate STRUCTURE runs, with K from 1 to 5
qsub $scripts/execute_structure.sh $input/Ash_farm_172_pacbio_contigs_unmasked_3_filtered_thinned_1000.recode.struc 1 1 5
qsub $scripts/execute_structure.sh $input/Ash_farm_172_pacbio_contigs_unmasked_3_filtered_thinned_1000.recode.struc 1 2 5
qsub $scripts/execute_structure.sh $input/Ash_farm_172_pacbio_contigs_unmasked_3_filtered_thinned_1000.recode.struc 1 3 5
qsub $scripts/execute_structure.sh $input/Ash_farm_172_pacbio_contigs_unmasked_3_filtered_thinned_1000.recode.struc 1 4 5
qsub $scripts/execute_structure.sh $input/Ash_farm_172_pacbio_contigs_unmasked_3_filtered_thinned_1000.recode.struc 1 5 5
```


# Linkage Disequlibrium

Commands to run analysis of linkage disequilibrium (LD) in Ash Farm V. inaequalis isolates

```bash
input=/home/groups/harrisonlab/project_files/venturia/SNP_calling
scripts=/home/passet/git_repos/scripts/popgen/summary_stats/
vcftools=/home/sobczm/bin/vcftools/bin
```

Calculate D, D' and r^2 for SNPs sparated by 1 and 100kbp in Ash Farm (program calculates the stats using only individuals listed after "--indiv" switch) and plot D' and r2 versus SNP physical distance, histogram of D' values


```bash
$vcftools/vcftools --vcf SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_filtered_thinned_1000.recode.vcf \
--hap-r2 --ld-window-bp-min 1000 --ld-window-bp 100000 --max-missing 1 \
--indv 007 --indv 024 --indv 025 --indv 030 --indv 044 --indv 049 --indv 083 --indv 096 --indv 097 --indv 098 --indv 101 --indv 106 --indv 119 --indv 172 --indv 173 --indv 182 --indv 190 --indv 196 --indv 197 --indv 199 --indv 202
mv out.hap.ld ld.Ash_farm_all

qsub $scripts/sub_plot_ld.sh ld.Ash_farm_all
```

Repeated with the 14 isolates from Bramley or Worcester only
```bash
$vcftools/vcftools --vcf SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_filtered_thinned_1000.recode.vcf \
--hap-r2 --ld-window-bp-min 1000 --ld-window-bp 100000 --max-missing 1 \
--indv 007 --indv 024 --indv 025 --indv 030 --indv 044 --indv 049 --indv 172 --indv 173 --indv 182 --indv 190 --indv 196 --indv 197 --indv 199 --indv 202
mv out.hap.ld ld.Ash_farm_BvW

qsub $scripts/sub_plot_ld.sh ld.Ash_farm_BvW
```

<!--
Repeated with isolates from Bramley and Worcester as above but removed 049 and 199 as they group to the opposite population and 057 due to poor sequencing of isolate (Therefore only 5 Bramley isolates remain and 7 Worcester)

```bash
$vcftools/vcftools --vcf Ash_farm_172_pacbio_contigs_unmasked_2_filtered_thinned_1000.recode.vcf \
--hap-r2 --ld-window-bp-min 1000 --ld-window-bp 100000 --max-missing 1 \
--indv 007 --indv 024 --indv 025 --indv 030 --indv 044 --indv 172 --indv 173 --indv 182 --indv 190 --indv 196 --indv 197 --indv 202
mv out.hap.ld ld.Ash_farm_BvW_minus_rogues

qsub $scripts/summary_stats/sub_plot_ld.sh ld.Ash_farm_BvW_minus_rogues
```
-->

Repeated with (6) Bramley isolates only

```bash
$vcftools/vcftools --vcf SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_filtered_thinned_1000.recode.vcf \
--hap-r2 --ld-window-bp-min 1000 --ld-window-bp 100000 --max-missing 1 \
--indv 007 --indv 024 --indv 025 --indv 030 --indv 044 --indv 199
mv out.hap.ld ld.Ash_farm_Bramley

qsub $scripts/sub_plot_ld.sh ld.Ash_farm_Bramley
```

Repeated with (7) Cox isolates only

```bash
$vcftools/vcftools --vcf SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_2_filtered_thinned_1000.recode.vcf \
--hap-r2 --ld-window-bp-min 1000 --ld-window-bp 100000 --max-missing 1 \
--indv 083 --indv 096 --indv 097 --indv 098 --indv 101 --indv 106 --indv 119
mv out.hap.ld ld.Ash_farm_Cox

qsub $scripts/sub_plot_ld.sh ld.Ash_farm_Cox
```

Repeated with (8) Worcsester isolates only

```bash
$vcftools/vcftools --vcf SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_2_filtered_thinned_1000.recode.vcf \
--hap-r2 --ld-window-bp-min 1000 --ld-window-bp 100000 --max-missing 1 \
--indv 049 --indv 172 --indv 173 --indv 182 --indv 190 --indv 196 --indv 197 --indv 202
mv out.hap.ld ld.Ash_farm_Worcester

qsub $scripts/sub_plot_ld.sh ld.Ash_farm_Worcester
```


Tidy LD work into a sub_directory
```bash
mkdir LD_analysis
mv ld.* LD_analysis/
```

Using R version 3.2.2 installed locally:
```bash
export PATH=/home/armita/prog/R/R-3.2.2/bin:${PATH}
```
And libraries stored in
```bash
export R_LIBS=/home/sobczm/R/x86_64-pc-linux-gnu-library/3.2:$R_LIBS
```

Appended libraries
```bash
R
.libPaths( c( .libPaths(), "/home/sobczm/R/x86_64-pc-linux-gnu-library/3.2") )
.libPaths( c( .libPaths(), "/home/armita/prog/R/R-3.2.2/library") )
```

LD plot (heatmap) for r2 values per contig

```bash
cd /SNP_calling/LD_analysis
scripts=/home/sobczm/bin/popgen/summary_stats

qsub $scripts/sub_ld_plot.sh ld.Ash_farm_all

qsub $scripts/sub_ld_plot.sh ld.Ash_farm_BvW

qsub $scripts/sub_ld_plot.sh ld.Ash_farm_BvW_minus_rogues

qsub $scripts/sub_ld_plot.sh ld.Bramley

qsub $scripts/sub_ld_plot.sh ld.Cox

qsub $scripts/sub_ld_plot.sh ld.Worcester
```

Could not get R libraries to work so Maria ran from her profile. Moved files into Venturia project file

```bash
mkdir SNP_calling/LD_analysis/maria
cp -r /home/sobczm/other/LD_analysis/ SNP_calling/LD_analysis/maria/
```
-->

# Randomness of fixed SNPs

## Need the numbers of total number of SNPs within the whole orchard population on each contig

First removed unwanted lines from filtered vcf file

```bash
sed '1,269d' SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_filtered.vcf > SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_filtered_onlysnps.vcf
```

And from unfiltered vcf file

```bash
sed '1,267d' SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3.vcf > SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_onlysnps.vcf
```

From filtered bw population = 625550

```bash
sed '1,267d' SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_bw_filtered.recode.vcf > SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_bw_filtered_onlysnps.vcf
```

From filtered bc population = 584854

```bash
sed '1,267d' SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_bc_filtered.recode.vcf > SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_bc_filtered_onlysnps.vcf
```

From filtered cw population = 605764

```bash
sed '1,267d' SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_cw_filtered.recode.vcf > SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_cw_filtered_onlysnps.vcf
```

unfiltered bw population

```bash
sed '1,267d' SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_bw.vcf > SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_bw_onlysnps.vcf
```

And fixed SNPs between b/w = 7168
```bash
sed '1,267d' /home/sobczm/popgen/other/passey/Maria/vcf_files/Ash_farm_172_pacbio_contigs_unmasked_3_bw_filtered_fixed.vcf > SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_bw_filtered_fixed_onlysnps.vcf
```

And fixed SNPs between b/c = 160
```bash
sed '1,267d' SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_bc_filtered_fixed.vcf > SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_bc_filtered_fixed_onlysnps.vcf
```

No fixed SNPs between c/w



```bash
grep -w 'contig_1' SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_filtered_onlysnps.vcf | wc -l
grep -w 'contig_1' SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_onlysnps.vcf | wc -l
grep -w 'contig_1' SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_bw_filtered_onlysnps.vcf | wc -l
grep -w 'contig_1' SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_bw_onlysnps.vcf | wc -l
grep -w 'contig_1' SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_bw_filtered_fixed_onlysnps.vcf | wc -l
```


## Need to re-run those files involoving SNPS in genes and non-synonymous

```bash
sed '1,275d' SNP_calling/2018/vcf_files/Ash_farm_172_pacbio_contigs_unmasked_3_bw_filtered_fixed_gene.vcf > SNP_calling/2018/vcf_files/Ash_farm_172_pacbio_contigs_unmasked_3_bw_filtered_fixed_gene_onlysnps.vcf
```

```bash
sed '1,275d' SNP_calling/2018/vcf_files/Ash_farm_172_pacbio_contigs_unmasked_3_bw_filtered_fixed_nonsyn.vcf > SNP_calling/2018/vcf_files/Ash_farm_172_pacbio_contigs_unmasked_3_bw_filtered_fixed_nonsyn_onlysnps.vcf
```

## Code into binary where a SNP is present compared to another file, i.e. when a fixed SNP occurs at a position code 1 at all other positions code 0, for use in Run test

For all fixed SNPs in all SNPs in Bramley and Worcester isolates
```Bash
ProgDir=/home/armita/git_repos/emr_repos/scripts/venturia_inaequalis
$ProgDir/vcf2randomness_test.py --fixed_SNPs /home/groups/harrisonlab/project_files/venturia/SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_bw_filtered_fixed.vcf --all_SNPs /home/groups/harrisonlab/project_files/venturia/SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_bw_filtered.recode.vcf > SNP_calling/Vi_all_fixed_snps_input.txt
```

For fixed SNPs in genes in all SNPs in Bramley and Worcester isolates
```Bash
ProgDir=/home/armita/git_repos/emr_repos/scripts/venturia_inaequalis
$ProgDir/vcf2randomness_test.py --fixed_SNPs /home/groups/harrisonlab/project_files/venturia/SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_bw_filtered_fixed_gene.vcf --all_SNPs /home/groups/harrisonlab/project_files/venturia/SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_bw_filtered.recode.vcf > SNP_calling/Vi_genes_in_all_snps_input.txt
```

For fixed SNPs in genes in all fixed SNPs in Bramley and Worcester isolates
```Bash
ProgDir=/home/armita/git_repos/emr_repos/scripts/venturia_inaequalis
$ProgDir/vcf2randomness_test.py --fixed_SNPs /home/groups/harrisonlab/project_files/venturia/SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_bw_filtered_fixed_gene.vcf --all_SNPs /home/groups/harrisonlab/project_files/venturia/SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_bw_filtered_fixed.vcf > SNP_calling/Vi_genes_in_fixed_snps_input.txt
```

For fixed nonsynonymous SNPs in all SNPs in Bramley and Worcester isolates
```Bash
ProgDir=/home/armita/git_repos/emr_repos/scripts/venturia_inaequalis
$ProgDir/vcf2randomness_test.py --fixed_SNPs /home/groups/harrisonlab/project_files/venturia/SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_bw_filtered_fixed_nonsyn.vcf --all_SNPs /home/groups/harrisonlab/project_files/venturia/SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_bw_filtered.recode.vcf > SNP_calling/Vi_nonsyn_in_all_snps_input.txt
```

For fixed nonsynonymous SNPs in all fixed SNPs in Bramley and Worcester isolates
```Bash
ProgDir=/home/armita/git_repos/emr_repos/scripts/venturia_inaequalis
$ProgDir/vcf2randomness_test.py --fixed_SNPs /home/groups/harrisonlab/project_files/venturia/SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_bw_filtered_fixed_nonsyn.vcf --all_SNPs /home/groups/harrisonlab/project_files/venturia/SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_bw_filtered_fixed.vcf > SNP_calling/Vi_nonsyn_in_fixed_snps_input.txt
```

For fixed nonsynonymous SNPs in fixed SNPs in genes in Bramley and Worcester isolates
```Bash
ProgDir=/home/armita/git_repos/emr_repos/scripts/venturia_inaequalis
$ProgDir/vcf2randomness_test.py --fixed_SNPs /home/groups/harrisonlab/project_files/venturia/SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_bw_filtered_fixed_nonsyn.vcf --all_SNPs /home/groups/harrisonlab/project_files/venturia/SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_bw_filtered_fixed_gene.vcf > SNP_calling/Vi_nonsyn_in_fixed_snps_gene_input.txt
```
