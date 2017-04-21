# SNP comparison of V. inaequalis isolates from different cultivars in same mixed orchard

Work performed in directory:
```
/home/groups/harrisonlab/project_files/venturia
```

##Alignment of trimmed MiSeq reads of isolates to Assembled PacBio genome

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
##Pre SNP calling cleanup

Setting up correct formatting for SNP analysis

```
input=alignment/bowtie/
scripts=home/passet/git_repos/scripts/popgen/snp
```

###Rename input mapping files in each folder by prefixing with the strain ID

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
##SNP calling
Run a SNP calling script developed by Maria

To change in each analysis:

input=/home/groups/harrisonlab/project_files/venturia/alignment/bowtie
reference=/home/groups/harrisonlab/project_files/venturia/repeat_masked/v.inaequalis/172_pacbio/filtered_contigs_repmask/172_pacbio_contigs_unmasked.fa


filename=$(basename "$reference")
output="${filename%.*}.dict"

###Prepare genome reference indexes required by GATK
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

###Start SNP calling with GATK
The submission script required needs to be custom-prepared for each analysis, depending on what samples are being analysed.
See inside the submission script below.

```bash
scripts=/home/passet/git_repos/scripts/venturia_inaequalis
qsub $scripts/sub_SNP_calling_multithreaded.sh
```

##Determining genentic structure


###Removal of Isolates for analysis
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


###Only retain biallelic high-quality SNPS with no missing data for genetic analyses.
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


###In some organisms, may want to thin (subsample) SNPs in high linkage diseqilibrium down to
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


###General VCF stats (remember that vcftools needs to have the PERL library exported)
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


###Calculate the index for percentage of shared SNP alleles between the individuals
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

<!--
Rscript --vanilla $scripts/distance_matrix.R SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_filtered_distance.log
-->
```


###Carry out PCA and plot the results
```bash
scripts=/home/passet/git_repos/scripts/popgen/snp

Rscript --vanilla $scripts/pca.R SNP_calling/172_pacbio_contigs_unmasked_filtered.vcf

Rscript --vanilla $scripts/pca.R SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_filtered.vcf

Rscript --vanilla $scripts/pca.R SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_2_filtered.vcf

<!--
Rscript --vanilla $scripts/pca.R SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_3_filtered.vcf
-->
```


###Calculate an NJ tree based on all the SNPs. Outputs a basic diplay of the tree, plus a Newick file to be used for displaying the tree in FigTree and beautifying it.

```bash
scripts=/home/passet/git_repos/scripts/popgen/snp
$scripts/nj_tree.sh SNP_calling/172_pacbio_contigs_unmasked_filtered.vcf

$scripts/nj_tree.sh SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_filtered.vcf
```
No tree produced above as no missing data allowed


Need to filter SNPs to retain those with no missing data in any individual
```bash
scripts=/home/sobczm/bin/popgen/snp
qsub $scripts/sub_vcf_parser.sh SNP_calling/Ash_farm_172_pacbio_contigs_unmasked_2_filtered.vcf 40 30 10 30 1 Y
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

###AMOVA analysis
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

<!--
Re-ran structure with 21 isolate group using structure_analysis_3.sh; Changed iterations to 1M (from 10M) and burin to 100k (from 1M)
-->

#Linkage Disequlibrium

Commands to run analysis of linkage disequilibrium (LD) in Ash Farm V. inaequalis isolates

```bash
input=/home/groups/harrisonlab/project_files/venturia/SNP_calling
scripts=/home/passet/git_repos/scripts/popgen
vcftools=/home/sobczm/bin/vcftools/bin
```

Calculate D, D' and r^2 for SNPs sparated by 1 and 100kbp in Ash Farm (program calculates the stats using only individuals listed after "--indiv" switch)

```bash
$vcftools/vcftools --vcf Ash_farm_172_pacbio_contigs_unmasked_2_filtered_thinned_1000.recode.vcf \
--hap-r2 --ld-window-bp-min 1000 --ld-window-bp 100000 \
--indv 007 --indv 024 --indv 025 --indv 030 --indv 044 --indv 049 --indv 057 --indv 083 --indv 096 --indv 097 --indv 098 --indv 101 --indv 106 --indv 119 --indv 172 --indv 173 --indv 182 --indv 190 --indv 196 --indv 197 --indv 199 --indv 202
mv out.hap.ld ld.Ash_farm_all

qsub $scripts/summary_stats/sub_plot_ld.sh ld.Ash_farm_all
```

Repeated with isolates from Bramley and Worcester only
```bash
$vcftools/vcftools --vcf Ash_farm_172_pacbio_contigs_unmasked_2_filtered_thinned_1000.recode.vcf \
--hap-r2 --ld-window-bp-min 1000 --ld-window-bp 100000 \
--indv 007 --indv 024 --indv 025 --indv 030 --indv 044 --indv 049 --indv 057 --indv 172 --indv 173 --indv 182 --indv 190 --indv 196 --indv 197 --indv 199 --indv 202
mv out.hap.ld ld.Ash_farm_BvW

qsub $scripts/summary_stats/sub_plot_ld.sh ld.Ash_farm_BvW
```
Repeated with isolates from Bramley and Worcester as above but removed 049 and 199 as they group to the opposite population and 057 due to poor sequencing of isolate (Therefore only 5 Bramley isolates remain and 7 Worcester)
```bash
$vcftools/vcftools --vcf Ash_farm_172_pacbio_contigs_unmasked_2_filtered_thinned_1000.recode.vcf \
--hap-r2 --ld-window-bp-min 1000 --ld-window-bp 100000 \
--indv 007 --indv 024 --indv 025 --indv 030 --indv 044 --indv 172 --indv 173 --indv 182 --indv 190 --indv 196 --indv 197 --indv 202
mv out.hap.ld ld.Ash_farm_BvW_minus_rogues

qsub $scripts/summary_stats/sub_plot_ld.sh ld.Ash_farm_BvW_minus_rogues
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

Plot D' and r2 versus SNP physical distance, histogram of D' values

```bash
cd LD_analysis
scripts=/home/sobczm/bin/popgen/summary_stats
vcftools=/home/sobczm/bin/vcftools/bin
qsub $scripts/sub_plot_ld.sh ld.Ash_farm_all
'''


LD plot (heatmap) for r2 values per contig

```bash
scripts=/home/sobczm/bin/popgen/summary_stats
qsub $scripts/sub_ld_plot.sh LD_analysis/ld.Ash_farm_all
```