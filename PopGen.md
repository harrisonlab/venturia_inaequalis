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

