# SNP comparison of V. inaequalis isolates from different cultivars in same mixed orchard

##Alignment of MiSeq reads of isolates to Assembled PacBio genome

```bash
	for Assembly in $(ls repeat_masked/v.inaequalis/172_pacbio/filtered_contigs_repmask/172_pacbio_contigs_unmasked.fa); do
		Jobs=$(qstat | grep 'sub_bo' | grep 'qw' | wc -l)
		while [ $Jobs -gt 1 ]; do
		sleep 10
		printf "."
		Jobs=$(qstat | grep 'sub_bo' | grep 'qw' | wc -l)
	done
	Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
	Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
	Paired=$(echo $Assembly | rev | cut -d '/' -f5 | rev)
	echo "$Organism - $Strain"
	for IsolateDir in $(ls raw_dna/paired/v.inaequalis/*); do
		Read_F=$(ls $IsolateDir/F/*.fastq.gz)
		Read_R=$(ls $IsolateDir/R/*.fastq.gz)
		OutDir=alignment/$Paired/$Organism/$Strain/
		ProgDir=/home/passet/git_repos/tools/seq_tools/genome_alignment/bowtie/
		qsub $ProgDir/sub_bowtie.sh $Assembly $Read_F $Read_R $OutDir
		done
	done
```