# SNP comparison of V. inaequalis isolates from different cultivars in same mixed orchard

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

