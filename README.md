# venturia_inaequalis
command for analysis of venturia genome

==========

Scripts used for the analysis of venturia genomes
Note - all this work was performed in the directory:
/home/groups/harrisonlab/project_files/venturia

The following is a summary of the work presented in this Readme.

The following processes were applied to venturia genomes prior to analysis:
Data qc
Genome assembly
Repeatmasking
Gene prediction
Functional annotation

Analyses performed on these genomes involved BLAST searching for:




#Building of directory structure

Project directory created for venturia work

```bash
ProjectDir=/home/groups/harrisonlab/project_files/venturia
mkdir -p $ProjectDir
cd $ProjectDir
```

Data was copied into this directory

```bash
	mkdir -p raw_dna/paired/v.inaequalis/049/F
mkdir -p raw_dna/paired/v.inaequalis/049/R
mkdir -p raw_dna/paired/v.inaequalis/057/F
mkdir -p raw_dna/paired/v.inaequalis/057/R
mkdir -p raw_dna/paired/v.inaequalis/098/F
mkdir -p raw_dna/paired/v.inaequalis/098/R
mkdir -p raw_dna/paired/v.inaequalis/119/F
mkdir -p raw_dna/paired/v.inaequalis/119/R
mkdir -p raw_dna/paired/v.inaequalis/199/F
mkdir -p raw_dna/paired/v.inaequalis/199/R
mkdir -p raw_dna/paired/v.inaequalis/007/F
mkdir -p raw_dna/paired/v.inaequalis/007/R
mkdir -p raw_dna/paired/v.inaequalis/044/F
mkdir -p raw_dna/paired/v.inaequalis/044/R
mkdir -p raw_dna/paired/v.inaequalis/083/F
mkdir -p raw_dna/paired/v.inaequalis/083/R
mkdir -p raw_dna/paired/v.inaequalis/172/F
mkdir -p raw_dna/paired/v.inaequalis/172/R
mkdir -p raw_dna/paired/v.inaequalis/182/F
mkdir -p raw_dna/paired/v.inaequalis/182/R
mkdir -p raw_dna/paired/v.inaequalis/025/F
mkdir -p raw_dna/paired/v.inaequalis/025/R
mkdir -p raw_dna/paired/v.inaequalis/096/F
mkdir -p raw_dna/paired/v.inaequalis/096/R
mkdir -p raw_dna/paired/v.inaequalis/118/F
mkdir -p raw_dna/paired/v.inaequalis/118/R
mkdir -p raw_dna/paired/v.inaequalis/196/F
mkdir -p raw_dna/paired/v.inaequalis/196/R
mkdir -p raw_dna/paired/v.inaequalis/202/F
mkdir -p raw_dna/paired/v.inaequalis/202/R
mkdir -p raw_dna/paired/v.inaequalis/024/F
mkdir -p raw_dna/paired/v.inaequalis/024/R
mkdir -p raw_dna/paired/v.inaequalis/030/F
mkdir -p raw_dna/paired/v.inaequalis/030/R
mkdir -p raw_dna/paired/v.inaequalis/101/F
mkdir -p raw_dna/paired/v.inaequalis/101/R
mkdir -p raw_dna/paired/v.inaequalis/106/F
mkdir -p raw_dna/paired/v.inaequalis/106/R
mkdir -p raw_dna/paired/v.inaequalis/197/F
mkdir -p raw_dna/paired/v.inaequalis/197/R
mkdir -p raw_dna/paired/v.inaequalis/036/F
mkdir -p raw_dna/paired/v.inaequalis/036/R
mkdir -p raw_dna/paired/v.inaequalis/097/F
mkdir -p raw_dna/paired/v.inaequalis/097/R
mkdir -p raw_dna/paired/v.inaequalis/173/F
mkdir -p raw_dna/paired/v.inaequalis/173/R
mkdir -p raw_dna/paired/v.inaequalis/190/F
mkdir -p raw_dna/paired/v.inaequalis/190/R
mkdir -p raw_dna/paired/v.inaequalis/saturn/F
mkdir -p raw_dna/paired/v.inaequalis/saturn/R

```

#Data qc

programs:
  fastqc
  fastq-mcf
  kmc

Data quality was visualised using fastqc:
```bash
	for RawData in $(ls raw_dna/paired/*/*/*/*.fastq.gz); do
		echo $RawData;
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		qsub $ProgDir/run_fastqc.sh $RawData
	done
```



Trimming was performed on data to trim adapters from
sequences and remove poor quality data. This was done with fastq-mcf




Trimming was performed on all isolates:

```bash
	for StrainPath in $(ls -d raw_dna/paired/*/*); 
		do
		ProgDir=/home/passet/git_repos/tools/seq_tools/rna_qc
		IlluminaAdapters=/home/passet/git_repos/tools/seq_tools/ncbi_adapters.fa
		ReadsF=$(ls $StrainPath/F/*.fastq*)
		ReadsR=$(ls $StrainPath/R/*.fastq*)
		echo $ReadsF
		echo $ReadsR
		qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters DNA
	done
```




Data quality was visualised once again following trimming:
```bash
	for RawData in $(ls qc_dna/paired/*/*/*/*.fq.gz); do
		ProgDir=/home/passet/git_repos/tools/seq_tools/dna_qc
		echo $RawData;
		qsub $ProgDir/run_fastqc.sh $RawData
	done
```


kmer counting was performed using kmc
This allowed estimation of sequencing depth and total genome size


```bash
for TrimPath in $(ls -d raw_dna/paired/*/*); do
ProgDir=/home/passet/git_repos/tools/seq_tools/dna_qc
TrimF=$(ls $TrimPath/F/*.fastq*)
TrimR=$(ls $TrimPath/R/*.fastq*)
echo $TrimF
echo $TrimR
qsub $ProgDir/kmc_kmer_counting.sh $TrimF $TrimR
done
```


mode kmer abundance prior to error correction was reported using the following
commands:

```bash
	for File in $(ls qc_dna/kmc/*/*/*_true_kmer_summary.txt); do
		basename $File;
		cat $File | grep -e 'abundance' -e 'size'
	done
```
Results of kmer counting (numbers in brackets indicate estimated kmer coverage from histogram)
'''
The estimated genome size is:  64158121
024_true_kmer_summary.txt
The mode kmer abundance is:  30
The estimated genome size is:  59324376
025_true_kmer_summary.txt
The mode kmer abundance is:  29
The estimated genome size is:  69351314
030_true_kmer_summary.txt
The mode kmer abundance is:  24
The estimated genome size is:  64520872
036_true_kmer_summary.txt
The mode kmer abundance is:  5 
The estimated genome size is:  35589866
044_true_kmer_summary.txt
The mode kmer abundance is:  5 (30)
The estimated genome size is:  452258529
049_true_kmer_summary.txt
The mode kmer abundance is:  34
The estimated genome size is:  56424510
057_true_kmer_summary.txt
The mode kmer abundance is:  5
The estimated genome size is:  265385077
083_true_kmer_summary.txt
The mode kmer abundance is:  25
The estimated genome size is:  61502122
096_true_kmer_summary.txt
The mode kmer abundance is:  5 (20)
The estimated genome size is:  348646667
097_true_kmer_summary.txt
The mode kmer abundance is:  22
The estimated genome size is:  63868207
098_true_kmer_summary.txt
The mode kmer abundance is:  5 (40)
The estimated genome size is:  495964379
101_true_kmer_summary.txt
The mode kmer abundance is:  33
The estimated genome size is:  57884967
106_true_kmer_summary.txt
The mode kmer abundance is:  20
The estimated genome size is:  63717914
118_true_kmer_summary.txt
The mode kmer abundance is:  7
The estimated genome size is:  262481354
119_true_kmer_summary.txt
The mode kmer abundance is:  33
The estimated genome size is:  58479942
172_true_kmer_summary.txt
The mode kmer abundance is:  24
The estimated genome size is:  63294175
173_true_kmer_summary.txt
The mode kmer abundance is:  5 (35)
The estimated genome size is:  437136955
182_true_kmer_summary.txt
The mode kmer abundance is:  32
The estimated genome size is:  62170328
190_true_kmer_summary.txt
The mode kmer abundance is:  5 (20)
The estimated genome size is:  309277674
196_true_kmer_summary.txt
The mode kmer abundance is:  28
The estimated genome size is:  59732382
197_true_kmer_summary.txt
The mode kmer abundance is:  28
The estimated genome size is:  61229552
199_true_kmer_summary.txt
The mode kmer abundance is:  5 (40)
The estimated genome size is:  520820360
202_true_kmer_summary.txt
The mode kmer abundance is:  32
The estimated genome size is:  58471759
saturn_true_kmer_summary.txt
The mode kmer abundance is:  5 (40)
The estimated genome size is:  471789289
'''


#Assembly

Assembly was performed with:
* Spades

## Spades Assembly
Isolates with mode kmer abundance 20 or above (i.e all except 036, 057 and 118)


```bash
	for StrainPath in $(ls -d qc_dna/paired/*/* | grep -v -w -e "036" -e "057" -e "118"); do
		Jobs=$(qstat | grep 'submit_S' | grep 'qw' | wc -l)
		while [ $Jobs -gt 1 ]; do
		sleep 10
		printf "."
		Jobs=$(qstat | grep 'submit_S' | grep 'qw' | wc -l)
		done
		ProgDir=/home/passet/git_repos/tools/seq_tools/assemblers/spades
		Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
		Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
		F_Read=$(ls $StrainPath/F/*.fq.gz)
		R_Read=$(ls $StrainPath/R/*.fq.gz)
		OutDir=assembly/spades/$Organism/$Strain
		echo $F_Read
		echo $R_Read
		qsub $ProgDir/submit_SPAdes_HiMem.sh $F_Read $R_Read $OutDir correct 15
	done
```
Spades assembly with the 3 isolates with mode kmer below 10

```bash
	for StrainPath in $(ls -d qc_dna/paired/*/* | grep -w -e "036" -e "057" -e "118"); do
		Jobs=$(qstat | grep 'submit_S' | grep 'qw' | wc -l)
		while [ $Jobs -gt 1 ]; do
		sleep 10
		printf "."
		Jobs=$(qstat | grep 'submit_S' | grep 'qw' | wc -l)
		done
		ProgDir=/home/passet/git_repos/tools/seq_tools/assemblers/spades
		Strain=$(echo $StrainPath | rev | cut -f1 -d '/' | rev)
		Organism=$(echo $StrainPath | rev | cut -f2 -d '/' | rev)
		F_Read=$(ls $StrainPath/F/*.fq.gz)
		R_Read=$(ls $StrainPath/R/*.fq.gz)
		OutDir=assembly/spades/$Organism/$Strain
		echo $F_Read
		echo $R_Read
		qsub $ProgDir/submit_SPAdes_HiMem.sh $F_Read $R_Read $OutDir correct 10
	done
```


Quast

```bash
	ProgDir=/home/passet/git_repos/tools/seq_tools/assemblers/assembly_qc/quast
		for Assembly in $(ls assembly/spades/*/*/filtered_contigs/contigs_min_500bp.fasta); do
		Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
		Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
		OutDir=assembly/spades/$Organism/$Strain/filtered_contigs
		qsub $ProgDir/sub_quast.sh $Assembly $OutDir
	done
```



Contigs were renamed in accordance with ncbi recomendations.

```bash
	ProgDir=~/git_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
	touch tmp.csv
		for Assembly in $(ls assembly/spades/*/*/filtered_contigs/contigs_min_500bp.fasta); do
		Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
		Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
		OutDir=assembly/spades/$Organism/$Strain/filtered_contigs
		$ProgDir/remove_contaminants.py --inp $Assembly --out $OutDir/contigs_min_500bp_renamed.fasta --coord_file tmp.csv
	done
	rm tmp.csv
```

# Summary of assemblies

```bash
	for File in $(ls assembly/spades/*/*/filtered_contigs/report.tsv); do 
		Organism=$(echo $File | rev |cut -f4 -d '/' | rev) 
		Strain=$(echo $File | rev |cut -f3 -d '/' | rev)
		echo $Organism > tmp_"$Strain".txt
		echo $Strain >> tmp_"$Strain".txt 
		cat $File | tail -n+2 | cut -f2 >> tmp_"$Strain".txt
	done
	paste tmp*.txt > assembly/spades/assembly_summary.tsv
	rm tmp*.txt
```

# Repeatmasking

Repeat masking was performed and used the following programs:
	Repeatmasker
	Repeatmodeler

The best assemblies were used to perform repeatmasking

```bash
	ProgDir=/home/passet/git_repos/tools/seq_tools/repeat_masking
		for BestAss in $(ls assembly/spades/*/*/*/contigs_min_500bp_renamed.fasta); do
		qsub $ProgDir/rep_modeling.sh $BestAss
		qsub $ProgDir/transposonPSI.sh $BestAss
	done
```

The number of bases masked by transposonPSI and Repeatmasker were summarised
using the following commands:

```bash
	for RepDir in $(ls -d repeat_masked/v.*/*/*); do
		Strain=$(echo $RepDir | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $RepDir | rev | cut -f3 -d '/' | rev)  
		RepMaskGff=$(ls $RepDir/*_contigs_hardmasked.gff)
		TransPSIGff=$(ls $RepDir/*_contigs_unmasked.fa.TPSI.allHits)
		printf "$Organism\t$Strain\n"
		printf "The number of bases masked by RepeatMasker:\t"
		sortBed -i $RepMaskGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
		printf "The number of bases masked by TransposonPSI:\t"
		sortBed -i $TransPSIGff | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
		printf "The total number of masked bases are:\t"
		cat $RepMaskGff $TransPSIGff | sortBed | bedtools merge | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
		]echo
	]done
```
Results of bases masked by repeatmasker
'''
v.inaequalis    007
The number of bases masked by RepeatMasker:     34407509
The number of bases masked by TransposonPSI:    0
The total number of masked bases are:   34407509

v.inaequalis    024
The number of bases masked by RepeatMasker:     33510544
The number of bases masked by TransposonPSI:    0
The total number of masked bases are:   33510544

v.inaequalis    025
The number of bases masked by RepeatMasker:     34645619
The number of bases masked by TransposonPSI:    0
The total number of masked bases are:   34645619

v.inaequalis    030
The number of bases masked by RepeatMasker:     33066245
The number of bases masked by TransposonPSI:    0
The total number of masked bases are:   33066245

v.inaequalis    036
The number of bases masked by RepeatMasker:     71808
The number of bases masked by TransposonPSI:    0
The total number of masked bases are:   71808

v.inaequalis    044
The number of bases masked by RepeatMasker:     33705055
The number of bases masked by TransposonPSI:    0
The total number of masked bases are:   33705055

v.inaequalis    049
The number of bases masked by RepeatMasker:     31089941
The number of bases masked by TransposonPSI:    0
The total number of masked bases are:   31089941

v.inaequalis    057
The number of bases masked by RepeatMasker:     7009434
The number of bases masked by TransposonPSI:    0
The total number of masked bases are:   7009434

v.inaequalis    083
The number of bases masked by RepeatMasker:     30796550
The number of bases masked by TransposonPSI:    0
The total number of masked bases are:   30796550

v.inaequalis    096
The number of bases masked by RepeatMasker:     33861179
The number of bases masked by TransposonPSI:    0
The total number of masked bases are:   33861179

v.inaequalis    097
The number of bases masked by RepeatMasker:     33777696
The number of bases masked by TransposonPSI:    0
The total number of masked bases are:   33777696

v.inaequalis    098
The number of bases masked by RepeatMasker:     32129773
The number of bases masked by TransposonPSI:    0
The total number of masked bases are:   32129773

v.inaequalis    101
The number of bases masked by RepeatMasker:     30249028
The number of bases masked by TransposonPSI:    0
The total number of masked bases are:   30249028

v.inaequalis    106
The number of bases masked by RepeatMasker:     29372196
The number of bases masked by TransposonPSI:    0
The total number of masked bases are:   29372196

v.inaequalis    118
The number of bases masked by RepeatMasker:     9267255
The number of bases masked by TransposonPSI:    0
The total number of masked bases are:   9267255

v.inaequalis    119
The number of bases masked by RepeatMasker:     33887998
The number of bases masked by TransposonPSI:    0
The total number of masked bases are:   33887998

v.inaequalis    172
The number of bases masked by RepeatMasker:     32246642
The number of bases masked by TransposonPSI:    0
The total number of masked bases are:   32246642

v.inaequalis    173
The number of bases masked by RepeatMasker:     29861271
The number of bases masked by TransposonPSI:    0
The total number of masked bases are:   29861271

v.inaequalis    182
The number of bases masked by RepeatMasker:     33169962
The number of bases masked by TransposonPSI:    0
The total number of masked bases are:   33169962

v.inaequalis    190
The number of bases masked by RepeatMasker:     0
The number of bases masked by TransposonPSI:    0
The total number of masked bases are:   0

v.inaequalis    196
The number of bases masked by RepeatMasker:     31586869
The number of bases masked by TransposonPSI:    0
The total number of masked bases are:   31586869

v.inaequalis    197
The number of bases masked by RepeatMasker:     31393780
The number of bases masked by TransposonPSI:    0
The total number of masked bases are:   31393780

v.inaequalis    199
The number of bases masked by RepeatMasker:     33016953
The number of bases masked by TransposonPSI:    0
The total number of masked bases are:   33016953

v.inaequalis    202
The number of bases masked by RepeatMasker:     31903088
The number of bases masked by TransposonPSI:    0
The total number of masked bases are:   31903088

v.inaequalis    saturn
The number of bases masked by RepeatMasker:     31550846
The number of bases masked by TransposonPSI:    0
The total number of masked bases are:   31550846
'''

# RNA-seq data download

Dowloaded Thakur et al RNA-seq data from NCBI 

prefetch -o raw_rna/unpaired/v.inaequalis SRR2164202
fastq-dump -O raw_rna/unpaired/v.inaequalis SRR2164202
prefetch -o raw_rna/unpaired/v.inaequalis SRR2164317
fastq-dump -O raw_rna/unpaired/v.inaequalis SRR2164317
prefetch -o raw_rna/unpaired/v.inaequalis SRR2164320
fastq-dump -O raw_rna/unpaired/v.inaequalis SRR2164320
prefetch -o raw_rna/paired/v.inaequalis SRR2164233
fastq-dump -O raw_rna/paired/v.inaequalis SRR2164233
prefetch -o raw_rna/paired/v.inaequalis SRR2164324
fastq-dump -O raw_rna/paired/v.inaequalis SRR2164324
prefetch -o raw_rna/paired/v.inaequalis SRR2164325
fastq-dump -O raw_rna/paired/v.inaequalis SRR2164325


# Gene Prediction


Gene prediction followed three steps:
	Pre-gene prediction
		- Quality of genome assemblies were assessed using Cegma to see how many core eukaryotic genes can be identified.
	Gene model training
		- Gene models were trained using assembled RNAseq data as part of the Braker1 pipeline
	Gene prediction
		- Gene models were used to predict genes in genomes as part of the the Braker1 pipeline. This used RNAseq data as hints for gene models.

# Pre-gene prediction

Quality of genome assemblies was assessed by looking for the gene space in the assemblies.
```bash
	ProgDir=/home/passet/git_repos/tools/gene_prediction/cegma
		cd /home/groups/harrisonlab/project_files/venturia
		for Genome in $(ls repeat_masked/v.*/*/*/*_contigs_unmasked.fa); do
		echo $Genome;
		qsub $ProgDir/sub_cegma.sh $Genome dna;
	done
```
<!--
Outputs were summarised using the commands:
```bash
	for File in $(ls gene_pred/cegma/F*/FOP1/*_dna_cegma.completeness_report); do
		Strain=$(echo $File | rev | cut -f2 -d '/' | rev);
		Species=$(echo $File | rev | cut -f3 -d '/' | rev);
		printf "$Species\t$Strain\n";
		cat $File | head -n18 | tail -n+4;printf "\n";
	done > gene_pred/cegma/cegma_results_dna_summary.txt

	less gene_pred/cegma/cegma_results_dna_summary.txt
```
for num in $(qstat | cut -f1 -d ‘ '); do 
echo $num; 
qdel $num; 
done

<!--
#Gene prediction

Gene prediction was performed for Fusarium genomes. Two gene prediction
approaches were used:

Gene prediction using Braker1
Prediction of all putative ORFs in the genome using the ORF finder (atg.pl)
approach.


## Gene prediction 1 - Braker1 gene model training and prediction

Gene prediction was performed using Braker1.

First, RNAseq data was aligned to Fusarium genomes.
* qc of RNA seq data is detailed below:

```bash
 for Folder in $(ls -d raw_rna/paired/F.oxysporum_fsp_cepae/*); do
	 FolderName=$(echo $Folder | rev | cut -f1 -d '/' | rev);
	 echo $FolderName;
	 ls $Folder/F;
	 ls $Folder/R;
	done
```
This contained the following data:
```
	55_72hrs_rep1
	sample013_1.combined.fastq.gz
	sample013_2.combined.fastq.gz
	55_72hrs_rep2
	sample014_1.combined.fastq.gz
	sample014_2.combined.fastq.gz
	55_72hrs_rep3
	sample015_1.combined.fastq.gz
	sample015_2.combined.fastq.gz
	control_72hrs_rep1
	sample002_1.combined.fastq.gz
	sample002_2.combined.fastq.gz
	control_72hrs_rep2
	sample004_1.combined.fastq.gz
	sample004_2.combined.fastq.gz
	control_72hrs_rep3
	sample005_1.combined.fastq.gz
	sample005_2.combined.fastq.gz
	FO47_72hrs_rep1
	sample006_1.combined.fastq.gz
	sample006_2.combined.fastq.gz
	FO47_72hrs_rep2
	sample007_1.combined.fastq.gz
	sample007_2.combined.fastq.gz
	FO47_72hrs_rep3
	sample012_1.combined.fastq.gz
	sample012_2.combined.fastq.gz
	Fus2_0hrs_prelim
	1_S1_L001_R1_001.fastq.gz
	1_S1_L001_R2_001.fastq.gz
	Fus2_16hrs_prelim
	3_S2_L001_R1_001.fastq.gz
	3_S2_L001_R2_001.fastq.gz
	Fus2_24hrs_prelim_rep1
	4_S3_L001_R1_001.fastq.gz
	4_S3_L001_R2_001.fastq.gz
	Fus2_24hrs_prelim_rep2
	Fus2_36hrs_prelim
	36hr-root_S10_L001_R1_001.fastq.gz
	36hr-root_S10_L001_R2_001.fastq.gz
	Fus2_48hrs_prelim
	6_S4_L001_R1_001.fastq.gz
	6_S4_L001_R2_001.fastq.gz
	Fus2_4hrs_prelim
	4hr-root_S7_L001_R1_001.fastq.gz
	4hr-root_S7_L001_R2_001.fastq.gz
	Fus2_72hrs_prelim
	7_S5_L001_R1_001.fastq.gz
	7_S5_L001_R2_001.fastq.gz
	Fus2_72hrs_rep1
	sample016_1.combined.fastq.gz
	sample016_2.combined.fastq.gz
	Fus2_72hrs_rep2
	sample018_1.combined.fastq.gz
	sample018_2.combined.fastq.gz
	Fus2_72hrs_rep3
	sample019_1.combined.fastq.gz
	sample019_2.combined.fastq.gz
	Fus2_8hrs_prelim
	8hr-root_S8_L001_R1_001.fastq.gz
	8hr-root_S8_L001_R2_001.fastq.gz
	Fus2_96hrs_prelim
	8_S6_L001_R1_001.fastq.gz
	8_S6_L001_R2_001.fastq.gz
	Fus2_CzapekDox
	6_S2_L001_R1_001_fastqc  6_S2_L001_R1_001.fastq.gz
	6_S2_L001_R2_001_fastqc  6_S2_L001_R2_001.fastq.gz
	Fus2_GlucosePeptone
	7_S3_L001_R1_001_fastqc  7_S3_L001_R1_001.fastq.gz
	7_S3_L001_R2_001_fastqc  7_S3_L001_R2_001.fastq.gz
	Fus2_PDA
	9_S4_L001_R1_001_fastqc  9_S4_L001_R1_001.fastq.gz
	9_S4_L001_R2_001_fastqc  9_S4_L001_R2_001.fastq.gz
	Fus2_PDB
	4_S1_L001_R1_001_fastqc  4_S1_L001_R1_001.fastq.gz
	4_S1_L001_R2_001_fastqc  4_S1_L001_R2_001.fastq.gz
```

Perform qc of RNAseq timecourse data
```bash
	for num in $(qstat | grep 'rna_qc' | cut -f1 -d ' '); do
		Treatment=$(head -n1 rna_qc_fastq-mcf.sh.o$num | cut -f10 -d '/')
		for FilePath in $(raw_rna/paired/F.oxysporum_fsp_cepae/55_72hrs_rep1); do
			echo $FilePath
			FileF=$(ls $FilePath/F/*.gz)
			FileR=$(ls $FilePath/R/*.gz)
			IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
			ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
			qsub $ProgDir/rna_qc_fastq-mcf.sh $FileF $FileR $IlluminaAdapters RNA
			sleep 10
		done
	done
```
```
raw_rna/paired/F.oxysporum_fsp_cepae/55_72hrs_rep1
Your job 6436212 ("rna_qc_fastq-mcf.sh") has been submitted
raw_rna/paired/F.oxysporum_fsp_cepae/55_72hrs_rep2
Your job 6436213 ("rna_qc_fastq-mcf.sh") has been submitted
raw_rna/paired/F.oxysporum_fsp_cepae/55_72hrs_rep3
Your job 6436214 ("rna_qc_fastq-mcf.sh") has been submitted
raw_rna/paired/F.oxysporum_fsp_cepae/control_72hrs_rep2
Your job 6436215 ("rna_qc_fastq-mcf.sh") has been submitted
raw_rna/paired/F.oxysporum_fsp_cepae/control_72hrs_rep3
Your job 6436216 ("rna_qc_fastq-mcf.sh") has been submitted
raw_rna/paired/F.oxysporum_fsp_cepae/FO47_72hrs_rep1
Your job 6436217 ("rna_qc_fastq-mcf.sh") has been submitted
raw_rna/paired/F.oxysporum_fsp_cepae/FO47_72hrs_rep2
Your job 6436218 ("rna_qc_fastq-mcf.sh") has been submitted
raw_rna/paired/F.oxysporum_fsp_cepae/FO47_72hrs_rep3
Your job 6436219 ("rna_qc_fastq-mcf.sh") has been submitted
```

Data quality was visualised using fastqc:
```bash
	for RawData in $(ls qc_rna/paired/F.oxysporum_fsp_cepae/*/*/*.fq.gz); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
		echo $RawData;
		qsub $ProgDir/run_fastqc.sh $RawData
	done
```

#### Aligning

Insert sizes of the RNA seq library were unknown until a draft alignment could
be made. To do this tophat and cufflinks were run, aligning the reads against a
single genome. The fragment length and stdev were printed to stdout while
cufflinks was running.

```bash
	for Assembly in $(ls repeat_masked/*/Fo47/*/*_contigs_unmasked.fa); do
		Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
		Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
		echo "$Organism - $Strain"
		for RNADir in $(ls -d qc_rna/paired/F.oxysporum_fsp_cepae/*); do
			Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
			echo "$Timepoint"
			FileF=$(ls $RNADir/F/*_trim.fq.gz)
			FileR=$(ls $RNADir/R/*_trim.fq.gz)
			OutDir=alignment/$Organism/$Strain/$Timepoint
			ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
			qsub $ProgDir/tophat_alignment.sh $Assembly $FileF $FileR $OutDir
		done
	done
```
Alignments were concatenated prior to running cufflinks:
Cufflinks was run to produce the fragment length and stdev statistics:

```bash
	for Assembly in $(ls repeat_masked/*/Fo47/*/*_contigs_softmasked.fa); do
		Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
		Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
		AcceptedHits=alignment/$Organism/$Strain/concatenated/concatenated.bam
		OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated_prelim
		echo "$Organism - $Strain"
		mkdir -p $OutDir
		samtools merge -f $AcceptedHits \
		alignment/$Organism/$Strain/55_72hrs_rep1/accepted_hits.bam \
		alignment/$Organism/$Strain/55_72hrs_rep2/accepted_hits.bam \
		alignment/$Organism/$Strain/55_72hrs_rep3/accepted_hits.bam \
		alignment/$Organism/$Strain/FO47_72hrs_rep1/accepted_hits.bam \
		alignment/$Organism/$Strain/FO47_72hrs_rep2/accepted_hits.bam \
		alignment/$Organism/$Strain/FO47_72hrs_rep3/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_0hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_16hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_24hrs_prelim_rep1/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_36hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_48hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_4hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_72hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_72hrs_rep1/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_72hrs_rep2/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_72hrs_rep3/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_8hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_96hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_CzapekDox/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_GlucosePeptone/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_PDA/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_PDB/accepted_hits.bam
		cufflinks -o $OutDir/cufflinks -p 8 --max-intron-length 4000 $AcceptedHits 2>&1 | tee $OutDir/cufflinks/cufflinks.log
	done
```

Output from stdout included:
```
	Processed 22484 loci.                        [*************************] 100%
	Map Properties:
	     Normalized Map Mass: 50507412.55
	     Raw Map Mass: 50507412.55
	     Fragment Length Distribution: Empirical (learned)
	                   Estimated Mean: 181.98
	                Estimated Std Dev: 78.39
	[13:02:48] Assembling transcripts and estimating abundances.
	Processed 22506 loci.                        [*************************] 100%
```

The Estimated Mean: 181.98 allowed calculation of of the mean insert gap to be
-20bp 182-(97*2) where 97 was the mean read length. This was provided to tophat
on a second run (as the -r option) along with the fragment length stdev to
increase the accuracy of mapping.


Then Rnaseq data was aligned to each genome assembly:

```bash
	for Assembly in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa | grep -v 'FOP1'); do
	# for Assembly in $(ls assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_supercontigs.fasta); do
		Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
		Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
		echo "$Organism - $Strain"
		for RNADir in $(ls -d qc_rna/paired/F.oxysporum_fsp_cepae/* | grep -e '_rep'); do
			Timepoint=$(echo $RNADir | rev | cut -f1 -d '/' | rev)
			echo "$Timepoint"
			FileF=$(ls $RNADir/F/*_trim.fq.gz)
			FileR=$(ls $RNADir/R/*_trim.fq.gz)
			OutDir=alignment/$Organism/$Strain/$Timepoint
			InsertGap='-20'
			InsertStdDev='78'
			Jobs=$(qstat | grep 'tophat' | grep 'qw' | wc -l)
			while [ $Jobs -gt 1 ]; do
				sleep 10
				printf "."
				Jobs=$(qstat | grep 'tophat' | grep 'qw' | wc -l)
			done
			printf "\n"
			ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
			qsub $ProgDir/tophat_alignment.sh $Assembly $FileF $FileR $OutDir $InsertGap $InsertStdDev
		done
	done
```


#### Braker prediction

```bash
	for Assembly in $(ls repeat_masked/*/*/*/*_contigs_softmasked.fa | grep -e 'FOP1'); do
		Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
		Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
		echo "$Organism - $Strain"
		mkdir -p alignment/$Organism/$Strain/concatenated
		samtools merge -f alignment/$Organism/$Strain/concatenated/concatenated.bam \
			alignment/$Organism/$Strain/55_72hrs_rep1/accepted_hits.bam \
			alignment/$Organism/$Strain/55_72hrs_rep2/accepted_hits.bam \
			alignment/$Organism/$Strain/55_72hrs_rep3/accepted_hits.bam \
			alignment/$Organism/$Strain/FO47_72hrs_rep1/accepted_hits.bam \
			alignment/$Organism/$Strain/FO47_72hrs_rep2/accepted_hits.bam \
			alignment/$Organism/$Strain/FO47_72hrs_rep3/accepted_hits.bam \
			alignment/$Organism/$Strain/Fus2_0hrs_prelim/accepted_hits.bam \
			alignment/$Organism/$Strain/Fus2_16hrs_prelim/accepted_hits.bam \
			alignment/$Organism/$Strain/Fus2_24hrs_prelim_rep1/accepted_hits.bam \
			alignment/$Organism/$Strain/Fus2_36hrs_prelim/accepted_hits.bam \
			alignment/$Organism/$Strain/Fus2_48hrs_prelim/accepted_hits.bam \
			alignment/$Organism/$Strain/Fus2_4hrs_prelim/accepted_hits.bam \
			alignment/$Organism/$Strain/Fus2_72hrs_prelim/accepted_hits.bam \
			alignment/$Organism/$Strain/Fus2_72hrs_rep1/accepted_hits.bam \
			alignment/$Organism/$Strain/Fus2_72hrs_rep2/accepted_hits.bam \
			alignment/$Organism/$Strain/Fus2_72hrs_rep3/accepted_hits.bam \
			alignment/$Organism/$Strain/Fus2_8hrs_prelim/accepted_hits.bam \
			alignment/$Organism/$Strain/Fus2_96hrs_prelim/accepted_hits.bam \
			alignment/$Organism/$Strain/Fus2_CzapekDox/accepted_hits.bam \
			alignment/$Organism/$Strain/Fus2_GlucosePeptone/accepted_hits.bam \
			alignment/$Organism/$Strain/Fus2_PDA/accepted_hits.bam \
			alignment/$Organism/$Strain/Fus2_PDB/accepted_hits.bam
		OutDir=gene_pred/braker/$Organism/"$Strain"_braker_new
		AcceptedHits=alignment/$Organism/$Strain/concatenated/concatenated.bam
		GeneModelName="$Organism"_"$Strain"_braker_new
		rm -r /home/armita/prog/augustus-3.1/config/species/"$Organism"_"$Strain"_braker_new
		ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/braker1
		qsub $ProgDir/sub_braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
	done
```

Fasta and gff files were extracted from Braker1 output.

```bash
	for File in $(ls gene_pred/braker/F.*/*_braker_new/*/augustus.gff); do
		getAnnoFasta.pl $File
		OutDir=$(dirname $File)
		echo "##gff-version 3" > $OutDir/augustus_extracted.gff
		cat $File | grep -v '#' >> $OutDir/augustus_extracted.gff
	done
```

The relationship between gene models and aligned reads was investigated. To do
this aligned reads needed to be sorted and indexed:

Note - IGV was used to view aligned reads against the Fus2 genome on my local
machine.

```bash
	InBam=alignment/F.oxysporum_fsp_cepae/Fus2_edited_v2/concatenated/concatenated.bam
	ViewBam=alignment/F.oxysporum_fsp_cepae/Fus2_edited_v2/concatenated/concatenated_view.bam
	SortBam=alignment/F.oxysporum_fsp_cepae/Fus2_edited_v2/concatenated/concatenated_sorted
	samtools view -b $InBam > $ViewBam
	samtools sort $ViewBam $SortBam
	samtools index $SortBam.bam
```


Cufflinks was run to compare the predicted genes to assembled transcripts:

```bash
	for Assembly in $(ls repeat_masked/*/FOP1/*/*_contigs_softmasked.fa); do
		Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
		Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
		AcceptedHits=alignment/$Organism/$Strain/concatenated/concatenated.bam
		OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated
		echo "$Organism - $Strain"
		mkdir -p $OutDir
		samtools merge -f $AcceptedHits \
		alignment/$Organism/$Strain/55_72hrs_rep1/accepted_hits.bam \
		alignment/$Organism/$Strain/55_72hrs_rep2/accepted_hits.bam \
		alignment/$Organism/$Strain/55_72hrs_rep3/accepted_hits.bam \
		alignment/$Organism/$Strain/FO47_72hrs_rep1/accepted_hits.bam \
		alignment/$Organism/$Strain/FO47_72hrs_rep2/accepted_hits.bam \
		alignment/$Organism/$Strain/FO47_72hrs_rep3/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_0hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_16hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_24hrs_prelim_rep1/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_36hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_48hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_4hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_72hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_72hrs_rep1/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_72hrs_rep2/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_72hrs_rep3/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_8hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_96hrs_prelim/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_CzapekDox/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_GlucosePeptone/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_PDA/accepted_hits.bam \
		alignment/$Organism/$Strain/Fus2_PDB/accepted_hits.bam
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
		qsub $ProgDir/sub_cufflinks.sh $AcceptedHits $OutDir/cuflfinks
	# cufflinks -o $OutDir/cufflinks -p 8 --max-intron-length 4000 $AcceptedHits 2>&1 | tee $OutDir/cufflinks/cufflinks.log
	done
```



## Supplimenting Braker gene models with CodingQuary genes

Additional genes were added to Braker gene predictions, using CodingQuary in
pathogen mode to predict additional regions.

Fistly, aligned RNAseq data was assembled into transcripts using Cufflinks.

Note - cufflinks doesn't always predict direction of a transcript and
therefore features can not be restricted by strand when they are intersected.

```bash
	for Assembly in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa | grep -e 'FOP1'); do
		Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
		Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
		echo "$Organism - $Strain"
		OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated
		mkdir -p $OutDir
		AcceptedHits=alignment/$Organism/$Strain/concatenated/concatenated.bam
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
		qsub $ProgDir/sub_cufflinks.sh $AcceptedHits $OutDir
	done
```

Secondly, genes were predicted using CodingQuary:

```bash
	for Assembly in $(ls repeat_masked/*/*/*/*_contigs_softmasked.fa | grep -e 'FOP1'); do
		Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
		Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
		echo "$Organism - $Strain"
		OutDir=gene_pred/codingquary/$Organism/$Strain
		CufflinksGTF=gene_pred/cufflinks/$Organism/$Strain/concatenated/transcripts.gtf
		ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
		qsub $ProgDir/sub_CodingQuary.sh $Assembly $CufflinksGTF $OutDir
	done
```

Then, additional transcripts were added to Braker gene models, when CodingQuary
genes were predicted in regions of the genome, not containing Braker gene
models:

```bash
	for BrakerGff in $(ls gene_pred/braker/F.*/*_braker_new/*/augustus.gff3 | grep -e 'FOP1'); do
		Strain=$(echo $BrakerGff| rev | cut -d '/' -f3 | rev | sed 's/_braker_new//g')
		Organism=$(echo $BrakerGff | rev | cut -d '/' -f4 | rev)
		echo "$Organism - $Strain"
		# BrakerGff=gene_pred/braker/$Organism/$Strain/F.oxysporum_fsp_cepae_Fus2_braker/augustus_extracted.gff
		Assembly=$(ls repeat_masked/$Organism/$Strain/*/"$Strain"_contigs_softmasked.fa)
		CodingQuaryGff=gene_pred/codingquary/$Organism/$Strain/out/PredictedPass.gff3
		PGNGff=gene_pred/codingquary/$Organism/$Strain/out/PGN_predictedPass.gff3
		AddDir=gene_pred/codingquary/$Organism/$Strain/additional
		FinalDir=gene_pred/codingquary/$Organism/$Strain/final
		AddGenesList=$AddDir/additional_genes.txt
		AddGenesGff=$AddDir/additional_genes.gff
		FinalGff=$AddDir/combined_genes.gff
		mkdir -p $AddDir
		mkdir -p $FinalDir

		bedtools intersect -v -a $CodingQuaryGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' > $AddGenesList
		bedtools intersect -v -a $PGNGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' >> $AddGenesList
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
		$ProgDir/gene_list_to_gff.pl $AddGenesList $CodingQuaryGff CodingQuarry_v2.0 ID CodingQuary > $AddGenesGff
		$ProgDir/gene_list_to_gff.pl $AddGenesList $PGNGff PGNCodingQuarry_v2.0 ID CodingQuary >> $AddGenesGff
		ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
		# GffFile=gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2_edited_v2/additional/additional_genes.gff
		# GffFile=gene_pred/codingquary/F.oxysporum_fsp_cepae/Fus2_edited_v2/out/PredictedPass.gff3

		$ProgDir/add_CodingQuary_features.pl $AddGenesGff $Assembly > $FinalDir/final_genes_CodingQuary.gff3
		$ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_CodingQuary.gff3 $FinalDir/final_genes_CodingQuary
		cp $BrakerGff $FinalDir/final_genes_Braker.gff3
		$ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_Braker.gff3 $FinalDir/final_genes_Braker
		cat $FinalDir/final_genes_Braker.pep.fasta $FinalDir/final_genes_CodingQuary.pep.fasta | sed -r 's/\*/X/g' > $FinalDir/final_genes_combined.pep.fasta
		cat $FinalDir/final_genes_Braker.cdna.fasta $FinalDir/final_genes_CodingQuary.cdna.fasta > $FinalDir/final_genes_combined.cdna.fasta
		cat $FinalDir/final_genes_Braker.gene.fasta $FinalDir/final_genes_CodingQuary.gene.fasta > $FinalDir/final_genes_combined.gene.fasta
		cat $FinalDir/final_genes_Braker.upstream3000.fasta $FinalDir/final_genes_CodingQuary.upstream3000.fasta > $FinalDir/final_genes_combined.upstream3000.fasta

		GffBraker=$FinalDir/final_genes_CodingQuary.gff3
		GffQuary=$FinalDir/final_genes_Braker.gff3
		GffAppended=$FinalDir/final_genes_appended.gff3
		cat $GffBraker $GffQuary > $GffAppended

		# cat $BrakerGff $AddDir/additional_gene_parsed.gff3 | bedtools sort > $FinalGff
	done
```

The final number of genes per isolate was observed using:
```bash
for DirPath in $(ls -d gene_pred/codingquary/F.*/*/final); do
echo $DirPath;
cat $DirPath/final_genes_Braker.pep.fasta | grep '>' | wc -l;
cat $DirPath/final_genes_CodingQuary.pep.fasta | grep '>' | wc -l;
cat $DirPath/final_genes_combined.pep.fasta | grep '>' | wc -l;
echo "";
done
```


## ORF finder

The genome was searched in six reading frames for any start codon and following
translated identification of a start codon translating sequence until a stop
codon was found. This is based upon the atg.pl script used in paper describing
the P. infestans genome. Additional functionality was added to this script by
also printing ORFs in .gff format.


```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
	for Genome in $(ls repeat_masked/F.*/*/*/*_contigs_unmasked.fa); do
		qsub $ProgDir/run_ORF_finder.sh $Genome
	done
```

The Gff files from the the ORF finder are not in true Gff3 format. These were
corrected using the following commands:

```bash
	ProgDir=~/git_repos/emr_repos/tools/seq_tools/feature_annotation
	for ORF_Gff in $(ls gene_pred/ORF_finder/*/*/*_ORF.gff | grep -v '_F_atg_' | grep -v '_R_atg_'); do
		ORF_Gff_mod=$(echo $ORF_Gff | sed 's/_ORF.gff/_ORF_corrected.gff3/g')
		echo ""
		echo "Correcting the following file:"
		echo $ORF_Gff
		echo "Redirecting to:"
		echo $ORF_Gff_mod
		$ProgDir/gff_corrector.pl $ORF_Gff > $ORF_Gff_mod
	done
```

#Functional annotation

## A) Interproscan

Interproscan was used to give gene models functional annotations.
Annotation was run using the commands below:

Note: This is a long-running script. As such, these commands were run using
'screen' to allow jobs to be submitted and monitored in the background.
This allows the session to be disconnected and reconnected over time.

Screen ouput detailing the progress of submission of interporscan jobs
was redirected to a temporary output file named interproscan_submission.log .

```bash
	screen -a
	cd /home/groups/harrisonlab/project_files/fusarium
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
	for Genes in $(ls gene_pred/codingquary/F.*/*/*/final_genes_combined.pep.fasta | grep -e 'FOP1'); do
		echo $Genes
		$ProgDir/sub_interproscan.sh $Genes
	done 2>&1 | tee -a interproscan_submisison.log
```

Following interproscan annotation split files were combined using the following
commands:

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
	for Proteins in $(ls gene_pred/codingquary/F.*/*/*/final_genes_combined.pep.fasta | grep 'FOP1'); do
		Strain=$(echo $Proteins | rev | cut -d '/' -f3 | rev)
		Organism=$(echo $Proteins | rev | cut -d '/' -f4 | rev)
		echo "$Organism - $Strain"
		echo $Strain
		InterProRaw=gene_pred/interproscan/$Organism/$Strain/raw
		$ProgDir/append_interpro.sh $Proteins $InterProRaw
	done
```


## B) SwissProt



```bash
	for Proteome in $(ls gene_pred/codingquary/F.*/*/*/final_genes_combined.pep.fasta | grep -e 'FOP1'); do
		Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
		OutDir=gene_pred/swissprot/$Organism/$Strain
		SwissDbDir=../../uniprot/swissprot
		SwissDbName=uniprot_sprot
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/swissprot
		qsub $ProgDir/sub_swissprot.sh $Proteome $OutDir $SwissDbDir $SwissDbName
	done
```

```bash
	for SwissTable in $(ls gene_pred/swissprot/*/*/swissprot_v2015_10_hits.tbl | grep -e 'FOP1'); do
		# SwissTable=gene_pred/swissprot/Fus2/swissprot_v2015_10_hits.tbl
		Strain=$(echo $SwissTable | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $SwissTable | rev | cut -f3 -d '/' | rev)
		echo "$Organism - $Strain"
		OutTable=gene_pred/swissprot/$Organism/$Strain/swissprot_v2015_tophit_parsed.tbl
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/swissprot
		$ProgDir/swissprot_parser.py --blast_tbl $SwissTable --blast_db_fasta ../../uniprot/swissprot/uniprot_sprot.fasta > $OutTable
	done
```


#Genomic analysis

## Comparison to FoL 4287

BLast searches were used to identify which genes had homologs on which
chromosomes of the Fusarium lycopersici genome.

```bash
	FoLGenomeFa=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum.FO2.31.dna.chromosome.fa
	for Proteome in $(ls gene_pred/codingquary/F.*/*/*/final_genes_combined.pep.fasta | grep -e 'FOP1'); do
	# for Proteome in $(ls assembly/external_group/F.oxysporum/fo47/broad/fusarium_oxysporum_fo47_1_proteins.fasta); do
	# for Proteome in $(ls gene_pred/external_group/F.oxysporum_fsp_lycopersici/4287/Fusox1/Fusox1_GeneCatalog_proteins_20110522_parsed.fa); do
		Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
		echo "$Organism - $Strain"
		OutDir=analysis/blast_homology/$Organism/$Strain
		ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
		qsub $ProgDir/run_blast2csv.sh $Proteome protein $FoLGenomeFa $OutDir
	done
```

Convert top blast hits into gff annotations

```bash
	for BlastHitsCsv in $(ls analysis/blast_homology/F.*/*/4287_chromosomal_final_genes_combined.pep.fasta_hits.csv | grep -e 'FOP1'); do
		Organism=$(echo $BlastHitsCsv | rev | cut -f3 -d '/' | rev)
		Strain=$(echo $BlastHitsCsv | rev | cut -f2 -d '/' | rev)
		echo "$Organism - $Strain"
		HitsGff=$(echo $BlastHitsCsv | sed  's/.csv/.gff/g')
		Column2="$Strain"_gene_homolog
		NumHits=1
		ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
		$ProgDir/blast2gff.pl $Column2 $NumHits $BlastHitsCsv > $HitsGff
	done
```

#### 2.6.3) Intersecting blast hits with genes from FoL

```bash
	for HitsGff in $(ls analysis/blast_homology/F.*/*/4287_chromosomal_final_genes_combined.pep.fasta_hits.gff); do
		Organism=$(echo $HitsGff | rev | cut -f3 -d '/' | rev)
		Strain=$(echo $HitsGff| rev | cut -f2 -d '/' | rev)
		echo "$Organism - $Strain"
		HitsDir=$(dirname $HitsGff)
		FoLGenes=assembly/external_group/F.oxysporum_fsp_lycopersici/4287_chromosomal/ensembl/Fusarium_oxysporum.FO2.31.chr.gff3
		FoLIntersect=$HitsDir/4287_chromosomal_final_genes_combined_intersect.bed
		bedtools intersect -wo -a $HitsGff -b $FoLGenes > $FoLIntersect
	done
```

## Effector genes

Putative pathogenicity and effector related genes were identified within Braker
gene models using a number of approaches:

 * A) From Augustus gene models - Identifying secreted proteins
 * B) From Augustus gene models - Effector identification using EffectorP
 * D) From ORF fragments - Signal peptide & RxLR motif  
 * E) From ORF fragments - Hmm evidence of WY domains  
 * F) From ORF fragments - Hmm evidence of RxLR effectors  


### A) From Augustus gene models - Identifying secreted proteins

Required programs:
 * SignalP-4.1
 * TMHMM

Proteins that were predicted to contain signal peptides were identified using
the following commands:

```bash
	SplitfileDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
	CurPath=$PWD
	for Proteome in $(ls gene_pred/codingquary/F.*/*/*/final_genes_combined.pep.fasta | grep -e 'FOP1'); do
		Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
		SplitDir=gene_pred/final_genes_split/$Organism/$Strain
		mkdir -p $SplitDir
		BaseName="$Organism""_$Strain"_final_preds
		$SplitfileDir/splitfile_500.py --inp_fasta $Proteome --out_dir $SplitDir --out_base $BaseName
		for File in $(ls $SplitDir/*_final_preds_*); do
			Jobs=$(qstat | grep 'pred_sigP' | grep 'qw' | wc -l)
			while [ $Jobs -gt 1 ]; do
				sleep 10
				printf "."
				Jobs=$(qstat | grep 'pred_sigP' | grep 'qw' | wc -l)
			done
			printf "\n"
			echo $File
			qsub $ProgDir/pred_sigP.sh $File signalp-4.1
		done
	done
```


The batch files of predicted secreted proteins needed to be combined into a
single file for each strain. This was done with the following commands:
```bash
	for SplitDir in $(ls -d gene_pred/final_genes_split/*/FOP1); do
		Strain=$(echo $SplitDir | rev |cut -d '/' -f1 | rev)
		Organism=$(echo $SplitDir | rev |cut -d '/' -f2 | rev)
		InStringAA=''
		InStringNeg=''
		InStringTab=''
		InStringTxt=''
		SigpDir=final_genes_signalp-4.1
		for GRP in $(ls -l $SplitDir/*_final_preds_*.fa | rev | cut -d '_' -f1 | rev | sort -n); do  
			InStringAA="$InStringAA gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp.aa";  
			InStringNeg="$InStringNeg gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp_neg.aa";  
			InStringTab="$InStringTab gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp.tab";
			InStringTxt="$InStringTxt gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_final_preds_$GRP""_sp.txt";  
		done
		cat $InStringAA > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_sp.aa
		cat $InStringNeg > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_neg_sp.aa
		tail -n +2 -q $InStringTab > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_sp.tab
		cat $InStringTxt > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_final_sp.txt
	done
```

Some proteins that are incorporated into the cell membrane require secretion.
Therefore proteins with a transmembrane domain are not likely to represent
cytoplasmic or apoplastic effectors.

Proteins containing a transmembrane domain were identified:

```bash
	for Proteome in $(ls gene_pred/codingquary/F.*/*/*/final_genes_combined.pep.fasta | grep -e 'FOP1'); do
		Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/transmembrane_helices
		qsub $ProgDir/submit_TMHMM.sh $Proteome
	done
```


### B) From Augustus gene models - Effector identification using EffectorP

Required programs:
 * EffectorP.py

```bash
	for Proteome in $(ls gene_pred/codingquary/F.*/*/*/final_genes_combined.pep.fasta | grep -e 'FOP1'); do
		Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
		Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
		BaseName="$Organism"_"$Strain"_EffectorP
		OutDir=analysis/effectorP/$Organism/$Strain
		ProgDir=~/git_repos/emr_repos/tools/seq_tools/feature_annotation/fungal_effectors
		qsub $ProgDir/pred_effectorP.sh $Proteome $BaseName $OutDir
	done

```

### C) Identification of MIMP-flanking genes

```bash
	for Genome in $(ls repeat_masked/F.*/*/*/*_contigs_unmasked.fa | grep -e 'FOP1'); do
		Organism=$(echo "$Genome" | rev | cut -d '/' -f4 | rev)
		Strain=$(echo "$Genome" | rev | cut -d '/' -f3 | rev)
		BrakerGff=$(ls gene_pred/codingquary/$Organism/"$Strain"/final/final_genes_CodingQuary.gff3)
		QuaryGff=$(ls gene_pred/codingquary/$Organism/"$Strain"/final/final_genes_Braker.gff3)
		OutDir=analysis/mimps/$Organism/$Strain
		mkdir -p "$OutDir"
		echo "$Organism - $Strain"
		ProgDir="/home/armita/git_repos/emr_repos/tools/pathogen/mimp_finder"
		$ProgDir/mimp_finder.pl $Genome $OutDir/"$Strain"_mimps.fa $OutDir/"$Strain"_mimps.gff > $OutDir/"$Strain"_mimps.log
		$ProgDir/gffexpander.pl +- 2000 $OutDir/"$Strain"_mimps.gff > $OutDir/"$Strain"_mimps_exp.gff
		echo "The number of mimps identified:"
		cat $OutDir/"$Strain"_mimps.fa | grep '>' | wc -l
		bedtools intersect -u -a $BrakerGff -b $OutDir/"$Strain"_mimps_exp.gff > $OutDir/"$Strain"_genes_in_2kb_mimp.gff
		bedtools intersect -u -a $QuaryGff -b $OutDir/"$Strain"_mimps_exp.gff >> $OutDir/"$Strain"_genes_in_2kb_mimp.gff
		echo "The following transcripts intersect mimps:"
		MimpGenesTxt=$OutDir/"$Strain"_genes_in_2kb_mimp.txt
		cat $OutDir/"$Strain"_genes_in_2kb_mimp.gff | grep -w 'mRNA' | cut -f9 | cut -f1 -d';' | cut -f2 -d'=' > $MimpGenesTxt
		cat $MimpGenesTxt | wc -l
		echo ""
	done
```

```bash
	cat gene_pred/final_genes_signalp-4.1/F.oxysporum_fsp_cepae/Fus2_edited_v2/Fus2_edited_v2_final_sp.aa | grep '>' | sed 's/>//g' | cut -f1 -d '-' | sed -r 's/\s//g' > gene_pred/final_genes_signalp-4.1/F.oxysporum_fsp_cepae/Fus2_edited_v2/Fus2_edited_v2_final_sp_headers.txt
	cat analysis/effectorP/F.oxysporum_fsp_cepae/Fus2/F.oxysporum_fsp_cepae_Fus2_EffectorP.txt | grep 'Effector' | cut -f1 | sed 's/|//g' > analysis/effectorP/F.oxysporum_fsp_cepae/Fus2/F.oxysporum_fsp_cepae_Fus2_EffectorP_headers.txt
	cat gene_pred/trans_mem/F.oxysporum_fsp_cepae/Fus2_edited_v2/Fus2_edited_v2_TM_genes_pos.txt | cut -f1 > gene_pred/trans_mem/F.oxysporum_fsp_cepae/Fus2_edited_v2/Fus2_edited_v2_TM_genes_pos_headers.txt
```



# 4. Genomic analysis



## 4.2 Orthology

Orthomcl was used to identify orthologous groups between Fusarium spp. genomes

Genomes were grouped by subspecies and orthology was determined within each
subspecies group. Orthology was also determined between subspecies groups.

| Pathogenic | non-pathogenic | Intermediate |
| ---------- | -------------- | -------------|
| 125        | A28            | 55           |
| A23        | D2             |              |
| Fus2       | PG             |              |



### 4.2.c) Orthology between pathogenic and non-pathogenic isolates

The Commands used to run this analysis are shown in
pathogen/orthology/F.oxysporum_fsp.cepae_pathogen_vs_non-pathogen_orthology.md



## 5. BLAST Searches

## 5.1.A) Identifying SIX genes

Protein sequence of previously characterised SIX genes used to BLAST against
assemlies.

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
	for Assembly in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa | grep 'FOP2'); do
		echo $Assembly
		Query=analysis/blast_homology/Fo_path_genes/Fo_path_genes_CRX.fa
		qsub $ProgDir/blast_pipe.sh $Query dna $Assembly
	done
```

Once blast searches had completed, the BLAST hits were converted to GFF
annotations:

```bash
	for BlastHits in $(ls analysis/blast_homology/*/*/*Fo_path_genes_CRX.fa_homologs.csv); do
		Strain=$(echo $BlastHits | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $BlastHits | rev | cut -f3 -d '/' | rev)
		ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
		HitsGff=$(echo $BlastHits | sed  's/.csv/.gff/g')
		Column2=SIX_homolog
		NumHits=1
		$ProgDir/blast2gff.pl $Column2 $NumHits $BlastHits > $HitsGff
	done
```

	The blast hits were summarised in a single table for all the genomes. The top
	identity of the top blast hit in relation to the enquire query sequence was
	presented for each blast hit.

```bash
	OutFile=analysis/blast_homology/Fo_path_genes_CRX_summary.tab
	cat analysis/blast_homology/F.proliferatum/A8/A8_Fo_path_genes_CRX.fa_homologs.csv | cut -f1 > tmp2.tab
	for BlastHits in $(ls analysis/blast_homology/*/*/*Fo_path_genes_CRX.fa_homologs.csv); do
		Strain=$(echo $BlastHits | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $BlastHits | rev | cut -f3 -d '/' | rev)
		echo "$Organism" > tmp.tab
		echo "$Strain" >> tmp.tab
		cat $BlastHits | cut -f10 >> tmp.tab
		paste tmp2.tab tmp.tab > $OutFile
		cp $OutFile tmp2.tab
	done
	rm tmp.tab
	rm tmp2.tab
```

```bash
	for HitsGff in $(ls analysis/blast_homology/*/*/*Fo_path_genes_CRX.fa_homologs.gff | grep -v 'trinity'); do
		Strain=$(echo $HitsGff | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $HitsGff | rev | cut -f3 -d '/' | rev)
		echo "$Organism - $Strain"
		GffBraker=gene_pred/codingquary/$Organism/$Strain/final/final_genes_Braker.gff3
		GffQuary=gene_pred/codingquary/$Organism/$Strain/final/final_genes_CodingQuary.gff3
		OutDir=$(dirname $HitsGff)
		SixIntersect=$OutDir/"$Strain"_Fo_path_genes_CRX.fa_hit_genes.bed
		bedtools intersect -wo -a $HitsGff -b $GffBraker > $SixIntersect
		bedtools intersect -wo -a $HitsGff -b $GffQuary >> $SixIntersect
		cat $SixIntersect | grep -w 'gene' | cut -f9,18
		echo ""
	done > analysis/blast_homology/Fo_path_genes/Fo_path_genes_CRX_hit_genes_summary.tab
```


## 5.1.B) Identifying FTF genes

```bash
	ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
	for Assembly in $(ls repeat_masked/*/Fus2*/*/*_contigs_unmasked.fa); do
		echo $Assembly
		Query=analysis/blast_homology/Fo_path_genes/FTF_cds_Sanchez_et_al_2016.fasta
		qsub $ProgDir/blast_pipe.sh $Query dna $Assembly
	done
```

```bash
	for BlastHits in $(ls analysis/blast_homology/*/*/*_FTF_cds_Sanchez_et_al_2016.fasta_homologs.csv); do
		Strain=$(echo $BlastHits | rev | cut -f2 -d '/' | rev)
		Organism=$(echo $BlastHits | rev | cut -f3 -d '/' | rev)
		ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
		HitsGff=$(echo $BlastHits | sed  's/.csv/.gff/g')
		Column2=FTF_homolog
		NumHits=1
		$ProgDir/blast2gff.pl $Column2 $NumHits $BlastHits > $HitsGff
	done
```



## 5.2 Identifying PHIbase homologs

The PHIbase database was searched against the assembled genomes using tBLASTx.

```bash
	for Assembly in $(ls repeat_masked/F.*/*/*/*_contigs_unmasked.fa); do
		ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
		qsub $ProgDir/blast_pipe.sh analysis/blast_homology/PHIbase/PHI_36_accessions.fa protein $Assembly
	done
```

following blasting PHIbase to the genome, the hits were filtered by effect on
virulence.

First the a tab seperated file was made in the clusters core directory containing
PHIbase. These commands were run as part of previous projects but have been
included here for completeness.


### Expressed Genes

As a preliminary analysis of the RNAseq data, highly expressed genes at 72hrs
post infection were identified in Fus2.

Note - cufflinks doesn't always predict direction of a transcript and
therefore features can not be restricted by strand when they are intersected.

```bash
	# samtools merge -f alignment/$Organism/$Strain/concatenated/Fus2_72hpi.bam alignment/$Organism/$Strain/Fus2_72hrs_rep1/accepted_hits.bam alignment/$Organism/$Strain/Fus2_72hrs_rep2/accepted_hits.bam alignment/$Organism/$Strain/Fus2_72hrs_rep3/accepted_hits.bam
	for AcceptedHits in $(ls alignment/*/*/concatenated/concatenated.bam | grep -v -e 'Fus2_edited_v2'); do
		Strain=$(echo $AcceptedHits | rev | cut -f3 -d '/' | rev)
		Organism=$(echo $AcceptedHits | rev | cut -f4 -d '/' | rev)
		OutDir=$(dirname $AcceptedHits)
		echo "$Organism - $Strain"
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
		qsub $ProgDir/sub_cufflinks.sh $AcceptedHits $OutDir
	done
```

```bash
	for Transcripts in $(ls alignment/F.*/*/concatenated/transcripts.gtf); do
		Strain=$(echo $Transcripts | rev | cut -f3 -d '/' | rev)
		Organism=$(echo $Transcripts | rev | cut -f4 -d '/' | rev)
		echo "$Organism - $Strain"
		GeneGff=gene_pred/codingquary/$Organism/$Strain/final/final_genes_appended.gff3
		ExpressedGenes=alignment/$Organism/$Strain/concatenated/expressed_genes.bed
		bedtools intersect -wao -a $Transcripts -b $GeneGff > $ExpressedGenes
	done
```


## 6. Summarising the Fusarium Proteome

```bash

```
 -->