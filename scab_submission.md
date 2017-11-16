
#Submission of V.inaequalis genomes to ncbi 


A Bioproject and biosample number was prepared for the genome submission at:
https://submit.ncbi.nlm.nih.gov


## Submission of sequence data to SRA

Reads were submitted to the SRA at https://submit.ncbi.nlm.nih.gov/subs/sra/ .
To do this, a metadata file was provided detailing each of the files in the
bioproject. The file was downloaded in excel format and edited manually. A copy
of the edited file and the final .tsv file is present at:

```bash
  ls genome_submission/PRJNA354841_SRA_metadata_acc.txt genome_submission/PRJNA354841_SRA_metadata_acc.xlsx
  ```

As these files included a file > 500 Mb, a presubmission folder was requested.
This aids submission of large data files. This file was created on the ftp server
at ftp-private.ncbi.nlm.nih.gov, with a private folder named
uploads/tom.passey@emr.ac.uk_GHO2Umdl Ncbi provided a username a password.
Files were uploaded into a folder created within my preload folder using ftp.

For genomic reads:
```bash
	# Bioproject="PRJNA354841"
	SubFolder="Vi_sra_PRJNA354841"
	mkdir $SubFolder
	for Read in $(ls raw_dna/paired/*/*/*/*.fastq.gz); do
		echo $Read;
		cp $Read $SubFolder/.
	done
	cp raw_dna/pacbio/*/*/extracted/concatenated_pacbio.fastq $SubFolder/.
	cd $SubFolder
	gzip concatenated_pacbio.fastq
	ftp ftp-private.ncbi.nlm.nih.gov
	cd uploads/tom.passey@emr.ac.uk_GHO2Umdl
	mkdir Vi_sra_PRJNA354841
	cd Vi_sra_PRJNA354841
	put Vi_sra_PRJNA354841
	prompt
	mput *
	bye
	cd ../
	#rm -r $SubFolder
```

# Submission of merged assembled genome of isolate 172 

Genome coverage required so ran following (Genome size based on assembled genome)

```bash
 count_nucl.pl -i raw_dna/pacbio/v.inaequalis/172_pacbio/extracted/concatenated_pacbio.fastq  -g 72
```

Following output:

```
The estimated genome size is: 72000000 bp

The input file is: raw_dna/pacbio/v.inaequalis/172_pacbio/extracted/concatenated_pacbio.fastq

Results for: raw_dna/pacbio/v.inaequalis/172_pacbio/extracted/concatenated_pacbio.fastq
 Within this file of 6996972081 bp there were 944907 fastq sequences
 of these 0 lines were empty.


Total results:

 There are a total of 6996972081 nucleotides in this file.

 This equates to an estimated genome coverage of 97.18 .
```

The following note was provided in the WGS submission page on NCBI in the box
labeled "Private comments to NCBI staff":

```
I have been advised to submit my assemblies to NCBI early in my submission process to ensure that my contigs pass the contamination screen. This assembly will be revised as appropriate, including renaming of contigs where needed. Please allow me to modify this submission at a later date, including upload of the final gene models.

'For future submissions, you could send us the fasta files early
in the submission process so we can run them through our foreign
contamination screen. We will let you know if we find any
sequences to exclude or trim before you generate your final
WGS submission.'...'*IMPORTANT* Include a comment that you are submitting
the fasta files to be screened by the contamination screen
prior to creating your final annotated submission.'
```

Uploaded file:

```
home/groups/harrisonlab/project_files/venturia/assembly/merged_canu_spades/v.inaequalis/172_pacbio/filtered_contigs/contigs_min_500bp_renamed.fasta
```

#NCBI response

NCBI responded that the file looked good so no edits required

# SbtFile

The genbank submission template tool was used at: http://www.ncbi.nlm.nih.gov/WebSub/template.cgi This produce a template file detailing the submission.

# Final Submission

These commands were used in the final submission of the 172_pacbio genome:


```bash
for Assembly in $(ls repeat_masked/v.*/*/filtered_contigs_repmask/*_contigs_unmasked.fa | grep -w "172_pacbio"); do
# tbl2asn options:
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
#
ProjDir=/home/groups/harrisonlab/project_files/venturia
cd $ProjDir
OutDir="genome_submission/$Organism/$Strain"
mkdir -p $OutDir

# Program locations:
AnnieDir="/home/armita/prog/annie/genomeannotation-annie-c1e848b"
ProgDir="/home/passet/git_repos/tools/genbank_submission"
# File locations:
# Assembly=$(ls repeat_masked/$Organism/$Strain/filtered_contigs_repmask/*_contigs_unmasked.fa)
InterProTab=$(ls gene_pred/interproscan/$Organism/$Strain/"$Strain"_interproscan.tsv)
SwissProtBlast=$(ls gene_pred/swissprot/$Organism/$Strain/swissprot_v2015_tophit_parsed.tbl)
SwissProtFasta=$(ls /home/groups/harrisonlab/uniprot/swissprot/uniprot_sprot.fasta)
GffFile=$(ls gene_pred/codingquary/v.inaequalis/172_pacbio/final/final_genes_appended.gff3)
SbtFile=genome_submission/v.inaequalis/172_pacbio/template.sbt

#SRA_metadata=$(ls genome_submission/PRJNA354841_SRA_metadata_acc.txt)
#BioProject=$(cat $SRA_metadata | sed 's/PRJNA/\nPRJNA/g' | grep "$StrainOfficial" | cut -f1 | head -n1)
#BioSample=$(cat $SRA_metadata | sed 's/PRJNA/\nPRJNA/g' | grep "$StrainOfficial" | cut -f2 | head -n1)


# ncbi_tbl_corrector script options:
SubmissionID="SUB2310658"
LocusTag="Vi05172"
LabID="harrisonlab"
# Final submisison file name:
FinalName="$Organism"_"$Strain"_Passey_2017

python3 $AnnieDir/annie.py -ipr $InterProTab -g $GffFile -b $SwissProtBlast -db $SwissProtFasta -o $OutDir/annie_output.csv --fix_bad_products
$ProgDir/edit_tbl_file/annie_corrector.py --inp_csv $OutDir/annie_output.csv --out_csv $OutDir/annie_corrected_output.csv

mkdir -p $OutDir/gag/round1
gag.py -f $Assembly -g $GffFile -a $OutDir/annie_corrected_output.csv --fix_start_stop -o $OutDir/gag/round1 2>&1 | tee $OutDir/gag_log1.txt
sed -i 's/Dbxref/db_xref/g' $OutDir/gag/round1/genome.tbl

cp $Assembly $OutDir/gag/round1/genome.fsa  
cp $SbtFile $OutDir/gag/round1/genome.sbt
mkdir -p $OutDir/tbl2asn/round1
tbl2asn -p $OutDir/gag/round1/. -t $OutDir/gag/round1/genome.sbt -r $OutDir/tbl2asn/round1 -M n -X E -Z $OutDir/gag/round1/discrep.txt -j "[organism=$OrganismOfficial] [strain=$StrainOfficial]"

mkdir -p $OutDir/gag/edited
$ProgDir/edit_tbl_file/ncbi_tbl_corrector.py --inp_tbl $OutDir/gag/round1/genome.tbl --inp_val $OutDir/tbl2asn/round1/genome.val --locus_tag $LocusTag --lab_id $LabID --gene_id "remove" --edits stop pseudo unknown_UTR correct_partial --rename_genes "g" --remove_product_locus_tags "True" --out_tbl $OutDir/gag/edited/genome.tbl

printf "StructuredCommentPrefix\t##Genome-Annotation-Data-START##
Annotation Provider\tHarrison Lab NIAB-EMR
Annotation Date\tSEP-2016
Annotation Version\tRelease 1.01
Annotation Method\tAb initio gene prediction: Braker 1.9 and CodingQuary 2.0; Functional annotation: Swissprot (2015 release) and Interproscan 5.18-57.0" \
> $OutDir/gag/edited/annotation_methods.strcmt.txt

sed -i 's/_pilon//g' $OutDir/gag/edited/genome.tbl
sed -i 's/\. subunit/kDa subunit/g' $OutDir/gag/edited/genome.tbl
sed -i 's/, mitochondrial//g' $OutDir/gag/edited/genome.tbl

cp $Assembly $OutDir/gag/edited/genome.fsa
cp $SbtFile $OutDir/gag/edited/genome.sbt
mkdir $OutDir/tbl2asn/final
tbl2asn -p $OutDir/gag/edited/. -t $OutDir/gag/edited/genome.sbt -r $OutDir/tbl2asn/final -M n -X E -Z $OutDir/tbl2asn/final/discrep.txt -j "[organism=$OrganismOfficial] [strain=$StrainOfficial]" -l paired-ends -a r10k -w $OutDir/gag/edited/annotation_methods.strcmt.txt
cat $OutDir/tbl2asn/final/genome.sqn > $OutDir/tbl2asn/final/$FinalName.sqn
done
```

```bash
	# Bioproject="PRJNA354841"
	SubFolder="Vi_annotated_PRJNA354841"
	mkdir $SubFolder
	for Read in $(ls genome_submission/v.inaequalis/172_pacbio/tbl2asn/final/v.inaequalis_172_pacbio_Passey_2017.sqn); do
		echo $Read;
		cp $Read $SubFolder/.
	done
	cp genome_submission/v.inaequalis/172_pacbio/tbl2asn/final/v.inaequalis_172_pacbio_Passey_2017.sqn $SubFolder/.
	cd $SubFolder
	gzip v.inaequalis_172_pacbio_Passey_2017.sqn
	ftp ftp-private.ncbi.nlm.nih.gov
	cd uploads/tom.passey@emr.ac.uk_GHO2Umdl
	mkdir Vi_annotated_PRJNA354841
	cd Vi_annotated_PRJNA354841
	put v.inaequalis_172_pacbio_Passey_2017.sqn.gz
	mput *
	bye
	cd ../
	#rm -r $SubFolder
```

NCBI reponse:
```
>
>[1] 2 genes completely overlapped by other genes. The corresponding CDS 
>spans overlap another CDS span which is very unusual. Looks like these 
>features should be annotated as alternatively spliced genes as 
>described at 
>https://www.ncbi.nlm.nih.gov/genbank/eukaryotic_genome_submission_annot
>ati
>on/#Alternativelysplicedgenes
>If not, some of these features should be deleted.
>
>Gene   Vi05172_g2590   lcl|contig_4:c>169213-<168934   Vi05172_g2590
>Gene   Vi05172_g12556  lcl|contig_60:c>210937-<210550  Vi05172_g12556
>
>
>[2] Some product names should be improved (see file posted on the portal).
>Some product names contain database identifier more appropriate in note.
>Remove them from the product names. Here are some examples:
>C584.13, C417.12, YGL114W, YPL109C, UM03490, AO090003000058, ...
>Please remove any identifiers with similar format. Note that this list 
>is not exhaustive. There might be additional identifiers formats in 
>your product names
```

# Final sumission

The following script, based on that used to submit Alternaria genome, was used to re-submit 172_pacbio genome after removal of duplicates:

## Setting varibales
Vairables containing locations of files and options for scripts were set:

```bash
# Program locations:
AnnieDir="/home/armita/prog/annie/genomeannotation-annie-c1e848b"
ProgDir="/home/passet/git_repos/tools/genbank_submission"
# File locations:
SbtFile="genome_submission/template.sbt"
LabID="harrisonlab"
```

## Generating .tbl file (GAG)

The Genome Annotation Generator (GAG.py) can be used to convert gff files into
.tbl format, for use by tbl2asn.

It can also add annotations to features as provided by Annie the Annotation
extractor.

### Extracting annotations (Annie)

Interproscan and Swissprot annotations were extracted using annie, the
ANNotation Information Extractor. The output of Annie was filtered to
keep only annotations with references to ncbi approved databases.
Note - It is important that transcripts have been re-labelled as mRNA by this
point.

```bash
for Assembly in $(ls repeat_masked/v.*/*/filtered_contigs_repmask/*_contigs_unmasked.fa | grep -w "172_pacbio"); do
  Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
  Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
  echo "$Organism - $Strain"
  OutDir="genome_submission/$Organism/$Strain"
  GffFile=$(ls gene_pred/final/$Organism/"$Strain"*/final/final_genes_appended_renamed.gff3)

  InterProTab=$(ls gene_pred/interproscan/$Organism/"$Strain"*/"$Strain"*_interproscan.tsv)
  SwissProtBlast=$(ls gene_pred/swissprot/$Organism/"$Strain"*/swissprot_vJul2016_tophit_parsed.tbl)
  SwissProtFasta=$(ls /home/groups/harrisonlab/uniprot/swissprot/uniprot_sprot.fasta)
  python3 $AnnieDir/annie.py -ipr $InterProTab -g $GffFile -b $SwissProtBlast -db $SwissProtFasta -o $OutDir/annie_output.csv --fix_bad_products
  ProgDir=/home/passet/git_repos/tools/genbank_submission
  $ProgDir/edit_tbl_file/annie_corrector.py --inp_csv $OutDir/annie_output.csv --out_csv $OutDir/annie_corrected_output.csv
done
```

### Running GAG

Gag was run using the modified gff file as well as the annie annotation file.
Gag was noted to output database references incorrectly, so these were modified.

```bash
for Assembly in $(ls repeat_masked/v.*/*/filtered_contigs_repmask/*_contigs_unmasked.fa | grep -w "172_pacbio"); do
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
echo "$Organism - $Strain"
OutDir="genome_submission/$Organism/$Strain"
GffFile=$(ls gene_pred/final/$Organism/"$Strain"*/final/final_genes_appended_renamed.gff3)
mkdir -p $OutDir/gag/round1
gag.py -f $Assembly -g $GffFile -a $OutDir/annie_corrected_output.csv --fix_start_stop -o $OutDir/gag/round1 2>&1 | tee $OutDir/gag_log1.txt
sed -i 's/Dbxref/db_xref/g' $OutDir/gag/round1/genome.tbl
done
```

<!-- ## manual edits

The gene NS_04463 was found to use the same start and stop codon as predicted
gene CUFF_4598_1_205. Both of these genes were predicted by codingquary. Neither
of these genes were predicted as having alternative splicing. As such the gene
NS_04463 was removed. The same was found for genes CUFF_11067_2_85 and
CUFF_11065_1_82 and as a result CUFF_11067_2_85 was removed.

```bash
  nano $OutDir/gag/round1/genome.tbl
``` -->

## tbl2asn round 1

tbl2asn was run an initial time to collect error reports on the current
formatting of the .tbl file.
Note - all input files for tbl2asn need to be in the same directory and have the
same basename.

```bash
	for Assembly in $(ls repeat_masked/v.*/*/filtered_contigs_repmask/*_contigs_unmasked.fa | grep -w "172_pacbio"); do
		Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
		Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
		echo "$Organism - $Strain"
		OutDir="genome_submission/$Organism/$Strain"

		cp $Assembly $OutDir/gag/round1/genome.fsa
		SbtFile=$(ls genome_submission/v.*/*/template.sbt)
		cp $SbtFile $OutDir/gag/round1/genome.sbt
		mkdir -p $OutDir/tbl2asn/round1
		tbl2asn -p $OutDir/gag/round1/. -t $OutDir/gag/round1/genome.sbt -r $OutDir/tbl2asn/round1 -M n -X E -Z $OutDir/gag/round1/discrep.txt -j "[organism=$Organism] [strain=$Strain]"
	done
```

## Editing .tbl file

The tbl2asn .val output files were observed and errors corrected. This was done
with an in house script. The .val file indicated that some cds had premature
stops, so these were marked as pseudogenes ('pseudo' - SEQ_FEAT.InternalStop)
and that some genes had cds coordinates that did not match the end of the gene
if the protein was hanging off a contig ('stop' - SEQ_FEAT.NoStop).
Furthermore a number of other edits were made to bring the .tbl file in line
with ncbi guidelines. This included: Marking the source of gene
predictions and annotations ('add_inference'); Correcting locus_tags to use the
given ncbi_id ('locus_tag'); Correcting the protein and transcript_ids to
include the locus_tag and reference to submitter/lab id ('lab_id'), removal of
annotated names of genes if you don't have high confidence in their validity
(--gene_id 'remove'). If 5'-UTR and 3'-UTR were not predicted during gene
annotation then genes, mRNA and exon features need to reflect this by marking
them as incomplete ('unknown_UTR').

```bash
	for Assembly in $(ls repeat_masked/v.*/*/filtered_contigs_repmask/*_contigs_unmasked.fa | grep -w "172_pacbio"); do
		Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
		Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
		echo "$Organism - $Strain"
		OutDir="genome_submission/$Organism/$Strain"
		SubmissionID="SUB2310658"
		LocusTag="Vi05172"
		LabID="harrisonlab"
		echo $SubmissionID
		mkdir -p $OutDir/gag/edited
		ProgDir=/home/passet/git_repos/tools/genbank_submission
		$ProgDir/edit_tbl_file/ncbi_tbl_corrector.py --inp_tbl $OutDir/gag/round1/genome.tbl --inp_val $OutDir/tbl2asn/round1/genome.val --locus_tag $SubmissionID --lab_id $LabID --gene_id "remove" --edits stop pseudo unknown_UTR correct_partial --remove_product_locus_tags "True" --del_name_from_prod "True" --out_tbl $OutDir/gag/edited/genome.tbl
	done
```


## Generating a structured comment detailing annotation methods

```bash
  for Assembly in $(ls repeat_masked/v.*/*/filtered_contigs_repmask/*_contigs_unmasked.fa | grep -w "172_pacbio"); do
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    echo "$Organism - $Strain"
    OutDir="genome_submission/$Organism/$Strain"
    printf "StructuredCommentPrefix\t##Genome-Annotation-Data-START##
    Annotation Provider\tHarrison Lab NIAB-EMR
    Annotation Date\tNOV-2017
    Annotation Version\tRelease 1.01
    Annotation Method\tAb initio gene prediction: Braker 1.9 and CodingQuary 2.0; Functional annotation: Swissprot (July 2016 release) and Interproscan 5.18-57.0" \
    > $OutDir/gag/edited/annotation_methods.strcmt.txt
  done
```

## Final run of tbl2asn

Following correction of the GAG .tbl file, tbl2asn was re-run to provide the
final genbank submission file.

The options -l paired-ends -a r10k inform how to handle runs of Ns in the
sequence, these options show that paired-ends have been used to estimate gaps
and that runs of N's longer than 10 bp should be labelled as gaps.

```bash
	for Assembly in $(ls repeat_masked/v.*/*/filtered_contigs_repmask/*_contigs_unmasked.fa | grep -w "172_pacbio"); do
		Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
		Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
		echo "$Organism - $Strain"
		OutDir="genome_submission/$Organism/$Strain"
		FinalName="$Organism"_"$Strain"_Nov_2017_Passey
		cp $Assembly $OutDir/gag/edited/genome.fsa
		cp $SbtFile $OutDir/gag/edited/genome.sbt
		mkdir $OutDir/tbl2asn/final
		tbl2asn -p $OutDir/gag/edited/. -t $OutDir/gag/edited/genome.sbt -r $OutDir/tbl2asn/final -M n -X E -Z $OutDir/tbl2asn/final/discrep.txt -j "[organism=$Organism] [strain=$Strain]" -l paired-ends -a r10k -w $OutDir/gag/edited/annotation_methods.strcmt.txt
		cat $OutDir/tbl2asn/final/genome.sqn | sed 's/_pilon//g' | sed 's/title "Saccharopine dehydrogenase.*/title "Saccharopine dehydrogenase/g' | sed 's/"Saccharopine dehydrogenase.*"/"Saccharopine dehydrogenase"/g' > $OutDir/tbl2asn/final/$FinalName.sqn
	done
```

```bash
# Bioproject="PRJNA354841"
SubFolder="Vi_annotated_PRJNA354841"
mkdir $SubFolder
	for Read in $(ls genome_submission/v.inaequalis/172_pacbio/tbl2asn/final/v.inaequalis_172_pacbio_Nov_2017_Passey.sqn); do
		echo $Read;
		cp $Read $SubFolder/.
	done
cp genome_submission/v.inaequalis/172_pacbio/tbl2asn/final/v.inaequalis_172_pacbio_Nov_2017_Passey.sqn $SubFolder/.
cd $SubFolder
gzip v.inaequalis_172_pacbio_Nov_2017_Passey.sqn
ftp ftp-private.ncbi.nlm.nih.gov
cd uploads/tom.passey@emr.ac.uk_GHO2Umdl
mkdir Vi_annotated_PRJNA354841
cd Vi_annotated_PRJNA354841
put v.inaequalis_172_pacbio_Nov_2017_Passey.sqn.gz
mput v.inaequalis_172_pacbio_Nov_2017_Passey.sqn.gz
bye
cd ../
#rm -r $SubFolder
```
<!--
```bash
for File in $(ls genome_submission/v.*/*_ncbi/tbl2asn/final/errorsummary.val); do
Organism=$(echo $File | rev | cut -f5 -d '/' | rev);
Strain=$(echo $File | rev | cut -f4 -d '/' | rev);
echo "$Organism - $Strain";
cat $File;
echo "Duplicated genes:"
cat genome_submission/$Organism/$Strain/tbl2asn/round1/genome.val | grep 'DuplicateFeat' | cut -f4 -d ':' | cut -f2 -d' '
echo "";
done > genome_submission/172_pacbio_isolate_errors.txt
```


The final error report contained the following warnings. These were judged to be
legitimate concerns but biologically explainable.
```
67 WARNING: SEQ_FEAT.PartialProblem
 5 WARNING: SEQ_FEAT.ProteinNameEndsInBracket
211 WARNING: SEQ_FEAT.ShortExon
18 WARNING: SEQ_FEAT.SuspiciousFrame
 5 INFO:    SEQ_FEAT.PartialProblem

 Note -
 *SEQ_FEAT.partial problem. In this case, upon investigation these genes were hannging
 off the end of a contig but did not have an mRNA feature that went off of the
 end of the contig. This was occuring due to an intron being predicted hanging
 off the contig. An example on the ncbi guidelines here shows this to be
 acceptable:
 http://www.ncbi.nlm.nih.gov/genbank/eukaryotic_genome_submission_annotation#Partialcodingregionsinincompletegenomes
 *SEQ_FEAT.ProteinNameEndsInBracket. These gene names include brackets for good
 reason
 ```
```bash
 for Assembly in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa | grep 'N139_ncbi' | grep -v 'old'); do
 Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
 Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
 mkdir -p tmp_assembly/$Organism/$Strain
 cp $Assembly tmp_assembly/$Organism/$Strain/.
 GffFile=$(ls gene_pred/final_genes/$Organism/$Strain/final/final_genes_appended.gff3)
 cp $GffFile tmp_assembly/$Organism/$Strain/.
 GeneConversions=$(ls genome_submission/$Organism/$Strain/gag/edited/genome_gene_conversions.tsv)
 cp $GeneConversions tmp_assembly/$Organism/$Strain/.
 done
 ```

 -->