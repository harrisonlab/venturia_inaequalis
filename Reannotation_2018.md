# Venturia inaequalis Reannotation

This document sets out the reannotation and gene prediction for the V. inaequalis genome based on that used for Alternaria. This has been run seperately (see PacBio_assembly.mdown) due to not being able to work out why the original assembled genome would not pass initail validation on submission to NCBI after removal of duplicates identified in the first submission.

## Gene Prediction


#### Aligning



```bash
  for Assembly in $(ls repeat_masked/*/*/filtered_contigs_repmask/*_contigs_unmasked.fa | grep -w '172_pacbio'); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    for rna_file in $(ls qc_rna/*/*/*/*/*.gz | grep -w 'paired'); do
      Timepoint=$(echo $rna_file | rev | cut -f3 -d '/' | rev)
      echo "$Timepoint"
      OutDir=alignment/star/$Organism/$Strain/$Timepoint/$Prefix
      mkdir -p $OutDir
      ProgDir=/home/passet/git_repos/tools/seq_tools/RNAseq
      qsub $ProgDir/sub_star.sh $Assembly $rna_file $OutDir
    done
  done
```


Accepted hits .bam file were concatenated and indexed for use for gene model training:

```bash
  for OutDir in $(ls -d alignment/star/*/*); do
    Strain=$(echo $OutDir | rev | cut -d '/' -f1 | rev)
    Organism=$(echo $OutDir | rev | cut -d '/' -f2 | rev)
    echo "$Organism - $Strain"
    # For all alignments
    BamFiles=$(ls $OutDir/*/*.sortedByCoord.out.bam | tr -d '\n' | sed 's/.bam/.bam /g')
    mkdir -p $OutDir/concatenated
    samtools merge -@ 12 -f $OutDir/concatenated/concatenated.bam $BamFiles
  done
```

#### Braker prediction

Before braker predictiction was performed, I double checked that I had the
genemark key in my user area and copied it over from the genemark install
directory:

```bash
	ls ~/.gm_key
	cp /home/armita/prog/genemark/2017/gm_key_64 ~/.gm_key
```


```bash
  for Assembly in $(ls repeat_masked/*/*/filtered_contigs_*/*_contigs_softmasked.fa | grep '172_pacbio'); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/braker/$Organism/"$Strain"_braker/2018
    mkdir -p $OutDir
    AcceptedHits=$(ls alignment/star/$Organism/$Strain/concatenated/concatenated.bam)
    GeneModelName="$Organism"_"$Strain"_braker_2018
    rm -r /home/armita/prog/augustus-3.1/config/species/"$Organism"_"$Strain"_braker_2018
    ProgDir=/home/passet/git_repos/tools/gene_prediction/braker1
    qsub $ProgDir/sub_braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
  done
```


## Supplimenting Braker gene models with CodingQuary genes

Additional genes were added to Braker gene predictions, using CodingQuary in
pathogen mode to predict additional regions.

Firstly, aligned RNAseq data was assembled into transcripts using Cufflinks.

Note - cufflinks doesn't always predict direction of a transcript and
therefore features can not be restricted by strand when they are intersected.

```bash
  for Assembly in $(ls repeat_masked/*/*/filtered_contigs_*/*_contigs_softmasked.fa | grep "172_pacbio"); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated/2018
    mkdir -p $OutDir
    AcceptedHits=$(ls alignment/star/$Organism/$Strain/concatenated/concatenated.bam)
    ProgDir=/home/passet/git_repos/tools/seq_tools/RNAseq
    qsub $ProgDir/sub_cufflinks.sh $AcceptedHits $OutDir
  done
```


Secondly, genes were predicted using CodingQuary:

```bash
  for Assembly in $(ls repeat_masked/*/*/filtered_contigs_*/*_contigs_softmasked.fa | grep "172_pacbio"); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    OutDir=gene_pred/codingquary/$Organism/$Strain/2018
    mkdir -p $OutDir
    CufflinksGTF=$(ls gene_pred/cufflinks/$Organism/$Strain/concatenated/2018/transcripts.gtf)
    ProgDir=/home/passet/git_repos/tools/gene_prediction/codingquary
    qsub $ProgDir/sub_CodingQuary.sh $Assembly $CufflinksGTF $OutDir
  done
```


Then, additional transcripts were added to Braker gene models, when CodingQuary
genes were predicted in regions of the genome, not containing Braker gene
models:

```bash
for BrakerGff in $(ls gene_pred/braker/*/*_braker/2018/*/augustus.gff3); do
  Strain=$(echo $BrakerGff| rev | cut -d '/' -f4 | rev | sed 's/_braker_new//g' | sed 's/_braker_pacbio//g' | sed 's/_braker//g')
  Organism=$(echo $BrakerGff | rev | cut -d '/' -f5 | rev)
  echo "$Organism - $Strain"
  Assembly=$(ls repeat_masked/*/*/filtered_contigs_*/*_contigs_softmasked.fa | grep "172_pacbio")
  CodingQuaryGff=$(ls gene_pred/codingquary/$Organism/$Strain/2018/out/PredictedPass.gff3)
  PGNGff=$(ls gene_pred/codingquary/$Organism/$Strain/2018/out/PGN_predictedPass.gff3)
  AddDir=gene_pred/codingquary/$Organism/$Strain/2018/additional
  FinalDir=gene_pred/final/$Organism/$Strain/final_2018
  AddGenesList=$AddDir/additional_genes.txt
  AddGenesGff=$AddDir/additional_genes.gff
  FinalGff=$AddDir/combined_genes.gff
  mkdir -p $AddDir
  mkdir -p $FinalDir

    for x in $CodingQuaryGff $PGNGff; do
      bedtools intersect -v -a $x -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';'
    done > $AddGenesList

    for y in $CodingQuaryGff $PGNGff; do
      ProgDir=/home/passet/git_repos/tools/seq_tools/feature_annotation
      $ProgDir/gene_list_to_gff.pl $AddGenesList $y CodingQuarry_v2.0 ID CodingQuary
    done > $AddGenesGff
  ProgDir=/home/passet/git_repos/tools/gene_prediction/codingquary
  # -
  # This section is edited
  $ProgDir/add_CodingQuary_features.pl $AddGenesGff $Assembly > $AddDir/add_genes_CodingQuary_unspliced.gff3
  $ProgDir/correct_CodingQuary_splicing.py --inp_gff $AddDir/add_genes_CodingQuary_unspliced.gff3 > $FinalDir/final_genes_CodingQuary.gff3
  # -
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
done
```

<!--
The final number of genes per isolate was observed using:
```bash
  for DirPath in $(ls -d gene_pred/final/*/*/final_2018); do
    Strain=$(echo $DirPath| rev | cut -d '/' -f2 | rev)
    Organism=$(echo $DirPath | rev | cut -d '/' -f3 | rev)
    echo "$Organism - $Strain"
    cat $DirPath/final_genes_Braker.pep.fasta | grep '>' | wc -l;
    cat $DirPath/final_genes_CodingQuary.pep.fasta | grep '>' | wc -l;
    cat $DirPath/final_genes_combined.pep.fasta | grep '>' | wc -l;
    echo "";
  done
```

The number of genes predicted by Braker, supplimented by CodingQuary and in the
final combined dataset was shown:


In preperation for submission to ncbi, gene models were renamed and duplicate gene features were identified and removed.
 * no duplicate genes were identified


```bash
for GffAppended in $(ls gene_pred/final/*/*/final_2018/final_genes_appended_renamed.gff3); do
Strain=$(echo $GffAppended | rev | cut -d '/' -f3 | rev)
Organism=$(echo $GffAppended | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
FinalDir=gene_pred/final/$Organism/$Strain/final_2018
GffFiltered=$FinalDir/filtered_duplicates.gff
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
$ProgDir/remove_dup_features.py --inp_gff $GffAppended --out_gff $GffFiltered
GffRenamed=$FinalDir/final_genes_appended_renamed.gff3
LogFile=$FinalDir/final_genes_appended_renamed.log
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
$ProgDir/gff_rename_genes.py --inp_gff $GffFiltered --conversion_log $LogFile > $GffRenamed
rm $GffFiltered
Assembly=$(ls repeat_masked/$Organism/$Strain/*/*_softmasked_repeatmasker_TPSI_appended.fa)
$ProgDir/gff2fasta.pl $Assembly $GffRenamed gene_pred/final/$Organism/$Strain/final_2018/final_genes_appended_renamed
# The proteins fasta file contains * instead of Xs for stop codons, these should
# be changed
sed -i 's/\*/X/g' gene_pred/final/$Organism/$Strain/final_2018/final_genes_appended_renamed.pep.fasta
done
```
-->


