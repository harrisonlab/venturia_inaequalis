
#Submission of V.inaequalis genomes to ncbi 


A Bioproject and biosample number was prepared for the genome submission at:
https://submit.ncbi.nlm.nih.gov


## Submission of sequence data to SRA

Reads were submitted to the SRA at https://submit.ncbi.nlm.nih.gov/subs/sra/ .
To do this, a metadata file was provided detailing each of the files in the
bioproject. The file was downloaded in excel format and edited manually. A copy
of the edited file and the final .tsv file is present at:

```bash
  ls genome_submission/SRA_metadata_acc.txt genome_submission/SRA_metadata_acc.xlsx
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
# put Vi_sra_PRJNA354841
prompt
mput *
bye
cd ../
rm -r $SubFolder
```


