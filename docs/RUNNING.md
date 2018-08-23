# Running Guide

## Setup the environment

To build and run the docker environment:

```bash
sudo docker build --tag tnpm:latest --rm /PATH/TO/pipeline-TNPM
sudo docker run --name tnpm --volume /PATH/TO/WORKING_DIR:/data -it tnpm
```

## Prepare references

To get the reference genomes prepared:

```bash
# Decompress data files
gzip --decompress --stdout GCF_000001405.38_GRCh38.p12_genomic.fna.gz > GRCh38.p12.fa
gzip --decompress --stdout GCF_000001405.38_GRCh38.p12_genomic.gff.gz > GRCh38.p12.gff

# Prepare references for STAR
gffread GRCh38.p12.gff -T -o GRCh38.p12.gtf
mkdir GRCh38.p12-STAR
STAR --runThreadN 8 --runMode genomeGenerate \
  --genomeDir ./GRCh38.p12-STAR --genomeFastaFiles GRCh38.p12.fa --sjdbGTFfile GRCh38.p12.gtf

# Prepare references for GATK
samtools faidx GRCh38.p12.fa
gatk CreateSequenceDictionary -R GRCh38.p12.fa -O GRCh38.p12.dict
```

## Align reads

To map reads and mark duplicates:

```bash
# Align reads
STAR --runThreadN 8 --runMode alignReads \
  --outSAMtype BAM SortedByCoordinate --outSAMmapqUnique 60 \
  --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 \
  --genomeDir genomes/GRCh38.p12-STAR \
  --outFileNamePrefix DAxxx_ --readFilesIn DAxxx_1.fastq DAxxx_2.fastq

# Mark duplicates
STAR --runThreadN 8 --runMode inputAlignmentsFromBAM \
  --bamRemoveDuplicatesType UniqueIdentical \
  --outFileNamePrefix DAxxx_ --inputBAMfile DAxxx_Aligned.sortedByCoord.out.bam

# Prepare BAM for GATK
gatk AddOrReplaceReadGroups -LB lib -PL illumina -PU unit \
  -SM DAxxx -I DAxxx_Processed.out.bam -O DAxxx.bam
samtools index DAxxx.bam
```

## Call variants

To use GATK from the Broad Institute to call somatic variants:

```bash
gatk Mutect2 -R genomes/GRCh38.p12.fa \
  -I DAxxx.bam -tumor DAxxx -I DAyyy.bam -normal DAyyy -O PRzzz.vcf.gz
```
