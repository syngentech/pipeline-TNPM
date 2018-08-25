# Running Guide

## Setup the environment

```bash
# Build the image
sudo docker build --tag tnpm:latest --rm /PATH/TO/pipeline-TNPM

# Run in container
sudo docker run --name tnpm --volume /PATH/TO/WORKING_DIR:/data -it tnpm
```

## Prepare references

```bash
# Store in /data/genomes
mkdir -p /data/genomes && cd /data/genomes

# Download data files
wget \
  ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/Homo_sapiens_assembly38.fasta.gz \
  ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
  ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/dbsnp_146.hg38.vcf.gz \
  ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_22/gencode.v22.chr_patch_hapl_scaff.annotation.gtf.gz

# Decompress data files
gzip --decompress --stdout Homo_sapiens_assembly38.fasta.gz \
  > Homo_sapiens_assembly38.fasta
gzip --decompress --stdout 1000G_phase1.snps.high_confidence.hg38.vcf.gz \
  > 1000G_phase1.snps.high_confidence.hg38.vcf
gzip --decompress --stdout dbsnp_146.hg38.vcf.gz \
  > dbsnp_146.hg38.vcf
gzip --decompress --stdout gencode.v22.chr_patch_hapl_scaff.annotation.gtf.gz \
  > gencode.v22.chr_patch_hapl_scaff.annotation.gtf

# Convert genome names
~/utils/map_genomes.rb < gencode.v22.chr_patch_hapl_scaff.annotation.gtf > Homo_sapiens_assembly38.gtf

# Prepare references for STAR
mkdir -p /data/genomes/Homo_sapiens_assembly38.STAR
STAR --runThreadN 8 --runMode genomeGenerate \
  --genomeDir Homo_sapiens_assembly38.STAR \
  --genomeFastaFiles Homo_sapiens_assembly38.fasta \
  --sjdbGTFfile Homo_sapiens_assembly38.gtf

# Prepare references for GATK
samtools faidx Homo_sapiens_assembly38.fasta
gatk CreateSequenceDictionary -R Homo_sapiens_assembly38.fasta -O Homo_sapiens_assembly38.dict
gatk IndexFeatureFile -F 1000G_phase1.snps.high_confidence.hg38.vcf
gatk IndexFeatureFile -F dbsnp_146.hg38.vcf
```

## Align reads

```bash
# Working under /data
cd /data

# Align reads
STAR --runThreadN 8 --runMode alignReads \
  --outSAMtype BAM SortedByCoordinate --outSAMmapqUnique 60 \
  --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 \
  --genomeDir genomes/Homo_sapiens_assembly38.STAR \
  --outFileNamePrefix DAxxx_ --readFilesIn DAxxx_1.fastq DAxxx_2.fastq

# Mark duplicates
STAR --runThreadN 8 --runMode inputAlignmentsFromBAM \
  --bamRemoveDuplicatesType UniqueIdentical \
  --outFileNamePrefix DAxxx_ --inputBAMfile DAxxx_Aligned.sortedByCoord.out.bam

# Add @RG info
gatk AddOrReplaceReadGroups -LB lib -PL illumina -PU unit \
  -SM DAxxx -I DAxxx_Processed.out.bam -O DAxxx.star.bam

# Recalibrate base quality scores
gatk BaseRecalibrator -R genomes/Homo_sapiens_assembly38.fasta \
  --known-sites genomes/dbsnp_146.hg38.vcf \
  --known-sites genomes/1000G_phase1.snps.high_confidence.hg38.vcf \
  -I DAxxx.star.bam -O DAxxx.star.table
gatk ApplyBQSR -R genomes/Homo_sapiens_assembly38.fasta \
  -I DAxxx.star.bam --bqsr-recal-file DAxxx.star.table -O DAxxx.bam
```

## Call variants

```bash
# Call somatic variants
gatk Mutect2 -R genomes/Homo_sapiens_assembly38.fasta \
  -I DAxxx.bam -tumor DAxxx -I DAyyy.bam -normal DAyyy \
  -O PRzzz.vcf.gz -bamout PRzzz.bam
```
