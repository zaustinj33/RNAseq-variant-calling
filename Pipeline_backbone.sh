#!/bin/bash

## Backbone of RNAseq variant calling pipeline based on GATK best practices: 
#https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-

# Setup dirs
mkdir -p working_data
mkdir -p result

# Required software:
#SAMTools
#FastQC
#trim_galore!
#STAR
#Picard
#GATK

# Pre-requisite files (human GRChv38):
GTF=$PWD/Annotation/Homo_sapiens.GRCh38.105.gtf
IDX=$PWD/Annotation/STAR_hg38_humanIndex
FA=$PWD/Annotation/Homo_sapiens.GRCh38.dna.primary_assembly.fa
SNPS=$PWD/Annotation/human_variants/Homo_sapiens_assembly38.dbsnp138.vcf
ANNO=$PWD/Annotation/human_variants/Homo_sapiens_assembly38.known_indels.vcf.gz

## Steps:
#1 FastQC
#2 trim_galore!
#3 STAR
#4 Mark Duplicates (Picard)
#5 Add read-groups
#6 Split Reads (GATK)
#7 Recalibrate Reads
#8 Call Variants
#9 Generate GVCF file
#10 Filter variants
#11 Annotate Variants
#12 MultiQC check



#1 FastQC
fastqc $2/raw_data/$1/$1_1.fq.gz
fastqc $2/raw_data/$1/$1_2.fq.gz

#2 trim_galore!
trim_galore --paired --phred33 --fastqc --illumina \
--clip_R1 6 --clip_R2 6 --dont_gzip -q 30 --length 30 \
$2/raw_data/$1/$1_1.fq.gz $2/raw_data/$1/$1_2.fq.gz \
--output_dir $2/working_data/$1

#3 STAR - also outputs expression matrix
STAR --runThreadN 32 --sjdbGTFfile $GTF \
--readFilesCommand gunzip \
-c --outSAMstrandField intronMotif \
--genomeDir $IDX --limitBAMsortRAM=40000000000 \
--readFilesIn $2/working_data/$1/$1_1_val_1.fq.gz $2/working_data/$1/$1_1_val_1.fq.gz \
--quantMode GeneCounts --outFileNamePrefix  $2/result/$1/$1_STARout

#4 Mark duplicates
picard MarkDuplicates \
INPUT=$2/result/$1/$1_STARout OUTPUT=$1_markDups.bam \
METRICS_FILE=$1_markDups_metrics.txt \
REMOVE_DUPLICATES=false ASSUME_SORTED=true PROGRAM_RECORD_ID='null' \
VALIDATION_STRINGENCY=LENIENT \
CREATE_INDEX=TRUE

#5 Add read-groups
picard AddOrReplaceReadGroups \
INPUT=$1_markDups.bam OUTPUT=$1_markDups_groups.bam \
RGLB=LB RGPL=ILLUMINA RGPU=PU RGSM=$1 \
CREATE_INDEX=TRUE

#6 Split Reads
gatk SplitNCigarReads \
-I $1_markDups.bam \
-R $FA \
-O $1_split.bam -OBI

#7 Recalibrate Reads
gatk BaseRecalibrator \
-I $1_split.bam \
-R $FA \
-O $1_table.recal \
--known-sites $SNPS \
--verbosity INFO \
--java-options -Xmx 50g

gatk ApplyBQSR \
-R $FA \
-I $1_split.bam \
--bqsr-recal-file $1_table.recal \
-O $1_recal.bam \
--create-output-bam-index true \
--java-options -Xmx 50g

#8 Call Variants
gatk HaplotypeCaller \
-I $1_recal_bam \
-R $FA \
-O $1_raw_variants.vcf \
--dont-use-soft-clipped-bases \
-stand-call-conf 20.0 -ERC GVCF -OBI \
--annotation MappingQualityRankSumTest \
        --annotation QualByDepth \
        --annotation ReadPosRankSumTest \
        --annotation RMSMappingQuality \
        --annotation FisherStrand \
        --annotation Coverage \
        --dbsnp $SNPS \
        --verbosity INFO \
        --java-options -Xmx 50g

#9 Generate GVCF file
gatk GenotypeGVCFs \
-R $FA \
--dbsnp $SNPS \
-V $1_raw_variants.vcf \
-O $1_gvcf.vcf

#10 Filter variants
gatk VariantFiltration \
-R $FA \
-V $1_raw_variants.vcf \
-O $1_filtered.vcf \
-window 35 -cluster 3 \
--filter-name FS \
--filter-expression "FS > 30.0" \
--filter-name QD \
--filter-expression "QD < 2.0"

#11 Annotate Variants
table_annovar.pl $1_filtered.vcf $ANNO \
--outfile $1 \
--buildver $FA \
--protocol refGene,cosmic87_coding,cosmic87_noncoding,clinvar_20180603,avsnp150,1000g2015aug_all,gnomad_genome,dbnsfp35a,dbscsnv11 \
--operation g,f,f,f,f,f,f,f,f \
--vcfinput \
--polish \
--remove

# MultiQC check
multiqc .

echo 'DONE'

