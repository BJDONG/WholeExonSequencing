#!/bin/bash
# pipeline for WES

sample=$1
workdir=/home/Bingjie/opt/02-working/test

# path for tools
cutbarcode=/home/Bingjie/opt/02-working/test/cutbarcode.py
BWA=/home/Bingjie/opt/tools/bwa/bwa-0.7.15/bwa
SAMTOOLS=/home/Bingjie/opt/tools/samtools/samtools-1.4.1/samtools
PICARD=/home/Bingjie/opt/tools/picard/picard-2.9.2/picard.jar
GATK=/home/Bingjie/opt/tools/GATK/GenomeAnalysisTK-3.7.0.jar


# path for reference 
hg19=/home/Bingjie/opt/data/human/hg19
known1000G_indels=/home/Bingjie/opt/data/human/gatk/1000G_phase1.indels.hg19.vcf
GoldStandard_indels=/home/Bingjie/opt/data/human/gatk/Mills_and_1000G_gold_standard.indels.hg19.vcf
dbSNP=/home/Bingjie/opt/data/human/gatk/dbsnp_138.hg19.vcf
IntervalList=/home/Bingjie/opt/data/human/Agilent/V6_S07604514/S07604514.interval_list

ID=${sample}
SM=${sample}
LB=`less ${workdir}/01-cleandata/testfile.1.clean.fq.gz | head -1 | cut -d ':' -f 3`
PU=${LB}'.'`less ${workdir}/01-cleandata/testfile.1.clean.fq.gz | head -1 | cut -d ':' -f 4`'.'${ID}
PL='Illumina'

python ${cutbarcode} ${workdir}/01-cleandata/${sample}.2.clean.fq.gz ${workdir}/02-cutbarcode/${sample}.2.clean.fq.gz

${BWA} mem -R @RG\\tID:$ID\\tSM:$SM\\tPU:$PU\\tLB:$LB\\tPL:$PL -M ${hg19}/hg19 ${workdir}/01-cleandata/${sample}.1.clean.fq.gz ${workdir}/02-cutbarcode/${sample}.2.clean.fq.gz | ${SAMTOOLS} view -bS -o ${workdir}/03-mapping/${sample}.bam

${SAMTOOLS} sort ${workdir}/03-mapping/${sample}.bam -T ${workdir}/03-mapping/${sample} -o ${workdir}/03-mapping/${sample}.sort.bam

java -jar ${PICARD} MarkDuplicates VALIDATION_STRINGENCY=SILENT INPUT=${workdir}/03-mapping/${sample}.sort.bam OUTPUT=${workdir}/04-picard/${sample}.marked.bam METRICS_FILE=${workdir}/04-picard/${sample}.Mkdup.metrics

${SAMTOOLS} index ${workdir}/04-picard/${sample}.marked.bam

java -jar ${GATK} -T RealignerTargetCreator -U ALLOW_N_CIGAR_READS -R ${hg19}/hg19.fa -I ${workdir}/04-picard/${sample}.marked.bam -o ${workdir}/05-GATK/${sample}.marked.intervals -known ${known1000G_indels} -known ${GoldStandard_indels}

java  -jar ${GATK} -T IndelRealigner -R ${hg19}/hg19.fa -I ${workdir}/04-picard/${sample}.marked.bam -targetIntervals ${workdir}/05-GATK/${sample}.marked.intervals -o ${workdir}/05-GATK/${sample}.marked.realn.bam -known $known1000G_indels -known $GoldStandard_indels 

java -jar ${GATK} -T BaseRecalibrator -R ${hg19}/hg19.fa -I ${workdir}/05-GATK/${sample}.marked.realn.bam -o ${workdir}/05-GATK/${sample}.marked.realn.recal -knownSites ${known1000G_indels} -knownSites ${GoldStandard_indels} -knownSites ${dbSNP}

java -jar ${GATK}  -T PrintReads -R ${hg19}/hg19.fa -I ${workdir}/05-GATK/${sample}.marked.realn.bam -o ${workdir}/05-GATK/${sample}.marked.realn.recal.bam --BQSR ${workdir}/05-GATK/${sample}.marked.realn.recal

java -jar ${GATK} -T HaplotypeCaller -R ${hg19}/hg19.fa -L ${IntervalList} --dbsnp ${dbSNP} --emitRefConfidence GVCF -I ${workdir}/05-GATK/${sample}.marked.realn.recal.bam  -o ${workdir}/06-germline/${sample}.marked.realn.recal.bam.g.vcf

java -jar ${PICARD} CollectHsMetrics I=${workdir}/05-GATK/${sample}.marked.realn.bam O=${workdir}/07-HsMetrics/${sample}.marked.realn.recal.bam.metrics R=${hg19}/hg19.fa BAIT_INTERVALS=${IntervalList} TARGET_INTERVALS=${IntervalList}

#https://broadinstitute.github.io/picard/picard-metric-definitions.html#HsMetrics
