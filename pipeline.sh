# Evolution.Inbreeding
Evolution Submission

# Pipeline developed for Smith et al. (2018) - Strikingly high levels of heterozygosity despite 
# 20 years of inbreeding in a clonal honey bee

# You will need the following in your PATH
# bwa
# samtools
# picard
# GATK
# VCFtools
# bedtools
# blastn
# SnpEff

# BWA index (more information @ http://bio-bwa.sourceforge.net/bwa.shtml)
# Generates the following files:
# .amb
# .ann
# .bwt
# .pac
# .sa

bwa index -a bwtsw /path/to.fasta

# SAMtools index (more information @ http://www.htslib.org/doc/samtools.html):
# Generates the following files:
# .fai

samtools faidx /path/to/.fasta

# Picard index dictionary (more information @ https://broadinstitute.github.io/picard/command-line-overview.html#Overview):
# Generates the following files:
# .dict

java -jar /path/to/CreateSequenceDictionary.jar R=/path/to/.fasta
O=/path/to/.dict

# BWA mem (more information @ http://bio-bwa.sourceforge.net/bwa.shtml)
# Generates the following file:
# .bam 

bwa mem -M -t -R '@RG\tID:\tPL:\tPU:\tSM:\tCN:' /path/to.fasta
/path/toR1.fastq.gz /path/to/R2.fastq.gz | samtools view - bSho
/path/to.bam -

# SAMtools sort (more information @ http://www.htslib.org/doc/samtools.html):
# sorts .bam file created above 

samtools sort /path/to.bam /path/to/_sort

# Picard mark duplicates (more information @ https://broadinstitute.github.io/picard/command-line-overview.html#Overview):
# Generates the following files:
# _rmdup.bam
# _markCloneDupMetrics.txt

java -jar /path/to/MarkDuplicates.jar INPUT=/path/to/_sort.bam
OUTPUT=/path/to/_rmdup.bam METRICS_FILE=/path/to_markCloneDupMetrics.txt

# SAMtools index (more information @ http://www.htslib.org/doc/samtools.html):
# Generates the following file:
# _rmdup.bam.bai

samtools index /path/to_rmdup.bam

# Realign around indels with GATK (more information @ https://software.broadinstitute.org/gatk/gatkdocs/3.7-0/):
# Generates the following files:
# _realigner.intervals
# _realigned.bam

java -jar /path/to/GenomeAnalysisTK.jar -T RealignerTargetCreator -R
/path/to.fasta -I /path/to/_rmdup.bam -o /path/to/_realigner.intervals

java -jar /path/to/GenomeAnalysisTK.jar -T IndelRealigner -R
/path/to/.fasta -I /path/to/_rmdup.bam -targetIntervals
/path/to/_realigner.intervals -o /path/to/_realigned.bam

# Recalibrate base qualities with GATK (more information @ https://software.broadinstitute.org/gatk/gatkdocs/3.7-0/):
# Generates the following files:
# _recal_data.table
# _final.bam

java -jar /path/to/GenomeAnalysisTK.jar -T BaseRecalibrator –R
/path/to/.fasta -I /path/to/_realigned.bam -knownSites /path/to/.vcf -o
/path/to/_recal_data.table

java -jar /path/to/GenomeAnalysisTK.jar -T PrintReads -R /path/to/.fasta
-I /path/to/_realigned.bam -BQSR /path/to/_recal_data.table –o
/path/to_final.bam

# Call SNPs and indels with GATK Haplotype Caller (more information @ https://software.broadinstitute.org/gatk/gatkdocs/3.7-0/):
# Generates the following file:
# output.raw.snps.indels.vcf

java -jar /path/to/GenomeAnalysisTK.jar -T HaplotypeCaller -R
/path/to/.fasta -I /path/to/_final.bam --dbsnp /path/to/.vcf
-stand_call_conf 30 -stand_emit_conf 10 -o
/path/output.raw.snps.indels.vcf

# create file with only SNPs

java -jar /path/to/GenomeAnalysisTK.jar -T SelectVariants -R 
/path/to/.fasta -V /path/to/.vcf -selectType SNP -o /path/to/.vcf

# create file with only INDELs

java -jar /path/to/GenomeAnalysisTK.jar -T SelectVariants -R 
/path/to/.fasta -V /path/to/.vcf -selectType INDEL -o /path/to/.vcf

# remove data from scaffolds 17 & 18 

sed '/^17./ d' /path/to/.vcf > /path/to/.vcf 
sed '/^18./ d' /path/to/.vcf > /path/to/.vcf 

# remove low quality SNPs

gawk '$0~/^#/ || $7!~/LowQual/' /path/to/.vcf  > /path/to/.vcf

# identify SNPs with a Fisher's strand bias > 60

java -jar /path/to/GenomeAnalysisTK.jar -T VariantFiltration -R 
/path/to/.fasta -V /path/to/.vcf --filterExpression "FS > 60.0" --filterName "rFS" -o /path/to/.vcf

# remove SNPs with a strand bias

grep -v 'rFS' path/to/.vcf > path/to/.vcf


# identify SNPs with a mapping quality < 40

java -jar /path/to/GenomeAnalysisTK.jar -T VariantFiltration -R
/path/to/.fasta -V /path/to/.vcf --filterExpression "MQ < 40.0" --filterName "rMQ" 
-o /path/to/

# remove SNPs with a mapping quality < 40

grep -v 'rMQ' /path/to/.vcf > /path/to/.vcf

# identify SNPs with a quality of depth value < 2

java -jar /path/to/GenomeAnalysisTK.jar -T VariantFiltration -R
/path/to/.fasta -V /path/to/.vcf --filterExpression "QD < 2.0" --filterName "rQD" 
-o /path/to/

# remove SNPs with a quality of depth value < 2

grep -v 'rQD' /path/to/.vcf > /path/to/.vcf

