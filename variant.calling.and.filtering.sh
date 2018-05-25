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
# bcftools
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

# identify SNPs with a genotype quality < 20

vcftools --vcf /path/to/.vcf --minGQ 20 --recode --out /path/to/.vcf

# remove SNPs with a genotype quality < 20

sed '/\.\/\./d' /path/to/.vcf > /path/to/.vcf

# identify SNPs with a read depth < 10

vcftools --vcf /path/to/.vcf --minDP 10 --recode --out /path/to/.vcf

# remove SNPs with a read depth < 10

sed '/\.\/\./d' /path/to/.vcf > /path/to/.vcf

# identify SNPs within 10 bps of indels

java -jar /path/to/GenomeAnalysisTK.jar -T VariantFiltration -R
/path/to/.fasta -V /path/to/.vcf --mask /path/to/.vcf --maskExtension 10 
--maskName "InDel" -o /path/to/.vcf

# remove SNPs within 10 bps of indels

grep -v 'InDel' /path/to/.vcf > /path/to/.vcf

# remove SNPs within repeat regions

vcftools --vcf /path/to/.vcf --exclude-bed repeats.bed --recode --out /path/to/.vcf

# Now need to find centromere positions on each chromosome... Solignac et al. (2007) to the rescue.
# Create database from reference genome 

makeblastdb -in /path/to/.fasta -parse_seqids -dbtype nucl

# blastn microsatellite sequences in centromeric regions from Solignac et a. (2007) 
# to A. mellifera reference genome for each chromosome

blastn -query /path/to/solignac_sequences.blast -db /path/to/.fasta -task

blastn -outfmt 7 -max_target_seqs 10 -evalue 0.5 -perc_identity 95 >
path/to/.out

# Now that we know where the centromere positions are we can remove variants within them

vcftools --vcf /path/to/.vcf --exclude-bed centromere.positions.bed --recode --out /path/to/.vcf

# remove sites with more than two alleles

vcftools --vcf /path/to/.vcf --recode --remove-filtered-all --max-alleles 2 --out /path/to/.vcf

# identify heterozygous SNPs with an allele ratio < 0.15 and SNPs within potentially duplicated regions

bcftools annotate -x INFO,^FORMAT/AD,^FORMAT/GT /path/to/.vcf > /path/to/.vcf

sed '/^#/d' /path/to/.vcf > /path/to/.vcf

cat /path/to/.vcf.vcf | awk '{print $1 "\t" $2 "\t" $10}' > /path/to/.bed

bedtools subtract -header -a /path/to/.vcf -b /path/to/.bed > /path/to/.vcf

bedtools subtract -header -a /path/to/.vcf -b /path/to/.bed > /path/to/.vcf

# identify and remove heterozygous and homozygous SNPs

grep -v '1/1:' /path/to/.vcf > /path/to/.vcf

grep -v '0/1:' /path/to/.vcf > /path/to/.vcf

# annotate SNPs with SnpEff

snpeff amel.gff3 -s /path/to/.html -no-downstream -no-upstream /path/to/.vcf > /path/to/.vcf

