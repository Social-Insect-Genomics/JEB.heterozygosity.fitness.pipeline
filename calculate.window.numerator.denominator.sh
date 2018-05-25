# Create 10kb, 100kb and 250kb windows across each scaffold and chromosome, 
# note need to change denominator values

bedtools makewindows -g .genome -w 10000 > amel.10.kb.windows.bed
bedtools makewindows -g .genome -w 100000 > amel.100.kb.windows.bed
bedtools makewindows -g .genome -w 250000 > amel.250.kb.windows.bed

# Calculate number of het snps per 10kb, 100kb, 250kb window
# Plot number of het snps per 10kb, 100kb and 250kb window Figure 1, Supplementary Figure 4 & 5

bedtools coverage -a amel.10.kb.windows.bed -b .vcf -counts > .bed
bedtools coverage -a amel.100.kb.windows.bed -b .vcf -counts > .bed
bedtools coverage -a amel.250.kb.windows.bed -b .vcf -counts > .bed

# Now create correct denominator values
# Calculate number of callable sites per 10kb, 100kb, 250kb window

bedtools coverage -a amel.10.kb.windows.bed -b .bed -counts > .bed
bedtools coverage -a amel.100.kb.windows.bed -b .bed -counts > .bed
bedtools coverage -a amel.250.kb.windows.bed -b .bed -counts > .bed
