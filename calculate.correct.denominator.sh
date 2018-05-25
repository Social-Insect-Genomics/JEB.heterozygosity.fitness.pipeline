#find callable sites within windows

#get every single position of amel4.5 genome and the read depth

 bedtools genomecov -ibam .bam -d > .bedgraph

# convert to bed and print 1,2,2,3 first four columns for .bed format

 awk '{print $1 "\t" $2 "\t" $2 "\t" $3}' .bedgraph > .bed

# remove scaffolds 17 and 18

sed '/^17./ d' .bed > .bed
sed '/^18./ d' .bed > .bed

# remove centromeres and repeats

bedtools subtract -a .bed -b centromere.positions.bed > .bed
bedtools subtract -a .bed -b repeats.bed > .bed

# remove low and high depth sites

awk '$4 < 10'  .bed > .bed
awk '$4 > twice.average.depth'  .bed > .bed
