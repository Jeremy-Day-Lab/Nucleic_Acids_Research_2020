#Load bedtools
module load BEDTools

cd /data/user/rphill3/Enhancer_Ident_ATAC/Final_PeakCalls

for i in *.narrowPeak.bed
do
#sort the bed
sort -k1,1 -k2,2n $i > Sorted_beds/$i.sorted.bed
#merge peaks within 1000bps of each other 
bedtools merge -i Sorted_beds/$i.sorted.bed -d 1000 > Merged_beds/$i.sorted.merged.bed
done

