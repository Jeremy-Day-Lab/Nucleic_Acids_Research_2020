module load SAMtools

cd /data/user/rphill3/Enhancer_Ident_ATAC

samtools merge Str_ATAC_merged.bam AS0008_K1_S11.bam AS0008_K2_S12.bam AS0008_K3_S13.bam AS0008_V1_S14.bam AS0008_V2_S15.bam AS0008_V4_S16.bam 

samtools sort Str_ATAC_merged.bam -o Str_ATAC_merged.sorted.bam

samtools index Str_ATAC_merged.sorted.bam Str_ATAC_merged.sorted.bam.bai



