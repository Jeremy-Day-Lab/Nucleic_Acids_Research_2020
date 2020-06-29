module load SAMtools

cd /data/user/rphill3/Enhancer_Ident_ATAC

samtools merge Hpc_ATAC_merged.bam AS0011_K1_S2.bam AS0011_K2_S4.bam AS0011_K4_S6.bam AS0011_V1_S1.bam AS0011_V2_S3.bam AS0011_V4_S5.bam

samtools sort Hpc_ATAC_merged.bam -o Hpc_ATAC_merged.sorted.bam

samtools index Hpc_ATAC_merged.sorted.bam Hpc_ATAC_merged.sorted.bam.bai
