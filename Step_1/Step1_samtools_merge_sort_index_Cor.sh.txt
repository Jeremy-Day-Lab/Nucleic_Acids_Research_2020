module load SAMtools

cd /data/user/rphill3/Enhancer_Ident_ATAC

samtools merge Cor_ATAC_merged.bam  ATAC_K2_S1.bam ATAC_K3_S2.bam ATAC_K4_S3.bam ATAC_V1_S4.bam ATAC_V2_S5.bam ATAC_V3_S6.bam

samtools sort Cor_ATAC_merged.bam -o Cor_ATAC_merged.sorted.bam

samtools index Cor_ATAC_merged.sorted.bam Cor_ATAC_merged.sorted.bam.bai
