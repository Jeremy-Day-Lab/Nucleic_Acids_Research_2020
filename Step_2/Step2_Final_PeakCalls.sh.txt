module load MACS2

macs2 callpeak \
        --treatment /data/user/rphill3/Enhancer_Ident_ATAC/Str_ATAC_merged.bam \
        --qvalue 0.00001 \
        --gsize 2729862089 \
        --format BAMPE \
        --outdir /data/user/rphill3/Enhancer_Ident_ATAC/Final_PeakCalls \
        --name Str_Merged_BAMPE_narrow_0.00001_FinalPeaks

macs2 callpeak \
        --treatment /data/user/rphill3/Enhancer_Ident_ATAC/Cor_ATAC_merged.bam \
        --qvalue 0.00001 \
        --gsize 2729862089 \
        --format BAMPE \
        --outdir /data/user/rphill3/Enhancer_Ident_ATAC/Final_PeakCalls \
        --name Cor_Merged_BAMPE_narrow_0.00001_FinalPeaks


macs2 callpeak \
        --treatment /data/user/rphill3/Enhancer_Ident_ATAC/Hpc_ATAC_merged.bam \
        --qvalue 0.00001 \
        --gsize 2729862089 \
        --format BAMPE \
        --outdir /data/user/rphill3/Enhancer_Ident_ATAC/Final_PeakCalls \
        --name Hpc_Merged_BAMPE_narrow_0.00001_FinalPeaks
