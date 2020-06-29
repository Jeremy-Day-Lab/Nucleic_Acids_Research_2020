#BAM files for individual samples were merged using samtools. Then, MACS2 was used to callpeaks. Peaks <=1000bps apart were merged within region. Here, I will read in the merged peak files, create a length column, and 
#subset peaks that are >146bps, or the length of the DNA wrapped around a single nucleosome. 
#Str
Str_Merged_1Kb_Peaks <- read.delim(file = "/Volumes/JDLab$/RobertPhillips/Bioinformatics/JD0032_EnhancerIdentification/PeakCalls_Merged_WithinRegion/Str_Merged_BAMPE_narrow_0.00001_FinalPeaks_peaks.narrowPeak.bed.sorted.merged.bed",header = FALSE)
Str_Merged_1Kb_Peaks$Length <- Str_Merged_1Kb_Peaks$V3 - Str_Merged_1Kb_Peaks$V2
Str_Merged_1Kb_Peaks <- subset(Str_Merged_1Kb_Peaks,subset = (Length >= 146))
write.table(x = Str_Merged_1Kb_Peaks,
            file = "/Volumes/JDLab$/RobertPhillips/Bioinformatics/JD0032_EnhancerIdentification/PeakCalls_Merged_WithinRegion/SizeSelected/Str_Merged_BAMPE_narrow_0.00001_FinalPeaks_peaks.narrowPeak.bed.sorted.merged.sizeselected.bed",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            sep = "\t")
#Cor
Cor_Merged_1Kb_Peaks <- read.delim(file = "/Volumes/JDLab$/RobertPhillips/Bioinformatics/JD0032_EnhancerIdentification/PeakCalls_Merged_WithinRegion/Cor_Merged_BAMPE_narrow_0.00001_FinalPeaks_peaks.narrowPeak.bed.sorted.merged.bed",header = FALSE)
Cor_Merged_1Kb_Peaks$Length <- Cor_Merged_1Kb_Peaks$V3 - Cor_Merged_1Kb_Peaks$V2
Cor_Merged_1Kb_Peaks <- subset(Cor_Merged_1Kb_Peaks,subset = (Length >= 146))
write.table(x = Cor_Merged_1Kb_Peaks,
            file = "/Volumes/JDLab$/RobertPhillips/Bioinformatics/JD0032_EnhancerIdentification/PeakCalls_Merged_WithinRegion/SizeSelected/Cor_Merged_BAMPE_narrow_0.00001_FinalPeaks_peaks.narrowPeak.bed.sorted.merged.sizeselected.bed",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            sep = "\t")
#Hpc 
Hpc_merged_1Kb_Peaks <- read.delim(file = "/Volumes/JDLab$/RobertPhillips/Bioinformatics/JD0032_EnhancerIdentification/PeakCalls_Merged_WithinRegion/Hpc_Merged_BAMPE_narrow_0.00001_FinalPeaks_peaks.narrowPeak.bed.sorted.merged.bed",header = FALSE)
Hpc_merged_1Kb_Peaks$Length <- Hpc_merged_1Kb_Peaks$V3 - Hpc_merged_1Kb_Peaks$V2
Hpc_merged_1Kb_Peaks <- subset(Hpc_merged_1Kb_Peaks,subset = (Length >= 146))
write.table(x = Hpc_merged_1Kb_Peaks,
            file = "/Volumes/JDLab$/RobertPhillips/Bioinformatics/JD0032_EnhancerIdentification/PeakCalls_Merged_WithinRegion/SizeSelected/Hpc_Merged_BAMPE_narrow_0.00001_FinalPeaks_peaks.narrowPeak.bed.sorted.merged.sizeselected.bed",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            sep = "\t")

#Now rbind all Peaks
All_Regions <- rbind(Str_Merged_1Kb_Peaks,Cor_Merged_1Kb_Peaks,Hpc_merged_1Kb_Peaks)
write.table(x = All_Regions,
            file = "/Volumes/JDLab$/RobertPhillips/Bioinformatics/JD0032_EnhancerIdentification/PeakCalls_Merged_WithinRegion/SizeSelected/AllRegions_BAMPE_narrow_0.00001_FinalPeaks_peaks.narrowPeak.bed.sorted.merged.sizeselected.bed",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            sep = "\t")

#Using bedtools the All_Regions file generated above was then 