# **Enhancer RNAs predict enhancer-gene regulatory links and are critical for enhancer function in neuronal systems**


This page contains bash and R code, as well as SeqMonk workflow information, for the enhancer identification pipeline used in "Enhancer RNAs predict enhancer-gene regulatory links and are critical for enhancer function in neuronal systems"  


## **General Enhancer Identification Workflow**

The diagram below outlines a general workflow for the identifcation of enhancers. Code corresponding to each step can be found within the repository. 

```mermaid
graph TD;
A[1. For each brain region separately, merge BAM files with SAMtools] --> B[2. For each brain region separately, call peaks with MACS2];
  B --> C[3. For each brain region separately, merge peaks within 1kb using bedtools];
  C --> D[4. For each brain region separately, remove any peaks smaller than 146 bp or the legnth of DNA wrapped around a nucleosome];
  D --> E[5. Combine all peaks  with rbind in R. These peaks are the Regions of Open Chromatin or ROCs];
  E --> F[6. Keep only intergenic ROCs or iROCs];
  F --> G[7. Quantify transcription within iROCs];
  G --> H[8. Identify all genes within 1MB up and downstream from TAPE center];
  H --> I[9. Correlate eRNA and mRNA transcriptional abundance];
	style A fill:#98EFFF,stroke:#000,stroke-width:2px;	
	style B fill:#98EFFF,stroke:#000,stroke-width:2px;	
	style C fill:#98EFFF,stroke:#000,stroke-width:2px;	
	style D fill:#98EFFF,stroke:#000,stroke-width:2px;	
	style E fill:#98EFFF,stroke:#000,stroke-width:2px;	
	style F fill:#98EFFF,stroke:#000,stroke-width:2px;	
	style G fill:#98EFFF,stroke:#000,stroke-width:2px;	
	style H fill:#98EFFF,stroke:#000,stroke-width:2px;	
	style I fill:#98EFFF,stroke:#000,stroke-width:2px;	
```

## **NGS Experimental Details**

RNA-Seq and ATAC-Seq datasets were generated from striatal, cortical, and hippocampal primary neuron cultures treated with 10mM KCl or a vehicle solution for one hour. Library preparation details can be found in the methods section of the manuscript. 

## **Citation**

Nancy V.N. Carullo, Robert A. Phillips III, Rhiana C. Simon, Salomon A. Roman Soto, Jenna E. Hinds, Aaron J. Salisbury, Jasmin S. Revanna, Kendra D. Bunner, Lara Ianov, Faraz A. Sultan, Katherine E. Savell, Charles A. Gersbach, Jeremy J. Day. Enhancer RNAs predict enhancer-gene regulatory links and are critical for enhancer function in neuronal systems. 2020. *Nucleic Acids Research*. 


## **Links**

All Day lab resources may be found at the [Day Lab website](http://day-lab.org/resources)  
[BioRxiv preprint](https://www.biorxiv.org/content/10.1101/270967v3)  


## **Raw data**

Raw data will be added once the GEO submission is finalized. 

