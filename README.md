# **Enhancer RNAs predict enhancer-gene regulatory links and are critical for enhancer function in neuronal systems**


This page contains bash and R code, as well as SeqMonk workflow information, for the enhancer identification pipeline used in "Enhancer RNAs predict enhancer-gene regulatory links and are critical for enhancer function in neuronal systems"  


```mermaid
graph TD;
  A-->B;
  A-->C;
  B-->D;
  C-->D;
```

## **Study design**

snRNA-seq datasets generated here were obtained from two experimental preparations:

1. The rat nucleus accumbens (NAc), with samples taken from four experimental groups:  
	*Male adult rats, 1hr following saline injection  
	*Male adult rats, 1hr following cocaine injection (20mg/kg, I.P.)  
	*Female adult rats, 1hr following saline injection  
	*Female adult rats, 1hr following cocaine injection (20mg/kg, I.P.)  

2. Primary striatal neurons (mixed from male and female E18 rat brains and cultured to DIV11) from four experimental groups:  
	*Vehicle (media alone, 1hr)  
	*Dopamine (50µM, 1hr)  
	*SKF-38393 (1µM, 1hr)  
	*Potassium chloride (25mM, 1hr)  

Datasets in this manuscript used the 10X Genomics Chromium Single Cell 3’ library construction kit to generate barcoded libraries for Illumina sequencing. For each dataset, we provide R Markdown and Word documents containing analysis parameters and code, as well as Python notebooks for generation of enhanced experimental signal (EES) values from the MELD package (Krishnaswamy Lab).  


## **Citation**

Savell, K.E.*, Tuscher, J.J.*, Zipperly, M.E.*, Duke, C.G.*, Phillips III, R.A.*, Bauman, A.J., Thukral, S., Sultan, F.A., Goska, N.A., Ianov, L., & Day, J.J. (2020). A dopamine-induced gene expression signature regulates neuronal function and cocaine response. *Science Advances* 6(26): eaba4221.
[Link](https://advances.sciencemag.org/content/6/26/eaba4221)


## **Links**

All Day lab resources may be found at the [Day Lab website](http://day-lab.org/resources)  
[Ratlas Shiny App](https://day-lab.shinyapps.io/ratlas/)  
[Published manuscript](https://advances.sciencemag.org/content/6/26/eaba4221)  
[BioRxiv preprint](https://www.biorxiv.org/content/10.1101/781872v1)  


## **Raw data**

Primary striatal neuron culture snRNA-seq - GEO accession [GSE137759](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137759)  
Adult NAc snRNA-seq - GEO accession [GSE137763](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137763)  


## **Source code**

[Manifold Enhancement of Latent Dimensions (MELD)](https://github.com/KrishnaswamyLab/MELD)  
[Seurat v3](https://github.com/satijalab/seurat)  