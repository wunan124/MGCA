# MGCA
Welcome to the GitHub repository associated with our research  
**Involvement of *CDH11* in ASD suggested by matched-gene co-expression analysis and mouse behavioral studies**  
Nan Wu, Yue Wang, Jing-Yan Jia, Yi-Hsuan Pan, Xiao-Bing Yuan

## Abstract
A large number of putative risk genes of ASD have been reported. The functions of most of these susceptibility genes in developing brains remain unknown, and a causal relationship between their variations and autism traits has not been established. The aim of this study is to predict the functional importance of putative risk genes at the whole-genome level based on the analysis of gene co-expression with a group of high confidence ASD risk genes (hcASDs). We found that three gene features, including gene size, mRNA abundance, and guanine-cytosine content, affect the genome-wide co-expression profiles of hcASDs. To circumvent the interference of these gene features on gene co-expression analysis (GCA), we developed a method to determine whether a gene is significantly co-expressed with hcASDs independent of confounding gene features by statistically comparing the co-expression of this gene with hcASDs to that of this gene with permuted gene sets of feature-matched genes. This method is referred to as "matched-gene co-expression analysis" (MGCA). MGCA demonstrated the convergence in the developmental expression profiles of hcASDs and improved the efficacy of risk gene prediction. Further analysis of two recently reported ASD candidate genes *CDH11* and *CDH9* suggested the involvement of *CDH11*, but not *CDH9*, in ASD. Consistent with this prediction, behavioral studies showed that *Cdh11*-null mice, but not *Cdh9*-null mice, have multiple autistic-like behavioral alterations. This study highlighted the power of MGCA in revealing functionally important ASD-relevant genes and suggested an important role of *CDH11* in ASD. 

## What you'll find here
* `Cross-Correlation.pl`: This script is designed to calculate the cross-correlation of gene expression.
* `Matrix.pl`: This script is designed to get the matrix of correlation coefficient between two interested gene sets.
* `200_mRand-Gene.pl`: This script is designed to get 200 random genesets of matched gene expression level, gene length or GC content relative to the geneset of interest and calculate their internal CC, CC to the input geneset and CC to whole genes. 
* `200_Rand-Gene.pl`: This script is designed to get 200 random genesets and calculate their internal CC and CC to the input geneset. The input files of this script are AllGene list file, the input gene list file, and the correlation matrix dataset. 
* `100000_mRand-Gene.pl`: This script is designed to get 100000 random genesets of matched gene expression level, gene length or GC content relative to the geneset of interest and calculate their internal CC and CC to the input geneset. 
* `100000_Rand-Gene.pl`: This script is designed to get 100000 random genesets and calculate their internal CC and CC to the input geneset. 
* `MGCA_Range50.pl`: This script is designed to calculate P-value of gene by MGCA under gene expression level, gene length or GC content ranked.
* `FDR.pl`: This script is designed to calculate FDR of whole genes by MGCA.
* `Heatmap_Cluster.R`:This script is designed to get the heatmap of the matrix of correlation coefficient between two interested gene sets.
* `hcASD_Gene`: high confidence ASD risk gene set.
* `mRNA-Abundance_Ranked_Gene`: Whole gene set ranked by gene expression.
* `gDNA-Size_Ranked_Gene`: Whole gene set ranked by gene size.
* `GC-Content_Ranked_Gene`: Whole gene set ranked by GC content.
* `mRNA-Abundance_Ranked_Data`: The human brain transcriptome dataset ranked by gene expression.
* `gDNA-Size_Ranked_Data`: The human brain transcriptome dataset ranked by gene size.
* `GC-Content_Ranked_Data`: The human brain transcriptome dataset ranked by GC content.

## Questions/Issues
Please direct any questions or issues with the code published here to `xbyuan@brain.ecnu.edu.cn`
