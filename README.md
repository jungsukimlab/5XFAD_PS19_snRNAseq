# 5XFAD_PS19_snRNAseq

Code associated with Park et al., 2026, *Alzheimer's & Dementia*.

**"Unique transcriptomic alterations in 5XFAD;PS19 mouse model identify glial lipid dysregulation and coordinated microglial–oligodendrocyte responses."**
https://alz-journals.onlinelibrary.wiley.com/doi/10.1002/alz.71742 

## Study overview
<img width="300" height="330" alt="image" src="https://github.com/user-attachments/assets/06ea3da4-a4de-4088-b2a4-e1b0de2474f4" />
<img width="300" height="330" alt="image" src="https://github.com/user-attachments/assets/67b69a16-f250-48bc-839a-2ae594cd41c8" />


We performed sex-balanced single-nucleus RNA sequencing (snRNA-seq) of cortical tissue from four mouse genotypes: wild-type (WT), 5XFAD, PS19, and combined 5XFAD;PS19 mice.

The analysis evaluates cell type-specific transcriptional changes associated with amyloid beta (Aβ) pathology, tau pathology, and combined Aβ-tau pathology. Downstream analyses include cell composition, gene co-expression network analysis, differential gene expression, cell–cell communication analysis, and cross-species alignment with human Alzheimer's disease datasets from AMP-AD.

## Analysis workflow

The scripts are organized in the order of the analysis workflow.

| Script                      | Analysis                                                                |
| --------------------------- | ----------------------------------------------------------------------- |
| `Step01_preprocessing.R`    | Ambient RNA correction, quality control, and doublet detection          |
| `Step02_integration.R`      | Integration of sample-level snRNA-seq datasets                          |
| `Step03_clustering.R`       | Clustering and cluster-resolution selection                             |
| `Step04_cluster_markers.R`  | Identification of cluster marker genes                                  |
| `Step05_annotation.R`       | Cell-type and cluster annotation                                   |
| `Step06_cell_composition.R` | Cell composition analysis across genotypes                              |
| `Step07_DEG.R`              | Cluster–specific differential gene expression analysis                |
| `Step08_DEG_downstream.R`   | Downstream analysis of differentially expressed genes, including enrichment analysis                   |
| `Step09_hdWGCNA.R`          | High-dimensional weighted gene co-expression network analysis (hdWGCNA) |
| `Step10_CellChat.R`         | Cell–cell communication analysis using CellChat                         |
| `Step11_DEG_pseudobulk.R` | Pseudobulk differential expression analysis across genotypes |
| `Step12_cross_species.R` | Cross-species transcriptomic alignment with AMP-AD human AD datasets |

## Data

The snRNA-seq data generated in this study are available through the Gene Expression Omnibus (GEO) under accession **GSE304499**.

Additional datasets used for cross-species analyses are available through Synapse:

- AMP-AD consensus co-expression modules: **syn11932957** 
- AMP-AD human RNA-seq differential expression data: **syn14237651**
- Mouse-human ortholog mappings: **syn17010253**
- AD biodomain definitions and labels: **syn25428992**, **syn26856828**

## Citation

If you use this code or dataset, please cite:

Park JH et al. Unique transcriptomic alterations in 5XFAD;PS19 mouse model identify glial lipid dysregulation and coordinated microglial–oligodendrocyte responses. *Alzheimer's & Dementia*. 2026;22(8):e71742. doi:10.1002/alz.71742
