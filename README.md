# Human Cerebrospinal fluid affects chemoradiotherapy sensitivities in tumour cells from glioblastoma patients.

### Data Access and Information

The unprocessed (fastq) and processed (counts, features, barcodes) transcriptomic data of this paper can be found at GSE243501 on GEO. 

### Data Description
### **Rscripts:**
This folder contains all the R scripts used to process transcriptomic data for this project. This includes;

1. ****GBM_CSF_QualityControl:****

   Quality Control parameters of single cell transcript data.

    ***NOTE: This script only contains quality control parameters for a single patient line. Each patient-line_condition file was processed individually before merging into a single Seurat object.***

    Standard quality control parameters include:
   - mitochondrial reads (exclude cells with > 20%)
   - housekeeping genes (exclude cells with < 70)
   - genes per cell (depending on distribution across each dataset. keep cells > 2000 genes per cell for all datasets except 159 (> 800 genes) and 448 (> 500 genes).
   - reads per cell (depending on distribution across each dataset. keep cells with < 55000 reads per cell for all datasets except MN1 (< 80000 reads).
   
3. ****GBM_CSF_Analysis:**** Analysis of transcriptomic data and code required to generate figures of the paper. 
   
4. ****SANTB00134_GM_NM_CSF:**** Comparisons between GM, NM and CSF for patient SANTB00134 (for Supplementary figure 15m). This script also contains quality control parameters for SANTB00134. 

### **Data:** 
1. Metadata: contains CSV files of transcriptomic data before and after pre-processing.
2. Reference sheets: contains files used to generate figures.

  ***Gene marker files:***
   - *housekeepers.txt:* Housekeeping gene lists used for quality control.
   
   -  *cellcycle_tirosh.tsv:* Cell cycle genes used to calculate proliferation scores.
   
   -    *top50quiescent_atkins.txt:* Quiescent gene lists used to calculate quiescence scores.
   
   -  *IDHwt.GBM.MetaModules.tsv:* MES, AC, NPC and OPC marker lists from Neftel et al, 2019.
   
   -  *Stemcell.Markers.csv:* Stem cell marker lists for Supplementary Figure 6e. Only Lathia gene lists were used for this figure. 

  ***Other files:*** 

   *TME_CSF_Treatment.csv:* Average cell viability data of twenty five cell lines following treatment (temzolomide (TMZ), irradiation (Irr) and trifluoperazine (TFP).  

### Code and Software

Sequenced BCL files were processed using standard CellRanger 6.1.1 pipelines. 

This dataset contains 10 patient cell lines (111, 134, 159, 448, 468, 469, 497, BAH1, MN1 and HW1). Cell lines were cultured in either GM (TME) or CSF for three days prior to sequencing. Transcriptomic data for this paper is named patient-line_condition (e.g. 111_CSF, 111_TME).

Requirements for downstream data analysis: 
- R
- RStudio
- Seurat

*All other required R packages are listed within Rscripts.* 

***NOTE: 'TME' is synonymous with 'GM' mentioned in this paper'.*** 
    
