# Human Cerebrospinal fluid affects chemoradiotherapy sensitivities in tumour cells from glioblastoma patients.

### Data Access and Information

The unprocessed (fastq) and processed (counts, features, barcodes) transcriptomic data of this paper can be found at GSE243501 on GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243501.

This dataset contains 10 patient cell lines (111, 134, 159, 448, 468, 469, 497, BAH1, MN1 and HW1). Cell lines were cultured in either GM (TME) or CSF for three days prior to sequencing. Transcriptomic data for this paper is named patient-line_condition (e.g. 111_CSF, 111_TME). 'TME' is synonymous with 'GM' mentioned in the paper.

**Dates of data collection:** 2019-2021

**Location:** Adelaide, South Australia

### Sample provision
**South Australian Neurological Tumour Bank (SANTB):** Patient GBM tumour samples (111, 134, 159, 448, 468, 469, 497) and human cerebrospinal fluid.

**Queensland Institute of Medical Research:** Patient GBM tumour samples (BAH1, MN1 and HW1).

**RNA sequencing:** The South Australian Genomics Facility and Australian Genome Research Facility.

### List of files within this repository
1. Rscripts
  - GBM_CSF_QualityControl.Rmd
  - GBM_CSF_Paper_Analysis.Rmd
  - SANTB00134_GM_NM_CSF.Rmd
2. data
   - RData/metadata sheets
        - GBM_MetaData.csv
        - GBM_MetaData_Downsampled.csv
        - Prefilter_GBM_MetaData.csv
   - reference sheets
      - housekeepers.txt
      - cellcycle_tirosh.tsv
      - top50quiescent_atkins.txt
      - IDHwt.GBM.MetaModules.tsv
      - Stemcell.Markers.csv
      - TCGA.Verhaak.GBM.txt

### Data Description ###
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
   
3. ****GBM_CSF_Paper_Analysis:**** Analysis of transcriptomic data and code required to generate figures of the paper. 
   
4. ****SANTB00134_GM_NM_CSF:**** Comparisons between GM, NM and CSF for patient SANTB00134 (for Supplementary figure 15m). This script also contains quality control parameters for SANTB00134. 

### **data:** 
1. RData/metadata sheets: contains CSV files of transcriptomic data before and after pre-processing.

   **Information regarding individual files:** 
   - *Prefilter_GBM_MetaData.csv:* This file contains the information for ALL cells prior to quality control. Please note: This file was compiled for illustrative purposes to generate figures for Supplementary Fig 5. NO processing was performed on this dataset. Quality control was performed separately for each *patient-line_condition* before being merged for downstream analysis.  
        - Number of variables: 10
        - Number of cases/rows: 67,817
        - Variable List: see below
        - Linked Rscript: GBM_CSF_QualityControl
          
   - GBM_MetaData.csv: This file contains the information for cells post-quality control. Quality control was performed separately for each *patient-line_condition* and also downsampled to contain equal numbers of cells in each condition (TME and CSF). 
        - Number of variables: 56
        - Number of cases/rows: 37,466
        - Variable List: see below
        - Linked Rscript: GBM_CSF_QualityControl, GBM_CSF_Paper_Analysis

   - GBM_MetaData_Downsampled.csv: This file contains the information of cells randomly downsampled from those in 'GBM_MetaData.csv'. The Seurat object was downsampled to contain 700 cells across each *patient-line_condition* (e.g. 700 cells in TME and CSF per patient or a total of 1400 cells per patient). This downsampled data was used to conduct differential expression analyses. 
        - Number of variables: 56
        - Number of cases/rows: 14,000
        - Variable List: see below
        - Linked RScript: GBM_CSF_Paper_Analysis

***Information regarding variables in metadata: These metadata files contain information for each cell (rows) and contain common variables (columns) which are listed below. 'TME' is synonymous with 'GM' mentioned in the paper. 'Linked files' are files used to calculate content within the columns.*** 


***NOTE: The metadata presented here has been presented in the associated paper. If the data is reanalysed using provided fastqs or processed data please note that values in the metadata maybe subject to change as;***
***(1) Some variables within these metadata objects are automatically generated when creating a Seurat Object and are indicated by '(Seurat)'.***
***(2) The content within some columns is also subject to change depending on the analysis conducted and the parameters used. These columns are denoted by '^'.***

  **Variables across ALL metadata files:**
  - orig.ident (Seurat): the identity of each cell
  - nCount_RNA (Seurat): total number of reads per cell
  - nFeature_RNA (Seurat): total number of genes per cell
  - PatientID: patient ID. The prefix 'SANTB' was opted for lines obtained from SANTB. All other lines are from the QIMR.
  - Condition: culture condition for each cell. Either 'TME' or 'CSF'
  - Condition_PatientID: Combined PatientID and culture condition. 
  - percent.mt: percentage of mitochondrial reads
  - scaled.percent.mt: scaled percentage of mitochondrial reads
  - log.percent.mt: log transformed percentage of mitochondrial reads
  - n.exp.hkgenes: number of housekeeping genes per cell
       - Linked files: housekeepers.txt (see folder: reference sheets)

 **Variables across GBM_MetaData and GBM_MetaData_Downsampled files:**
- RNA_snn_res.0.05 (Seurat)^: shared neighbourhoods calculated through seurat 'FindNeighbours' at a resolution of 0.05. See RScripts for more details.
- seurat_clusters (Seurat)^: Seurat clusters calculated through Seurat 'FindClusters'. 
- MESlike.Score^: Mesenchymal (MES) score calculated using Seurat 'AddModuleScore'. 
     - Linked files: IDHwt.GBM.MetaModules.tsv (see folder: reference sheets).
- AClike.Score^: Astrocyte (AC) score calculated using Seurat 'AddModuleScore'. 
     - Linked files: IDHwt.GBM.MetaModules.tsv (see folder: reference sheets).
- OPClike.Score^: Oligodendrocyte precursor (OPC) score calculated using Seurat 'AddModuleScore'. 
     - Linked files: IDHwt.GBM.MetaModules.tsv (see folder: reference sheets)
- NPClike.Score^: Neural progenitor (NPC) score calculated using Seurat 'AddModuleScore'. 
     - Linked files: IDHwt.GBM.MetaModules.tsv (see folder: reference sheets)
- MaxScore_OPC.NPC: For each cell, the maximum score between the OPC and NPC scores. 
- MaxScore_AC.MES: For each cell, the maximum score between the AC and MES scores. 
- D1: For each cell, difference between MaxScore_AC.MES and MaxScore_OPC.NPC. 
- combGroup1:  Either AC_MES or OPC_NPC based on D1 > 0 or < 0
- MaxScore_AC.OPC: For each cell, the maximum score between the AC and OPC scores. 
- MaxScore_MES.NPC: For each cell, the maximum score between the MES and NPC scores. 
- D2: For each cell, difference between MaxScore_AC.OPC and MaxScore_MES.NPC.
- combGroup2: Either AC_OPC or MES_NPC based on D2 > 0 or < 0
- rel.metamodule.y.score: for AC_MES and OPC_NPC cells, -log2(AClike.Score - MESlike.Score) or log2(OPClike.Score - NPClike.Score) values. See RScripts for full details on analysis. 
- rel.metamodule.x.score: for AC_OPC and MES_NPC cells, -log2(AClike.Score - OPClike.Score) or log2(MESlike.Score - NPClike.Score) values. See RScripts for full details on analysis. 
- Quadrant: Classification of cells as MESlike, AClike, NPClike or OPClike based on rel.metamodule.y.score and rel.metamodule.x.score.
- MESlike1.Score^: MESlike1 score calculated using seurat 'AddModuleScore'. 
     - Linked files: IDHwt.GBM.MetaModules.tsv (see folder: reference sheets)
- MESlike2.Score^: MESlike2 score calculated using seurat 'AddModuleScore'. 
     - Linked files: IDHwt.GBM.MetaModules.tsv (see folder: reference sheets)
- NPClike1.Score^: NPClike1 score calculated using seurat 'AddModuleScore'. 
     - Linked files: IDHwt.GBM.MetaModules.tsv (see folder: reference sheets)
- NPClike2.Score^: NPClike2 score calculated using seurat 'AddModuleScore'.
     - Linked files: IDHwt.GBM.MetaModules.tsv (see folder: reference sheets)
- MESlike: classification of each MESlike cell as MESlike1, MESlike2 or notMESlike for those cells previously classified as AC, NPC or OPClike. 
- NPClike: classification of each NPClike cell as NPClike1, NPClike2 or notNPClike for those cells previously classified as AC, MES or OPClike. 
- G1S.Score^: G1S score calculated using seurat 'AddModuleScore'. 
     - Linked files: cellcycle_tirosh.tsv (see folder: reference sheets)
- G2M.Score^: G2M score calculated using seurat 'AddModuleScore'. 
     - Linked files: cellcycle_tirosh.tsv (see folder: reference sheets)
- Proliferation.Score^: Average of G1S and G2M scores.
- Cell_Cycle: Classification of cells as cycling or non-cycling based on proliferation scores.
- CellCycle_Condition: Combination of Cell_Cycle and Condition. 
- CellCycleDifference.Score^: Difference between G1S and G2M scores.
- Quiescent.Score: Quiescent score calculated using seurat 'AddModuleScore'. 
     - Linked files: top50quiescent_atkins.txt (see folder: reference sheets)
- Quiescence: Classification of cells as quiescent or non-quiescent based on quiescence scores.
- TCGA_Mesenchymal^: MES score calculated using Seurat 'AddModuleScore'. 
     - Linked files: TCGA.Verhaak.GBM.txt (see folder: reference sheets)
- TCGA_Classical^: Classical (CL) score calculated using Seurat 'AddModuleScore'. 
     - Linked files: TCGA.Verhaak.GBM.txt (see folder: reference sheets)
- TCGA_Proneural^: Proneural (PN) score calculated using Seurat 'AddModuleScore'. 
     - Linked files: TCGA.Verhaak.GBM.txt (see folder: reference sheets)
- Condition_Quadrant: Combination of Condition and Quadrant. 
- MaxScore_MES.CL: For each cell, the maximum score between MES and CL scores. 
- MaxScore_CL.PN: For each cell, the maximum score between CL and PN scores. 
- MaxScore_MES.PN: For each cell, the maximum score between MES and PN scores. 
- TCGA: Classification of cells according to the TCGA subtypes: Classical, Mesenchymal, Proneural.
- Lathia_GSC.Score^: Glioma stem cell (GSC) scores calculated using gene lists obtained from Lathia et al, 2015.
     - Linked files: Stemcell.Markers.csv (see folder: reference sheets)
- Richards_Dev_GSC.Score^: Developmental GSC scores calculated using gene lists obtained from Richards et al, 2021.
     - Linked files: Stemcell.Markers.csv (see folder: reference sheets)
- Richards_Inj_GSC.Score^: Injury response GSC scores calculated using gene lists obtained from Richards et al, 2021.
     - Linked files: Stemcell.Markers.csv (see folder: reference sheets)
- Lathia_GSC: Classification of cells as GSC or non-GSC based on Lathia_GSC.Score
- Richards_GSC: Classification of cells as Dev-GSC or Inj-GSC based on Richards_Dev_GSC.Score and Richards_Inj_GSC.Score
- Tirosh_GSC.Score^: GSC scores calculated using gene lists obtained from Tirosh et al, 2016.
     - Linked files: Stemcell.Markers.csv (see folder: reference sheets)
- Tirosh_GSC: Classification of cells as GSC or non-GSC based on Tirosh_GSC.Score.

***NOTE: Detailed analysis methods (score calculations, classifications etc.) can be found within GBM_CSF_Paper_Analysis.Rmd. All scores are arbitrary values and therefore have no units***

2. reference sheets: contains files used to generate figures.

  ***Gene marker files:***
   - *housekeepers.txt:* Housekeeping genes used for quality control.
      - Reference study: Tirsoh et al, 2016 (https://www.science.org/doi/abs/10.1126/science.aad0501)
      - Variables: 1
      - Number of cases/rows: 98
      - Linked RScript: GBM_CSF_QualityControl.Rmd
   
   - *cellcycle_tirosh.tsv:* Cell cycle genes used to calculate proliferation scores.
      - Reference study: Tirosh et al, 2016 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5465819/)
      - Variables: 2
      - Number of cases/rows: 55
      - Variable list:
           - G1/S (Gap 1/Synthesis): Genes defining the G1S cell cycle phase. 
           - G2/M (Gap 2/Mitosis): Genes defining the G2M cell cycle phase. 
      - Linked RScript: GBM_CSF_Paper_Analysis.Rmd
        
   - *top50quiescent_atkins.txt:* Quiescent gene lists used to calculate quiescence scores.
      - Reference study: Atkins et al, 2019 (https://www.sciencedirect.com/science/article/abs/pii/S0014482718311698)
      - Variables: 1
      - Number of cases/rows: 50
      - Linked RScript: GBM_CSF_Paper_Analysis.Rmd
   
   - *IDHwt.GBM.MetaModules.tsv:* MES, AC, NPC and OPC genes from Neftel et al, 2019.
      - Reference study: Neftel et al, 2019 ()
      - Variables: 8
      - Number of cases/rows: 50
      - Variable list:
           - MESlike2: Genes defining the MESlike2 cell state
           - MESlike1: Genes defining the MESlike1 cell state
           - AClike: Genes defining the AClike cell state
           - OPClike: Genes defining the OPClike cell state
           - NPClike1: Genes defining the NPClike1 cell state
           - NPClike2: Genes defining the NPClike2 cell state
           - G1/S: Genes defining the G1S cell state. This gene list is the same as in cellcycle_tirosh.tsv.
           - G2/M: Genes defining the G2M cell state. This gene list is the same as in cellcycle_tirosh.tsv.
      - Linked RScript: GBM_CSF_Paper_Analysis.Rmd
        
   - *Stemcell.Markers.csv:* Stem cell marker lists for Supplementary Figure 6e. Only Lathia gene lists were used for this figure.
     - Reference studies: Lathia et al, 2015 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4495393/), Richards et al, 2021 (https://www.nature.com/articles/s43018-020-00154-9), Tirosh et al, 2016 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5465819/)
     - Variables: 4
     - Number of cases/rows: 100
     - Variable list:
        - GSC_Lathia: GSC genes from Lathia et al, 2015
        - Richards_2021_DevelopmentalGSC: Developmental GSC genes from Richards et al, 2021
        - Richards_2021_InjuryResponseGSC: Injury response GSC genes from Richards et al, 2021
        - GSC_Tirosh: GSC genes from Tirosh et al, 2016
     - Linked RScript: GBM_CSF_Paper_Analysis.Rmd
       
   - *TCGA.Verhaak.GBM.txt:* Neural genes from Verhaak et al, 2010. 
     - Reference study: Verhaak et al, 2010 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2818769/)
     - Variables: 4
     - Number of cases/rows: 129
     - Variable list:
        - Neural: Genes defining the Neural subtype
        - Proneural: Genes defining the PN subtype
        - Mesenchymal: Genes defining the MES subtype
        - Classical: Genes defining the CL subtype
      - Linked RScript: GBM_CSF_Paper_Analysis.Rmd
        
  ***other files:*** 

   *TME_CSF_Treatment.csv:* Average cell viability data of 25 patient-derived cell lines following treatment.  
   - Variables: 4
   - Number of cases/rows: 25
   - Variable list:
       -  TME TMZ: Average cell survival of cells cultured in TME (GM) for three days before treatment with 100 µM temozolomide (TMZ) for 7 days (percent). 
       -  TME Irradiation: Average cell survival of cells cultured in TME for three days before treatment with 5 fractions of 2 Gy irradiation (Irr) for 5 days (percent). 
       -  TME TFP	TME TMZ + Irr: Average cell survival of cells cultured in TME for three days before treatment with 10 µM TFP for 24 hours (percent). 
       -  TME TMZ + Irr +TFP: Average cell survival of cells cultured in TME for three days before combined treatment with 100 µM temozolomide (TMZ) for 7 days, 5 fractions of 2 Gy irradiation (Irr) for 5 days and 10 µM TFP for 24 hours (percent). 
       -  CSF TMZ: Average cell survival of cells cultured in CSF for three days before treatment with 100 µM temozolomide (TMZ) for 7 days (percent). 
       -  CSF Irradiation: Average cell survival of cells cultured in CSF for three days before treatment with 5 fractions of 2 Gy irradiation (Irr) for 5 days (percent). 
       -  CSF TFP	TME TMZ + Irr: Average cell survival of cells cultured in CSF for three days before treatment with 10 µM TFP for 24 hours (percent). 
       -  CSF TMZ + Irr +TFP: Average cell survival of cells cultured in CSF for three days before combined treatment with 100 µM temozolomide (TMZ) for 7 days, 5 fractions of 2 Gy irradiation (Irr) for 5 days and 10 µM TFP for 24 hours (percent). 
   - Linked RScript: GBM_CSF_Paper_Analysis.Rmd

### Code and Software

Sequenced BCL files were processed using standard CellRanger 6.1.1 pipelines. 

Requirements for downstream data analysis: 
- R
- RStudio
- Seurat

*All other required R packages are listed within RScripts.* 




    
