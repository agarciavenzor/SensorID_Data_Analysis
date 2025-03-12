# SPARK: mapping the Sensor Proteome Alliance for Repair Kinetics during double-strand break (DSB) repair

## Alfredo Garcia-Venzor<sup>1,2</sup>, Uri Grupel<sup>1,2</sup>, Shai-Kaluski-Kopach<sup>1,2</sup>, Miguel Lurgi<sup>3</sup>, and <ins>Debra Toiber<sup>2,1*</sup></ins> 

1. Department of Life Sciences, Ben-Gurion University of the Negev, Beer Sheva, 84105, Israel. 
2. The School of Brain Sciences and Cognition, Ben Gurion University of the Negev, Beer Sheva, 84105, Israel.
3. Department of Biosciences, Swansea University, Singleton Park, Swansea SA2 8PP, UK.

* Corresponding author: **Debra Toiber**, <ins>toiber@bgu.ac.il</ins>

## Scientific Abstract.

DNA double-strand breaks (DSB) are critical genomic lesions. DSBs are repaired by two main pathways: homologous recombination (HR) or non-homologous end joining (NHEJ). The activity of the sensor proteins at the DSB site has a major role in influencing the repair outcome. However, little is known about how each DSB sensor proteins may simultaneously integrate different nuclear processes to efficiently cope with the DNA lesions. In this study, we generated SPARK, the dynamic Protein-protein interaction (PPI) networks of three DSB-sensors (Sirt6, Ku80, and Mre11) during DSB-repair. Using proximity-labeling assays, we unravel the DSB-sensor dynamic-interactome. Our results show that each sensor has a unique interactome, which reacts to DSB induction and progressively expands to include new specific nuclear functions. Although the sensors rely on a common core of nuclear processes, sensor-specific interactors allow functional differences and confer robustness during DDR. Moreover, we discovered that PPI networks contain connector nodes, which coordinate the function of large nuclear processes occurring simultaneously. Altogether, our study shows how three DSB-sensor proteins coordinate and regulate broad nuclear processes through deep changes in the structure and composition of their interactome, relying on a small set of connector proteins to coordinate big macromolecular complexes.

## Summary of the SensorID data analysis repository.

The present Github repository is intended to serve as a tool for the public to access the computational analysis done to understand and generate the results published in the SPARK article. The information and R scripts contained in this repository will follow the same order as presented in the results section of the SPARK article, and it will include the original files as well as the code used to generate relevant results and analysis. Therefore, to understand the computational procedures contained in this repository, it is necessary to follow the results discussion of the main SPARK article. However, operations that were done in other software, such as Excel or Cytoscape, to manually curate or annotate the information contained in the manuscript are not contained in this repository.

Note that if the user needs to consult or review the analyzed results, the networks, or the node-level metrics calculations explained in the SPARK manuscript, this information is summarized in the [SPARK website](https://sparkid.bgu.ac.il/). The SPARK project website was developed as a tool for easily navigating the interactomic results in a protein or functional group-based manner and retrieving networks and metrics depending on the interests of the user.

Furthermore, the complete database containing all the obtained results and metrics for each protein-interactor in each one of the DSB-sensor interactomes is encompassed in the supplementary tables of the SPARK (SensorID) scientific article. Therefore, for a deeper query of the interactomic results, the authors recommend the public to refer to the original paper and use it as a guide for navigating the [SPARK website](https://sparkid.bgu.ac.il/) and the Github repository.

## Index of the data analysis scripts.

### [Mass Spectrometry Raw Data Analysis by LIMMA multiple T-test.](/SensorID_Data_Analysis/Limma_t_Test_Analysis/)
The Mass Spectrometry intensity values of the entire SensorID experiment were split into each DSB-Sensor proteins (Ku80, Mre11, and Sirt6) and each file contain the three independent biological replicates of the two negative controls (Cneg_IR5, Cneg_IR24) and the 7 experimental conditions (nonIR5, IR5, IR30, IR2h, IR8h, IR24h, and nonIR24h). The raw data tables were already manually filtered to remove those proteins that were detected by less than 2 peptides.
The raw data and the R script used to perform the LIMMA multiple T-test analysis against the negative controls are contained in the (place here the folder) section of this repository.
The Mass spectrometry intensity values of the SensorID shNucleolin experiments are contained in a single dataset, and each column is named by the DSB-sensor and the time point. The raw data tables were already manually filtered to remove those proteins that were detected by less than 2 peptides. The raw data and the R script used to perform the LIMMA multiple T-test analysis against the negative controls are contained in the (place here the folder) section of this repository.

**SensorID LIMMA multiple T-Test analysis**
The section Limma T-Test contains all the R scripts and files to analyze the SensorID experiments' 
raw data for mass spectrometry. The section is divided into two parts with independent set of files:
 - 1. **SensorID dynamic interactome analysis**, containes the next files:
    - + **SensorID_Limma_T_Test.R**: R Script for analyzing the SensorID dynamic interactomes mass spectrometry raw data.
    - + **Sirt6ID_RawData.csv**: Sirt6ID mass spectrometry Intensity raw data, each sample contains 3 biological replicates.
    - + **Ku80ID_RawData.csv**: Ku80ID mass spectrometry Intensity raw data, each sample contains 3 biological replicates.
    - + **Mre11ID_RawData.csv**: Mre11ID mass spectrometry Intensity raw data, each sample contains 3 biological replicates.
    - + **Mre11ID_RawData.csv**: Mre11ID mass spectrometry Intensity raw data, each sample contains 3 biological replicates.
    - + **ColAnnoDF.csv**: SensorID dynamic interactomes Sample metadata.
    - + **SensorID_FoldChange_Matrix.csv**: SensorID dynamic interactome Fold Change values versus the negative controls. Mannually filtered.
 - 2. **SensorID shNucleolin Experiments Interactome analysis**, contains the next files:
    - + **SensorID_shNclExp_Limma_T_Test.R**: R Script for analyzing the raw data of SensorID shNucleolin Experiments Interactome mass spectrometry.
    - + **SensorID_shNclExp_RawData_All.csv**: SensorID shNucleolin Experiments Interactome mass spectrometry Intensity raw data, each sample contains 3 biological replicates..
    - + **SensorID_shNclList_Diff_Prots.csv**: List of the Differentially bound proteins in the SensorID shNcl experiments and their annotations. Mannually curated.
    - + **ColAnnoDF_shNclExp.csv**: SensorID shNcl experiments Sample metadata.
 
###



