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

### [Mass Spectrometry Raw Data Analysis by LIMMA multiple T-test.](./Limma_t_Test_Analysis)
The files in this folder contain the R scripts and raw data files required to identify the DSB-sensor interactions and quantify its Fold Change against the background signals. The Mass Spectrometry intensity values of the entire SensorID experiment were split into each DSB-Sensor proteins (Ku80, Mre11, and Sirt6) and each file contain the three independent biological replicates of the two negative controls (Cneg_IR5, Cneg_IR24) and the 7 experimental conditions (nonIR5, IR5, IR30, IR2h, IR8h, IR24h, and nonIR24h). The raw data tables were already manually filtered to remove those proteins that were detected by less than 2 peptides.
The raw data and the R script used to perform the LIMMA multiple T-test analysis against the negative controls are contained in the (place here the folder) section of this repository.
The Mass spectrometry intensity values of the SensorID shNucleolin experiments are contained in a single dataset, and each column is named by the DSB-sensor and the time point. The raw data tables were already manually filtered to remove those proteins that were detected by less than 2 peptides. The raw data and the R script used to perform the LIMMA multiple T-test analysis against the negative controls are contained in the (place here the folder) section of this repository.

**SensorID LIMMA multiple T-Test analysis**
The section Limma T-Test contains all the R scripts and files to analyze the SensorID experiments' 
raw data for mass spectrometry. The section is divided into two parts with independent sets of files:
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
 
### [Dynamic clustering by Fuzzy C Means Analysis.](./Fuzzy_C_Means_Analysis)
The files in this folder contain the R-script and the files required to cluster the DSB-sensor interactors into dynamic Fuzzy C means clusters using the previously calculated Fold Changes as a normalized measurement of its binding intensities. The protein interactors are clustered depending on their interaction dynamics; however, since Fuzzy C means clustering gives different results on each run, it is possible that the obtained results do not follow the numbering system used in the original paper, but the structure and behaviour of the clusters should be the same. This section of the repository contains the following files:
**Fuzzy_C_Means_Analysis**
The section Fuzzy_C_Means_Analysis contains all the R scripts and files to analyze the dynamic clustering by the Fuzzy C means algorithm:
 - 1. **Fuzzy_C_Means_Analysis**, containes the next files:
    - + **SensorID_FuzzyCMeans_Analysis.R**: R Script for analyzing the SensorID dynamic clusters using the mFuzz clustering package.
    - + **Sirt6_Data.csv**: Sirt6ID Fold Change values per time point, Fold Change calculated against the Negative controls.
    - + **Ku80_Data.csv**: Ku80ID Fold Change values per time point, Fold Change calculated against the Negative controls.
    - + **Mre11_Data.csv**: Mre11ID Fold Change values per time point, Fold Change calculated against the Negative controls.
    - + **Mre11ID_RawData.csv**: Mre11ID mass spectrometry Intensity raw data, each sample contains 3 biological replicates.
    - + **ColAnnoDF.csv**: SensorID dynamic interactomes Sample metadata.
    - + **Protein_Annotations.csv**: SensorID interactions annotations. Mannually filtered and curated based on Gene Ontology annotations.

### [SensorID dynamic interactome Network analysis with igraph](./SensorID_Network_Analysis)
This repository section encompasses the network analysis of the SensorID dynamic interactomes. The network analysis was done using the functions contained in the igraph R package. However, the visualization of the network was performed on Cytoscape due to its multiple layout software and easier mapping of variables in the network graphs. Notably, each network is uploaded as a separate edge list. Each list was retrieved from the String database by quering all the proteins contained in the specific time point interactome and filtering the physical interactions by a threshold of 0.4 confidence. The SensorID matrix file contains the node metadata, including the functional annotations and the fold change at specific time points.
This section of the repository contains the following files:
**SensorID_Network_Analysis**
The section SensorID_Network_Analysis contains all the R scripts and files for the network analysis of the SensorID dynamic interactomes:
 - 1. **SensorID_Network_Analysis**, containes the next files:
    - + **SensorID_Network_Analysis.R**: R Script for analyzing the SensorID dynamic interactomes using the functions of the igraph R package.
    - + **SensorID_Matrix.csv**: Matrix containing the metadata of the SensorID interactomes, it contains the name and annotation of each interaction, as well as its Fold Change versus the negative controls at each studied time point.
    - + **S6_"Time.Point"_Edges.csv**: Family of files containing the Edge list of the Sirt6ID interactomes at each specific time point.
    - + **Ku_"Time.Point"_Edges**: Family of files containing the Edge list of the Ku80ID interactomes at each specific time point.
    - + **Mre_"Time.Point"_Edges**: Family of files containing the Edge list of the Mre11ID interactomes at each specific time point.

### [SensorID shNcl Experiment Network analysis with igraph](./shNcl_Network_Analysis)
This repository section encompasses the network analysis of the SensorID IR2h interactomes after the silencing of Nucleolin using a specific shRNA (shNcl). The comparisons are done against the SensorID IR2h interactome treated with a control shRNA (shScr). The network analysis was done using the functions contained in the igraph R package. However, the visualization of the network was performed on Cytoscape due to its multiple layout software and easier mapping of variables in the network graphs. Notably, each network is uploaded as a separate edge list. Each list was retrieved from the String database by quering all the proteins contained in the specific time point interactome and filtering the physical interactions by a threshold of 0.4 confidence. The SensorID shNclExp Matrix file contains the node metadata, including the functional annotations and the fold change versus the negative controls of the experiments. It also contains the statistical test of the shNcl condition versus the shScr condition.
This section of the repository contains the following files:
**shNcl_Network_Analysis**
The section shNcl_Network_Analysis contains all the R scripts and files for the network analysis of the SensorID shNcl experiment interactomes:
 - 1. **shNcl_Network_Analysis**, containes the next files:
    - + **shNcl_Network_Analysis.R**: R Script for analyzing the SensorID shNcl experiment interactomes using the functions of the igraph R package.
    - + **SensorID_shNclExp_Matrix.csv**: Matrix containing the metadata of the SensorID shNcl and SensorID shScr interactomes. It contains the name and annotation of each interactor, as well as its Fold Change versus the negative controls.
    - + **SensorID_shNclExp_CompleteNetwork_EdgeList.tsv**: Edge list that includes all the interactions from the SensorID shNcl experiment. It is similar to the Full Nuclear Network (FNN) used in the SensorID dynamic interactome network analysis.
    - + **Sirt6ID_shScrlvsCneg_Network_EdgeList.tsv**: Edge list of the Sirt6ID IR2h treated with the control shRNA (shScr).
    - + **Sirt6ID_shNclvsCneg_Network_EdgeList.tsv**: Edge list of the Sirt6ID IR2h treated with the Nucleolin-targeting shRNA (shNcl).
    - + **Ku80ID_shScrlvsCneg_Network_EdgeList.tsv**: Edge list of the Ku80ID IR2h treated with the control shRNA (shScr).
    - + **Ku80ID_shNclvsCneg_Network_EdgeList.tsv**: Edge list of the Ku80ID IR2h treated with the Nucleolin-targeting shRNA (shNcl).
    - + **Mre11ID_shScrlvsCneg_Network_EsgeList.tsv**: Edge list of the Mre11ID IR2h treated with the control shRNA (shScr).
    - + **Mre11ID_shNclvsCneg_Network_EdgeList.tsv**: Edge list of the Mre11ID IR2h treated with the Nucleolin-targeting shRNA (shNcl).

### [SensorID dynamic interactome rNetcarto Modularity Analysis](./SensorID_Rnetcarto_Modularity_Analysis)
This repository section contains the SensorID networks edgelists, the node data frame, and the R Script required to perform the modularity analysis and node classification using the Rnetcarto algorithm. The modularity analysis of the PPI networks was done using the simulated annealing Rnetcarto modularity algorithm, and calculating the connectivity (within-module degree) and the participation coefficient, which are used to further distinguish network structural roles. The connectivity measures the number of connections of a particular node within its module. On the other hand, the participation coefficient depicts the proportions of a particular node's connections with other modules in the network.
The results were later manually verified and merged with previous results to generate a comprehensive data frame that contains all the analysis results. The files contained in this folder are:
**SensorID_Rnetcarto_Modularity_Analysis**
The section SensorID_Rnetcarto_Modularity_Analysis contains all the R scripts and files for the network modularity analysis of the SensorID interactomes:
 - 1. **SensorID_Rnetcarto_Modularity_Analysis**, containes the next files:
    - + **SensorID_rNetcarto_ModularityAnalysis.R**: R Script for analyzing the SensorID interactomes using the functions of the Rnetcarto and igraph R packages.
    - + **SensorID_Node_DataFrame.csv**: Matrix containing the metadata of the SensorID interacting proteins. It contains the name and annotation of each interactor, as well as its Fold Change versus the negative controls at each studied time point, and the modularity analysis results obtained.
    - + **S6_"Time.Point"_Edges.csv**: Edge lists of the Sirt6ID networks. Each file contains the edge lists of each time point; the full network refers to the edge list of all the interactors found in the Sirt6ID networks.
    - + **Ku_"Time.Point"_Edges.csv**: Edge lists of the Ku80ID networks. Each file contains the edge lists of each time point; the full network refers to the edge list of all the interactors found in the Ku80ID networks.
    - + **Mre_"Time.Point"_Edges.csv**: Edge lists of the Mre11ID networks. Each file contains the edge lists of each time point; the full network refers to the edge list of all the interactors found in the Mre11ID networks.

### [SensorID shNcl Experiment Rnetcarto modularity analysis](./SensorID_shNcl_Experiment_Rnetcarto_Modularity_Analysis)
This repository section is intended to contain the mass spectrometry results, the edge list, and the R script to identify the structural modules of the SensorID protein-protein interacting networks after Nucleolin silencing. The modularity analysis of the PPI networks was done using the simulated annealing Rnetcarto modularity algorithm, and calculating the connectivity (within-module degree) and the participation coefficient, which are used to further distinguish network structural roles.
The results from the Rnetcarto modularity analysis were verified and manually merged with previous analysis tables to generate a data frame that contains a comprehensive description of each node in the network. 
This section of the repository contains the following files:
**SensorID_shNcl_Experiment_Rnetcarto_Modularity_Analysis**
The section SensorID_shNcl_Experiment_Rnetcarto_Modularity_Analysis contains the R scripts and files for the network modularity analysis of the SensorID shNcl experiment interactomes:
 - 1. **SensorID_shNcl_Experiment_Rnetcarto_Modularity_Analysis**, containes the next files:
    - + **SensorID_shNclExp_rNetcarto_Modularity_Analysis.R**: R Script for analyzing the SensorID shNcl experiment interactomes using the functions of the Rnetcarto and igraph R packages.
    - + **Sirt6ID_shNclExp_DataFrame.csv**: Matrix containing the metadata of the Sirt6ID shNcl and SensorID shScr interactomes. It contains the name and annotation of each interactor, as well as its Fold Change versus the negative controls.
    - + **Ku80ID_shNclExp_DataFrame.csv**: Matrix containing the metadata of the Ku80ID shNcl and SensorID shScr interactomes. It contains the name and annotation of each interactor, as well as its Fold Change versus the negative controls.
    - + **Mre11ID_shNclExp_DataFrame.csv**: Matrix containing the metadata of the Mre11ID shNcl and SensorID shScr interactomes. It contains the name and annotation of each interactor, as well as its Fold Change versus the negative controls.
    - + **Sirt6ID_shNclExp_FullNetwork_EdgeList.tsv**: Edge list of the Sirt6ID Full network containing all the interactions found in the control shRNA (shScr) and the Nucleolin targeting shRNA (shNcl).
    - + **Sirt6ID_shScrlvsCneg_Network_EdgeList.tsv**: Edge list of the Sirt6ID IR2h treated with the control shRNA (shScr).
    - + **Sirt6ID_shNclvsCneg_Network_EdgeList.tsv**: Edge list of the Sirt6ID IR2h treated with the Nucleolin-targeting shRNA (shNcl).
    - + **Ku80ID_shNclExp_FullNetwork_EdgeList.tsv**: Edge list of the Sirt6ID Full network containing all the interactions found in the control shRNA (shScr) and the Nucleolin targeting shRNA (shNcl).
    - + **Ku80ID_shScrlvsCneg_Network_EdgeList.tsv**: Edge list of the Ku80ID IR2h treated with the control shRNA (shScr).
    - + **Ku80ID_shNclvsCneg_Network_EdgeList.tsv**: Edge list of the Ku80ID IR2h treated with the Nucleolin-targeting shRNA (shNcl).
    - + **Mre11ID_shNclExp_FullNetwork_EdgeList.tsv**: Edge list of the Sirt6ID Full network containing all the interactions found in the control shRNA (shScr) and the Nucleolin targeting shRNA (shNcl).
    - + **Mre11ID_shScrlvsCneg_Network_EsgeList.tsv**: Edge list of the Mre11ID IR2h treated with the control shRNA (shScr).
    - + **Mre11ID_shNclvsCneg_Network_EdgeList.tsv**: Edge list of the Mre11ID IR2h treated with the Nucleolin-targeting shRNA (shNcl).

### [Final Notes]
The files and script published in this repository were used to generate the results published in the scientific article "SPARK: mapping the Sensor Proteome Alliance for Repair Kinetics during double-strand break (DSB) repair", for further clarification or discussion on the analysis published here please refer to Dr. Debra Toiber (<ins>toiber@bgu.ac.il</ins>) from the Ben Gurion University of the Neguev, Israel. 




