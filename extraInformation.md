# Additional information
Additional information regarding the purpose of each of the Notebooks and the Folder structure and content.  
## Notebook purpose
**01_inputData**  
- <ins>01_retrievingBAMfiles.ipynb</ins>: Provides information about how to treat FASQ files to obtain bam files suitable for TADbit.  
  
- <ins>02_saveMatrixFiles.ipynb</ins>: Contains the code to both obtain the biases for the PRINT normalisation and create the interaction matrices files. 
  
**02_modelling**  
- <ins>01_modellingParametersGridSearch.ipynb</ins>: Contains the code to run the optimization step of the modelling in order to select the best combination of input parameters for the final models.   
  
- <ins>02_chooseBestParameters.ipynb</ins>: Contains the code to use the output from the previous Notebook to select the best combination of input parameters for the final models. Depends on data from <ins>01_modellingParametersGridSearch.ipynb</ins>.  
  
- <ins>03_modelling.ipynb</ins>: Contains the code to compute the final set of models. Depends on data from <ins>02_chooseBestParameters.ipynb</ins>.  
  
**03_modelAnalysis**  
- <ins>01_convertTADdynModels_toTADbitModels.ipynb</ins>: Contains the code to convert TADdyn format models to TADbit format.  
   
- <ins>02_clusterModelsEnsemble.ipynb</ins>: Contains the code to cluster the models and store this information to be used in the next Notebooks.  
  
- <ins>03_compareModelMatrices.ipynb</ins>: Contains the code to compute and plot the contact matrices that are inferred from the models.  
  
- <ins>04_distanceBetweenPairsOfParticles.ipynb</ins>: Contains the code to obtain distances between sets of particles of interest and then display them in boxplot (Figure 3C from the manuscript). Depends on data from <ins>02_clusterModelsEnsemble.ipynb</ins>.  
  
- <ins>05_linePlot.ipynb</ins>: Contains the code to measure distances from a particle of interest with the rest of the particles of the model and to generate the lineal plot that we display in Figure 4A from the manuscript. Depends on data from <ins>02_clusterModelsEnsemble.ipynb</ins>.  
  
- <ins>06_dRMSDclusteringOfModelEnsembles.ipynb</ins>: Contains the code to cluster by dRMSD all the models from the same region (comparing even models from different cells) and to display the resulting clustering tree (Figure 3B from the manuscript).   
  
- <ins>07_radialPlot.ipynb</ins>: Contains the code to get the radialPlots (Figure 3D from the manuscript). Depends on data from <ins>02_clusterModelsEnsemble.ipynb</ins>.  
  
- <ins>08_distanceBtwTranscribedBins.ipynb</ins>: Contains the code to get and display the distances between all the particles containing transcribed genes (Figure 4B from the manuscript). Depends on data from <ins>02_clusterModelsEnsemble.ipynb</ins>.  
  
- <ins>09_co-occurrenceMatrix.ipynb</ins>: Contains the code to get and display the co-occurrence matrices (Figure 4C-E from the manuscript). Depends on data from <ins>02_clusterModelsEnsemble.ipynb</ins>.   
  
- <ins>10_communityAnalysis.ipynb</ins>: Contains the code to analyse the communities of genes obtained in <ins>09_co-occurrenceMatrix.ipynb</ins> in terms of distances and expression (Figure 4F,G and Supplementary Figure 6 from the manuscript). Depends on data from <ins>02_clusterModelsEnsemble.ipynb</ins> and <ins>09_co-occurrenceMatrix.ipynb</ins>.   
  
## Folder structure and content
Note: Avoid using '_' in any folder inside or in the root of the location of this repository in your computer.  

- "additionalInput" folder contains text files with additional data to be loaded, like enhancer and promoter coordinates, and methylation or gene expression data.  
  
- "bamfiles" folder contains one level of subfolders stating the cell ID. Inside ...bamfiles/cellID/ we should store its correspondent sorted .bam file and the associated index .bai file. This repository does not provide any bam file due to the large size of the files.  
  
- "code" folder contains sets of python 2 functions that will be using in the Notebooks.  
  
- "codeScripts" folder contains the scripts involving the Notebooks from "01_inputData" and "02_modelling" so that these steps can be run in a cluster.  
  
- "containers" folder contains the recipes for building the Singularity environments of both TADdyn and TADbit.  
  
- "fastTest" folder contains a small interaction matrix to test that the optimisation and modelling steps work well.  
  
- "jupyterNotebooks" folder contains the Notebooks used for the modelling and analysis of the pcHi-C datasets.  
  
- "matrices" folder contains two levels of subfolders. The first one stating the cell ID, and the second one the ID of the region defined in the interaction matrix. Inside ...matrices/cellID/regionID/ we will find its correspondent interaction matrix. The interaction matrix file format used in this Notebooks is a tab-delimited text version of the matrix with a number of columns equal to the number of elements delimited by a tab, and a number of rows equal to the number of lines in the file. The naming format of the matrices must follow Matrix[extraIfWanted]_[cellID]_[regionID]_[Chromosome]-[startCoord]-[EndCoord]_[resolution]bp, where:  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*extraIfWanted: is any additional text we might want to add after Matrix. It can be left empty  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*Chromosome: is the chromosome name of the region as stated in the bam file  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*startCoord: is the starting coordinates of the region stored in the file  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*EndCoord: is the coordinate of the beginning of the last bin included in the matrix. I. E.: If the resolution is 5 kb, and EndCoord is 15000, the matrix would include the bin that goes from 15000 to 19999 (15000 + 4999)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*resolution: is the resolution of the matrix in bp  
  
- "misc" folder contains the Notebooks workflow image displayed in the Readme file.
  
- "models" folder contains two levels of subfolders. The first one stating the cell ID, and the second one the ID of the region defined in the .models file. Inside ...models/cellID/regionID/ we will find its correspondent ".TADdynDict" or ".models" file. Naming format:  
[cellID]_[regionID]_C[distanceCutoff]L[lowfreq]U[upfreq]M[maxdist]Res[resolution].[TADdynDict/models]  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*distanceCutoff: is the optimal distance cutoff value that was set in the optimization for the modelling  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*lowfreq: is the optimal lowfreq value that was set in the optimization for the modelling  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*upfreq: is the optimal upfreq value that was set in the optimization for the modelling  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*maxdist: is the optimal Maxdist value that was set in the optimization for the modelling  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*resolution: is the resolution of the experiment in bp.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*TADdynDict/models: is the format of the file. "TADdynDict" for TADdyn models and "models" for TADbit models.  
  
- "optimization" folder contains two levels of subfolders. The first one stating the cell ID, and the second one the ID of the region. Inside ...optimization/cellID/regionID/ we will find its correspondent optimization output files, that will be used in <ins>02_chooseBestParameters.ipynb</ins> to define the best modelling parameters per each cell and region combination.  
  
- "outData" folder contains the output files generated during the analysis.  
  
- "outPlot" folder contains the plot pdf files generated during the analysis.  

