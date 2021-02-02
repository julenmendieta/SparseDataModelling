# Sparse Data Modelling
We provide [IPython Notebooks](https://jupyter.org/) (in "jupyterNotebooks" folder) containing the code to perform each step of the analysis in:  
    
**3D reconstruction of genomic regions from sparse interaction data**.  
Julen Mendieta-Esteban, Marco Di Stefano, David Castillo, Irene Farabella, and Marc A Marti-Renom.  
https://www.biorxiv.org/content/10.1101/2020.10.11.334847v1.full.pdf  
  
from the extraction of the interaction matrix of interest from the BAM file, to matrix normalisation, modelling, and models analysis.  

Intermediate files and additional information are also included.

To start using this repository, first follow the [installation](#installation) steps, and then run with Python 2 each of the Notebooks located in the "jupyterNotebooks" folder in the [order](#notebooks-ordering-and-programs-used-for-them) stated below. 

## Installation
### Repository files:
First download or clone the repository files to your computer. Since some example files files are too heavy, we have upploaded them into Dropbox. In order to download them, you can follow two approaches:  
- Using the browser:  
Download the following two zip files to the "SparseDataModelling" folder, and unzip them by merging the output "models" and "outData" folders with the existent ones:  
https://www.dropbox.com/s/h718iid33lq41hx/models.zip?dl=0  
https://www.dropbox.com/s/xsf9g5l8fdw9907/outData.zip?dl=0  
  
- Using the terminal:  
Enter into de "SparseDataModelling" folder and run the following code:
```
wget -O models.zip https://www.dropbox.com/s/h718iid33lq41hx/models.zip?dl=0 ; unzip -o models.zip   
wget -O outData.zip https://www.dropbox.com/s/xsf9g5l8fdw9907/outData.zip?dl=0 ; unzip -o outData.zip   
```

Note: Avoid using '_' in any folder name inside or in the root of the location of this repository in your computer.  
  
### Programs  
This repository relies on some programs for the modelling and analysis steps. Users can choose between downloading a [Singularity](https://sylabs.io/) container with them or directly installing them into their computers.  
##### Using Singularity containers
Some users might want to skip the instalation process, specially in HPCs not directly accessible to the public internet. For those cases, the Singularity container recipes for TADdyn and TADbit are available in the "containers" folder. To use them, users will need to first [install Singularity](https://sylabs.io/guides/3.6/admin-guide/installation.html) (we do NOT recommend to install Singularity with Conda). The actual containers can be downloaded from:   
https://www.dropbox.com/sh/uz7iikid2w9wv0d/AADPVGm4dMIiv2OtROEFakEJa?dl=0  
    
To open the Notebooks users can use the following commands with each of them:  
TADdyn:  
`singularity exec singularity_TADdyn.sif jupyter notebook --port=8888 --notebook-dir=/PATH/TO/NOTEBOOKS`  
TADbit:  
`singularity exec singularity_TADbit.sif jupyter notebook --port=8888 --notebook-dir=/PATH/TO/NOTEBOOKS` 
  
Where /PATH/TO/NOTEBOOKS is the path where the "jupyterNotebooks" folder is located in the computer.

##### Installing TADdyn and TADbit in the computer
TADdyn can be installed following the steps here:  
https://github.com/3DGenomes/TADdyn  
TADbit can be installed following the steps here:  
https://github.com/3DGenomes/TADbit  
  
You will need to install the following additional python libraries for data analysis and plotting. These libraries are included in the provided singularity recipes and containers, but are absent in the installation instructions from 3Dgenomes.  
  
For TADdyn:  
- seaborn  
    
For TADbit:  
- seaborn  
- fastcluster  
- sklearn  

## Notebooks ordering and programs used for them
The ordering to run the Notebooks is stated in the numbers at the beginning of the folders in "jupyterNotebooks". The same rule is valid for the Notebooks located inside. In this way, we would run the Notebooks in the following order provided that we have a bam file as a starting point.  
First:  
- "01_inputData". All the Notebooks inside have to be executed with TADbit  

Then:  
- "02_modelling". All the Notebooks inside have to be executed with TADdyn  

Finally:  
- "03_modelAnalysis". All the Notebooks inside have to be executed with TADbit  

This repository is organised in a tree directory structure to facilitate the analysis of the users own data. In this way:  
- Starting from a bam file: Users can normalise and store the interaction matrices from their own data by adding their bam files into the “bamfiles” folder in the same tree directory structure as we state in the [Folder structure](#folder-structure-and-content) section below. Then, they would need to run the Notebooks inside "01_inputData". To ensure that the bam files have the right format, users can follow the instructions provided in <ins>01_inputData/01_retrievingBAMfiles.ipynb</ins>. To get the normalisation biases and the interaction matrices users have to run <ins>01_inputData/02_saveMatrixFiles.ipynb</ins>.  

- Starting from an interaction matrix file: Users can work with their own interaction matrices by emptying the "matrices" folder and including their own files in the same tree directory structure as we state in the [Folder structure](#folder-structure-and-content) section below. Then, they would need to run the Notebooks inside "02_modelling". By running those Notebooks in order, users will optimize, select best parameters, and model their matrices.

- Starting from a TADdyn or TADbit models file. Users can analyse their own TADdyn or TADbit format models by emptying the "models" folder and including their own files in the same tree directory structure as we state in the [Folder structure](#folder-structure-and-content) section below. Then, they would need to run the Notebooks inside "03_modelAnalysis". If they would add a TADdyn model they would need to first transform it to TADbit model format by running <ins>03_modelAnalysis/01_convertTADdynModels_toTADbitModels.ipynb</ins>. After this, they will have a TADbit format model file that will be analysed in the next Notebooks. 
  
## Parameters to be modified in the Notebooks
Every Notebook has a section named “Parameters to modify” that contains all the parameters that must be modified according to the users aim and data. The only Notebook that contains an additional section to be modified is <ins>02_chooseBestParameters.ipynb</ins>, that requires the user in section "Generate file with top parameters" to decide which are the dcutoff and maxdist values that best correlations have shown in the plots.
    
&nbsp;&nbsp;    
&nbsp;&nbsp;
&nbsp;&nbsp;
&nbsp;&nbsp;
&nbsp;&nbsp;
# Additional information

## Notebook purpose
**01_inputData**  
- <ins>01_retrievingBAMfiles.ipynb</ins>: Provides information about how to treat FASQ files to obtain bam files suitable for TADbit  
  
- <ins>02_saveMatrixFiles.ipynb</ins>: Contains the code to both obtain the biases for the PRINT normalisation and create the interaction matrices files 
  
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
  
## Testing that the modelling works well
Users can set as True the variable "runFastTest" in <ins>01_modellingParametersGridSearch.ipynb</ins> and <ins>03_modelling.ipynb</ins> to test in a short time that the code for the modelling is working well.  

## Using a cluster to run the modelling
The most time consuming steps in this repository are the ones involving 01_inputData and 02_modelling. For this reason, some users might want to run them in a cluster, so that more resources can be allocated for the task. Hence, the code needed for these steps has been stored in the "codeScripts" folder as a set of scripts that can be executed in the terminal and submitted to a job queue. Users still need to modify some [parameters](#parameters-to-be-modified-in-the-notebooks) present at the beginning of the scripts.

NOTE: For a matter of consistency, the code of the scripts in "codeScripts" is the same as the one found in the Notebooks (with minor differences). The only difference regards to the code referring to 02_chooseBestParameters.ipynb, which has been splited into two scripts (02a_chooseBestParameters.py and 02b_storeBestParameters.py) so that all the variables to be modified are easily found.

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

