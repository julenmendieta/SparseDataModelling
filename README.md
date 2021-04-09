# Sparse Data Modelling
We provide [IPython Notebooks](https://jupyter.org/) (in "jupyterNotebooks" folder) containing the code to perform each step of the analysis in:  
    
**3D reconstruction of genomic regions from sparse interaction data**.  
Julen Mendieta-Esteban, Marco Di Stefano, David Castillo, Irene Farabella, and Marc A Marti-Renom.  
https://www.biorxiv.org/content/10.1101/2020.10.11.334847v1.full.pdf  
  
from the extraction of the interaction matrix of interest from the BAM file, to matrix normalisation, modelling, and models analysis.  

Intermediate files and additional information are also included.

To start using this repository, first follow the [installation](#installation) steps, and then run with Python 2 each of the Notebooks located in the "jupyterNotebooks" folder in the [order](#notebooks-workflow) stated below. 

## Installation
### Repository files:
First download or clone the repository files to your computer. Since some example files are too heavy, we have upploaded them to another repository. In order to download them, you can follow two approaches:  
- Using the browser:  
Download the whole additional data file to the "SparseDataModelling" folder, unzip it, and unzip the models.zip and outData.zip by merging the output "models" and "outData" folders with the existent ones. Then move the .sif files to the containers folder:  
http://sgt.cnag.cat/3dg/datasets/files/20210201_Mendieta-Esteban_etal_CodeData.zip  
  
- Using the terminal:  
Enter into de "SparseDataModelling" folder and run the following code:
```
wget -O extra.zip http://sgt.cnag.cat/3dg/datasets/files/20210201_Mendieta-Esteban_etal_CodeData.zip ; unzip -o extra.zip
unzip -o models.zip   
unzip -o outData.zip   
mv *sif containers/
```

Note: Avoid using '_' in any folder name inside or in the root of the location of this repository in your computer.  
  
### Dependencies
This repository relies on some programs for the modelling and analysis steps. Users can choose between downloading a [Singularity](https://sylabs.io/) container with them or directly installing them into their computers.  
##### Using Singularity containers
Some users might want to skip the instalation process, specially in HPCs not directly accessible to the public internet. For those cases, the Singularity container recipes for TADdyn and TADbit are available in the "containers" folder. To use them, users will need to first [install Singularity](https://sylabs.io/guides/3.6/admin-guide/installation.html) (we do NOT recommend to install Singularity with Conda). The actual containers are available inside of the additional data file:   
http://sgt.cnag.cat/3dg/datasets/files/20210201_Mendieta-Esteban_etal_CodeData.zip  
    
To open the Notebooks users can use the following commands with each of them:  
TADdyn:  
`singularity exec singularity_TADdyn.sif jupyter notebook --port=8888 --notebook-dir=/PATH/TO/NOTEBOOKS`  
TADbit:  
`singularity exec singularity_TADbit.sif jupyter notebook --port=8888 --notebook-dir=/PATH/TO/NOTEBOOKS` 
  
Where /PATH/TO/NOTEBOOKS is the path where the "jupyterNotebooks" folder is located in the computer.

##### Installing TADdyn and TADbit in the computer
TADdyn can be installed following the steps here:  
https://github.com/3DGenomes/TADdyn  
\*Please, remember to set the LAMMPS exception flag in the [LAMMPSfolder]/src/Makefile.package file to:  
&nbsp;&nbsp;&nbsp;&nbsp;LMP_INC = -DLAMMPS_EXCEPTIONS  

TADbit can be installed following the steps here:  
https://github.com/3DGenomes/TADbit  
  
You will need to install the following additional python libraries for data analysis and plotting. These libraries are included in the provided singularity recipes and containers, but are absent in the installation instructions from 3Dgenomes.  
  
For TADdyn:  
- seaborn  
    
For TADbit:  
- seaborn  
- fastcluster  
- sklearn  
 
## Using a cluster computer to run the modelling
The most time consuming steps in this repository are the ones involving 01_inputData and 02_modelling. For this reason, some users might want to run them in a cluster, so that more resources can be allocated for the task. Hence, the code needed for these steps has been stored in the "codeScripts" folder as a set of scripts that can be executed in the terminal and submitted to a job queue. Users still need to modify some [parameters](#parameters-to-be-modified-in-the-notebooks) present at the beginning of the scripts.

NOTE: For a matter of consistency, the code of the scripts in "codeScripts" is the same as the one found in the Notebooks (with minor differences). The only difference regards to the code referring to 02_chooseBestParameters.ipynb, which has been splited into two scripts (02a_chooseBestParameters.py and 02b_storeBestParameters.py) so that all the variables to be modified are easily found.

## Parameters to be modified in the Notebooks
This respository contains multiple Notebooks and [scripts](using-a-cluster-computer-to-run-the-modelling) to run the main analysis done in the work cited [before](#sparse-data-modelling). Every Notebook has a section named “Parameters to modify” that contains all the parameters that must be modified according to the users aim and data. 
  
  
## Notebooks workflow
The ordering to run the Notebooks is stated in the numbers at the beginning of the folders in "jupyterNotebooks". The same rule is valid for the Notebooks located inside. In this way, and with a BAM file as an starting point, we will run the Notebooks as stated in the workflow below:  

<img src="https://github.com/julenmendieta/SparseDataModelling/blob/main/misc/scriptsFlow.jpeg" width="500" height="400">

As we see from the workflow, the repository is organised in a way that facilitates the analysis of the users own data. In this way:  
1. Starting from a bam file: Users can normalise and store the interaction matrices from their own data by adding their bam files into the “bamfiles” folder. Then, they would need to run the Notebooks inside "01_inputData". To ensure that the bam files have the right format, users can follow the instructions provided in <ins>01_inputData/01_retrievingBAMfiles.ipynb</ins>. To get the normalisation biases and the interaction matrices users have to run <ins>01_inputData/02_saveMatrixFiles.ipynb</ins>.  

2. Starting from an interaction matrix file: Users can work with their own interaction matrices by emptying the "matrices" folder and including their own files there. Then, they would need to run the Notebooks inside "02_modelling". By running those Notebooks in order, users will optimize, select best parameters, and model their matrices.

3. Starting from a TADdyn or TADbit models file. Users can analyse their own TADdyn or TADbit format models by emptying the "models" folder and including their own files. Then, they would need to run the Notebooks inside "03_modelAnalysis". If they would add a TADdyn model they would need to first transform it to TADbit model format by running <ins>03_modelAnalysis/01_convertTADdynModels_toTADbitModels.ipynb</ins>. After this, they will have a TADbit format model file that will be analysed in the next Notebooks. 
   
At the time to add your own files, please follow the same tree directory structure as we state in the [Folder structure](https://github.com/julenmendieta/SparseDataModelling/blob/main/extraInformation.md#folder-structure-and-content) section below.  

NOTE: if you are using Singularity containers to run the the Notebooks, run the ones located in "01_inputData  and in "03_modelAnalysis" using the container with TADbit, and the ones located in "02_modelling" using the container with TADdyn.
    
&nbsp;&nbsp;    
&nbsp;&nbsp;
&nbsp;&nbsp;
&nbsp;&nbsp;
&nbsp;&nbsp;
# Additional information
Users can get additional information regarding the purpose of each of the Notebooks and the Folder structure and content in [here](https://github.com/julenmendieta/SparseDataModelling/blob/main/extraInformation.md).
