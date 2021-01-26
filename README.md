# Sparse Data Modelling
We provide IPython Notebooks (in "jupyterNotebooks" folder) containing the code to perform each step of the analysis in:  
**3D reconstruction of genomic regions from sparse interaction data**.  
Julen Mendieta-Esteban, Marco Di Stefano, David Castillo, Irene Farabella, and Marc A Marti-Renom.  
https://www.biorxiv.org/content/10.1101/2020.10.11.334847v1.full.pdf  
from the extraction of the interaction matrix of interest from the BAM file, to matrix normalisation, modelling, and models analysis.  

Intermediate files and additional information are also included.

To start using this Git, enter the "jupyterNotebooks" folder and run the Notebooks inside each sub-directory following the numerical order (Python 2).  

## Before running the Notebooks:
Download the model files and additional data that could not be uploaded to GitHub.  
-Using the browser:  
Download the following two zip files, and unzip them by merging the output "models" and "outData" folders with the existent ones:  
https://www.dropbox.com/s/h718iid33lq41hx/models.zip?dl=0  
https://www.dropbox.com/s/xsf9g5l8fdw9907/outData.zip?dl=0  
  
-Using the terminal:  
```
wget -O models.zip https://www.dropbox.com/s/h718iid33lq41hx/models.zip?dl=0 ; unzip -o models.zip   
wget -O outData.zip https://www.dropbox.com/s/xsf9g5l8fdw9907/outData.zip?dl=0 ; unzip -o outData.zip   
```



## Running the Notebooks:  
### Using Singularity containers
In case is needed, the Singularity container recipes for TADdyn and TADbit are available in the "containers" folder. The actual containers can be downloaded from:   
https://www.dropbox.com/sh/uz7iikid2w9wv0d/AADPVGm4dMIiv2OtROEFakEJa?dl=0

To open the Notebooks you can use the following commands with each of them:  
TADdyn:  
singularity exec singularity_TADdyn.sif jupyter notebook --port=8888 --notebook-dir=/PATH/TO/NOTEBOOKS  
TADbit:  
singularity exec singularity_TADbit.sif jupyter notebook --port=8888 --notebook-dir=/PATH/TO/NOTEBOOKS

### Installing TADdyn and TADbit in your computer
TADdyn can be installed following the steps here:  
https://github.com/3DGenomes/TADdyn  
TADbit can be installed following the steps here:  
https://github.com/3DGenomes/TADbit  

You will need to install the following additional python libraries for data analysis and plotting. These libraries are included in the provided singularity recipes and containers, but are absent in the installation instructions from 3Dgenomes.  

For TADdyn:  
-seaborn  
   
For TADbit:  
-seaborn  
-fastcluster  
-sklearn  

### Notebooks ordering and programs used for them
The ordering to run the Notebooks is stated in the numbers at the beginning of the folders in jupyterNotebooks. The same rule is valid for the Notebooks located inside. In this way we would run the Notebooks in the following order if we would have a bam file as an starting point.  
First:  
"01_inputData". Notebooks have to be executed with TADbit  
Then:  
"02_modelling". Notebooks have to be executed with TADdyn  
Finally:  
"03_modelAnalysis". Notebooks have to be executed with TADbit  
  
Depending on our starting data we can avoid some steps. In this way, if we would already have the interaction matrices files we would start from "02_modelling", and if we would already had the models and just wanted to analyse them we would start from "03_modelAnalysis".  

### Testing that the modelling works well
Users can set as True the variable "runFastTest" in "01_modellingParametersGridSearch" and "03_modelling.ipynb" to test in a short time that the code for the modelling is working well.  

## Using a cluster to run the modelling
The most time consuming steps in this repository are the ones involving 01_inputData and 02_modelling. For this reason some users might want to run them in a cluster so that more resources can be allocated for the task. Hence, the code needed for these steps has been stored in the "codeScripts" folder as a set of scripts that can be executed in the terminal and hence submited to a job queue.  

NOTE: For a matter of consistency the code of the scripts in "codeScripts" is the same as the one found in the Notebooks (with minor differences). The only change that users have to note is that the section with the variables to be modified has been moved to the beggining of the script.

## Folder structure
-"additionalInput" folder contains text files with additional data to be loaded, like enhancer and promoter coordinates, and methylation or gene expression data.  
-"bamfiles" folder contains one level of subfolders stating the cell ID. Inside ...bamfiles/cellID/ we should store its correspondent sorted .bam file and the associated index .bai file. This repository does not provide any bam file due to the large size of the files.  
-"code" folder contains sets of python 2 functions that will be using in the Notebooks.  
-"codeScripts" folder contains the scripts involving the Notebooks from "01_inputData" and "02_modelling" so that this steps can be ran in a cluster.  
-"containers" folder contains the recipes for building the Singularity environments of both TADdyn and TADbit.  
-"fastTest" folder contains an small interaction matrix to test that the optimisation and modelling steps work well.  
-"jupyterNotebooks" folder contains the Notebooks used for the modelling and analysis of the pcHi-C datasets.  
-"matrices" folder contains two levels of subfolders. The first one stating the cell ID, and the second one the ID of the region defined in the interaction matrix. Inside ...matrices/cellID/regionID/ we will find its correspondent interaction matrix.  
-"models" folder contains two levels of subfolders. The first one stating the cell ID, and the second one the ID of the region defined in the .models file. Inside ...models/cellID/regionID/ we will find its correspondent .models file.  
-"optimization" folder contains two levels of subfolders. The first one stating the cell ID, and the second one the ID of the region. Inside ...optimization/cellID/regionID/ we will find its correspondent optimization output files, that will be used in "02_chooseBestParameters.ipynb" to define the best modelling parameters per each cell and region combination.  
-"outData" folder contains the output files generated during the analysis.  
-"outPlot" folder contains the plot pdf files generated during the analysis.  

