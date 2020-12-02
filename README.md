# TADdyn_tutorial
We provide IPython Notebooks (in jupyterNotebooks) containing all the steps from the obtention of the interest interaction matrices from BAM files, to their normalisation, modelling, and further analysis. As it was done in:  
**3D reconstruction of genomic regions from sparse interaction data**.  
Julen Mendieta-Esteban, Marco Di Stefano, David Castillo, Irene Farabella, and Marc A Marti-Renom. 
https://www.biorxiv.org/content/10.1101/2020.10.11.334847v1.full.pdf

Intermediate files and additional information are also included.

To start using this Git, go to jupyterNotebooks folder and run the notebooks inside of the folders in numeric order with Python 2.

## Running the Notebooks:  
### Using Singularity containers
In case is needed, the Singularity container recipes for TADdyn and TADbit are available in the "containers" folder. The actual containers can be downloaded from:   
https://www.dropbox.com/sh/uz7iikid2w9wv0d/AADPVGm4dMIiv2OtROEFakEJa?dl=0

To open the Notebooks we could use the following commands with each of them:  
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
The ordering to run the Notebooks is stated in the numbers at the beginning of their names.  
First:  
01_inputData scripts are executed with TADbit  
Then:  
02_modelling scripts are executed with TADdyn  
Finally:  
03_modelAnalysis scripts are executed with TADbit  


## Folder structure
-"additionalInput" folder contains text files with additional data to be loaded, like enhancer and promoter coordinates, and methylation or gene expression data.  
-"code" folder contains sets of python 2 functions that will be using in the Notebooks.  
-"containers" folder contains the recipes for building the Singularity environments of both TADdyn and TADbit.  
-"jupyterNotebooks" folder contains the Notebooks used for the modelling and analysis of the pcHi-C datasets.  
-"matrices" folder contains two levels of subfolders. The first one stating the cell ID, and the second one the ID of the region defined in the interaction matrix. Inside ...matrices/cellID/regionID/ we will find its correspondent interaction matrix.  
-"models" folder contains two levels of subfolders. The first one stating the cell ID, and the second one the ID of the region defined in the .models file. Inside ...models/cellID/regionID/ we will find its correspondent .models file.  
-"optimization" folder contains two levels of subfolders. The first one stating the cell ID, and the second one the ID of the region. Inside ...optimization/cellID/regionID/ we will find its correspondent optimization output files, that will be used in "02_chooseBestParameters.ipynb" to define the best modelling parameters per each cell and region combination.  
-"outData" folder contains the output files generated during the analysis.  
-"outPlot" folder contains the plot pdf files generated during the analysis.  

