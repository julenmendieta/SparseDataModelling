#!/usr/bin/env python
# coding: utf-8


### PARAMETERS TO MODIFY ###

## Data Paths (Location of the base folder downloaded from GitHub)
basePath = '/home/julen/TADdyn/SparseDataModelling/'

## MODIFY THIS VALUES ACCORDING TO YOUR RESULTS
# gather top dcutoff and maxdist parameters
topDuctoff = 200
topMaxdist = 300

######




# In this Notebook we will use the estimated best parameters for the modelling 
# to reate a file with the best combinations of lowfreq and upfreq for the given
# best dcutof and maxdist
# 

# # Libraries and functions
import sys


# # Run 

# ## Define additional paths 


scriptsPath = basePath + 'code/modellingScripts/'
GeneralOptimOutPath = basePath + 'optimization/'


# ## Import additional libraries 


sys.path.append(basePath + 'code')
import fileHandling
import plotting
sys.path.append(basePath + 'code/modellingScripts')
import Check_optimization


# ## Get matrix paths 


matricesLength, regionsAll, matrices = fileHandling.getMatricesPaths(basePath, 
                                                                    starting='Matrix')


# ## Generate file with top parameters 

# We will create a file called modellinParams.txt in the optimization folder 
# before the folders of each cell and region. This file will contain the parameters 
# with the top correlation for the given dcutoff and maxdist. These files are 
# created to facilitate the modelling process, but the selected parameter combinations 
# should be checked prior to running the final modelling step.   
# 
# In the case of Ery, the best combination of paramteres for the final ensemble of 
# models was set to:  
# lowfreq = 0.0  
# upfreq = 0.0  
# maxdist = 300    
# dcutoff = 200 # NOTE: this is set to 200 as is the best parameter for the comparative 
# analyis with Mon and nCD4


## get file with top optimization paramteres to model
inputPaths = [[],[]]
for regi in regionsAll:
    for cell in matrices:
        matrixPaths = matrices[cell][regi]
        optimPaths = basePath + 'optimization/%s/%s/' %(cell, regi)
        inputPaths[0] += [matrixPaths]
        inputPaths[1] += [optimPaths]
        

        
topModels = '%s_%s' %(topDuctoff, topMaxdist)
topCorrelations=Check_optimization.checkAll(GeneralOptimOutPath, [inputPaths[0],
                                                                  inputPaths[1]], 
                                            show_dcut=False, dcut=False, topModels=topModels)



print topCorrelations




