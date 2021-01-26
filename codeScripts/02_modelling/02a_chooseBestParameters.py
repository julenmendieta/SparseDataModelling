#!/usr/bin/env python
# coding: utf-8


### PARAMETERS TO MODIFY ###

## Data Paths (Location of the base folder downloaded from GitHub)
basePath = '/home/julen/TADdyn/SparseDataModelling/'

# the plots that help in the selection of dcutoff and maxdist
# will be stored in basepath + outPlot/optimization

## IMPORTANT NOTE ##
# The best dcutoff and maxdist values found by this script will be
# used as input in 02b_storeBestParameters.py to generate the
# file for the modelling

######




# In this Notebook we will visually display the correlation data obtained 
# in "01_exploreBestParameters.ipynb". This will allow us to define the 
# best optimization parameters to compute the final ensemble of 3D models.
# 

# # Libraries and functions
import sys
import os
#import matplotlib
#matplotlib.use('Agg')
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
import optTXTparser2


# ## Get matrix paths 


matricesLength, regionsAll, matrices = fileHandling.getMatricesPaths(basePath, 
                                                                    starting='Matrix')


# ## Join output files 

# To be more clean and have the results more handy, we will combine the correlation 
# values splited in each parameter combination file into a single file per different 
# dcutoff values.

for regi in regionsAll:
    print('--- %s ---' %regi)
    for cell in matrices:
        print('## %s ##' %(cell))
        scriptp = basePath + 'code/modellingScripts/optTXTparser.py'
        inPath = '/'.join(matrices[cell][regi].split('/')[:-1]) + '/'
        outPath = basePath + 'optimization/%s/%s/' %(cell, regi)
        # create folder if absent
        if not os.path.exists(outPath):
            os.makedirs(outPath)
        # run script
        optTXTparser2.parseOptimOutput(inPath, outPath)
        

# ## Check best dcutoff 

# This step allows us to decide the optimal dcutoff for the modelling. In this 
# example we just prepared the optimization in one interaction matrix (b-globin of Ery). 
# See how the correlation values descend when increasing and decreasing the dcutoff 
# value. In this case a dcutoff of 300 would be the best option, since its correlation 
# values are higher than values above and below.   
# 
# NOTE: If we aim to compare matrices from different conditions we should select a 
# dcutoff that fits best all the matrices, even if it is not the optimal in all the cases.


for regi in regionsAll:
    print('--- %s ---' %regi)
    for cell in matrices:
        print('## %s ##' %(cell))
        savePath = basePath + 'outPlot/optimization/%s/' %regi
        if not os.path.exists(savePath):
            os.makedirs(savePath)

        optimPaths = [basePath + 'optimization/%s/%s/' %(cell, regi)]
        topCor = Check_optimization.checkAll(GeneralOptimOutPath, optimPaths, 
                                            show_dcut=True, 
                                             dcut=False, topModels=False,
                                             savePath=savePath, outOff=True)


# This step allows us to decide the optimal dcutoff for the modelling

# ## Check best maxdist 

# Here we find that the best correlationw is obtained using a maxdist of 300 nm

topDuctoff = 200
for regi in regionsAll:
    print('--- %s ---' %regi)
    for cell in matrices:
        print('## %s ##' %(cell))
        savePath = basePath + 'outPlot/optimization/%s/' %regi
        if not os.path.exists(savePath):
            os.makedirs(savePath)

        optimPaths = [basePath + 'optimization/%s/%s/' %(cell, regi)]
        topCor = Check_optimization.checkAll(GeneralOptimOutPath, optimPaths, 
                                            show_dcut=False, 
                                             dcut=topDuctoff, topModels=False,
                                             savePath=savePath, outOff=True)


# The code bellow will produce a plot with the correlation values associated with 
# the "topDuctoff" to assess the chosen top parameters. The best correlation should 
# usually be surrounded by, or nearby to, other top ranked correlation values in the 
# following plot


# Here we define the number of models we did in the optimization, and how many of 
# them we kept
nmodels = 100
nkeep = 100  

inputPaths = [[],[]]
for regi in regionsAll:
    for cell in matrices:
        matrixPaths = matrices[cell][regi]
        optimPaths = basePath + 'optimization/%s/%s/' %(cell, regi)
        inputPaths[0] += [matrixPaths]
        inputPaths[1] += [optimPaths]
        

for nm, matPath in enumerate(inputPaths[0]):
    plotting.optimPlot1(matPath, inputPaths[1][nm], nmodels = nmodels, 
                        nkeep = nkeep, ductoff=topDuctoff, outOff=True)