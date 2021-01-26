#!/usr/bin/env python
# coding: utf-8


### PARAMETERS TO MODIFY ###

## paths. Avoid paths with white spaces or _
# Location of the base folder downloaded from GitHub
basePath = '/home/julen/TADdyn/SparseDataModelling/'


## Execution parameters
# days-hours:minutes:seconds of maximum running time of the process
jobTime = '0-08:00:00' 
# how many models will we create
nmodels = 160
# Steps in which you what to divide the modelling (to divide the modelling
# in n steps and save the output)
steps = 20
# How many CPU's will we allow for the process
ncpu = 26

#If we set outputAppart to False the models will be stored inside of a folder 
# called finalModels inside of the same folder as the matrix. If set to True, 
# models will be stored in a folder called "models" at the same level of the 
# folder called "matrices".
outputAppart = True

## If you want to run a fast test with a small matrix set runFastTest to True
# This option will use a short matrix to compute the optimisation
# only in one  test cell.
# WARNING: It wont use any of the matrices stored in the matrices 
# subfolders, only the test one that is located in matrices/testCell
# output models will be stored at basepath + fastTest/models/
runFastTest = False


######





# # Libraries and functions 


import sys
import os
import shutil
import datetime



# # Run 

# ## Define additional paths 

scriptsPath = basePath + 'code/modellingScripts/'
GeneralOptimOutPath = basePath + 'optimization/'
tempOut = basePath + 'temporal/modelling'


# ## Import additional libraries

sys.path.append(basePath + 'code')
import metrics
sys.path.append(basePath + 'code/modellingScripts')
import combineModels2
import runmodelling2


# ## Recover best optimization parameters 

# Here we will recover the top parameters for the final ensemble of 
# models (1000 out of 1500 modelled), calculate the time needed to 
# compute all models and get a variable with the commands to run. 
# These parameters were defined in 02_chooseBestParameters.ipynb

# Load modelling info
if runFastTest == True:
    print 'NOTE: Running fast test instead of inputed matrices'
    testM = basePath + '''fastTest/matrices/testCell/testReg/MatrixNormFreq\
Filtrd_testCell_testReg_chr8-132755000-133555000_5000bp'''
    content = '%s\t0.0\t0.0\t200.0\t300.0' %testM
    with open(basePath + 'fastTest/optimization/modellinParams.txt', 'w') as f:
        f.write(content)
        
    nmodels2 = 15
    steps2 = 3
    ncpu2 = 5
    
    cmds, times = metrics.getModellingCommand(basePath + 'fastTest/optimization/', 
                                              tempOut, jobTime, nmodels2, ncpu2, 
                                              outputAppart, scriptsPath, step=steps2)
    
else:
    cmds, times = metrics.getModellingCommand(GeneralOptimOutPath, tempOut, jobTime,
                       nmodels, ncpu, outputAppart, scriptsPath, step=steps)
    ncpu2 = ncpu
   


# get expected modelling time
totalTime2 = str(datetime.timedelta(seconds=times))
print("Stimated time with assigned number of %s CPU's: %s" %(ncpu2,
                                    totalTime2))


# ## Reconstruct the 3D organisation of the chromatin

# This code will produce the models (.modelsTemp files) and merge them (in 
# case were done in steps) into a .TADdynDict file. If previous runs with 
# the same combinations of parameters are found in the folder, will also 
# be merged. WARNING: Merged files will be erased.


# Make models (results in .model files)
for nc, cmd in enumerate(cmds):
    print('Run %s of %s' %((nc + 1), len(cmds)))
    matPath = cmd.split(' ')[12]
    low = float(cmd.split(' ')[2])
    up = float(cmd.split(' ')[8])
    maxd = cmd.split(' ')[6]
    if '.' in maxd:
        maxd = int(float(maxd))
    else:
        maxd = int(maxd)
    c = float(cmd.split(' ')[4])
    lammpsOut = cmd.split(' ')[10]
    n_cpu = int(cmd.split(' ')[18])
    pathOut = cmd.split(' ')[20]
    runmodelling2.runModelling(matPath, low, up, maxd, c, lammpsOut, jobTime,
                            tempOut=None, nmodels=int(cmd.split(' ')[16]), 
                            n_cpu=n_cpu, pathOut=pathOut)
    
# Merge, filter, and delete .model files
# This step is important since files finishing in .models
# can be overwriten in future modelling runs
outpaths = set()
for cmd in cmds:
    outpaths.add(cmd.split()[-1])
print('--- Merging ---')
for out in outpaths:
    print(out)
    combineModels2.combineModels(out)


# ## Clean temporal folders 

# The content inside of the generated temporal folders inside 
# $tempOut is usually deleted after the modelling successfully 
# finishes, not the folders itself though. Besides, if the modelling 
# process breaks for any reason, the temporal files will remain there. 
# TADdyn modelling generates a lot of temporal files so, in order to ensure 
# that they don't accumulate, we will remove the container folders after 
# each optimisation process.

tempFolders = os.listdir(basePath + 'temporal/')


for t in tempFolders:
    print t
    shutil.rmtree(basePath + 'temporal/%s' %t)

