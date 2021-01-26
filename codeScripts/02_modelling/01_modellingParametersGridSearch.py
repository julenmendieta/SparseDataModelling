#!/usr/bin/env python
# coding: utf-8


### PARAMETERS TO MODIFY ###

# Now we will define a range of values per each parameter, and we will check how well the models represent the input data.  
# 
# NOTE: The interaction matrix should ideally be included inside a structure of folders without '\_'. Nested in a folder that indicates first the cell type and then the project or regions. Then it must start by Matrix, and contain the following elements separated by '\_':  
# .../Cell/Region/[Matrix]\_[Cell]\_[Region]\_[Chromosome]-[startCoord]-[EndCoord]\_[resolution]bp  
# 
# Example:  
# .../Ery/b-globin/MatrixNormFreqFiltrd_Ery_b-globin_chr11-4615000-6175000_5000bp
# In this example, all matrices will be in a tree of folders inside "matrices" folder

# You can find a guide to decide the range of parameters to check for your regions of interest in here:  
# https://github.com/3DGenomes/MethodsMolBiol/blob/master/Notebooks/Methods-10-Modeling%20parameters%20optimization.ipynb  
# 
# The following code explores lowfreq values ranging from -1 to 0.5 in steps of 0.5, upfreq values from -0.5 to 0.5 in steps of 0.5. For the shake of speed, we will only look for one value of maxdist (300) and dcutoff (200) in this example.


## Optimization parameters
# range of lowfreq values
lowfreq_arange = np.arange(-1.0,1.0,0.5)
# range of dcutoff values
dcutoff_range= [200,300,100]  # start, end + step, step 
# range of maxdist values
m_range= np.arange(300,400,100)
# upfreq range values
upfreq_range= np.arange(-0.5,1.0,0.5)
# number of CPU to use
cpus = 8 
# number of models to reconstruct in the optimization
nmodels = 100
# days-hours:minutes:seconds of maximum running time per model
jobTime = '0-08:00:00' 

## Data Paths (Location of the base folder downloaded from GitHub)
basePath = '/home/julen/TADdyn/SparseDataModelling/'

## If you want to run a fast test with a small matrix set to True
# This option will use a short matrix to compute the optimisation
# only in one  test cell.
# WARNING: It wont use any of the matrices stored in the matrices 
# subfolders, only the test one that is located in matrices/testCell
runFastTest = True


###


# In this Notebook we will explore different combinations of input parameters (lowfreq, upfreq, maxdist, and dcutoff) in order to find the optimal combination to compute reliable 3D models of our regions of interest \*. The code allows generating files with the correlation between the models' ensembles and the input interaction matrix that will be assessed in "02_chooseBestParameters.ipynb" in order to select the best parameters combination, or (if needed) define a different range of input parameters to continue the exploration of the parameters.  
#  
#  \* Serra, F., Baù, D., Goodstadt, M., Castillo, D. Filion, G., & Marti-Renom, M.A. (2017). Automatic analysis and 3D-modelling of Hi-C data using TADbit reveals structural features of the fly chromatin colors. PLOS Comp Bio 13(7) e1005665. doi:10.1371/journal.pcbi.1005665 
#  
#  \* Stefano, M.D., Stadhouders, R., Farabella, I., Castillo, D., Serra, F., Graf, T., Marti-Renom, M.A. (2020). Transcriptional activation during cell reprogramming correlates with the formation of 3D open chromatin hubs. Nature Communications, 11(1). ISSN 20411723.; doi: 10.1038/s41467-020-16396-1

# # Libraries and functions 



from taddyn import Chromosome, Experiment
import numpy as np
import os
import shutil
from taddyn.modelling.impoptimizer  import IMPoptimizer
import copy
import sys
import datetime


# # Run

# ## Define additional paths 



scriptsPath = basePath + 'code/modellingScripts/'
tempOut = basePath + 'temporal/'


# ## Import additional libraries  



sys.path.append(basePath + 'code')
import fileHandling
import metrics


# ## Get matrix paths 



addToBase = ''
if runFastTest:
    addToBase = 'fastTest/'
    # for a fast test we will only get 25 models.
    # WARNING: This number of models is only useful for a test, not for selecting the 
    # best comination of parameters in the optimisation
    nmodels = 25
    print 'Running FAST TEST modelling'
matricesLength, regionsAll, matrices = fileHandling.getMatricesPaths(basePath + addToBase, 
                                                                     starting='Matrix')


# ## Preparing optimization run

# This code will create a set of commands to run the optimization for the given combinations of parameters

# ### Get combinations of paramteres and run commands 



combinations = {}
for cell in matrices:
    combinations[cell] = {}
    for regi in matrices[cell]:
        matPath = matrices[cell][regi]
        c_rangeDot= metrics.dotToC('_'.join(str(j) for j in dcutoff_range))
        combinations[cell][regi] = metrics.getParamCombi(lowfreq_arange, m_range, c_rangeDot, upfreq_range, 
                         scriptsPath, dcutoff_range, matPath, jobTime, nmodels,
                         tempOut, cpu=cpus)


# ### Stimate total modelling time

# It is important to ensure that the number of models you want to build in the optimization will be generated in a feasible calculation time. The calculation time depends on both the number of particles of the models (which depends on the length of the region and the resolution of the experiment) and the number of CPUs available. Please, pay special attention to the output of the following sections and adjust the region to be modelled, its resolution, or the number of CPUs accordingly.
# 
# WARNING: We discorage users to work with models bigger than 1000 bins/particles due to the exponential increase in modelling timing.



for regi in regionsAll:
    print('--- %s ---' %regi)
    ncell = len(matrices.keys())
    print('Counting %s cells' %ncell)
    ncombi = len(combinations[cell][regi]) * ncell
    
    # each combination has n models, so we need to multiply
    totalModels = ncombi * nmodels
    
    # Finally we get the total time to do all the models in each region
    timePerModel = metrics.stimateTime(matricesLength[regi])
    totalTime = totalModels * timePerModel
    totalTime2 = str(datetime.timedelta(seconds=totalTime))
    
    print('%s models will be computed, in a median stimated time (with 1 CPU) of %s' %(
                                        totalModels, totalTime2))
    
    totalTime2 = str(datetime.timedelta(seconds=totalTime/cpus))
    print("Stimated time with assigned number of %s CPU's: %s" %(cpus,
                                        totalTime2))
    print('')


# ### Run optimization 

# Output files with the correlation of the models from each of the combinations will be stored in the same folder as each of the matrices.  
# 
# Example: opt_LF-1.0UF0.0C100.0-200.0-300.0Mdis200_5000bp.txt  
# 
# In the same directory, a folder called "lammpsSteps" will be created to store the models for the parameter sets that yet need to be started or whose run ended without generating all the required models.



for cell in combinations:
    print('## %s ##' %(cell))
    for regi in combinations[cell]:
        print('--- %s ---' %regi)
        for nc, combi in enumerate(combinations[cell][regi]):
            print('Combination %s' %(nc))
            get_ipython().system(u' python {combi}')


# ### Continuing optimization run

# NOTE: Sometimes the simulation run may not get to the end. In these cases, they are restarted from a different initial random positioning of the particles. By default this step will be repeated at most 10 times.



## we will rerun the models until we finish all them or we reach 
# 10 steps
combinations_t = copy.deepcopy(combinations)
for cell in combinations_t:
    print('## %s ##' %(cell))
    for regi in combinations_t[cell]:
        print('--- %s ---' %regi)
        matPath = matrices[cell][regi]
        nchecks = 0
        while len(combinations_t[cell][regi]) > 0 or nchecks < 10:
            combinations2 = copy.copy(combinations_t[cell][regi])
            combinations_t[cell][regi] = []
            for nc, combi in enumerate(combinations2):
                # get paths and check if the modelling finished
                path = '/'.join(matPath.split('/')[:-1]) + '/'
                jobName = 'LF%sUF%sMdis%s_%sbp' %(combi.split()[2], combi.split()[8], combi.split()[6], 
                                            matPath.split('_')[-1][:-2])
                keep_restart_out_dir = path + 'lammpsSteps/jobArray_%s/' %jobName
                if os.path.isdir(keep_restart_out_dir):
                    print('Combination %s' %(nc))
                    # then it didnt finished
                    combinations_t += [combi]
                    get_ipython().system(u' python {combi}')
            nchecks += 1
        
        


# ## Removing temporal folders

# The content inside of the generated temporal folders inside $tempOut is usually deleted after the modelling successfully finishes, not the folders itself though. Besides, if the modelling process breaks for any reason, the temporal files will remain there. TADdyn modelling generates a lot of temporal files so, in order to ensure that they don't accumulate, we will remove the container folders after each optimization process.



tempFolders = os.listdir(tempOut)




for t in tempFolders:
    print(t)
    shutil.rmtree(tempOut + t)

