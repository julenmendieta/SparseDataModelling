from math                           import sqrt
import numpy as np
import itertools
try:
    from pytadbit.modelling.structuralmodels  import load_structuralmodels,StructuralModels
    from pytadbit.modelling.structuralmodel import StructuralModel
    from pytadbit.modelling.impmodel import IMPmodel
    from pytadbit                     import HiC_data
    from pytadbit.parsers.hic_bam_parser import filters_to_bin
    from sklearn.metrics import calinski_harabaz_score
    from pysam                           import AlignmentFile
except:
    print('TADbit libraries not loaded, is ok if working with TADdyn alone')
import os, sys
import copy
from os import listdir
from scipy.cluster.hierarchy import fcluster
from scipy.cluster.hierarchy import linkage as linkage_sci
from collections import defaultdict
import warnings
# norm part
from collections                     import OrderedDict
import cPickle as pickle


## Normalisation part
# Obtain multiContacts from file just in a pairwise manner
def goThroughReads(section_pos, line, interPerBin,
                                resol, nRead=0):
    '''
    Function to obtain interaction frecuencies per bin from 
        normal TSV with no multiContact data spected
        
    :param section_pos: Dictionary with the chromosomes
        reference name as keys and a list or tuple of
        the range of bins in which they lie. Ej.:
        {'chr1':(0, 5000)}. Is 0 index, and last bin 
        from range, is not included inside, so next
        chromosome could be {'chr2':(5000, 10000)}
    :param line: list or tuple with:
        [chr1, startPos1, chr2, startPos2]
    :param interPerBin: defaultdict(int) with the number
        of times each bin interacts
    :param resol: Resolution at wich we are going to be 
        normalising our data
    :param 0 nRead: Integer indicating number of reads 
        counted
    
    '''
   
    # store each apparition of a fragment in a concatemer
    #we store the bin of the mapping start position
    fragment1 = (int(line[1]) / resol) + section_pos[line[0]][0]
    interPerBin[fragment1] += 1
    
    fragment2 = (int(line[3]) / resol) + section_pos[line[2]][0]
    interPerBin[fragment2] += 1
    
    # New concatemer seen
    nRead += 1
    
    return interPerBin, nRead

# Open tsv or bam and obtain contact frecuencies in a pairwise manner
def getInteractionsPerBin(infile, resol, locusCh=False,
                         regRange = False, returnNread = False,
                         filter_exclude=(1, 2, 3, 4, 6, 7, 8, 9, 10)):
    '''
    Function to get the number of concatemers were a bin of interest is 
        appearing (is bin based, so all fragments which start inside of
        a bin margin will be joined)
        
    :param infile: Path to the input file. If TADbit style TSV, be sure 
        it ends with .tsv. If usual BAM file, be sure it is sorted, the
        index is located in the same folder, and it ends with .bam
    :param resol: Resolution at wich we are going to be normalising our
        data
    :param False locusCh: Set to True if you want to return just data to
        normalise the bin between the smallest and biggest binned coordinate
        in regRange. Even if True the whole bam file must be cheked, so
        wont safe any time.
    :param False regRange: list or tuple with the first and last binned
        coordinates we wont to retrieve.
    :param False returnNread: wether you want or not nRead to be
        returned. It returns the number of reads check.
    :param (1, 2, 3, 4, 6, 7, 8, 9, 10) filter exclude: filters to define the
        set of valid pair of reads. Just valid for TADbit style BAMfiles. If
        working with already filtered non TADbit BAM set filter_exclude = ()
    '''
    
    # change filter exclude to binary
    if not isinstance(filter_exclude, int):
        filter_exclude = filters_to_bin(filter_exclude)
    # variable to store id
    prev = ''
    interPerBin = defaultdict(int)
    
    nRead = 0
    
    if infile.endswith('.bam'):
        # get section positions 
        bamfile = AlignmentFile(infile, 'rb')
        bam_refs = bamfile.references
        bam_lengths = bamfile.lengths

        sections = OrderedDict(list(zip(bam_refs,
                                   [x for x in bam_lengths])))
        total = 0
        section_pos = OrderedDict()
        for crm in sections:
            section_pos[crm] = (total, total + (sections[crm] / resol + 1))
            total += (sections[crm] / resol + 1)
            
            
        # check if this BAM file is not TADbit style
        if 'Hicup' in bamfile.text:
            print 'It seems this BAM file was produced outside TADbit, make \
sure it has already been filtered'
            if filter_exclude != ():
                print 'Consider changing filter_exclude so its value is () \
or you might get no reads'
                
        bamfile.close()

        bamfile = AlignmentFile(infile, 'rb')
        for r in bamfile.fetch():
            # Check if it is among positions to be filtered
            # BEWARE that it follows TADbit format
            if r.flag & filter_exclude:
                continue
            crm1 = r.reference_name
            pos1 = r.reference_start + 1
            crm2 = r.next_reference_name
            pos2 = r.next_reference_start + 1
            
            line = [crm1, pos1, crm2, pos2]
            interPerBin, nRead = goThroughReads(section_pos, line, interPerBin,
                                    resol, nRead=nRead)
        bamfile.close()
    

    elif infile.endswith('.tsv'):
        with open(infile, 'r') as f:
            # get chromosome lengths
            sections = OrderedDict()
            while True:
                line = f.readline()
                if line.startswith('#'):
                    line = line.split()
                    if line[1] == 'CRM':
                        sections[line[2]] = int(line[3])
                        
                elif line.startswith('@'):
                    pass
                else:
                    break
                    
            # create binned positioning for chromosomes
            total = 0
            section_pos = OrderedDict()
            for crm in sections:
                section_pos[crm] = (total, total + (sections[crm] / resol + 1))
                total += (sections[crm] / resol + 1)

            # iterate over reads
            line = line.split()
            line = [line[1], line[2], line[7], line[8]]

            # Run current line (first one)
            interPerBin, nRead = goThroughReads(section_pos, line, interPerBin,
                                    resol, nRead=nRead)

            # Go for next lines
            for line in f:
                line = line.split()
                line = [line[1], line[2], line[7], line[8]]
                
                interPerBin, nRead = goThroughReads(section_pos, line, interPerBin,
                                    resol, nRead=nRead)

    
     # Get genomic coordinates for region of interest
    if locusCh == True:
        regionStart = min(regRange)
        regionEnd = max(regRange) 

        ## modify if we want data from all genome
        # Get all the fragments that start inside this coordinates
        keys = [k for k in concatemers.keys() if (regionStart <= 
                                                            k <= 
                                                            regionEnd)]
        regInterPerBin = defaultdict(int)
        for k in keys:
            regInterPerBin[k] += interPerBin[k]
            
    # Or not
    else:
        regInterPerBin = interPerBin
            

    if returnNread == False:
        return regInterPerBin
    else:
        return regInterPerBin, nRead


# Normalise by frequencies given the presence of each interacting fragment                                           
def frequenciesNorm(hic_data, resol, regRange, concatemersBin, multResult=100, 
                    keep=False, mininter=0, positionAdjust=0):
    
    '''
    param 0 positionAdjust: In case the positions from concatemersBin are taking 
        into account bining from the whole genome, but we just load in hic_data
        one chromosome. Here concatemersBin will be substracted added to the 
        positions in regRange in order to compensate this
    
    '''
    # create HiC data for normalised interactions
    dict_sec = hic_data.sections
    genome_seq = hic_data.chromosomes
    size = sum(genome_seq[crm] for crm in genome_seq)
    norm_data = HiC_data((), size, genome_seq, dict_sec, resolution=resol)

    # Will remove from divider concatemers counted twice
    if keep == False:
        for nbin1, bin1 in enumerate(regRange):
            for bin2 in regRange[nbin1:]:
                # If diagonal or bellow sed interactions
                if bin1 == bin2 or hic_data[bin1, bin2] <= mininter:
                    pass  # Leave it as zero

                else:
                    # get divider
                    #if concatemersBin[bin1] == 0:
                    #    if concatemersBin[bin2] == 0:
                    #        divider = 1
                    #    else:
                    #        divider = concatemersBin[bin2]
                    #elif concatemersBin[bin2] == 0:
                    #    divider = concatemersBin[bin1]
                    #else:
                    divider = (concatemersBin[bin1 + positionAdjust] + 
                               concatemersBin[bin2 + positionAdjust] - 
                               hic_data[bin1, bin2])

                    #divider = float(concatemersBin[bin1] + concatemersBin[bin2])

                    #if divider == 0:
                    #    divider = 1
                    # if divider is 0 it means concatemersBin was taken with another index
                    #ie just checking a chromosome, or whole genome and here just 
                    #normalising a file were we loaded one chromosome
                    # if both are zero 
                    norm_data[bin1, bin2] = (hic_data[bin1, bin2] / float(divider)) * multResult
                    norm_data[bin2, bin1] = norm_data[bin1, bin2]

    else:
        for ke in keep:
            # If diagonal or bellow sed interactions
            if (ke[0] == ke[1] or hic_data[ke[0], ke[1]] <= mininter or
                (ke[0] not in regRange or ke[1] not in regRange)):
                pass  # Leave it as zero
            else:
                # get divider
                #if concatemersBin[ke[0]] == 0:
                #    if concatemersBin[ke[1]] == 0:
                #        divider = 1
                #    else:
                #        divider = concatemersBin[ke[1]]
                #elif concatemersBin[ke[1]] == 0:
                #    divider = concatemersBin[ke[0]]
                #else:
                divider = (concatemersBin[ke[0] + positionAdjust] + 
                               concatemersBin[ke[1] + positionAdjust] - 
                               hic_data[ke[0], ke[1]])

                #divider = float(concatemersBin[bin1] + concatemersBin[bin2])

                #if divider == 0:
                #    divider = 1
                # if both are zero 
                norm_data[ke[0], ke[1]] = hic_data[ke[0], ke[1]] / float(divider)
                norm_data[ke[1], ke[0]] = norm_data[ke[0], ke[1]]
        
    return norm_data


def getSectionPos(infile, resol):
    bamfile = AlignmentFile(infile, 'rb')
    bam_refs = bamfile.references
    bam_lengths = bamfile.lengths

    sections = OrderedDict(list(zip(bam_refs,
                               [x for x in bam_lengths])))
    total = 0
    section_pos = OrderedDict()
    for crm in sections:
        section_pos[crm] = (total, total + (sections[crm] / resol + 1))
        total += (sections[crm] / resol + 1)
        
    bamfile.close()
    return section_pos

def getCaptures(allPdir, chrom, regionStart, regionEnd, resol,
                       mlength):
    promList = set()
    with open(allPdir, 'r') as f:
        for line in f:
            if line.startswith(chrom):
                line2 = line.rstrip().split('\t')
                if (line2[0]) == chrom:
                    if (regionStart <= int(line2[1])/resol*resol <= regionEnd or 
                        regionStart <= int(line2[2])/resol*resol <= regionEnd):
                        pos1 = (int(line2[1]) - regionStart) / resol
                        pos2 = (int(line2[2]) - regionStart) / resol
                        for i in range(pos1, pos2 + 1):
                            # if capture located in the matrix
                            if i < mlength and i >= 0:
                                promList.add(i)
                            else:
                                print i
    return promList


## optimisation ##
def cToDot(text):
    newText = ''
    for t in text:
        if t == 'c':
            newText += '.'
        elif t == 'n':
            newText += '-'
        else:
            newText += t
    return newText

def dotToC(text):
    newText = ''
    for t in text:
        if t == '.':
            newText += 'c'
        elif t == '-':
            newText += 'n'
        else:
            newText += t
    return newText

def getParamCombi(lowfreq_arange, m_range, c_rangeDot, upfreq_range, 
                 scriptsPath, dcutoff_range, matPath, jobTime, nmodels,
                 tempOut, cpu=False):
    # Create a counter for the different lammps output directories
    allCombi = []
    ## Create the file with the commands to be run in the array
    #fout=open(runfile,'w')
    for x in lowfreq_arange:
        for m in m_range:
            # check if we have a dcutoff small enough compared to maxdist to allow running
            # we allow an overlap of 50nm between maxdist and dcutoff
            if m - float(cToDot(c_rangeDot).split('_')[0]) >= -50:
                for u in upfreq_range:
                    if u >= x:
                        cmd=''
                        cmd+='%s01_NR_optimisation.py -l %s '%(scriptsPath, x)
                        cmd+= '-d %s '%c_rangeDot
                        cmd+= '-m %s '%m
                        cmd+= '-u %s '%u
                        cmd+= '-p %s '%matPath
                        cmd+= '-t %s '%jobTime
                        cmd+= '-nm %s '%str(nmodels)
                        if cpu != False:
                            cmd+= '-cpu %s '%str(cpu)
                        cmd+= '-tp %s'%(tempOut)
                        #cmd+= '\n'
                        allCombi += [cmd]
                        #fout.write(cmd)
    #fout.close()
    return allCombi

def stimateTime(nparticle):
    return 0.001*nparticle**2 + (-0.135*nparticle) + 48.116


def PCHiC_filtering(exp, index=0):
    zeros = {}
    for i in range(exp.size):
        lineAdd = 0
        for j in xrange(exp.size):
            lineAdd += exp.norm[index]['matrix'].get(i * exp.size + j, 0)
        if lineAdd == 0:
            zeros[i] = 0
    exp._zeros = [zeros]


def getModellingCommand(GeneralOptimOutPath, tempOut, jobTime,
                       nmodels, ncpu, outputAppart, scriptsPath, step=1):
    
    cmds = []
    times = 0
    with open(GeneralOptimOutPath + 'modellinParams.txt', 'r') as f:
        for line in f:
            if len(line) != 0:
                line = line.split()
                for ste in range(step):
                    nmodels2 = nmodels / step
                    if nmodels2 == 0:
                        sys.exit('There are more steps than models')

                    if ste == step - 1:
                        nmodels2 += nmodels % step                    
                    
                    matPath = line[0]
                    lowfreq = line[1]
                    upfreq = line[2]
                    d_cutoff= line[3]
                    maxdist= line[4]
                    
                    # modify tempout accordingly
                    try:
                        iniStep = int(tempOut.split('_')[-1])
                        tempOut = '_'.join(tempOut.split('_')[:-1]) + '_%s' %(ste + iniStep + 1)
                    except:
                        tempOut = tempOut + '_%s' %(ste + 1)

                    # get the parameters
                    cmd = '-l %s -d %s -m %s -u %s -lf %s -p %s -t %s -nm %s -ncpu %s' %(
                                                lowfreq,
                                                d_cutoff,
                                                maxdist,
                                                upfreq,
                                                tempOut,
                                                matPath,
                                                jobTime,
                                                nmodels2,
                                                ncpu)
                    if outputAppart == True:
                        pathOut = '/'.join(matPath.split('/')[:-4] + 
                                           ['models'] + 
                                           matPath.split('/')[-3:-1]) + '/'
                        cmd += ' -po %s' %pathOut

                    # add the paths
                    cmd = '%s03_NR_runmodelling.py %s' %(scriptsPath, cmd)
                    cmds += [cmd]
                    
                    # get number of particles and models for timing
                    _, start, end = matPath.split('/')[-1].split('_')[3].split('-')
                    resol = int(matPath.split('/')[-1].split('_')[-1][:-2])
                    longi = ((int(end) - int(start)) / resol) + 1
                    time = (stimateTime(longi) * nmodels2) / min(ncpu, nmodels2)
                    times += time
    return cmds, times
    

def taddynToTadbit(fi):
    with open(fi, 'rb') as pickle_file:
        ensemble_of_models = pickle.load(pickle_file)
    if len(ensemble_of_models['stages']) == 0:
        models_stage = dict((i, StructuralModel(ensemble_of_models['models'][mod]))
                        for i, mod in enumerate(ensemble_of_models['models']))
        model1 =StructuralModels(
                ensemble_of_models['loci'], models_stage, {}, 
                resolution=ensemble_of_models['resolution'], original_data=ensemble_of_models['original_data'],
                zscores=ensemble_of_models['zscores'][0], config=ensemble_of_models['config'],
                zeros=ensemble_of_models['zeros'][0],restraints=ensemble_of_models['restraints'])
        
        # if naming is like in defined format will try to recover some info
        fi2 = fi.split('/')[-1]
        if len(fi2.split('_')) == 3 and '_C' in fi2 and 'Res' in fi2:
            cell = fi2.split('_')[0]
            resol = fi2.split('Res')[-1].split('.')[0]
            region = fi2.split('_')[1]
            identifier = '%s_%s' %(cell, region)
        else:
            resol = 0
            identifier = 'NA'
            cell = 'NA'
            region = 'NA'
    
        # Define description and set index values
        description = {'identifier'        : identifier,
                       'chromosome'        : ['NA'],
                       'start'             : [0],
                       'end'               : [0],
                       'species'           : 'NA',
                       'restriction enzyme': 'NA',
                       'cell type'         : cell,
                       'experiment type'   : region,
                       'resolution'        : resol,
                       'assembly'          : 'NA'}
        
        for nmo, mo in enumerate(model1):
            # change rand_init values for actual index and add descriptions
            mo['index'] = nmo
            mo['description'] = description
            
    
    return model1


##
def square_distance_to(mod, part1, part2):
        """
        :param part1: index of a particle in the model
        :param part2: external coordinate (dict format with x, y, z keys)
        :returns: square distance between one point of the model and an external
           coordinate
        """
        #print part1,part2
        return ((mod['x'][part1] - part2[0])**2 +
                (mod['y'][part1] - part2[1])**2 +
                (mod['z'][part1] - part2[2])**2)
        
def radius_of_gyration(mod):
        """
        Calculates the radius of gyration or gyradius of the model
        Defined as:
        .. math::
          \sqrt{\\frac{\sum_{i=1}^{N} (x_i-x_{com})^2+(y_i-y_{com})^2+(z_i-z_{com})^2}{N}}
        with:
        * :math:`N` the number of particles in the model
        * :math:`com` the center of mass
        :returns: the radius of gyration for the components of the tensor
        """
        longi = len(mod['x'])
        
        com = center_of_mass(mod)

        rog = sqrt(sum(square_distance_to(mod, i,
                                                (com['x'], com['y'], com['z']))
                       for i in xrange(longi)) / longi)
        return rog
    
def center_of_mass(model, selection=False):
        """
        Gives the center of mass of a model
        :returns: the center of mass of a given model
        """
        if selection == False:
            longi = len(model['x'])
            selection = xrange(longi)
        else:
            longi = len(selection)
            
        r_x = sum(model['x'][s] for s in selection)/longi
        r_y = sum(model['y'][s] for s in selection)/longi
        r_z = sum(model['z'][s] for s in selection)/longi
        return dict((('x', r_x), ('y', r_y), ('z', r_z)))


def square_distance_to2(pos1, pos2):
        """
        :param part1: index of a particle in the model
        :param part2: external coordinate (dict format with x, y, z keys)
        :returns: square distance between one point of the model and an external
           coordinate
        """
        #print part1,part2
        return ((pos1['x'] - pos2['x'])**2 +
                (pos1['y'] - pos2['y'])**2 +
                (pos1['z'] - pos2['z'])**2)
    
    
def median_absolute_deviation(x, axis=0, center=np.median, 
                              scale=1.0):
    r"""
    Copied from scipy
    """
    x = np.asarray(x)

    # Wrap the call to center() in expand_dims() so it acts like
    # keepdims=True was used.
    med = np.expand_dims(center(x, axis=axis), axis)
    mad = np.median(np.abs(x - med), axis=axis)

    return mad / scale


def load_model_from_xyz(f_name, rand_init, resolution, scale=0.01):

    model = IMPmodel((('x', []), ('y', []), ('z', []), ('rand_init', str(rand_init)),
                      ('index', 0), ('objfun', rand_init), ('radius', resolution*scale/2)))
    model['description'] = f_name.split('.')[-2]

    take_this = True #just to reduce the number of points
    for line in open(f_name):
        if not line.startswith('# '):

            line_val = line.split()
            #print line_val
            model['x'].append(float(line_val[2]))
            model['y'].append(float(line_val[3]))
            model['z'].append(float(line_val[4]))

    return model

def load_all_models(dir_models, resolution=5000):

    #max_models = 9999999
    #if nbr_models != 'all':
    #    max_models = nbr_models
    results = []
    n=0
    for fn in listdir(dir_models):
        #print fn
        if fn.endswith("xyz"):
        #    if n == 2:
         #       continue

            k= fn.split(".")[1]
            #print fn,"model",k
            results.append((str(k), load_model_from_xyz(dir_models + fn,int(k), resolution)))
            n += 1

    nloci = 0
    models = {}

    for i, (_, m) in enumerate(
        sorted(results, key=lambda x: x[1]['objfun'])):
        models[i] = m

    nloci = len(m['x'])

    return StructuralModels(
            nloci, models, [], resolution,description='EBF1', zeros=tuple([1 for i in xrange(nloci)]))


def getDistancesFromInterest(orderCell, models, regiones, 
                             allClusters,modelsKeep,
                            clusterCheck, interAll, promAll,
                            checkAllBin=True, originalRegionPos=None):
    '''
    A function to select the elements stored in interAll to calculate
    distances from them to other particlesin the model
    :param True checkAllBin: Wether to check the distance between elements
        in interAll with enhancers and promoters (False) or with
        all particles in the model (True)
        
    '''
    ## first we create an index of bins
    allBins = {}
    if checkAllBin == True:
        for regi in regiones:
            regiResol = regiones[regi][3]
            allBins[regi] = range(((regiones[regi][2] - regiones[regi][1]) / regiResol) + 1)

    else:
        for regi in regiones:
            allBins[regi] = set()
        for regi in enhAll:
            for k in enhAll[regi].keys():
                allBins[regi].add(k)

        for regi in promAll:
            for k in promAll[regi].keys():
                allBins[regi].add(k)
        
    ## Now we prepare for storing distances between points
    regionsAll = sorted(regiones)
    distancesBtw = {}
    for region in regionsAll:
        for cell in orderCell:
            fi = models[cell][regi]
            flag = '%s_%s' %(cell, region)
            print flag
            # Get region index
            # get regions start and end position (and chromosome)
            chrom = regiones[region][0]
            regionStart = regiones[region][1]
            regionEnd = regiones[region][2]

            if not region in distancesBtw:
                distancesBtw[region] = {}

            ## If we set coordinates to trim the plot
            if originalRegionPos != None:
                fromOri = (originalRegionPos[region][1]/resol) - (regionStart / resol)
                fromEnd = (regionEnd / resol) - (originalRegionPos[region][2] / resol)

                regionStart2 = regionStart + (fromOri * resol)
                regionEnd2 = regionEnd - (fromEnd * resol)
            else:
                regionStart2 = regionStart
                regionEnd2 = regionEnd


            # Check if we have promoters of interest
            # Load models and do clustering
            mods=load_structuralmodels(fi)
            # keep the amount of selected models
            mods.define_best_models(modelsKeep)
            mods.clusters = allClusters[flag]

            # check presence of cluster of interest
            distancesBtw[regi][cell] = {}
            for cluster in clusterCheck:
                if cluster in mods.clusters.keys():
                    distancesBtw[regi][cell][cluster] = {}

                    # First we look for distances from promoters
                    distancesBtw[regi][cell][cluster]['promoter'] = {}
                    for keypro in sorted(interAll[region]['promoter']):
                        # get all posible combinations between promoters
                        otherP = [p for p in promAll[region].keys() if p != keypro]
                        other = sorted(list(allBins[regi]))
                        combinationsBtw = list(itertools.product([keypro], set(otherP + other)))

                        # get point to point distances
                        distancesBtw[regi][cell][cluster]['promoter'][keypro] = {}

                        for comp in combinationsBtw:
                            values = mods.median_3d_dist(comp[0] + 1, comp[1] + 1, 
                                                         cluster = cluster, plot=False, 
                                                         median=False)
                            distancesBtw[regi][cell][cluster]['promoter'][keypro][comp] = values 

                    # Then distances from each enhancer to all promoters
                    distancesBtw[regi][cell][cluster]['enhancer'] = {}
                    for keypro in sorted(interAll[region]['enhancer']):
                        # get all posible combinations between promoters
                        otherP = [p for p in promAll[region].keys() if p != keypro]
                        other = sorted(list(allBins[regi]))
                        combinationsBtw = list(itertools.product([keypro], set(otherP + other)))

                        # get point to point distances
                        distancesBtw[regi][cell][cluster]['enhancer'][keypro] = {}

                        for comp in combinationsBtw:
                            values = mods.median_3d_dist(comp[0] + 1, comp[1] + 1, 
                                                         cluster = cluster, plot=False,
                                                         median=False)
                            distancesBtw[regi][cell][cluster]['enhancer'][keypro][comp] = values 

                else:
                    print 'No cluster %s in %s' %(cluster, fi)
                    print 'Measures not taken'
                    
    return distancesBtw



def writeModelCmm(regionsAll, outdata, orderCell, models,
                 modelsKeep, toCluster=False, allClusters=None,
                 maxCluster=2, minModel=50):
    '''
    :param False toCluster: Wether to split by clusters or not
    :param None allClusters: Dictionary with the information about the
        clusters.
    :param :
    :param :
    :param :
    :param :

    '''
    indexList = {}
    passFilter = {}

    for regi in regionsAll:
        outpath2 = outdata + regi + '/'
        if not os.path.exists(outpath2):
            os.makedirs(outpath2)
        indexList[regi] = {}
        pos0 = 0
        for c in orderCell:
            fi = models[c][regi]

            mods = load_structuralmodels(fi)
            # keep the amount of selected models
            mods.define_best_models(modelsKeep)

            passFilter['%s_%s' %(c, regi)] = []
            # create index
            indexList[regi][c] = [pos0, 0]

            if toCluster:
                mods.clusters = copy.copy(allClusters['%s_%s' %(c, regi)])
                # will assign each model to a cluster
                for clu in mods.clusters:
                    for rmo in mods.clusters[clu]:
                        mods[rmo]['cluster'] = clu

                for clu in sorted(mods.clusters):
                    # check just the first 10 clusters and among them the ones with more than 50 models
                    if clu <= maxCluster and len(mods.clusters[clu]) >= minModel:
                        #print clu
                        passFilter['%s_%s' %(c, regi)].append(clu)

                        # get assigned index of each model
                        modindexes = []
                        for m in mods.clusters[clu]:
                            modindexes.append(mods[m]['index'])
                        # change rand_init index of models (it gives name to .xyz files)
                        for nm, m in enumerate(modindexes, pos0):
                            mods[m]['rand_init'] = str(nm)

                        # update cluster index
                        mods.clusters[clu] = [str(r+pos0) for r in range(len(modindexes))]
                        mods.write_xyz(outpath2, cluster=clu, )
                        pos0 = nm + 1
            else:
                clu = 1
                # get assigned index of each model
                modindexes = []
                for m in range(len(mods)):
                    modindexes.append(mods[m]['index'])
                # change rand_init index of models (it gives name to .xyz files)
                for nm, m in enumerate(modindexes, pos0):
                    mods[m]['rand_init'] = str(nm)

                # update cluster index
                mods.clusters[clu] = [str(r+pos0) for r in range(len(modindexes))]
                mods.write_xyz(outpath2, cluster=clu, )
                pos0 = nm + 1


            indexList[regi][c][1] = nm
    return indexList, passFilter

def getBetweenBinDist(model, combinations, modelCheckCluster=False):
    '''
    :param model: TADbit StructuralModels object
    :param combinations: List with pairs of bin locations whose
        distance will be obtained
    
    '''
    cluster_models = []
    if modelCheckCluster != False:
        cluster_models = model.clusters[modelCheckCluster]
    else:
        for cluster in model.clusters:
            cluster_models += model.clusters[cluster]
            
    
                
    distances = []
    for nmo, clum in enumerate(cluster_models):
        mo = model[clum]
        for bin1, bin2 in combinations:
            pos1 = {'x': mo['x'][bin1], 'y': mo['y'][bin1], 'z': mo['z'][bin1]}
            pos2 = {'x': mo['x'][bin2], 'y': mo['y'][bin2], 'z': mo['z'][bin2]}

            dist = sqrt((pos2['x'] - pos1['x'])**2 + (pos2['y'] - pos1['y'])**2 + (pos2['z'] - pos1['z'])**2)
            distances += [dist]
    return distances

def getExpressedDistances(selected, markers,regionsAll, models,
                            modelsKeep, cluster, allClusters):
    distancesAll = {}
    for genePos, geneName in selected:
        distancesAll[(genePos, geneName)] = {}
        for regi in regionsAll:
            distancesAll[(genePos, geneName)][regi] = {}
            for cell in models:
                combinations = itertools.product(markers[regi][cell], [genePos])

                mod = models[cell][regi]
                flag = '%s_%s' %(cell, regi)
                model = load_structuralmodels(mod)
                # keep the amount of selected models
                model.define_best_models(modelsKeep)
                model.clusters = allClusters[flag]

                distancesAll[(genePos, geneName)][regi][cell] = getBetweenBinDist(model, combinations, modelCheckCluster=cluster)
    return distancesAll

def getBetweenBinDist2(regionsAll, models, expresedBins, modelsKeep,
                      allClusters, cluster):
    distancesPerModel = {}
    for regi in regionsAll:
        distancesPerModel[regi] = {}
        for cell in models:
            combinations = list(itertools.combinations(expresedBins[regi][cell], 2))
            
            mod = models[cell][regi]
            flag = flag = '%s_%s' %(cell, regi)

            model = load_structuralmodels(mod)
            # keep the amount of selected models
            model.define_best_models(modelsKeep)
            model.clusters = allClusters[flag]

            distancesPerModel[regi][cell] = {}
            cluster_models = model.clusters[cluster]
            # we store distance in the sorted order of combinations between expressed
            # elements
            for nmo, clum in enumerate(cluster_models):
                mo = model[clum]
                distances = []
                for bin1, bin2 in combinations:
                    pos1 = {'x': mo['x'][bin1], 'y': mo['y'][bin1], 'z': mo['z'][bin1]}
                    pos2 = {'x': mo['x'][bin2], 'y': mo['y'][bin2], 'z': mo['z'][bin2]}

                    dist = sqrt((pos2['x'] - pos1['x'])**2 + 
                                (pos2['y'] - pos1['y'])**2 + 
                                (pos2['z'] - pos1['z'])**2)
                    distances += [dist]

                distancesPerModel[regi][cell][nmo] = distances
    return distancesPerModel

## co-occurrence matrices ##
def getCooccurrenceMatrices(regionsAll, models, distancesPerModel, 
                            expresedBins, clusterRange=(2,10)):
    posConverter = {}
    coOcurrMatrices = {}
    definedClusters = {}

    for regi in regionsAll:
        posConverter[regi] = {}
        coOcurrMatrices[regi] = {}
        definedClusters[regi] = {}
        for cell in models:
            maxis = []
            coOcurrMatrices[regi][cell] = {}

            combinations = list(itertools.combinations(expresedBins[regi][cell], 2))
            sortPos = sorted(expresedBins[regi][cell])
            longi = len(sortPos)
            posDict = {j:i for i,j in enumerate(sortPos)}
            posConverter[regi][cell] = {i:j for i,j in enumerate(sortPos)}

            # cluster distance matrix
            distanceMatrix_Coocurrence = [[0 for i in range(longi)] for j in range(longi)]

            for nmo in distancesPerModel[regi][cell]:
                templist = distancesPerModel[regi][cell][nmo]

                # get distance matrix for one model
                distanceMatrix_temp = [[0 for i in range(longi)] for j in range(longi)]
                for nc, comb in enumerate(combinations):
                    posi = posDict[comb[0]]
                    posj = posDict[comb[1]]
                    distanceMatrix_temp[posi][posj] = templist[nc]
                    distanceMatrix_temp[posj][posi] = templist[nc]

                # get clusters in that model
                matrix = np.array(distanceMatrix_temp)
                #matrixTest = euclidean_distances(matrix)

                # get preferred number of clusters
                ones = []
                twos = []
                for k in range(clusterRange[0], clusterRange[1]):
                    cluster = fcluster(linkage_sci(np.array(matrix), method='ward',
                                               metric='euclidean'),k,'maxclust')
                    #print k, metrics.calinski_harabaz_score(matrix, labels)
                    ones += [k]
                    twos += [calinski_harabaz_score(matrix, cluster)]
                m = max(twos)
                maxi = [i for i, j in enumerate(twos) if j == m]
                maxi = ones[maxi[0]]
                maxis += [maxi]

                # do this number of clusters
                cluster = fcluster(linkage_sci(np.array(matrix), method='ward',
                                               metric='euclidean'),maxi,'maxclust')

                # get bins in each cluster
                clusterBins = defaultdict(list)
                for nc, clu in enumerate(cluster):
                    clusterBins[clu] += [posConverter[regi][cell][nc]]

                # assig coocurrence to ensemble matrix
                for clu in clusterBins:
                    combi_clu = list(itertools.combinations(clusterBins[clu], 2))
                    for comb1, comb2 in combi_clu:
                        posi = posDict[comb1]
                        posj = posDict[comb2]
                        distanceMatrix_Coocurrence[posi][posj] += 1
                        distanceMatrix_Coocurrence[posj][posi] += 1

            # store coocurrence matrix
            coOcurrMatrices[regi][cell] = copy.copy(distanceMatrix_Coocurrence)
            # convert coocurrence in percentage
            nmodel = len(distancesPerModel[regi][cell])
            for i in range(len(coOcurrMatrices[regi][cell])):
                for j in range(i, len(coOcurrMatrices[regi][cell])):
                    coOcurrMatrices[regi][cell][i][j] = coOcurrMatrices[regi][cell][i][j] / float(nmodel) * 100
                    coOcurrMatrices[regi][cell][j][i] = coOcurrMatrices[regi][cell][i][j]

            definedClusters[regi][cell] = maxis

    return definedClusters, coOcurrMatrices, posConverter


def getClusterPositions(coOcurrMatrices, posConverter,
                        newSignalPos, toshow=True):
    clustersBinPos = {}
    for regi in coOcurrMatrices:
        clustersBinPos[regi] = {}
        for cell in coOcurrMatrices[regi]:
            matrix = np.array(coOcurrMatrices[regi][cell])

            # get preferred number of clusters
            ones = []
            twos = []
            maxis = []
            for k in range(2, 10):
                # we will hide some warnings taht appear in this part
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    cluster = fcluster(linkage_sci(np.array(matrix), method='ward',
                                               metric='euclidean'),k,'maxclust')
                #print k, metrics.calinski_harabaz_score(matrix, labels)
                ones += [k]
                twos += [calinski_harabaz_score(matrix, cluster)]
            m = max(twos)
            maxi = [i for i, j in enumerate(twos) if j == m]
            maxi = ones[maxi[0]]
            maxis += [maxi]

            # do this number of clusters
            # we will hide some warnings taht appear in this part
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                cluster = fcluster(linkage_sci(np.array(matrix), method='ward',
                                               metric='euclidean'),maxi,'maxclust')
            #print cell, min(cluster), max(cluster)

            clustersBinPos[regi][cell] = defaultdict(list)
            for nbin in range(len(matrix)):
                clustersBinPos[regi][cell][cluster[nbin]].append(nbin)

            for clu in clustersBinPos[regi][cell]:
                clustersBinPos[regi][cell][clu] = [posConverter[regi][cell][c] 
                                                   for c in clustersBinPos[regi][cell][clu]]

    clusterGeneNames = {}
    for regi in clustersBinPos:
        clusterGeneNames[regi] = {}
        print '#' * 30, regi, '#' * 30
        for cell in coOcurrMatrices[regi]:
            clusterGeneNames[regi][cell] = {}
            print ':' * 10, cell, ':' * 10
            for clu in clustersBinPos[regi][cell]:
                clusterGeneNames[regi][cell][clu] = []
                print '_' * 5, clu, '_' * 5
                ids = []
                for p in clustersBinPos[regi][cell][clu]:
                    ids += newSignalPos[regi][cell][p]
                    clusterGeneNames[regi][cell][clu] += newSignalPos[regi][cell][p]
                if toshow == True:
                    print '\n'.join(ids)
                    
    return clustersBinPos, clusterGeneNames


def clusterCOMdistances(regiones,models, allClusters, clustersBinPos,
                        modelCheckCluster, modelsKeep):

    ## First we get the center of mass of each cluster in each model
    centersOfMass = {}
    for regi in regiones:
        centersOfMass[regi] = {}
        for cell in models:
            centersOfMass[regi][cell] = {}
            mod = models[cell][regi]
            flag = '%s_%s' %(cell, regi)

            model = load_structuralmodels(mod)
            # keep the amount of selected models
            model.define_best_models(modelsKeep)
            model.clusters = allClusters[flag]
            model.align_models(in_place=True)

            for clu in clustersBinPos[regi][cell]:
                # get particle selection
                selection = clustersBinPos[regi][cell][clu]

                # will check just models in the selected clusters
                cluster_models = []
                if modelCheckCluster != False:
                    cluster_models = model.clusters[modelCheckCluster]
                else:
                    for cluster in model.clusters:
                        cluster_models += model.clusters[cluster]


                # get gyration radi of selected particles in selected models 
                mas_centers = []
                for nmo, clum in enumerate(cluster_models):
                    mo = model[clum]
                    # get center of mass
                    com = center_of_mass(mo, selection=selection)
                    #
                    mas_centers += [com]

                centersOfMass[regi][cell][clu] = copy.copy(mas_centers)

    # Then we get the distances between the COM of all the clusters
    # in each model
    distancesCluster = {}
    for regi in centersOfMass:
        distancesCluster[regi] = {}
        for cell in centersOfMass[regi]:
            distancesCluster[regi][cell] = defaultdict(list)
            cluster_combinations = list(itertools.combinations(centersOfMass[regi][cell], 2))
            for nmodel in range(len(centersOfMass[regi][cell][1])):
                for nc, (clu1, clu2) in enumerate(cluster_combinations):
                    pos1 = centersOfMass[regi][cell][clu1][nmodel]
                    pos2 = centersOfMass[regi][cell][clu2][nmodel]

                    dist = sqrt(square_distance_to2(pos1, pos2))
                    distancesCluster[regi][cell][cluster_combinations[nc]] += [dist]

    ## compute the mean
    meandistancesCluster = {}
    for regi in regiones:
        meandistancesCluster[regi] = {}
        for cell in distancesCluster[regi]:
            meandistancesCluster[regi][cell] = {}
            for compare in distancesCluster[regi][cell]:
                meandistancesCluster[regi][cell][compare] = np.mean(distancesCluster[regi][cell][compare])

    return meandistancesCluster, centersOfMass

def rankExpression(clustersBinPos, newSignalPos, newSignal):
    clusterMeanExp = {}
    for regi in clustersBinPos:
        clusterMeanExp[regi] = {}
        for cell in clustersBinPos[regi]:
            clusterMeanExp[regi][cell] = {}
            for clu in clustersBinPos[regi][cell]:
                cluExp = []
                for particle in clustersBinPos[regi][cell][clu]:
                    genes = newSignalPos[regi][cell][particle]
                    for ge in genes:
                        cluExp += [newSignal[regi][cell][ge]]
                clusterMeanExp[regi][cell][clu] = np.mean(cluExp)

    # now we convert sum into ranking
    from scipy.stats import rankdata
    rankedSum = {}
    for regi in clusterMeanExp:
        rankedSum[regi] = {}
        for cell in clusterMeanExp[regi]:
            rankedSum[regi][cell] = {}

            ks = clusterMeanExp[regi][cell].keys()
            rnaki = rankdata([clusterMeanExp[regi][cell][k] for k in ks])
            for nk, k in enumerate(ks):
                rankedSum[regi][cell][k] = int(rnaki[nk])

    return rankedSum


def getExpressionAndCOMd(centersOfMass, models, modelsKeep, clustersBinPos,
                        newSignalPos, newSignal, allClusters, modelCheckCluster):

    expeCOMdist = {}
    for regi in models[models.keys()[0]]:
        expeCOMdist[regi] = {}
        for cell in models:
            expeCOMdist[regi][cell] = {}
            mod = models[cell][regi]
            flag = '%s_%s' %(cell, regi)

            model = load_structuralmodels(mod)
            # keep the amount of selected models
            model.define_best_models(modelsKeep)
            model.clusters = allClusters[flag]
            model.align_models(in_place=True)

            for clu in clustersBinPos[regi][cell]:
                expeCOMdist[regi][cell][clu] = {'expression':[], 'distanceToCOM':[]}
                # get particle selection
                selection = clustersBinPos[regi][cell][clu]

                # will check just models in the selected clusters
                cluster_models = []
                if modelCheckCluster != False:
                    cluster_models = model.clusters[modelCheckCluster]
                else:
                    for cluster in model.clusters:
                        cluster_models += model.clusters[cluster]


                # get distance from mass center to each element from cluster
                clu_distances = defaultdict(list)
                for nmo, clum in enumerate(cluster_models):
                    mass_center = centersOfMass[regi][cell][clu][nmo]
                    mo = model[clum]

                    # get distance from each element
                    for sel1 in selection:
                        pos1 = mass_center
                        pos2 = {'x':mo['x'][sel1], 'y':mo['y'][sel1], 'z':mo['z'][sel1]}

                        dist = sqrt(square_distance_to2(pos1, pos2))

                        # expression
                        genes = newSignalPos[regi][cell][sel1]
                        expe = 0
                        for ge in genes:
                            expe += newSignal[regi][cell][ge]
                        if expe != 0:
                            clu_distances[sel1] += [dist]

                # now we store the median distance per each particle
                for sel1 in selection:
                    expeCOMdist[regi][cell][clu]['distanceToCOM'] += [np.median(clu_distances[sel1])]

                    genes = newSignalPos[regi][cell][sel1]
                    expe = 0

                    for ge in genes:
                        expe += newSignal[regi][cell][ge]
                    expeCOMdist[regi][cell][clu]['expression'] += [expe]

    return expeCOMdist


def getCOMdistFocusAll(models, orderCell, interAll, modelsKeep):
    distances = {}
    for regi in models[models.keys()[0]]:
        focusPos = sum([interAll[regi][k].keys() 
                        for k in interAll[regi]], [])

        for cell in orderCell:
            distances[cell] = []
            mo = models[cell][regi]
            mods = load_structuralmodels(mo)
            # keep the amount of selected models
            mods.define_best_models(modelsKeep)

            for model in mods:
                COMselection = center_of_mass(model, selection=focusPos)
                COMfull = center_of_mass(model)
                dist = sqrt(square_distance_to2(COMselection, COMfull))
                distances[cell] += [dist]

    return distances
