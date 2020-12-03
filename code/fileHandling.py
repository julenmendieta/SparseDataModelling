import os
import copy
from collections import defaultdict

def listdirs(path):
    return [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]


def bamPaths(basePath, ending='bam'):
    bamfilesPath = '%sbamfiles/' %basePath
    bamfiles = {}
    cells = sorted(os.listdir(bamfilesPath))
    for cell in cells:
        bamfiles[cell] = []
        cellpath = bamfilesPath + '%s/' %cell
        bamfiles_ = os.listdir(cellpath)
        for ma in bamfiles_:
            if ma.endswith(ending):
                bamfiles[cell] += [cellpath + ma]     
        # check number of files
        if len(bamfiles[cell]) == 0:
            print('No bamfiles in %s' %(cell))
            del bamfiles[cell]
        elif len(bamfiles[cell]) > 1:
            print('There should only be a bam file in %s folder' %(cell))
            print('FOLLOWING STEPS WILL CRASH')
        else:
            bamfiles[cell] = bamfiles[cell][0]

    return bamfiles


def getMatricesPaths(basePath, starting='Matrix'):
    matricesPath = '%smatrices/' %basePath
    matrices = {}
    cells = sorted(os.listdir(matricesPath))
    regionsAll = set()
    for cell in cells:
        matrices[cell] = {} 
        cellpath = matricesPath + '%s/' %cell
        regions = listdirs(cellpath)
        for regi in regions:
            regionsAll.add(regi)
            matrices[cell][regi] = []
            regipath = cellpath + '%s/' %regi
            matrices_ = os.listdir(regipath)
            for ma in matrices_:
                if ma.startswith(starting):
                    matrices[cell][regi] += [regipath + ma]     
            # check number of files
            if len(matrices[cell][regi]) == 0:
                print('No matrices in %s, %s' %(cell, regi))
                del matrices[cell][regi]
            elif len(matrices[cell][regi]) > 1:
                print('There should only be a top ensemble in %s, %s folder' %(cell, regi))
                print('FOLLOWING STEPS WILL CRASH')
            else:
                matrices[cell][regi] = matrices[cell][regi][0]

    # get length of matrices
    matricesLength = {}
    for regi in regions:
        ma = matrices[cell][regi]
        _, start, end = ma.split('/')[-1].split('_')[3].split('-')
        resol = int(ma.split('/')[-1].split('_')[-1][:-2])
        matricesLength[regi] = ((int(end) - int(start))/ resol) + 1

    return matricesLength, sorted(list(regionsAll)), matrices


def getModelsPaths(basePath, ending='models'):
    modelsPath = '%smodels/' %basePath
    models = {}
    cells = sorted(os.listdir(modelsPath))
    regionsAll = set()
    for cell in cells:
        models[cell] = {} 
        cellpath = modelsPath + '%s/' %cell
        regions = listdirs(cellpath)
        for regi in regions:
            regionsAll.add(regi)
            models[cell][regi] = []
            regipath = cellpath + '%s/' %regi
            models_ = os.listdir(regipath)
            for mo in models_:
                if mo.endswith(ending):
                    models[cell][regi] += [regipath + mo]     
            # check number of files
            if len(models[cell][regi]) == 0:
                print('No models in %s, %s' %(cell, regi))
                del models[cell][regi]
            elif len(models[cell][regi]) > 1:
                print('There should only be a top ensemble in %s, %s folder' %(cell, regi))
                print('FOLLOWING STEPS WILL CRASH')
            else:
                models[cell][regi] = models[cell][regi][0]
    return cells, sorted(list(regionsAll)), models


def getREgiInfo(basePath, regi, cell):
    cellMatrixPath = basePath + 'matrices/%s/%s/' %(cell, regi)
    matrix = [m for m in os.listdir(cellMatrixPath) if m.startswith('Matrix')]
    if len(matrix) != 1:
        print('There should only be one matrix in this folder:')
        print(cellMatrixPath)
        return 'More matrix files than supposed in folder'

    else:
        matrix = matrix[0]
        resol = int(matrix.split('_')[-1][:-2])
        chromosome, start, end = matrix.split('_')[-2].split('-')
        return [chromosome, int(start), int(end), resol]


def getEnhPos(enhanFile, chrom, regionStart, regionStart2, regionEnd2,
                 resol):
    Enhinter = {}
    with open(enhanFile, 'r') as f:
        for line in f:
            if line.split(None, 1)[0] == '%s' %chrom:
                line = line.rstrip().split()
                # in enhancer we take midpoint between start and end for the bining
                midPoint = ((int(line[1]) + int(line[2])) / 2)
                if regionStart2 <= midPoint <= regionEnd2:
                    # turn to bin
                    midPointBin = (midPoint - regionStart) / resol
                    # Start storing data
                    if midPointBin not in Enhinter: 
                        Enhinter[midPointBin] = {}
                        Enhinter[midPointBin]['name'] = line[3]

                    else:

                        Enhinter[midPointBin]['name'] = '%s,%s' %(Enhinter[midPointBin]['name'], line[3])
    return Enhinter

def getPromPos(promFile, chrom, regionStart, regionStart2, regionEnd2,
                  resol):
    PPinter = {}
    with open(promFile, 'r') as f:
        for line in f:
            if line.split(None, 1)[0] == '%s' %chrom:
                line = line.rstrip().split()
                # in promoters we take starting point of the region for bining
                startPoint = int(line[1])
                if regionStart2 <= startPoint <= regionEnd2:
                    # turn to bin
                    startPointBin = (startPoint - regionStart) / resol
                    # Start storing data
                    genes = line[3].split(';')
                    if startPointBin not in PPinter: 
                        PPinter[startPointBin] = {}
                        PPinter[startPointBin]['name'] = ','.join(genes)
                    else:
                        genes2 = PPinter[startPointBin]['name'].split(',')
                        genes2 = ','.join(list(set(genes + genes2)))
                        PPinter[startPointBin]['name'] = genes2
    return PPinter

def getInterestPos(interestFile, promAll, enhAll, regiones):
    interestOne = []
    with open(interestFile, 'r') as f:
            for line in f:
                interestOne.append(line.split()[3])

    interAll = {}
    for regi in regiones:
        interAll[regi] = {'enhancer':{}, 'promoter':{}}
        for inte in interestOne:
            for k in promAll[regi]:
                if ',' in promAll[regi][k]['name']:
                    names = promAll[regi][k]['name'].split(',')
                    for na in names:
                        if inte == na:
                            interAll[regi]['promoter'][k] = {'name':promAll[regi][k]['name']}                                       
                else:
                    if inte == promAll[regi][k]['name']:
                        interAll[regi]['promoter'][k] = {'name':promAll[regi][k]['name']}

        for inte in interestOne:
            for k in enhAll[regi]:
                if ',' in enhAll[regi][k]['name']:
                    names = enhAll[regi][k]['name'].split(',')
                    for na in names:
                        if inte == na:
                            interAll[regi]['enhancer'][k] = {'name':enhAll[regi][k]['name']}                                       
                else:
                    if inte == enhAll[regi][k]['name']:
                        interAll[regi]['enhancer'][k] = {'name':enhAll[regi][k]['name']}

    return interAll

def getElementCoordinates(regionsAll, regiones, enhanFile,
                          promFile, interestFile,
                          originalRegionPos = None):

    promAll = {}
    enhAll = {}
    #geneAll = {}
    for region in regionsAll:
        reg = regiones[region]
        # get regions chromosome, start and end position
        chrom = reg[0]
        regionStart = reg[1]
        regionEnd = reg[2]
        # get resolution
        resol = reg[3]

        # In case we gave other coordinates to trim the plots coordiantes
        if originalRegionPos != None:
            fromOri = (originalRegionPos[region][1]/resol) - (reg[1] / resol)
            fromEnd = (reg[2] / resol) - (originalRegionPos[region][2] / resol)

            regionStart2 = reg[1] + (fromOri * resol)
            regionEnd2 = reg[2] - (fromEnd * resol)
        else:
            regionStart2 = regionStart
            regionEnd2 = regionEnd

        ## get region enhancers locations
        Enhinter = getEnhPos(enhanFile, chrom, regionStart, regionStart2, 
                            regionEnd2, resol)

        # get region promoters locations
        PPinter = getPromPos(promFile, chrom, regionStart, regionStart2, 
                             regionEnd2, resol)

        # Note overlap between enhancers and promoters
        for en in Enhinter.keys():
            if en in PPinter:
                print('Enhancer %s was in same bin as a promoter (%s). \
    Bin position %s' %(Enhinter[en]['name'],
                    PPinter[en]['name'], en))
                #del Enhinter[en]

        promAll[region] = PPinter
        enhAll[region] = Enhinter

    # save the elements in which we are interested in
    interAll = getInterestPos(interestFile, promAll, enhAll, regiones)
    
    return enhAll, promAll, interAll


# Function to get markers for linealPlot
def createMarkers(enhAll, promAll, regiones, orderCell,
                 signalData, threshold=0, getFullExp=False):

    ## We add enhancers and promoters in a marker list
    markers = {}
    if len(enhAll) != 0:
        for nameVar, var in [('Enhancers', enhAll), ('Promoters', promAll)]:
            markers[nameVar] = {}
            for regi in regiones:
                markers[nameVar][regi] = {}
                for c in orderCell:
                    markers[nameVar][regi][c] = var[regi]

    ## Then we get a list of elements, coordinates and signal from signalData
    signalValues = {}
    for regi in sorted(regiones):
        signalValues[regi] = {}
    # load signal data
    with open(signalData, 'r') as f:
        header = f.next()
        # first we check dimensions have sense (three columns of info + one per cell)
        header = header.split()
        if len(header) == (len(orderCell) + 3):
            # get the position of each cell
            cellPos = {}
            for cell in orderCell:
                for nh, h in enumerate(header):
                    if h == cell:
                        cellPos[cell] = nh
            for line in f:
                line = line.split()
                element = line[0]
                chrom_ = line[1]
                pos = int(line[2])
                # place expressed genes in each regions
                for regi in regiones:
                    chrom = regiones[regi][0]
                    regionStart = regiones[regi][1]
                    regionEnd = regiones[regi][2]
                    resol = regiones[regi][3]
                    if chrom_ == chrom:
                        if regionStart <= pos < regionEnd + resol:
                            signalValues[regi][element] = {}
                            signalValues[regi][element]['position'] = pos
                            for cell in orderCell:
                                signalValues[regi][element][cell] = float(line[cellPos[cell]])
        else:
            print('File number of columns dont match with expected dimensions')
            print('Expected: three columns of info plus one column per cell')

    ## Get a marker for elements with signal above thresshold
    ## In this case expressed genes
    markers['Signal'] = {}
    for regi in regiones:
        regionStart = regiones[regi][1]
        resol = regiones[regi][3]
        markers['Signal'][regi] = {}
        for cell in orderCell:
            markers['Signal'][regi][cell] = set()
        for gene in signalValues[regi]:
            geneBin = (signalValues[regi][gene]['position'] - regionStart)/resol
            for cell in orderCell:
                if signalValues[regi][gene][cell] > thresshold:
                    markers['Signal'][regi][cell].add(geneBin)

    if getFullExp:
        return markers, signalValues
    else:
        return markers


def readf2(line, signal=1, thresshold=None):
    return line[0], int(line[1]), signal


def readf3(line, thresshold=None):
    signal = float(line[2])
    if thresshold != None:
        signal = 0 if signal < thresshold else signal
    return line[0], int(line[1]), signal


def loadSignalFile(signalDataFolder, regiones, orderCell,
                  thresshold, binariseSignal=True):
    ## get list of files with signal of interest
    signalFiles = os.listdir(signalDataFolder)
    ## Then we get a list of elements, coordinates and signal from signalData
    signalValues = {}
    for regi in sorted(regiones):
        # Gather region info
        chrom = regiones[regi][0]
        regionStart = regiones[regi][1]
        regionEnd = regiones[regi][2]
        resol = regiones[regi][3]
        # get list with same length as bins in region
        binsList = [0 for i in range(regionStart, regionEnd + resol, resol)]

        signalValues[regi] = {}
        for cell in orderCell:
            signalValues[regi][cell] = {}
            # get different signal info
            for fi in signalFiles:
                if fi.startswith('%s_' %cell):
                    # get signal type name
                    markName = fi.split('_')[1].split('.')[0]
                    signalValues[regi][cell][markName] = copy.copy(binsList)
                    with open(signalDataFolder + fi, 'r') as fi2:
                        line = fi2.next().split()
                        if len(line) == 2:
                            readf = readf2
                            thresshold2 = None
                        elif len(line) == 3:
                            readf = readf3
                            thresshold2 = thresshold
                        else:
                            print('%s has %s columns, should be 2 or 3' %(fi,
                                                                          len(line)))
                            ERROR
                        # store first line
                        line = readf(line, thresshold = thresshold2)
                        if line[0] == chrom:
                            if regionStart <= line[1] < regionEnd + resol:
                                # Assign value to binned possition
                                peakPosBin = (line[1] - regionStart) / resol
                                signalValues[regi][cell][markName][peakPosBin] += line[2]

                        # next lines of file
                        for line in fi2:
                            line = readf(line.split(), thresshold = thresshold2)
                            if line[0] == chrom:
                                if regionStart <= line[1] < regionEnd + resol:
                                    # Assign value to binned possition
                                    peakPosBin = (line[1] - regionStart) / resol
                                    signalValues[regi][cell][markName][peakPosBin] += line[2]

    # If we want to make contingency tables comparing ares with and without
    # signal irrespective of the level of it
    if binariseSignal == True:
        for regi in signalValues:
            for cell in signalValues[regi]:
                for markName in signalValues[regi][cell]:
                    for pos, val in enumerate(signalValues[regi][cell][markName]):
                        if val != 0:
                            signalValues[regi][cell][markName][pos] = 1
    return signalValues


def getSignalAndPos(regiones, orderCell, signalData, signalThreshold=0):
    newSignal = {}
    newSignalPos = {}
    for regi in regiones:
        newSignal[regi] = {}
        newSignalPos[regi] = {}
        for cell in orderCell:
            newSignal[regi][cell] = {}
            newSignalPos[regi][cell] = defaultdict(list)

    with open(signalData, 'r') as f:
        header = f.next().split()
        for line in f:
            line = line.split()
            element = line[0]
            chrom_ = line[1]
            pos = int(line[2])
            # place expressed genes in each regions

            for regi in regiones:
                # region coordinates
                chrom = regiones[regi][0]
                regionStart = regiones[regi][1]
                regionEnd = regiones[regi][2]
                resol = regiones[regi][3]
                # Check if inside region
                if chrom_ == chrom:
                    if regionStart <= pos < regionEnd + resol:
                        # check if each cell has signal
                        posBin = (pos - regionStart) / resol
                        for ncell in range(3, len(line)):
                            if float(line[ncell]) > signalThresshold:
                                newSignalPos[regi][header[ncell]][posBin] += [element]
                                newSignal[regi][header[ncell]][element] = float(line[ncell])
    return newSignalPos, newSignal
                        
