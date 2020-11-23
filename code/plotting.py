import matplotlib.patches as mpatches
import seaborn as sns
import pandas as pd
from matplotlib import pyplot as plt
import copy
import numpy as np
import matplotlib.patches as patch
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LinearSegmentedColormap
from scipy.spatial import distance
import fastcluster
from fastcluster import linkage
import matplotlib.backends.backend_pdf
from scipy.cluster.hierarchy import dendrogram
import scipy.stats as stats
import scipy.cluster.hierarchy as sch
from scipy.stats.mstats import zscore
from pytadbit.modelling.structuralmodels  import load_structuralmodels
from taddyn import Chromosome
from taddyn.modelling.impoptimizer  import IMPoptimizer
from math import sqrt
import os
import sys
sys.path.append('.')
import metrics
from scipy.stats import linregress
from scipy.stats import spearmanr


def optimPlot1(matPath, optimPath, nmodels = 100, nkeep = 100):
    res = int(matPath.split('_')[-1][:-2])


    ## get all final Optimisation files in matrix folder
    optFpath = []
    optFiles = os.listdir(optimPath)
    for opt in optFiles:
            if opt.startswith('ValOptimisation') and not opt.endswith('pdf'):
                    optFpath.append(optimPath + opt)

    ## create experiment object
    test_chr = Chromosome(name='Test',centromere_search=False,
                          species='Homo sapiens', assembly='na')#, max_tad_size=260000)
    # If interaction index and not matrix
    test_chr.add_experiment('test',exp_type='Hi-C', resolution=res,
                             norm_data=matPath)
    exp = test_chr.experiments[0]

    ## Filter data
    metrics.PCHiC_filtering(exp)

    ## Plot
    for opt in optFpath:
        print opt
        optim=IMPoptimizer(exp,start=1, end=exp.size, close_bins=1, n_models=nmodels, n_keep=nkeep)
        optim.load_from_file(opt)
        print optim.get_best_parameters_dict()
        optim.plot_2d(savefig=opt[:-3]+'pdf',show_best=10)#[0.2,0.8])#skip={'maxdist':2000}

        plt.show()


# Plot HiC matrix in a way we can store it in subplots
def plotHiC(matrix1, bad_color=None, axe=None, transform=np.log2, 
            rescale_zeros=True, title = None, **kwargs):
    """
    Plot HiC matrix

    :param matrix: list of lists with values to be plotted
    :param None bad_color: plots NaNs in a given color
    :param kwargs: extra parameters for the imshow function of matplotlib

    :returns: Nothing but pain and despair
    """

    matrix = copy.deepcopy(matrix1)

    if bad_color is not None:
        kwargs['cmap'] = plt.get_cmap(kwargs.get('cmap', None))
        kwargs['cmap'].set_bad(bad_color, 1.)

    if not isinstance(matrix, (np.ndarray, np.generic)):
        matrix = np.asarray(matrix)

    # remove zeroes from the matrix in order to avoid -inf with log transform
    if rescale_zeros:
        try:
            mini = min(matrix[np.nonzero(matrix)]) / 2.
        except ValueError:
            mini = 0.
        matrix[matrix==0] = mini

    matrix = np.ma.masked_where(np.isnan(matrix), transform(matrix))
    if axe == None:
        im = plt.imshow(matrix, interpolation='None', origin='lower', **kwargs)

        plt.xlim(0 - 0.5, len(matrix[0]) - 0.5)
        plt.ylim(-0.5, len(matrix) - 0.5)
        
        if title != None:
            plt.title(title)
    else:
        im = axe.imshow(matrix, interpolation='None', origin='lower', **kwargs)

        axe.set_xlim(0 - 0.5, len(matrix[0]) - 0.5)
        axe.set_ylim(-0.5, len(matrix) - 0.5)

        if title != None:
            axe.set_title(title)

    return im
    

def plotBoxplot(distancesTemp, saveFig=False, title=False, meanDistance=False,
               order=['pVSp', 'enhVSp', 'enhVSenh'], orderCell=False,
               colors2 = [(69/255.,185/255.,235/255.), (204/255.,204/255.,204/255.)], pdf=[]):
    ''' 
    Function to plot boxplots of different distances arranged by comparison and 
        cell
    '''
            
    colors = copy.copy(colors2)
    # colorsp = []

    # for i in range(6):
    #     colorsp.append('b')
    #     colorsp.append(colors[i])
    #     colorsp.append(colors[i])
    # sns.set_palette(colorsp, 18)
    cells = []
    enhs = []
    allValues = []
    listOfValues = {}
    if orderCell == False:
        sortKeys = sorted(distancesTemp.keys())
    else:
        sortKeys = orderCell
    for cell in sortKeys:
        for enh in sorted(distancesTemp[cell].keys()):
            if len(distancesTemp[cell][enh]) != 0:
                if meanDistance == True:
                    values = np.mean(distancesTemp[cell][enh]) 

                elif meanDistance == False:
                    values = distancesTemp[cell][enh]
                    # Add the values for the p-test list
                listOfValues[(cell, enh)] = values
                for va in values:
                    cells.append(cell)
                    enhs.append(enh)
                    allValues.append(va)
            else:
                cells.append(cell)
                #keypros.append(keypro)
                enhs.append(enh)
                allValues.append(None)

    print sortKeys
    pd_dat = pd.DataFrame(
            {'cell': cells,
            'Gene':enhs,
            'Distance (nm)': allValues})
    f = plt.figure(figsize=(10,10))
   
    ax = sns.boxplot(x='Gene', hue="cell", y="Distance (nm)", data=pd_dat, 
                   order=order, notch=False)
       
    region = distancesTemp.keys()[0].split('reg')[-1]

    if title == False:
        plt.title("Bin to bin distances distribution in region %s\n" \
                    %(region))
    else:
        plt.title(title)

    cellLabel = '\n'.join(sortKeys)
    plt.xticks(range(6), [cellLabel] * 6, rotation=90)
    ax.legend_.remove()
    # color boxes
    # first we need to know in which columns we have data
    enhancers = order
    box = 0
    col = -1
    # store separated data frames by cell
    cell_dic = {}
    for c in set(sortKeys):
        pd_temp = pd_dat.loc[pd_dat['cell'] == c]
        cell_dic[c] = pd_temp
    # create a blacklist to skip where we dont have data
    blacklist = []
    for enh_ in enhancers:
        data = True
        for c in set(sortKeys):
            # selected color
            col += 1
            # If we have data
            if len(cell_dic[c].loc[cell_dic[c]['Gene'] == enh_]['Distance (nm)']) > 1:
                mybox = ax.artists[box]
                # Change the appearance of that box
                mybox.set_facecolor(colors[col])
                # we have passed one box
                box += 1
                # data checker
                data = data * True
            else:
                val = float(cell_dic[c].loc[cell_dic[c]['Gene'] == enh_]['Distance (nm)'])
                # if its NaN
                if float(val) != float(val):
                    data = data * False
                    blacklist.append(col)
                else:
                    mybox = ax.artists[box]
                    # Change the appearance of that box
                    mybox.set_facecolor(colors[col])
                    # we have passed one box
                    box += 1
                    # data checker
                    data = data * True
   
    ## Modify legend
    enh_patch = []
    for no in range(len(order)):
        enh_patch += [mpatches.Patch(color=colors[no * len(orderCell)], label=order[no])]
    # enh1_patch = mpatches.Patch(color=colors[0], label=order[0])
    # enh2_patch = mpatches.Patch(color=colors[3], label=order[1])
    # enh3_patch = mpatches.Patch(color=colors[6], label=order[2])
    # enh4_patch = mpatches.Patch(color=colors[9], label=order[3])
    # enh5_patch = mpatches.Patch(color=colors[12], label=order[4])
    # enh6_patch = mpatches.Patch(color=colors[15], label=order[5])
    # leg = plt.legend(handles=[enh1_patch, enh2_patch, enh3_patch,
    #                          enh4_patch, enh5_patch, enh6_patch])

    leg = plt.legend(handles=enh_patch)
    
    # Get the bounding box of the original legend
    bb = leg.get_bbox_to_anchor().inverse_transformed(ax.transAxes)

    # Change to location of the legend. 
    xOffset = 0.34
    bb.x0 += xOffset
    bb.x1 += xOffset
    leg.set_bbox_to_anchor(bb, transform = ax.transAxes)


    if saveFig == True:
        pdf.savefig( f , bbox_inches='tight')


def plotMultiTrack(toPlotBox, toPlot2, regi, regiones, promoter,
                   promName, markers, minMax, xlab, ylab, idlab,
                   cellTypes, markerOrder=None, featureColors=None,
                    figsize=(20,7), plotType='linePlot', focusColor=None,
                  ylim=False, dcutoff=False,
                  palette='Reds', padlist = [False, False, False],
                  sizelist=[False, False, False],
                  boxplot=False):

    '''
    :param toPlotBox: a Dictionary with three keys, indicating cell, y_axis, and x_axis.
        It contains the distribution of distances found in the ensemble of models between
        the selected particles o interest and the focus particle to measure the distance 
        against.
    :param toplot2: Dictionary with cells as keys, containing as items another dictionary
        wth keys named as each of the axis and the associated data on it. I.E.:
        toPlot2 = {cell1: {y_axis:[value1, value2], x_axis:[value1, value2], ...}
    param markers: dictionaries with the markers positions that will be shown above
        the line plot. The needed structure of each of the inner dicitonaries:
        {regi:{cell1:{pos1:{content1}, pos2:{content2},...}, cell2:{...} ...}}
        Structure of markers:
        markers = {'innerDic name': innerDict}
    param None markerOrder: List with the orther in which the markers will be displaced
    param None featureColors: Dictionary indicating HEX color code for the markers
    param False dcutoff: Distance boundary from which we set the interaction
    param 'Reds' palette: color list for each of the cell lines (same order as cellTypes)
    param [False, False, False] padlist: List with a value in case want to twick the spacing
        of the upper tracks (sometimes they go inside the plot). You can try by starting to
        change the first one to -0.5
    param None focusColor: RGB color style for marking the focus point
        I.E.: [67/255.,147/255.,195/255.]
    :param False boxplot: wether to add points as boxplot or not

    Parameter combinations per length
        941 particles: padlist=[0.0001,0.0001,0.002], sizelist=['2%', '2%', '2%'] figsize(30,10)
    '''

    if palette == 'Reds':
        palette = sns.color_palette('Reds')

    if focusColor == None:
        focusColor = [67/255.,147/255.,195/255.]  # blue
        #focusColor = [215/255.,48/255.,39/255.]  # red

    fig,ax=plt.subplots(1, 1, figsize=figsize)
    cellTypes = copy.copy(cellTypes)

    # right now i use this to avoid the plot from moving due to label space
    #but its clear that there is a problem with the proporiton the base plot
    #takes in comparison with the others. Maybe the added axis percentaje
    #should be updated depending on plot size
    left = max(1/float(minMax[1] - minMax[0])**2 * 250, 0.05)
    left = min(left, 0.05)
    bottom = max(1/float(minMax[1] - minMax[0])**2 * 300, 0.05)
    bottom = min(bottom, 0.05)
    plt.subplots_adjust(left=left, right=1-left, top=1- bottom, bottom=bottom)

    ##
    ## First plot, lines connecting distances with no distribution
    ##
    # plot the tracks to connect all the points
    # pointplot trets data as categorical, so no floats in x
    for nk, k in enumerate(cellTypes):
        if boxplot == True:
            xs = [x + 0.5 for x in toPlot2[k][xlab]]
        else:
            xs = toPlot2[k][xlab]
        ax.plot(xs,
                toPlot2[k][ylab], 'o-', color=palette[nk])

    # Add interaction boundary if provided
    if dcutoff != False:
        ax.axhline(dcutoff, color='grey', alpha=0.3, linestyle='--')


    ##
    ## Second plot, distributions
    ##
    # we adjust positions of boxes and spheres in the middle of
    #the coordinates
    if plotType == 'boxplot':
        promoterPos = promoter
        toPlotBox[xlab] = [x + 0.5 for x in toPlotBox[xlab]]
        for nc, cell in enumerate(cellTypes):
            pos = pd.DataFrame(toPlotBox)[idlab] == cell
            df = pd.DataFrame(toPlotBox)[pos]
            sns.boxplot(x=xlab, y=ylab, data=df, ax=ax, notch=True,
                       showfliers=False, width=2, color=palette[nc],
                       linewidth=0.5)
    elif plotType == 'violinplot':
        promoterPos = promoter + 0.5
        toPlotBox[xlab] = [x + 0.5 for x in toPlotBox[xlab]]
        for nc, cell in enumerate(cellTypes):
            pos = pd.DataFrame(toPlotBox)[idlab] == cell
            df = pd.DataFrame(toPlotBox)[pos]
            sns.violinplot(x=xlab, y=ylab, data=df, ax=ax,
                           color=palette[nc],
                       linewidth=0.5, inner='quartile', cut=0, scale='width')
    else:
        promoterPos = promoter
        # plot pointplot with std
        sns.pointplot(data=pd.DataFrame(toPlotBox), x=xlab, y=ylab, hue=idlab,
                      ax=ax, estimator=np.mean, ci='sd', palette=palette,
                      #ax=ax, estimator=np.median, ci='sd', palette=palette,
                     hue_order=cellTypes)


    # get the number of divisions where labels would be equally spread
    start = minMax[0] - 10
    end = minMax[1] + 10
    #start = 0
    #end = ((regionEnd - regionStart) / resol) + 1

    resol = regiones[regi][3]
    regionStart = regiones[regi][1] +  (start * resol)
    regionEnd = regiones[regi][1] + (end * resol)

    plt.xlim(start, end)
    print regionStart, regionEnd
    #plt.xlim(minMax[0] - 10, minMax[1] + 10)

    #start, end = allAxis[-3].get_xlim()
    #print start, end
    binrange = range(int(start), int(end) + 1)

    posrange = range(regionStart, regionEnd + resol, resol)

    divisor = 10
    # use range to select the values in bin (value used at the time to plo) and
    #genomic position (value that we actually want to see)

    binx = [binrange[i] for i in range(0, len(binrange), len(binrange)/divisor)]
    posx = [posrange[i] for i in range(0, len(posrange), len(binrange)/divisor)]

    plt.xticks(binx,
               posx,
               rotation=45)


    if ylim != False:
        plt.ylim(ylim[0], ylim[1])

    # store ylims
    ylims = ax.get_ylim()




    # Add the patch for the promoter location
    rect2 = patch.Rectangle((promoterPos - 0.5,ylims[0]), 1, ylims[1] + abs(ylims[0]),
                            color=focusColor, alpha=0.3)
    ax.add_patch(rect2)


    ##
    ## Thirth plot, feature markers from the upper part
    ## Is quite hard to adjust it for different sizes
    ## Should work on a way to automatically adjust tosize and topad
    ##
    if len(markers) != 0:
        if markerOrder == None:
            markerOrder = sorted(markers.keys())
        if featureColors == None:
            featureColors = {}
            featureColList = ['#d73027', '#542788', '#1a9850', '#4393c3']
            for i in range(len(markerOrder)):
                # if we ran out of colors repeat
                remv = (i / len(featureColList)) * len(featureColList)
                featureColors[markerOrder[i]] = featureColList[i - remv]

        # Add patches for feature location
        # reverse list so the orther is the same as in the labels
        cellTypes.reverse()
        divider = make_axes_locatable(ax)
        for cell in cellTypes:
            features = {}
            for nma, ma in enumerate(markerOrder):
                features[nma] = markers[ma][regi][cell]

            for fe in sorted(features):
                # choose separation between plots

                ## This part is not automatised and deppends on how many features you have
                if fe == 0:
                    if cell != cellTypes[0]:
                        if padlist[0] != False:
                            topad = padlist[0]
                        else:
                            topad = 0.05
                        if sizelist[0] != False:
                            tosize = sizelist[0]
                        else:
                            tosize = "10%"

                        sizelist
                        ax2 = divider.append_axes("top", size=tosize, pad=topad)
                    else:
                        if padlist[1] != False:
                            topad = padlist[1]
                        else:
                            topad = 0
                        if sizelist[1] != False:
                            tosize = sizelist[1]
                        else:
                            tosize = "5%"
                        ax2 = divider.append_axes("top", size=tosize, pad=topad)
                else:
                    if padlist[2] != False:
                        topad = padlist[2]
                    else:
                        topad = -float(minMax[1] - minMax[0])**2 / 120000
                    if sizelist[2] != False:
                        tosize = sizelist[2]
                    else:
                        tosize = "10%"
                    ax2 = divider.append_axes("top", size=tosize, pad=topad)

                ##

                # select color depending on the comparison
                cmap = LinearSegmentedColormap.from_list(
                        'tempCmap', colors=['#f7f7f7',
                                            featureColors[markerOrder[fe]]], N=2)

                toplot = []
                if plotType != 'lineplot':
                    for b in binrange:
                        if b in features[fe]:
                            toplot += [1]
                        else:
                            toplot += [0]
                else:
                    # we skeew the data in here
                    #not sure why I have to reduce binrange by one, but works
                    for nb, b in enumerate(binrange[:-1]):
                        if 0 < nb < len(binrange) - 2:
                            if b in features[fe]:
                                toplot += [1]
                                toplot += [1]
                            else:
                                toplot += [0]
                                toplot += [0]
                        elif nb == 0:
                            if b in features[fe]:
                                toplot += [1]
                            else:
                                toplot += [0]
                        elif nb == len(binrange) - 2:
                            if b in features[fe]:
                                toplot += [1]
                                toplot += [1]
                                toplot += [0] # is ok, just to not add false peaks
                            else:
                                toplot += [0]
                                toplot += [0]
                                toplot += [0]

                ax2.imshow([toplot], interpolation='nearest', origin='lower',
                           extent=[0,len(binrange),0,1], cmap = cmap,
                          vmin=0, vmax=1)

                ax2.set_yticks([])
                ax2.set_xticks([])
                #ax2.get_yaxis().set_visible(False)
                ax2.get_xaxis().set_visible(False)

                if fe == 1:
                    ax2.set_ylabel(cell, rotation=0)
                    ax2.yaxis.set_label_coords(-0.05,0.3)

                #ax2.set_position([0.2, 0.2, 0.4, 0.4])
                #ax2.set_axis_bgcolor('none')
                for axis in ['top','bottom','left','right']:
                    ax2.spines[axis].set_linewidth(0)



    # set figure tittle
    title = 'Region: %s, promoter: %s' %(regi,
                                        promName)

    fig.suptitle(title, fontsize=20)
    fig.subplots_adjust(top=0.9)


    plt.show()
    return fig


def setLinealPlot(models, orderCell, regiones, regi, distancesBtw, selected,
                     selected2, colors, pdf=None, markerOrder=None, markers={},
                     plotType='lineplot', padlist=[False, False, False],
                     sizelist=[False, False, False], figsize=(20,5),
                     boxplot=False, featureColors=None):

        # get dcutoff of this region
        dcutoff = float(models[orderCell[0]][regi].split('C')[-1].split('L')[0])
        # We fo for each cluster of the selected ones
        for clu in distancesBtw[regi][orderCell[0]]:
            # Separate focus by promoters and enhancers
            for focus in distancesBtw[regi][orderCell[0]][clu]:
                if len(distancesBtw[regi][orderCell[0]][clu][focus]) != 0:
                    # If its a promoter will paint it in blue
                    if focus == 'promoter':
                        focusColor = [67/255.,147/255.,195/255.]  # blue
                    # Ohterwise red
                    else:
                        focusColor = [215/255.,48/255.,39/255.]  # red

                    ## Define label info
                    idlab = 'Cell'
                    chrom = regiones[regi][0]
                    ylab = '3D distance in the models cluster %s' %clu
                    xlab = 'Genomic coordinates (%s)' %chrom

                    # go per each focus element in interAll
                    for promoter, promName in selected:
                        # go per each cell line
                        toPlot = {idlab:[], xlab:[], ylab:[]}
                        toPlot2 = {}
                        toPlotBox = {idlab:[], xlab:[], ylab:[]}
                        ## Prepare the information to be ploted
                        for nk, cell in enumerate(orderCell):
                            toPlot2[cell] = {xlab:[], ylab:[]}
                            # get median positions to be ploted
                            others = sorted(distancesBtw[regi][cell][clu][focus][promoter].keys() + [(promoter, promoter)])
                            others = sorted(list(set([o[1] for o in others])))
                            minMax = min(others), max(others)
                            longi = (regiones[regi][2] - regiones[regi][1]) / regiones[regi][3]

                            for i in range(longi):
                                if i == promoter:
                                    val = 0
                                    toPlot2[cell][xlab].append(i)
                                    toPlot2[cell][ylab].append(val)
                                elif i in others:
                                    val = np.median(distancesBtw[regi][cell][clu][focus][promoter][(promoter, i)])
                                    toPlot2[cell][xlab].append(i)
                                    toPlot2[cell][ylab].append(val)
                                else:
                                    val = float('Nan')

                                toPlot[idlab].append(cell)
                                toPlot[xlab].append(i)
                                toPlot[ylab].append(val)


                            ##toPlotBox = {'Id':[], 'x':[], 'y':[]}
                            ## Here we define the bins that will have the std bar or a boxplot
                            ## basically the ones in selected2
                            for i in range(longi):
                                if i == promoter:
                                    val = [float('Nan')]
                                elif i in selected2:
                                    val = distancesBtw[regi][cell][clu][focus][promoter][(promoter, i)]
                                else:
                                    val = [float('Nan')]
                                for d in val:
                                    toPlotBox[idlab].append(cell)
                                    toPlotBox[xlab].append(i)
                                    toPlotBox[ylab].append(d)

                        ## Call ploting function
                        f = plotMultiTrack(toPlotBox, toPlot2, regi, regiones, promoter,
                                        promName, markers, minMax, xlab, ylab, idlab,
                                        cellTypes = orderCell, figsize=figsize, plotType=plotType,
                                        palette=colors, dcutoff=dcutoff, markerOrder=markerOrder,
                                        focusColor=focusColor, boxplot=boxplot, padlist=padlist,
                                        sizelist=sizelist, featureColors=featureColors)#, focus='enhancer')

                        if pdf != None:
                            pdf.savefig(f, bbox_inches='tight')


# Function to make the dendogram of dRMSD values between cell models
def modelDistanceDendo(outdata, regi, index, orderCell,
                      colorsList, outplot, saveFig,
                      select=u'dRMSD'):
    indexReg = index.loc[regi, ]

    # now we will take the column of interest to compute distance matrix
    #f = open("%s%s/model_distances.txt" %(outpath, regi))
    f = open("%s/model_distances_%s.txt" %(outdata, regi))

    header = f.next()
    f.close()
    header = [h.split('#')[-1] for h in header.split()]
    for ni, i in enumerate(header):
        if i == select:
            index_select = ni

    # convert to matrix
    matrix = [[0.0 for i in range(indexReg[orderCell[-1]] + 1)]
                      for i in range(indexReg[orderCell[-1]] + 1)]
    matrix = np.array(matrix)
    #with open("%s%s/model_distances.txt" %(outpath, regi)) as f:
    with open("%s/model_distances_%s.txt" %(outdata, regi)) as f:
        # skip header
        f.next()
        for line in f:
            lsplit = line.split()
            i = int(lsplit[0])
            j = int(lsplit[1])
            val = float(lsplit[index_select])
            matrix[i,j] = val

    # compute euclidean distance matrix
    distanceMatrix = matrix
    # make it simetric
    #distanceMatrix = np.maximum( distanceMatrix, distanceMatrix.transpose() )

    # Compute average linkage
    A_dist = distance.squareform(distanceMatrix)
    Z = linkage(A_dist,method="ward", metric='euclidean')

    D_leaf_colors = {}
    for nc, cell in enumerate(orderCell):
        if nc == 0:
            for r in range(indexReg[cell] + 1):
                D_leaf_colors[r] = colorsList[nc]   # green
        else:
            for r in range(indexReg[orderCell[nc - 1]] + 1, indexReg[cell] + 1):
                D_leaf_colors[r] = colorsList[nc]   # red

    dflt_col = "#808080"   # Unclustered gray
    link_cols = {}
    for i, i12 in enumerate(Z[:,:2].astype(int)):
        c1, c2 = (link_cols[x] if x > len(Z) else D_leaf_colors[x]
        for x in i12)
        link_cols[i+1+len(Z)] = c1 if c1 == c2 else dflt_col


    # Dendrogram
    if saveFig:
        pdf = matplotlib.backends.backend_pdf.PdfPages(outplot + 'clusteringTree_%s.pdf' %regi)

    fig, ax = plt.subplots(1,1, figsize=(15,15))
    D = dendrogram(Z=Z,color_threshold=None,
      leaf_font_size=12, leaf_rotation=45, link_color_func=lambda x: link_cols[x], ax=ax)
    plt.xticks([])
    plt.ylabel("Wards's sum of squares from dRMSD data")
    plt.show()

    if saveFig:
        pdf.savefig( fig , bbox_inches='tight')
        pdf.close()


## Functions for radialPlot ##
def getMatrixOrder(histdict, method='ward', metric='euclidean'):
    '''
    Function to set matrix order by similarity
    '''
    marks = sorted(histdict.keys())
    # get data in matrix format
    matrix = [[histdict[mark][k] for mark in marks ]for k in sorted(histdict[marks[0]].keys())]
    # normalize all columns by zscore values
    matrix = zscore(matrix, axis=0)
    #matrix = zip(*[[(v/float(sum(m))) if sum(m) else 0 for v in m] for m in zip(*matrix)])
    #matrix = np.array(matrix)
    
    # transpose matrix to get column clusters
    matrix2 = np.matrix.transpose(np.array(matrix))


    # Cluster it
    Y = sch.linkage(matrix2, method=method, metric=metric)
    Z2 = sch.dendrogram(Y, no_plot=True)
    markOrder = Z2['leaves']
    
    # return clustering order
    return markOrder

def squarePlot(histdict, region, title = '', prevOrder='', saveFig = True, minMax = "", zscore = False):
    # Generate some variables we will need
    marks = sorted(histdict.keys())
    nmark = len(marks)
    
    # get title
    if title == '':
        title = 'Marker intensities'
        
    # get data in matrix format
    matrix = [[histdict[mark][k] for mark in marks ]for k in sorted(histdict[marks[0]].keys())]
    if zscore == True:
        # normalize all columns by zscore values
        matrix = zscore(matrix, axis=0)
    #matrix = zip(*[[(v/float(sum(m))) if sum(m) else 0 for v in m] for m in zip(*matrix)])
    else:
        matrix = np.array(matrix)
    
    if prevOrder == '':
        markOrder = range(nmark)
        
    else:
        markOrder = prevOrder
    matrix = matrix[:,markOrder]

    # set colorbar limits
    mini = 100
    maxi = 0
    for i in matrix:
        for ii in i:
            mini = min(mini, ii)
            maxi = max(maxi, ii)
    colorLim = max(abs(mini), abs(maxi))
    colorLim = [-colorLim, colorLim]
    # check if there are given ones and are bigger
    if minMax != '':
        if (minMax[0] > colorLim[0]) or (minMax[1] < colorLim[1]):
            print 'There are values smaller/bigger than the provided range'
            print 'Default range will be used'
        
        # Just using the given values limit for asiggning color
        else:
            colorLim = minMax

    

    fig = plt.figure(figsize=(20, 10))
    plt.imshow(matrix, 
               cmap = "RdBu_r", aspect='auto', interpolation='none', origin='lower', 
               vmin=colorLim[0], vmax=colorLim[1])
    plt.xticks(range(nmark), [marks[i] for i in markOrder])
    plt.colorbar()

    
    # Add y labels tick
    ylabels = [int(i) for i in sorted(histdict[marks[0]].keys())]
    plt.yticks(range(len(ylabels)), ylabels)

    # Add labels
    plt.ylabel("Distance from center (nm)")
    plt.xlabel("ChIP mark Zscore")
    plt.title('%s to perifery in %s\n' %(title, region), size=20)
    
    
    #plt.colorbar(boundaries=np.linspace(-colorLim,colorLim,20))
    if saveFig == True:
        pdf.savefig( fig, bbox_inches='tight')
    plt.show()
    

    
# Newer version
def radialPlot(histdict, region, valRange, markOrder='', timesEach=7, nylabels=10, 
               title='', minMax = "", oneCheese=False, color='RdBu_r',
               unEqualRadi = False, colorTrend='bimodal', fixedMinMax=False,
               divisions=True):
    

    '''
    Plot radial matrix with given data from selected point of the model outwards it
    :param histdict: Dictionary with two leves. First the data to separate into plot portions,
        second the value asociated to each radius in the plot
    :param nmark: Number of keys we have in the first level of histdict
    :param bins: Number of distance chunks we have (number of keys in the second level of histdict)
    :param valRange: Value ranges analyzed from start point to end point
    :param fixedMinMax: Wether the function can look for a range that contains or values (False), or 
        the range given in minMax cannot be changed (True)
    :param True divisions: Wether to add or not translucent lines marking the radi of the
        spherical shells
    '''

    # Generate some variables we will need
    if oneCheese == False:
        marks = sorted(histdict.keys())
        nmark = len(marks)
        # Need to add one value more to bin so we see all data
        bins = len(histdict[marks[0]].keys()) + 1
    else:
        marks = ['']
        nmark = 1
        # Need to add one value more to bin so we see all data
        bins = len(histdict.keys()) + 1



    # get title
    if title == '':
        title = 'Marker intensities'
    # Generate some data...
    # Note that all of these are _2D_ arrays, so that we can use meshgrid
    # You'll need to "grid" your data to use pcolormesh if it's un-ordered points
    # theta es el angulo y r es la distancia en 'radio' a la que esta cada sub circunferencia
    portions = nmark * timesEach + 1
    theta, r = np.mgrid[0:2*np.pi:portions+0j, 0:1:bins+0j]
    z = np.random.random(theta.size).reshape(theta.shape)
    # Por algun motivo la ultima lista que se mete en pcolormesh no aparece, asi que hay q darle indices +1
    #aunque no existan
    # Lo mismo pasa cuando mira cuantos radios tiene cada quesito, siempre tendra en cuenta como que hay uno mas
    zs = []

    # If we are dealing with radius with different sizes
    if unEqualRadi == True:
        r = []
        for i in range(portions):
            r.append([v / valRange[-1] for v in valRange])
        r = np.array(r)

    # get new marks order
    if oneCheese == False:
        if markOrder != '':
            newMarks = [marks[i] for i in markOrder]
        else:
            newMarks = [marks[i] for i in range(nmark)]
        for n in range(nmark):
            values = [histdict[newMarks[n]][i] for i in sorted(histdict[newMarks[n]].keys())]
            for i in range(timesEach):
                zs.append(values)
    else:
        newMarks = ''
        values = [histdict[i] for i in sorted(histdict.keys())]
        for i in range(timesEach):
            zs.append(values)
    





    fig, ax2 = plt.subplots(ncols=1, subplot_kw=dict(projection='polar'), figsize=(10,10))
    #fig = plt.figure(figsize=(10,10))
    #plt.polar()
    
    infinites = []
    # Take color limit 
    if fixedMinMax == False:
        for n in range(0, (nmark * timesEach), timesEach):
            if sum(zs[n]) != 0:
                ## set colorbar limits
                # keep track of infinites
                for iz, z in enumerate(zs[n]):
                    if z ==  float('Inf') or z ==  float('-Inf'):
                        infinites.append([n, iz])
                # Find extrem values for color ranges
                mini = min([z for z in zs[n] if (z !=  float('Inf') and z !=  float('-Inf'))])
                maxi = max([z for z in zs[n] if (z !=  float('Inf') and z !=  float('-Inf'))])
                colorLim = max(abs(mini), abs(maxi))
                colorLim = [-colorLim, colorLim]
                # check if there are given ones and are bigger
                if minMax != '':
                    if (minMax[0] > colorLim[0]) or (minMax[1] < colorLim[1]):
                        print 'There are values smaller/bigger than the provided range'
                        print 'Default range will be used'
                        minMax = colorLim
                        #print colorLim
                    # Just using the given values limit for asiggning color
                    else:
                        colorLim = minMax
            else:
                colorLim = minMax
    else:
        colorLim = minMax
    if colorTrend == 'unimodal':
        colorLim[0] = 0     
    for n in range(0, (nmark * timesEach), timesEach):
        #if sum(zs[n]) != 0:
        ## Plot cheese
        plt.pcolormesh(theta[n:n+timesEach+1], r[n:n+timesEach+1], zs[n:n+timesEach], 
                       cmap=color, vmin=colorLim[0], vmax=colorLim[1], edgecolors='face')
                       #cmap=colors[n/timesEach], )
        ## Add vertical line separating portions
        plt.axvline(x=theta[n][0], color='black', alpha=0.3, linestyle='--')

    # If there is an infinite value we add an asterisk
    #print infinites
    #for infi in infinites:
    #    plt.plot((theta[infi[0]][0] + theta[infi[0] + timesEach][0])/2, 
    #             (r[0][infi[1]] + r[0][infi[1] + 1])/2, 
    #             '*', c = 'white', markersize = 20)
    #plt.plot(0.5178449428994164, 0.25, '*', c='white', markersize = 20)
    ## Add labels
    # get position for x labels
    angles = [i[0] for i in theta]
    # we should put the label more or less in the middle portion
    labpos = timesEach / 2
    angles = [ angles[n + labpos] for n in range(0, len(angles) - 1, timesEach)]
    # Add x labels  
    plt.xticks(angles, newMarks)


    #ax2.set_ylim([0, 1])

    ## Add y tick values
    # we remove the starting point from valRAnge
    #valRange = valRange[1:]
    # we need to create a range from 0 to one with same numbers as our bins
    if unEqualRadi == True:
        binrange = r[0][1:]
        # The plot will for sure have too many divisions
        fibonacci_numbers = [0, 1]
        for i in range(2,len(binrange) + 3):
            fibonacci_numbers.append(fibonacci_numbers[i-1]+fibonacci_numbers[i-2])
        # sacrilegious change
        fibonacci_numbers[3] = 0
        fibonacci_numbers[4] = 2
        binrangePos = [i for i in fibonacci_numbers[3:] if i <= len(binrange)]
        valRangeS = [int(valRange[1:][i]) if i in binrangePos else '' for i in range(0, len(valRange[1:]))]
        plt.yticks(binrange, valRangeS)
        alpha=0.3
    else:
        nytics = bins
        binrange = np.linspace(0,1,nytics,endpoint=True)[1:]
        # If the plot has to many divisions we need to show just a few axis
        steps = (bins / nylabels) + 1
        if bins > nylabels :
            binrangePos = [i for i in range(0, len(binrange), steps)]
            #binrangeS = [binrange[i] if i in binrangePos else '' for i in range(0, len(binrange))]
            valRangeS = [int(valRange[1:][i]) if i in binrangePos else '' for i in range(0, len(valRange[1:]))]
            plt.yticks(binrange, valRangeS)
        # If there are just a few divisions we just use the range values to delimitate y axis
        else:
            plt.yticks(binrange, valRange[1:])
        # transparency for radius lines
        alpha=0.3



    # Add LINES if we dont have to many divisions
    if len(valRange) < 100 and divisions == True:
        plt.grid(axis='y', linestyle='--', linewidth=1, alpha=alpha)

    # Add colorbar
    if colorTrend == 'bimodal':
        cbar = plt.colorbar(ticks=[colorLim[0], 0, colorLim[1]], orientation='vertical', fraction=0.040, pad=0.15)
    else:
        cbar = plt.colorbar(ticks=[colorLim[0], colorLim[1]], orientation='vertical', fraction=0.040, pad=0.15)
    #cbar.ax.set_yticklabels(['Low', 'Average', 'High'])  # horizontal colorbar

    
    # title etc
    plt.title('%s to perifery in %s\n\n' %(title, region), size=20)

    # move positions of y labels
    if oneCheese == True:
        ax2.set_rlabel_position(90)

    plt.show()

    return colorLim, infinites, fig
    
# Function to know if we are dealing with a NaN
# works because NaN isn't equal to anything, even itself
def isNaN(num):
    return num != num
    
    
# Function to load coverage files info and stats
def covLoading(covFiles, regiones, resol, discrete=False):
    '''
    param False discrete: False if you want actual values, or threshold if you want
        0 if smaller or 1 if greater or equal
    '''

    notPresent = set()
    # Load files coverage
    covDict = {}
    for regi in regiones:
        covDict[regi] = {}
        for cfi in covFiles:
            marker = cfi.split('/')[-1].split('_')[0]
            marker = marker.replace('.', '')
            covDict[regi][marker] = []
        
            with open(cfi, 'r') as f:
                for line in f:
                    line = line.rstrip().split('\t')
                    # in line[0] we would have region id
                    if line[0] in covDict[regi][marker]:
                        if discrete == False:
                            covDict[regi][marker].append(float(line[1]))

                        else:
                            if float(line[1]) < discrete:
                                covDict[regi][marker].append(0)
                            else:
                                covDict[regi][marker].append(1)
                    else:
                        # just in case we want too see regions we havent add to analysis
                        notPresent.add(line[0])
    marks = sorted(covDict[regi].keys())         
    nmark = len(marks)

    # Check if region length is ok
    for regi in covDict.keys():
        for m in marks:
            reg = regiones[regi]
            regionStart = reg[1]
            regionEnd = reg[2]
            longi = ((regionEnd - regionStart) / resol) + 1
            if len(covDict[regi][m]) != longi:
                difference = len(covDict[regi][m]) - longi
                print 'Region %s has %s more/less positions in file %s' \
                %(regi, difference, m)
                # If more in file, we remove from the end #### CHANGE WHEN CORRECT FILES ###
                exit()
    # get coverage total average and standar deviation
    statsCov = {}
    for regi in covDict.keys():
        statsCov[regi] = {}
        for m in marks:
            # get total average
            statsCov[regi][m] = {'tAvg': np.mean(covDict[regi][m]), 'tSD':np.std(covDict[regi][m])}

    return covDict, statsCov, marks, nmark

# Function to get distance from interest point 
def getDistancesFromPoint(mods, cluster, interest='center'):
    '''
    :param mods: Clustered ensemble of models in tadbit StructuralModels object
    :param cluster: Integer indicating the models from wich cluster will
        be measured.
    :param 'center' interest: Point from wich we measure distances. If
        not set model center of mass will be set as default. If an integer
        is probided, the distances will be measured from the bin this integer
        relates to
    :returns: cdistDict, a dictionary with model bins as keys (integers) and
                    the distance from this bin, in all the models from the 
                    ensemble belonging to the selected cluster, to the interest
                    bin
              cdistDictMean, a dictionary like cdistDict but just containing the 
                    mean distances
              mini, a float with the minimum distance detected
              maxi, a float with the maximum distance detected
              
    
    '''
    # Here we store distances from interest point 
    # get cluster ids
    models = [mods[str(m)]['index'] for m in mods.clusters[cluster]]
    models = [mods[mdl] for mdl in models]

    # variables to know values range
    mini = 100
    maxi = 0

    cdistDict = {}
    # create bin indexes
    for nbin in range(len(models[0]['x'])):
        cdistDict[nbin] = []
    for mo in models:
        if interest == 'center':
            # get center of mass
            mcen = metrics.center_of_mass(mo)
        else:
            # In this case interest should be a bin
            mcen = {'x': mo['x'][interest], 'y': mo['y'][interest], 'z': mo['z'][interest]}

        # go for each bin
        for nbin, x in enumerate(mo['x']):
            pos = {'x': x, 'y': mo['y'][nbin], 'z': mo['z'][nbin]}
            # get distance 
            dist = sqrt((mcen['x'] - pos['x'])**2 + (mcen['y'] - pos['y'])**2 + (mcen['z'] - pos['z'])**2)
            cdistDict[nbin].append(dist)
            # store range
            mini = min(mini, dist)
            maxi = max(maxi, dist)



    # get distance dict with mean
    cdistDictMean = {}
    for nbin in cdistDict.keys():
        cdistDictMean[nbin] = np.mean(cdistDict[nbin])

    return cdistDict, cdistDictMean, mini, maxi


# Funtion to get the distribution from interest point of all bins
# Funtion to get the distribution from interest point of all bins
def getMarkerDistribution(cdistDict, cdistDictMean, covDict, statsCov, marks, region, 
                          mini, modelRadius, resThres, maxRadi=False, groupBy = 'binNumber',
                         discrete=False,method = 'divVolume', pval=0.01):


    
    '''
    zscore bool True: wether to use zscores or raw data instead
    groupBy str 'binNumber': Options are 'binNumber' and 'density'. binNumber stands for
        a radius distribution of fixed distance whereas 'density' will do all spheres with 
        mantaining equal volumes
    para 'zscore' method: Which method to use at the time to normalize. Can choose between
        zscore, contingency, divVolume and percentage
    param False discrete: NOT False if you want to build a contingency table and obtain the 
        odds
    '''
    # Create list to tell when there is no data
    noData = set()
    # Get histogram input
    histdict = {}
    # dictionary for pseudo normalization
    histdict2 = {}
    histdictMean = {}
    #histdictContinuous = {}
    #histdictContinuousMean = {}

    # compute model radius
    maximumDistance = resThres * (int(modelRadius/resThres) + 1)
    # Need to create this one in tha Density case to avoid white areas at the end
    maximumDistanceDensi = round(modelRadius + 1)
    valRange = range(0, int(maximumDistance) + 1, int(resThres))
    # set the limit of our search diameter
    if maxRadi == False:
        maxRadi2 = maximumDistance
        maxRadi2Densi = maximumDistanceDensi
    else:
        maxRadi2 = maxRadi
        maxRadi2Densi = maxRadi

    # If groupBy != 'density' we set valRangeTrimmed, otherwise will ve overwriten
    valRangeTrimmed = [i for i in valRange if i <= maxRadi2]
    if groupBy == 'density':
        # First of all we modify maximumDistance so its real
        valRangeTrimmed = [0]
        # get list with all cut sites
        # Since distances are from center of bin, we need to ad a radius of a bin
        #modelDiameter += resThres/2 ### !!!!!!!!! NO LO VEO NECESARIO

        # Now we need to compute de volume of the first sphere (with r equal to 
        #bin radius)
        firstVolume = (4 * np.pi * (resThres)**3) / 3
        valRangeTrimmed.append(resThres)
        # With info of the first volume we can multiply it for each radius in the plot
        #to obatain the distance for the next radius
        n = 2
        nextRadi = round((((n * firstVolume) * 3) / (4 * np.pi))**(1./3))
        # ULTIMO CAMBIO DE < A <=
        while nextRadi <= maxRadi2Densi:
            valRangeTrimmed.append(nextRadi)
            # Find in next radius away
            n += 1
            nextRadi = round((((n * firstVolume) * 3) / (4 * np.pi))**(1./3))


        # create the valrange we use to store distance data
        # if we have not reached the end of the model because it lies in the middle
        #of two radii checked
        if (round(modelRadius) + 1) > max(valRangeTrimmed):
            valRange = list(sorted(set(valRangeTrimmed + [round(modelRadius) + 1])))
        else:
            valRange = valRangeTrimmed
        #print valRange
        # in each marker
        for k in covDict.keys():
            histdict[k] = {}
            histdict2[k] = {}
            histdictMean[k] = {}
            #histdictContinuous[k] = defaultdict(int)
            #histdictContinuousMean[k] = defaultdict(int)

            # we create the value ranges we have prepared
            for mrange in valRange[1:]:
                #mrange = np.mean([v, valRange[nv + 1]])
                histdict[k][mrange] = 0
                histdict2[k][mrange] = []
                histdictMean[k][mrange] = 0

            # for each bin in our model for all values
            for nbin in cdistDict.keys():
                for nb in cdistDict[nbin]:
                    # we check between which range positions should our value be located and add it
                    pos = [n+1 for n, va in enumerate(valRange) if va <= nb < valRange[n+1]][0]
                    mrange = valRange[pos]

                    histdict[k][mrange] += covDict[k][nbin]
                    histdict2[k][mrange].append(covDict[k][nbin])

                    # Add values to the continuous list (not binned)
                    # Is a defaultdict, so if value doesnt exist is the same, will add it
                    #histdictContinuous[k][nb] += covDict[k][nbin]

                # Same with the mean distance values for a Bin in all models
                pos = [n+1 for n, va in enumerate(valRange) if va <= cdistDictMean[nbin] < valRange[n+1]][0]
                mrange = valRange[pos]
                histdictMean[k][mrange] += covDict[k][nbin]

                # Add values to the continuous list (not binned)
                # Is a defaultdict, so if value doesnt exist is the same, will add it
                #pos = cdistDictMean[nbin]
                #histdictContinuousMean[k][pos] += covDict[k][nbin]

    # We adjust valRange and maxRadi so the last one points to the range value locating the limit
    #maxRadi2 = maxRadi2 - ((valRange[1] - valRange[0]) / 2)
    #print valRange
    else:
        # in each marker
        for k in covDict.keys():
            histdict[k] = {}
            histdict2[k] = {}
            histdictMean[k] = {}
            #histdictContinuous[k] = defaultdict(int)
            #histdictContinuousMean[k] = defaultdict(int)

            # we create the value ranges we have prepared
            for nv, v in enumerate(valRange[:-1]):
                #mrange = np.mean([v, valRange[nv + 1]])
                mrange = valRange[nv + 1]
                histdict[k][mrange] = 0
                histdict2[k][mrange] = []
                histdictMean[k][mrange] = 0

            # for each bin in our model for all values
            for nbin in cdistDict.keys():
                for nb in cdistDict[nbin]:
                    # we check between which range positions should our value be located and add it
                    # if we substract the ranging starting point to our value and divide it by the range
                    #length we get the position in valRange where our value is located
                    pos = int((nb - int(mini)) / (valRange[1] - valRange[0]))
                    mrange = valRange[pos + 1]
                    #mrange = np.mean([valRange[pos], valRange[pos + 1]])

                    histdict[k][mrange] += covDict[k][nbin]
                    histdict2[k][mrange].append(covDict[k][nbin])

                    # Add values to the continuous list (not binned)
                    # Is a defaultdict, so if value doesnt exist is the same, will add it
                    #histdictContinuous[k][nb] += covDict[k][nbin]

                # Same with the mean distance values for a Bin in all models
                pos = int((cdistDictMean[nbin] - int(mini)) / (valRange[1] - valRange[0]))
                #histdictMean[k][np.mean([valRange[pos], valRange[pos + 1]])] += covDict[k][nbin]
                mrange = valRange[pos + 1]
                histdictMean[k][mrange] += covDict[k][nbin]

                # Add values to the continuous list (not binned)
                # Is a defaultdict, so if value doesnt exist is the same, will add it
                #pos = cdistDictMean[nbin]
                #histdictContinuousMean[k][pos] += covDict[k][nbin]
        #return histdict2

    # Get values distribution 
    #listMark = []
    #for m in histdict2.keys():
    #    tempi = []
    #    for rang in histdict2[m].keys():
    #        # Always PI - TB
    #        tempi.append(histdict2[m][rang])  
    #    listMark += tempi

    # If there is a maximum distance set from the point of interest we remove the ranges above
    if maxRadi != False:
        for k in marks:
            for di in histdict2[k].keys():
                if di > maxRadi2:
                    del histdict2[k][di]
                    del histdict[k][di]
                    del histdictMean[k][di]

            #for di in histdictContinuous[k].keys():
            #    if di > maxRadi:
            #        del histdictContinuous[k][di]
            #        del histdictContinuousMean[k][di]

    
    ## Normalize data in histdict2
    # obtain numer of models been analysed
    nmodels = len(cdistDict[cdistDict.keys()[0]])
    # start iterating over particles
    for k in marks:
        # get number of particles with positive and negative mark
        # number of particles with positive signal
        regiPositive = np.nansum(covDict[k]) * nmodels
        # number of particles with negative signal
        regiNegative = (len(covDict[k]) - np.nansum(covDict[k])) * nmodels
        for piece in histdict2[k].keys():
            if method == 'zscore':
                #prev = histdict2[k][piece] 
                histdict2[k][piece] = ((np.mean([i for i in histdict2[k][piece]]) - statsCov[k]['tAvg'])\
                                       / statsCov[k]['tSD'])
            # if we want odds ratio from a contingency table
            if discrete != False:
                if method == 'contingency':
                    # Get positive presence and negative signal number in the spherical shell
                    positive = np.nansum(histdict2[k][piece])
                    negative = len(histdict2[k][piece]) - positive
                    # get positive and negatives out of the shell
                    restPos = regiPositive - positive
                    restNeg = regiNegative - negative
                    contingencyTab = [[positive,negative],[restPos,restNeg]]
                    # get odds
                    #print contingencyTab, k, piece
                    oddsratio, pvalue = stats.fisher_exact(contingencyTab)
                    # convert to logarithm if significant and not 0
                    #if k == 'NKX61':
                    #print contingencyTab
                    #print piece, oddsratio, pvalue
                    if pvalue <= pval and oddsratio != 0:
                        oddsratio = np.log(oddsratio)
                    else:
                        oddsratio=0
                    # assign log Odds ratio value
                    #if k == 'NKX61':
                    #    print oddsratio
                    histdict2[k][piece] = oddsratio
                elif method == 'percentage':
                    # obtain percentage of positives in our spherical shell
                    shellPositives = np.nansum([i for i in histdict2[k][piece]])
                    wholePositives = np.nansum(covDict[k]) * nmodels
                    if wholePositives != 0:
                        #print shellPositives, wholePositives
                        histdict2[k][piece] = round((shellPositives / float(wholePositives)) * 100)
                    else:
                        histdict2[k][piece] = 0
                        noData.add(k)
            elif method == 'divVolume':
                # divide by volume and multiply result by 1000 to get higher values
                histdict2[k][piece] = (np.nansum([i for i in histdict2[k][piece]]) / firstVolume) * 1000

            # if the result is a nan, means that there wasnt enough bins in this distance, so we change 
            # it for a 0
            if isNaN(histdict2[k][piece]) or histdict2[k][piece] == float('Inf'):
                #print k, piece, 'removed'
                histdict2[k][piece] = 0
    
    print 'YOY SHOULD CHANGE THE PART THAT COLLAPSES RADIUS SMALLER IN DIFFERENCE THAN 1NM'
    print 'AT LEAST SHOW A WARNING OR STOP IN THERE'
    if maxRadi == False:
        return histdict2, valRange, noData, (histdict, histdictMean) # ,histdictContinuous, histdictContinuousMean)
    else:
        return histdict2, list(sorted(set(valRangeTrimmed))), noData, (histdict, histdictMean)
    

def getRadialPlot(selected, signalValues, outplot, regi,  orderCell, models, modelsKeep,
                     allClusters, firstSphere, minMax, thres=2, maxRadi=600, 
                      saveFig=False, orderByCluster=False, cluster=1,
                      markOrder=None):
    '''
    :param False orderByCluster: Wether to order the marks by similarity in the first
        cell from orderCell and then maintain that order for the rest of the cells
    
    '''
    
    if not os.path.exists(outplot + regi):
        os.makedirs(outplot + regi)

    marks = sorted(signalValues[regi][orderCell[0]].keys())
    if markOrder != None:
        markOrder_ = []
        for mm in markOrder:
            for nm, m in enumerate(marks):
                if m == mm:
                    markOrder_ += [nm]
    else:
        markOrder_ = markOrder

    # For each position of interest
    for interest, interestName in selected:
        # for each model
        for n, cell in enumerate(orderCell):
            model = models[cell][regi]     
            # load cluster info in model
            flag = '%s_%s' %(cell, regi)
            model = load_structuralmodels(model)
            # keep the amount of selected models
            model.define_best_models(modelsKeep)
            model.clusters = allClusters[flag]

            covDictDiscrt = copy.copy(signalValues[regi][cell])

            # we get the distances from interest point to other bins in our model
            cdistDict, cdistDictMean, mini, modelRadius = getDistancesFromPoint(model, 
                                                                                cluster, interest)

            ## delete center bin data from positions bins to avoid using his coverage data
            #del cdistDict[interest]
            # remove the marker from central point setting it as NaN
            #for m1 in covDictDiscrt.keys():
            #    covDictDiscrt[m1][regi][interest] = float('NaN')
            # get marker distribution and bin distribution
            histdict2, valRange, noData, _ = getMarkerDistribution(cdistDict, cdistDictMean, covDictDiscrt, None, 
                                                            marks, regi, mini, modelRadius, firstSphere, 
                                                            maxRadi = maxRadi, groupBy='density', method='contingency',
                                                            discrete=thres)
            

            for say in noData:
                print 'There is a Mark with no data at all in region %s: ' %regi
                print say

            if saveFig == True:
                pdf = matplotlib.backends.backend_pdf.PdfPages(outplot +
                                                            regi+ '/radialPlots_%snm_contingen_%s_%s.pdf' %(firstSphere, 
                                                                                                                interestName, 
                                                                                                                cell))
            if n == 0:   
                # If we want to cluster mark order by similarity in the first cell from list
                if orderByCluster:
                    markOrder_ = getMatrixOrder(histdict2, method='ward', metric='euclidean')
            
                # Get marker columns order
                minMax2, infinites, fig = radialPlot(histdict2, flag, valRange, markOrder_, 
                    title = 'ChIP intensity distribution from %s' %interestName,
                    minMax = minMax, unEqualRadi = True, 
                    fixedMinMax=True, divisions=False)

                if saveFig:
                    pdf.savefig( fig, bbox_inches='tight')

            #print n, minMax


            # Radial plot with binned data for mark intensity distribution
            if n != 0:
                _, infinites, fig = radialPlot(histdict2, flag, valRange, markOrder_, 
                    title = 'ChIP intensity distribution from %s' %interestName,
                    minMax = minMax, unEqualRadi = True, 
                    fixedMinMax = True, divisions=False)

                if saveFig:
                    pdf.savefig( fig, bbox_inches='tight')

            print 'Positions of infinite values: '
            print infinites
            if saveFig:
                pdf.close()

## Expression boxplots ##
def expressionBoxplot(distancesAll, orderCell, colors, pdf=None):
    sns.set_palette(colors, 3)
    for genePos, geneName in distancesAll:
        for regi in distancesAll[(genePos, geneName)]:
            cells = []
            distances = []
            for cell in distancesAll[(genePos, geneName)][regi]:
                for dist in distancesAll[(genePos, geneName)][regi][cell]:
                    cells += [cell]
                    distances += [dist]
            dataframe = pd.DataFrame.from_dict({'cell':cells, 'Particle distances':distances})


            f = plt.figure(figsize=(10,10))
            sns.boxplot(x='cell', y='Particle distances', data = dataframe, order = orderCell)
            plt.title('Distance distribution from %s to all the bins with expression', fontsize=15)
            plt.xticks(fontsize=15)
            plt.yticks(fontsize=15)
            plt.xlabel('Cell', fontsize=20)
            plt.ylabel('Particle distance', fontsize=20)

            if pdf != None:
                pdf.savefig( f , bbox_inches='tight')

    if pdf != None:
        pdf.close()

## Co-occurrence matrix ##
def plotClusterPlot(mtr, experiments2, title, colorlab=False,
                   valRange=[0,1], cmap='viridis', clusterThres=-np.inf,
                   labelColoring=False, labelsSize=15,
                   titleSize=20):
    '''
    Function to generate the clustered heatmap of correlations
    :param mtr: list of lists with the correlation matrix
    :param experiments2: list with the labels (in same order) of
        the correlation matrix
    :param False colorlab: dictionary with keys equal to labels
        in experiments2 and HEX color codes for the labels
    
    '''
    
    
    if colorlab == False:
        colorlab = {}
        for e in experiments2:
            colorlab[e] = '#252525'  # black
            
    row_colors = [colorlab[e] for e in experiments2] 
    


    # Fist obtain the future order of the dendogram
    linkage = fastcluster.linkage_vector(mtr, method='ward', 
                                 metric='euclidean')


    dendo = sch.dendrogram(linkage, no_plot=True,
                                        color_threshold=clusterThres)
    dendoOrd = dendo['leaves']
    # reorder correlation matrix for annotation
    anotCoor = [[mtr[d][dd] for dd in dendoOrd] for d in dendoOrd]


    # Draw the heatmap with the mask and correct aspect ratio
    clumap = sns.clustermap(mtr, method='ward', cmap=cmap, vmin=valRange[0], vmax=valRange[1], #annot=np.asanyarray(anotCoor),
                  linewidths=.5, xticklabels=experiments2, yticklabels=[' ' for e in experiments2],
                           row_colors=row_colors, col_colors=row_colors)

    fig = clumap.fig
    fig.suptitle(title, size=titleSize)

    # get the new order for the labeels
    newOrd = clumap.dendrogram_row.reordered_ind
    # check it was as obtained before
    if dendoOrd != newOrd:
        asdasdasd
    newExps = [experiments2[n] for n in newOrd]


    # change label color in y axis
    clumap.ax_heatmap.yaxis.axes.set_yticklabels(newExps, #fontweight='bold',
                                                size=labelsSize, va='center', rotation=0)
    if labelColoring == True:
        for tick_label in clumap.ax_heatmap.axes.get_yticklabels():
            tick_text = tick_label.get_text()
            tick_label.set_color(colorlab[tick_text])


    # change label color in x axis
    clumap.ax_heatmap.xaxis.axes.set_xticklabels(newExps, #fontweight='bold',
                                                size=labelsSize)
    if labelColoring == True:
        for tick_label in clumap.ax_heatmap.axes.get_xticklabels():
            tick_text = tick_label.get_text()
            tick_label.set_color(colorlab[tick_text])
    

    plt.show()
    
    return fig


def getColorLabel(expPresent1, cell, binedCoor1, expressionLimit, maxVal=False,
                 cmap='viridis', negCol=(1., 1., 1.)):
    colorlab = {}
    #expressionLimit

    for bin1 in binedCoor1:
        # take the most expressed value from the bin
        value = max(expPresent1[cell][gene] for gene in binedCoor1[bin1])
        geneName = ','.join(binedCoor1[bin1])# + ' %s' %bin1
        if value < expressionLimit:
            value = 0
        colorlab[geneName] = value

    nonZeros = [c for c in colorlab.values() if c != 0]
    if maxVal == False:
        maxVal = round(max(nonZeros), 1)

    # prepare 10 ranges
    steps = (maxVal - expressionLimit) / 10.
    colorRange = np.arange(expressionLimit, maxVal + steps, steps)

    # create palette
    viridisPalette = sns.color_palette(cmap, 10)

    # assigs colors
    colorlab2 = {}
    prev = 0
    for nc, c in enumerate(colorRange[1:]):
        for gene in colorlab:

            if prev <= colorlab[gene] < c:
                colorlab2[gene] = viridisPalette[nc]
        prev = c
        
    # reasign negative values
    for gene in colorlab:
        if colorlab[gene] <= 0:
            colorlab2[gene] = negCol
    
    # last value from range
    for gene in colorlab:
        if c <= colorlab[gene]:
            colorlab2[gene] = viridisPalette[nc]
    return colorlab2


def getCoOccurrencePlot(coOcurrMatrices, expressionLimit, maxVal, cluster,
                       newSignalPos, posConverter, newSignal, saveFig,
                       outplot=None):

    for regi in coOcurrMatrices:
        if saveFig == True:
            pdf = matplotlib.backends.backend_pdf.PdfPages(outplot + 
                            'CooccurrenceM_minExp%s_maxEsp%s_cluster%s_khindex.pdf' %(expressionLimit, 
                                                                                   maxVal,
                                                                                  cluster))

        for cell in coOcurrMatrices[regi]:
            # first get labels
            labels = [','.join(newSignalPos[regi][cell][posConverter[regi][cell][i]]) 
                      for i in sorted(posConverter[regi][cell])]

            matrix = coOcurrMatrices[regi][cell]


            colorlab2 = getColorLabel(newSignal[regi], cell, newSignalPos[regi][cell], 
                                      expressionLimit, maxVal=maxVal, cmap='Reds', 
                                      negCol = (189/255.,189/255.,189/255.))


            title = 'Matrix of expressed genes co-clustering in models in %s_%s' %(regi, cell)
            fig = plotClusterPlot(matrix, labels, title, colorlab=colorlab2,
                               valRange=[0,100], cmap='viridis', labelsSize=12)

            if saveFig == True:
                pdf.savefig( fig , bbox_inches='tight')
        if saveFig == True:
            pdf.close()

## Plots for the clustering of COM of clusters
def plotClusterPlot(mtr, experiments2, title, colorlab=False,
                   valRange=[0,1], cmap='viridis', clusterThres=-np.inf,
                   labelColoring=False, labelsSize=15,
                   titleSize=20):
    '''
    Function to generate the clustered heatmap of correlations
    :param mtr: list of lists with the correlation matrix
    :param experiments2: list with the labels (in same order) of
        the correlation matrix
    :param False colorlab: dictionary with keys equal to labels
        in experiments2 and HEX color codes for the labels
    
    '''
    
    
    if colorlab == False:
        colorlab = {}
        for e in experiments2:
            colorlab[e] = '#252525'  # black
            
    row_colors = [colorlab[e] for e in experiments2] 
    


    # Fist obtain the future order of the dendogram
    linkage = fastcluster.linkage_vector(mtr, method='ward', 
                                 metric='euclidean')


    dendo = sch.dendrogram(linkage, no_plot=True,
                                        color_threshold=clusterThres)
    dendoOrd = dendo['leaves']
    # reorder correlation matrix for annotation
    anotCoor = [[mtr[d][dd] for dd in dendoOrd] for d in dendoOrd]


    # Draw the heatmap with the mask and correct aspect ratio
    clumap = sns.clustermap(mtr, method='ward', cmap=cmap, vmin=valRange[0], vmax=valRange[1], #annot=np.asanyarray(anotCoor),
                  linewidths=.5, xticklabels=experiments2, yticklabels=[' ' for e in experiments2],
                           row_colors=row_colors, col_colors=row_colors)

    fig = clumap.fig
    fig.suptitle(title, size=titleSize)

    # get the new order for the labeels
    newOrd = clumap.dendrogram_row.reordered_ind
    # check it was as obtained before
    if dendoOrd != newOrd:
        ERROR
    newExps = [experiments2[n] for n in newOrd]


    # change label color in y axis
    clumap.ax_heatmap.yaxis.axes.set_yticklabels(newExps, #fontweight='bold',
                                                size=labelsSize, va='center', rotation=0)
    if labelColoring == True:
        for tick_label in clumap.ax_heatmap.axes.get_yticklabels():
            tick_text = tick_label.get_text()
            tick_label.set_color(colorlab[tick_text])


    # change label color in x axis
    clumap.ax_heatmap.xaxis.axes.set_xticklabels(newExps, #fontweight='bold',
                                                size=labelsSize)
    if labelColoring == True:
        for tick_label in clumap.ax_heatmap.axes.get_xticklabels():
            tick_text = tick_label.get_text()
            tick_label.set_color(colorlab[tick_text])
    

    plt.show()
    
    return fig


def plotCOM(distancesCluster, clustersBinPos, outplot, rankedSum,
           modelCheckCluster, saveFig=False):
    for regi in distancesCluster:
        if saveFig == True:
            pdf = matplotlib.backends.backend_pdf.PdfPages(outplot + 
                            'ClusterCOMdistance_averageClusterExp_cluster%s.pdf' %(
                                                                                  modelCheckCluster))

        for cell in distancesCluster[regi]:
            # get matrix
            nclust = len(clustersBinPos[regi][cell].keys())
            matrix = [[0 for i in range(nclust)] for i in range(nclust)]
            for c1, c2 in distancesCluster[regi][cell]:
                matrix[c1-1][c2-1] = np.median(distancesCluster[regi][cell][(c1, c2)]) * -1
                matrix[c2-1][c1-1] =matrix[c1-1][c2-1]

            # first get labels
            labels = range(1, nclust+1)

            Palette = sns.color_palette('Reds', nclust)
            colorlab2 = {}
            for k in rankedSum[regi][cell]:
                colorlab2[k] = Palette[rankedSum[regi][cell][k] - 1]

            title = 'Matrix of COM distance between expression clusters in %s_%s' %(regi, cell)
            fig = plotClusterPlot(matrix, labels, title, colorlab=colorlab2,
                               valRange=[0, np.amin(matrix)], cmap='viridis',
                                 labelsSize=30)

            if saveFig == True:
                pdf.savefig( fig , bbox_inches='tight')
        if saveFig == True:
            pdf.close()


def plotExpressionAndCOMd(outplot, expeCOMdist, orderCell, colors,
                          saveFig=False):
    if saveFig == True:
        pdf = matplotlib.backends.backend_pdf.PdfPages(outplot + 
                        'expressionVSdistanceToClustersCOM.pdf')

    statistics = {}
    for regi in expeCOMdist:
        statistics[regi] = {}
        markers = {1:'o', 2:'v', 3:'s', 4:'X'}
        distances_all = []
        expes_all = []
        for nc, cell in enumerate(orderCell):
            statistics[regi][cell] = {}
            distances = []
            expes = []
            print '#' * 3, cell, '#' * 3
            fig, ax = plt.subplots(figsize=(8, 5))
            for clu in expeCOMdist[regi][cell]:
                dist = expeCOMdist[regi][cell][clu]['distanceToCOM']
                expe = expeCOMdist[regi][cell][clu]['expression']

                distances += dist
                expes += expe

                distances_all += dist
                expes_all += expe

                plt.plot(dist, expe, markers[clu], color=colors[nc])


            sns.regplot(distances, expes, '.', color=colors[nc],
                       scatter_kws={'s':1})
            plt.title("%s, %s" %(regi, cell))
            plt.ylabel('Gene expression (log(FPKM))', size=18)
            plt.xlabel('Distance to the center of mass of the cluster', size=18)
            ax.tick_params(axis='both', which='major', labelsize=15)
            ax.tick_params(axis='both', which='minor', labelsize=15)
            #plt.ylim(0, 13)
            plt.show()
            statistics[regi][cell]['lineregress'] = linregress(distances, expes)
            statistics[regi][cell]['spearmanr'] = spearmanr(distances, expes)


            if saveFig == True:
                pdf.savefig( fig , bbox_inches='tight')

        fig = plt.figure(figsize=(8, 5))
        sns.regplot(distances_all, expes_all, color='blue', 
                    )
        plt.title("All cells")
        plt.ylabel('Gene expression (log(FPKM))')
        plt.xlabel('Distance to the center of mass of the cluster')
        #plt.ylim(0, 13)
        plt.show()
        statistics[regi]['allCell'] = {'lineregress': linregress(distances_all, expes_all)}
        statistics[regi]['allCell']['spearmanr'] = spearmanr(distances_all, expes_all)

        if saveFig == True:
            pdf.savefig( fig , bbox_inches='tight')


    if saveFig == True:
        pdf.close()
        
    return statistics


def plot11(outplot, orderCell, distances, colors, saveFig=False):

    if saveFig == True:
        pdf = matplotlib.backends.backend_pdf.PdfPages(outplot + 'distanceFromGlobToCOM.pdf')

    x = []
    y = []
    for cell in orderCell:
        y += distances[cell]
        x += [cell] * len(distances[cell])

    pd_dat = pd.DataFrame(
                {'COM distances':y,
                'Cell': x})

    fig = plt.figure()
    ax = sns.boxplot(x='Cell', y='COM distances', data=pd_dat, order=orderCell,palette=colors)
    #plt.plot(x, y, 'o', color='red')
    plt.title('Distance between globin genes COM and model COM')
    #plt.ylim(145, 250)
    #plt.axhline(0)
    plt.ylabel('Distance (nm)', size=15)
    plt.xlabel('Cell', size=15)
    plt.xticks(rotation=0, size=13)
    plt.yticks(size=13)

    if saveFig == True:
        pdf.savefig( fig , bbox_inches='tight')
    plt.show()

    if saveFig == True:
        pdf.close()