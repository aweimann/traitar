#!/usr/bin/env python
#adapted from Nathan Salomonis: http://code.activestate.com/recipes/578175-hierarchical-clustering-heatmap-python/

import matplotlib as mpl
#pick non-x display
mpl.use('Agg')
import matplotlib.pyplot as pylab
import scipy
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as dist
import numpy
import string
import time
import sys, os
import getopt
import numpy as np
import pandas as ps
from .PhenotypeCollection import PhenotypeCollection
import warnings
warnings.filterwarnings("ignore", category=FutureWarning) 
#ignore these warnings
#/usr/lib/pymodules/python2.7/matplotlib/collections.py:548: FutureWarning: elementwise comparison failed; returning scalar instead, but in the future will perform elementwise comparison
#  if self._edgecolors == 'face':
#/usr/lib/pymodules/python2.7/matplotlib/backends/backend_pdf.py:2184: FutureWarning: comparison to `None` will result in an elementwise object comparison in the future.
#  different = bool(ours != theirs)

################# Perform the hierarchical clustering #################

def heatmap(x, row_header, column_header, primary_pt_models, color_f, row_method,
            column_method, row_metric, column_metric,
            filename, sample_f, secondary_pt_models):
    
    print "\nrunning hiearchical clustering using %s for columns and %s for rows" % (column_metric,row_metric)
        
    """
    This below code is based in large part on the protype methods:
    http://old.nabble.com/How-to-plot-heatmap-with-matplotlib--td32534593.html
    http://stackoverflow.com/questions/7664826/how-to-get-flat-clustering-corresponding-to-color-clusters-in-the-dendrogram-cre

    x is an m by n ndarray, m observations, n genes
    """
    
    ### Define the color gradient to use based on the provided name
    #if color_gradient == 'red_white_blue':
    #    cmap=pylab.cm.bwr
    #if color_gradient == 'red_black_sky':
    #    cmap=RedBlackSkyBlue()
    #if color_gradient == 'red_black_blue':
    #    cmap=RedBlackBlue()
    #if color_gradient == 'red_black_green':
    #    cmap=RedBlackGreen()
    #if color_gradient == 'yellow_black_blue':
    #    cmap=YellowBlackBlue()
    #if color_gradient == 'seismic':
    #    cmap=pylab.cm.seismic
    #if color_gradient == 'green_white_purple':
    #    cmap=pylab.cm.PiYG_r
    #if color_gradient == 'coolwarm':
    #    cmap=pylab.cm.coolwarm

    ### Scale the max and min colors so that 0 is white/black
    #vmin=x.min()
    #vmax=x.max()
    #vmax = max([vmax,abs(vmin)])
    #vmin = vmax*-1
    #norm = mpl.colors.Normalize(vmin/2, vmax/2) ### adjust the max and min to scale these colors

    ### Scale the Matplotlib window size
    default_window_hight = 10.5
    default_window_width = 10
    fig = pylab.figure(figsize=(default_window_width,default_window_hight)) ### could use m,n to scale here
    color_bar_w = 0.015 ### Sufficient size to show
    color_bar_w = 0.015 ### Sufficient size to show
        
    ## calculate positions for all elements
    # ax1, placement of dendrogram 1, on the left of the heatmap
    #if row_method != None: w1 = 
    [ax1_x, ax1_y, ax1_w, ax1_h] = [0.05,0.42,0.2,0.4]   ### The second value controls the position of the matrix relative to the bottom of the view
    width_between_ax1_axr = 0.004
    height_between_ax1_axc = 0.004 ### distance between the top color bar axis and the matrix
    
    # axr, placement of row side colorbar
    [axr_x, axr_y, axr_w, axr_h] = [0.31,0.1,color_bar_w,0.6] ### second to last controls the width of the side color bar - 0.015 when showing
    axr_x = ax1_x + ax1_w + width_between_ax1_axr
    axr_y = ax1_y; axr_h = ax1_h
    width_between_axr_axm = 0.004

    # axc, placement of column side colorbar
    [axc_x, axc_y, axc_w, axc_h] = [0.4,0.63,0.5,color_bar_w] ### last one controls the hight of the top color bar - 0.015 when showing
    axc_x = axr_x + axr_w + width_between_axr_axm
    axc_y = ax1_y + ax1_h + height_between_ax1_axc
    height_between_axc_ax2 = 0.004

    # axm, placement of heatmap for the data matrix
    [axm_x, axm_y, axm_w, axm_h] = [0.4,0.9,2.5,0.5]
    axm_x = axr_x + axr_w + width_between_axr_axm
    axm_y = ax1_y; axm_h = ax1_h
    axm_w = axc_w

    # ax2, placement of dendrogram 2, on the top of the heatmap
    [ax2_x, ax2_y, ax2_w, ax2_h] = [0.3,0.72,0.6,0.15] ### last one controls hight of the dendrogram
    ax2_x = axr_x + axr_w + width_between_axr_axm
    ax2_y = ax1_y + ax1_h + height_between_ax1_axc + axc_h + height_between_axc_ax2
    ax2_w = axc_w

    # placement of the phenotype legend
    [axpl_x, axpl_y, axpl_w, axpl_h] = [0.78,0.84,0.05,0.13]
    # placement of the sample legend

    # axcb - placement of the sample legend
    [axsl_x, axsl_y, axsl_w, axsl_h] = [0.05,0.29,0.05,0.09]

    # axcb - placement of the color legend
    [axcb_x, axcb_y, axcb_w, axcb_h] = [0.05, 0.88,0.05,0.09]


    # Compute and plot top dendrogram
    if not column_method is None and x.shape[1] > 1:
        start_time = time.time()
        d2 = dist.pdist(x.T)
        D2 = dist.squareform(d2)
        ax2 = fig.add_axes([ax2_x, ax2_y, ax2_w, ax2_h], frame_on=True)
        Y2 = sch.linkage(D2, method=column_method, metric=column_metric) ### array-clustering metric - 'average', 'single', 'centroid', 'complete'
        Z2 = sch.dendrogram(Y2)
        ind2 = sch.fcluster(Y2,0.7*max(Y2[:,2]),'distance') ### This is the default behavior of dendrogram
        time_diff = str(round(time.time()-start_time,1))
        ax2.set_xticks([]) ### Hides ticks
        ax2.set_yticks([])
        #print 'Column clustering completed in %s seconds' % time_diff
    else:
        ind2 = ['NA']*len(column_header) ### Used for exporting the flat cluster data
        
    # Compute and plot left dendrogram.
    if not row_method is None and x.shape[0] > 1:
        start_time = time.time()
        x_bin = x.copy()
        x_bin[x_bin > 0] = 1
        d1 = dist.pdist(x_bin)
        D1 = dist.squareform(d1)  # full matrix
        ax1 = fig.add_axes([ax1_x, ax1_y, ax1_w, ax1_h], frame_on=True) # frame_on may be False
        ax1.set_xticks([]) ### Hides ticks
        ax1.set_yticks([])
        Y1 = sch.linkage(D1, method=row_method, metric=row_metric) ### gene-clustering metric - 'average', 'single', 'centroid', 'complete'
        Z1 = sch.dendrogram(Y1, orientation='right')
        ind1 = sch.fcluster(Y1,0.7*max(Y1[:,2]),'distance') ### This is the default behavior of dendrogram
        time_diff = str(round(time.time()-start_time,1))
        #print 'Row clustering completed in %s seconds' % time_diff
    else:
        ind1 = ['NA']*len(row_header) ### Used for exporting the flat cluster data
        
    # Plot heatmap color legend
    n = len(x[0]); m = len(x)
    if secondary_pt_models is not None:
        cmaplist = np.array([[247,247,247],[166,206,227],[178,223,138],[31,120,180]])/256.0
    else:
        cmaplist = np.array([[247,247,247],[31,120,180]])/256.0
    cmap = mpl.colors.ListedColormap(cmaplist)
    axcb = fig.add_axes([axcb_x, axcb_y, axcb_w, axcb_h], frame_on=False)  # axes for colorbar
    #cb = mpl.colorbar.ColorbarBase(axcb, cmap=cmap, orientation='horizontal')
    bounds = numpy.linspace(0, len(cmaplist), len(cmaplist) + 1) 
    norm = mpl.colors.BoundaryNorm(bounds, len(cmaplist))
    cb = mpl.colorbar.ColorbarBase(axcb, cmap=cmap, norm=norm, spacing='proportional', ticks=bounds, boundaries=bounds)
    if secondary_pt_models is not None:
        axcb.set_yticklabels(["negative", "%s positive" % primary_pt_models.get_name(), "%s positive" % secondary_pt_models.get_name(), "both predictors positive"], fontsize = 8)
        axcb.yaxis.set_ticks([0.125, 0.375, 0.625, 0.875])
    else:
        axcb.set_yticklabels(["%s negative" % primary_pt_models.get_name(), "%s positive" % primary_pt_models.get_name()], fontsize = 8)
        axcb.yaxis.set_ticks([0.25, 0.75])
    axcb.set_title("Heatmap colorkey", fontsize = 10, loc = "left")
    
    # Plot distance matrix.
    axm = fig.add_axes([axm_x, axm_y, axm_w, axm_h])  # axes for the data matrix
    xt = x
 
    if not column_method is None and x.shape[1] > 1:
        idx2 = Z2['leaves'] ### apply the clustering for the array-dendrograms to the actual matrix data
        xt = xt[:,idx2]
        ind2 = ind2[idx2] ### reorder the flat cluster to match the order of the leaves the dendrogram
        pass
    if not row_method is None and x.shape[0] > 1 :
        idx1 = Z1['leaves'] ### apply the clustering for the gene-dendrograms to the actual matrix data
        xt = xt[idx1,:]   # xt is transformed x
        ind1 = ind1[idx1] ### reorder the flat cluster to match the order of the leaves the dendrogram
    ### taken from http://stackoverflow.com/questions/2982929/plotting-results-of-hierarchical-clustering-ontop-of-a-matrix-of-data-in-python/3011894#3011894
    im = axm.matshow(xt, aspect='auto', origin='lower', cmap=cmap, norm=norm) ### norm=norm added to scale coloring of expression with zero = white or black

    axm.set_xticks([]) ### Hides x-ticks
    axm.set_yticks([])

    # Add text
    new_row_header=[]
    new_column_header=[]
    for i in range(x.shape[0]):
        margin = 0
        if len(row_header) > 0 :
            fontdict = {'fontsize': 7}
        if len(row_header) > 30 :
            fontdict = {'fontsize': 7}
            margin = 0.5
        if len(row_header) > 50 :
            fontdict = {'fontsize': 4}
        if len(row_header) > 100 :
            fontdict = {'fontsize': 2}
        if len(row_header) > 200:
            fontdict = {'fontsize': 1}
        #if len(row_header)<100: ### Don't visualize gene associations when more than 100 rows
        axm.plot([-0.5, len(column_header)], [i - 0.5, i - 0.5], color = 'black', ls = '-')
        if x.shape[0] > 1 and row_method is not None:
            label = row_header[idx1[i]]
        else: 
            label = row_header[i]
        fontdict.items
        axm.text(x.shape[1] + 0.2, i - margin , '  ' + label, fontdict = fontdict)
        new_row_header.append(label)
            
    for i in range(x.shape[1]):
        if not column_method is None and x.shape[1] > 1:
            axm.plot([i-0.5, i-0.5], [-0.5, len(row_header) - 0.5], color = 'black', ls = '-')
            axm.text(i-0.5, -0.5, ' '+ column_header[idx2[i]], fontdict = {'fontsize': 7}, rotation=270, verticalalignment="top") # rotation could also be degrees
            new_column_header.append(column_header[idx2[i]])
        else: ### When not clustering columns
            axm.plot([i-0.5, i-0.5], [-0.5, len(row_header) - 0.5], color = 'black', ls = '-')
            axm.text(i-0.5, -0.8, ' '+column_header[i], fontdict = {'fontsize': 7}, rotation=270, verticalalignment="top")
            new_column_header.append(column_header[i])
    
    pt2acc = primary_pt_models.get_pt2acc()
    pt2acc.index = pt2acc.loc[:, "accession"]
    #colors
    colors = ps.read_csv(color_f, index_col = None, sep = "\t")
    if "category" in pt2acc.columns:
        #assign categories to colors
        import sets
        #get unique categories in the order they appear in the pt mapping table
        cats = sorted(set(pt2acc.loc[:, "category"].tolist()), key=lambda x: pt2acc.loc[:, "category"].tolist().index(x))
        if not colors.shape[0] < len(cats):
            # Plot phenotype legend
            axpl = fig.add_axes([axpl_x, axpl_y, axpl_w, axpl_h], frame_on=False)  # axes for colorbar
            #for i in pt2cat2col.index:
            #    if pt2cat2col.loc[i,"Category"] not in cat2col: 
            #        cat2col[pt2cat2col.loc[i,"Category"]] = pt2cat2col.loc[i, ["r", "g", "b"]]
            #        col2id[pt2cat2col.loc[i,"Category"]] = j 
            #        j += 1
            pt2cat = dict([(pt2acc.loc[i, "accession"], pt2acc.loc[i, "category"]) for i in pt2acc.index])
            cat2id = dict([(cats[i - 1], i) for i in range(1, len(cats) + 1)])
            cmaplist = ps.DataFrame(colors.iloc[:len(cats),])
            cmaplist.index = cats
            cmaplist = cmaplist / 256.0
            cmap_p = mpl.colors.ListedColormap(cmaplist.values)
            bounds = numpy.linspace(0, cmaplist.shape[0], cmaplist.shape[0] + 1) 
            norm = mpl.colors.BoundaryNorm(bounds, cmaplist.shape[0])
            cb = mpl.colorbar.ColorbarBase(axpl, cmap=cmap_p, norm=norm, spacing='proportional', ticks=bounds, boundaries=bounds)
            axpl.set_yticklabels([i for i in cats], fontsize = 6)
            axpl.yaxis.set_ticks(np.arange(1.0 / len(cats) / 2, 1,  1.0 / len(cats)))
            axpl.set_title("Phenotype colorkey", fontsize =  10, loc = "left")
            # Plot colside colors
            # axc --> axes for column side colorbar
            axc = fig.add_axes([axc_x, axc_y, axc_w, axc_h])  # axes for column side colorbar
            dc = numpy.array([cat2id[pt2cat[i]]  for i in column_header]).T
            if x.shape[1] > 1 and column_method is not None:
                dc = dc[idx2]
            dc.shape = (1, x.shape[1])
            im_c = axc.matshow(dc, aspect='auto', origin='lower', cmap=cmap_p)
            axc.set_xticks([]) ### Hides ticks
            axc.set_yticks([])
    
    # Plot rowside colors
    if sample_f is not None and x.shape[0] > 1:
        samples = ps.read_csv(sample_f, sep = "\t", index_col = "sample_name")
        if "category" in samples.columns:
            #get unique sample categories and sort according to the order they appear in the sampling file
            sample_cats = sorted(set(samples.loc[:, "category"].tolist()), key = lambda x: samples.loc[:, "category"].tolist().index(x))
            cat2col = dict([(sample_cats[i - 1], i) for i in range(1, len(sample_cats) + 1)])
            cmaplist = ps.DataFrame(colors.iloc[:len(sample_cats),]) / 256.0
            cmap_p = mpl.colors.ListedColormap(cmaplist.values)
            axr = fig.add_axes([axr_x, axr_y, axr_w, axr_h])  # axes for row side colorbar
            dr = numpy.array([cat2col[samples.loc[i, "category"]]  for i in row_header]).T
            if row_method is not None:
                dr = dr[idx1]
            dr.shape = (samples.shape[0], 1)
            #cmap_r = mpl.colors.ListedColormap(['r', 'g', 'b', 'y', 'w', 'k', 'm'])
            im_r = axr.matshow(dr, aspect='auto', origin='lower', cmap=cmap_p)
            axr.set_xticks([]) ### Hides ticks
            axr.set_yticks([])
            # Plot sample legend
            axsl = fig.add_axes([axsl_x, axsl_y, axsl_w, axsl_h], frame_on=False)  # axes for colorbar
            bounds = numpy.linspace(0, len(sample_cats), len(sample_cats) + 1) 
            norm = mpl.colors.BoundaryNorm(bounds, len(sample_cats))
            cb = mpl.colorbar.ColorbarBase(axsl, cmap=cmap_p, norm=norm, spacing='proportional', ticks=bounds, boundaries=bounds)
            axsl.yaxis.set_ticks(np.arange(1.0 / len(sample_cats) / 2, 1,  1.0 / len(sample_cats)))
            axsl.set_yticklabels([i for i in sample_cats], fontsize = 6)
            axsl.set_title("Sample colorkey", loc = "left", fontsize = 10)
    
    
    #exportFlatClusterData(filename, new_row_header,new_column_header,xt,ind1,ind2)

    ### Render the graphic
    if len(row_header)>50 or len(column_header)>50:
        pylab.rcParams['font.size'] = 6
    else:
        pylab.rcParams['font.size'] = 8

    pylab.savefig(filename)
    pylab.savefig(filename, dpi=300) #,dpi=200
    #pylab.show()

def getColorRange(x):
    """ Determines the range of colors, centered at zero, for normalizing cmap """
    vmax=x.max()
    vmin=x.min()
    if vmax<0 and vmin<0: direction = 'negative'
    elif vmax>0 and vmin>0: direction = 'positive'
    else: direction = 'both'
    if direction == 'both':
        vmax = max([vmax,abs(vmin)])
        vmin = -1*vmax
        return vmax,vmin
    else:
        return vmax,vmin
    
################# Export the flat cluster data #################

def exportFlatClusterData(filename, new_row_header,new_column_header,xt,ind1,ind2):
    """ Export the clustered results as a text file, only indicating the flat-clusters rather than the tree """
    
    filename = string.replace(filename,'.pdf','.txt')
    export_text = open(filename,'w')
    column_header = string.join(['UID','row_clusters-flat']+new_column_header,'\t')+'\n' ### format column-names for export
    export_text.write(column_header)
    column_clusters = string.join(['column_clusters-flat','']+ map(str, ind2),'\t')+'\n' ### format column-flat-clusters for export
    export_text.write(column_clusters)
    
    ### The clusters, dendrogram and flat clusters are drawn bottom-up, so we need to reverse the order to match
    new_row_header = new_row_header[::-1]
    xt = xt[::-1]
    
    ### Export each row in the clustered data matrix xt
    i=0
    for row in xt:
        export_text.write(string.join([new_row_header[i],str(ind1[i])]+map(str, row),'\t')+'\n')
        i+=1
    export_text.close()
    
    ### Export as CDT file
    filename = string.replace(filename,'.txt','.cdt')
    export_cdt = open(filename,'w')
    column_header = string.join(['UNIQID','NAME','GWEIGHT']+new_column_header,'\t')+'\n' ### format column-names for export
    export_cdt.write(column_header)
    eweight = string.join(['EWEIGHT','','']+ ['1']*len(new_column_header),'\t')+'\n' ### format column-flat-clusters for export
    export_cdt.write(eweight)
    
    ### Export each row in the clustered data matrix xt
    i=0
    for row in xt:
        export_cdt.write(string.join([new_row_header[i]]*2+['1']+map(str, row),'\t')+'\n')
        i+=1
    export_cdt.close()

################# Create Custom Color Gradients #################
#http://matplotlib.sourceforge.net/examples/pylab_examples/custom_cmap.html

def RedBlackSkyBlue():
    cdict = {'red':   ((0.0, 0.0, 0.0),
                       (0.5, 0.0, 0.1),
                       (1.0, 1.0, 1.0)),
    
             'green': ((0.0, 0.0, 0.9),
                       (0.5, 0.1, 0.0),
                       (1.0, 0.0, 0.0)),
    
             'blue':  ((0.0, 0.0, 1.0),
                       (0.5, 0.1, 0.0),
                       (1.0, 0.0, 0.0))
            }

    my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return my_cmap

def RedBlackBlue():
    cdict = {'red':   ((0.0, 0.0, 0.0),
                       (0.5, 0.0, 0.1),
                       (1.0, 1.0, 1.0)),

             'green': ((0.0, 0.0, 0.0),
                       (1.0, 0.0, 0.0)),
    
             'blue':  ((0.0, 0.0, 1.0),
                       (0.5, 0.1, 0.0),
                       (1.0, 0.0, 0.0))
            }

    my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return my_cmap

def RedBlackGreen():
    cdict = {'red':   ((0.0, 0.0, 0.0),
                       (0.5, 0.0, 0.1),
                       (1.0, 1.0, 1.0)),
    
             'blue': ((0.0, 0.0, 0.0),
                       (1.0, 0.0, 0.0)),
    
             'green':  ((0.0, 0.0, 1.0),
                       (0.5, 0.1, 0.0),
                       (1.0, 0.0, 0.0))
            }
    
    my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return my_cmap

def YellowBlackBlue():
    cdict = {'red':   ((0.0, 0.0, 0.0),
                       (0.5, 0.0, 0.1),
                       (1.0, 1.0, 1.0)),
    
             'green': ((0.0, 0.0, 0.8),
                       (0.5, 0.1, 0.0),
                       (1.0, 1.0, 1.0)),
    
             'blue':  ((0.0, 0.0, 1.0),
                       (0.5, 0.1, 0.0),
                       (1.0, 0.0, 0.0))
            }
    ### yellow is created by adding y = 1 to RedBlackSkyBlue green last tuple
    ### modulate between blue and cyan using the last y var in the first green tuple
    my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return my_cmap

  
