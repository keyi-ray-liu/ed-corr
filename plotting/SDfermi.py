import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import sys
import os

from TNplottingutils import dataloader
from TNplottingutils import random_utils
from TNplottingutils import searchkey
from SDplotter.SDutils import labels_determine


def compare_fermi(groups, verbose =False, filetitle = ''):
    
    groupnames = groups.keys()
    fontsize = 10
    fig, axes = plt.subplots(2, 1, figsize=(15, 7))
    ax_s : plt.Axes = axes[0]
    ax_d : plt.Axes = axes[1]
    

    for i, groupname in enumerate(groupnames):


        group = groups[groupname]

        files = list(filter(lambda x: 'mixed' in x and 'biasAinit1.0' not in x, group['files']))
        _, _,  _, _, occs, *vals  = dataloader([], files)

        
        for j, file in enumerate(files):

            occ = occs[j]
            LR = np.loadtxt(file + '/LR')


            try:
                res = int(searchkey('Ls', file))
            except:
                res = 64

            arr =int(np.sqrt(occ.shape[-1] - res * 2))

            left = np.argwhere( LR>0).flatten()
            right = np.argwhere( LR<0).flatten()

            leftind = np.where(left >= res, left + arr ** 2, left)
            rightind = np.where( right>= res, right + arr ** 2, right)

            print(leftind)

            #marker, _, color, label = labels_determine(file, -1, arr=None, plottag=r'$\mu_{{init}} = $' +searchkey('biasAinit', file), verbose=verbose, colordetermine='init')

            marker, _, color, label = labels_determine(file, -1, arr=None, plottag=r'$v_{sd} = $' +searchkey('scoupling', file), verbose=verbose, colordetermine='coupling')

            plotax = ax_s 
            plotfont = fontsize
                
            plotax.plot(occ[0][leftind], label=label, marker=marker, color=color)

            plotax.set_ylabel('Max $S_{vN}$ on bipartile cut', fontsize = plotfont)
            plotax.set_xlabel('time step', fontsize = plotfont)
            plotax.set_title('Source levels', fontsize = plotfont)

            plotax = ax_d
            plotfont = fontsize
                
            plotax.plot(occ[0][rightind], label=label, marker=marker, color=color)

            plotax.set_ylabel('Max $S_{vN}$ on bipartile cut', fontsize = plotfont)
            plotax.set_xlabel('time step', fontsize = plotfont)
            plotax.set_title('Drain levels', fontsize = plotfont)

    
    random_utils.no_duplicate_label(ax_s, fontsize, ncol=1)
    random_utils.no_duplicate_label(ax_d, fontsize, ncol=1)
    fig.tight_layout()
    fig.subplots_adjust(wspace=0, hspace=0)
    fig.savefig('plots/{}.pdf'.format(filetitle))


