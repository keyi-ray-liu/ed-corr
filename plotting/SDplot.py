import numpy as np
import matplotlib.pyplot as plt
import sys
import matplotlib.cm as colormap
import json
from matplotlib import rc
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os
import glob

sys.path.append(os.path.abspath("/Users/knl20/Desktop/Code/TN"))
from TNplottingutils import dataloader
from TNplottingutils import random_utils
from TNplottingutils import searchkey
from SDutils import labels_determine




def plotone(occ, time, fig, ax : plt.Axes, paras):

    rescm = colormap.hot
    arrcm = colormap.cool
    total = occ.shape[1]
    L = paras['L']
    paras['t'] = -1.0
    timestep = paras['timestep']

    for i in range(occ.shape[1]):

        if i < L :
            cm = rescm
            linestyle = 'solid'
            selftotal = L * 2
            j = i

        elif i >= total - L:
            cm = rescm
            linestyle = 'solid'
            selftotal = L * 2    
            j = i - 9   
        
        else:
            cm = arrcm
            linestyle = 'dashed'
            selftotal = total - L * 2
            j = i - L

        ax.plot( time, occ[:, i], label='site' + str(i), c= cm(j/selftotal), linestyle=linestyle)

    

    # if abs(paras['scoupling']) <= 0.01:
    #     ax.set_ylim(0, 5e-3)

    # elif abs(paras['scoupling']) <= 0.1:
    #     ax.set_ylim(0, 0.1)

    #ax.set_ylim(0, 0.005)
    ax.set_xlim(0, 5)
    ax.legend(loc='center right', bbox_to_anchor=(1.1, 0.6))
    ax.set_title(' t : {}, s coupling : {:.5f}, d coupling : {:.5f},  res size = {}, timestep = {}'.format( paras["t"], paras['scoupling'],  paras['dcoupling'], L, timestep))




    

def quick():
    fs = sys.argv[1:]
    fig , axes = plt.subplots(len(fs), 1, figsize=(20, 10 * len(fs)))
    
    if len(fs) == 1:
        axes = [axes] 

    for i, f in enumerate(fs):
        try:
            occ = np.loadtxt( f + '/occs')

        except OSError as e:

            print(type(e))
            occup = np.loadtxt(f + '/occup')
            occdn = np.loadtxt(f + '/occdn')
            occ = occup+ occdn

        with open( f + '/sdpara.json', 'r') as io:
            paras = json.load( io)

        time = np.loadtxt(f + '/times')
        ax : plt.Axes = axes[i]
        plotone(occ, time, fig, ax, paras)
        
    fig.savefig('plots/{}.pdf'.format('_'.join(fs)))




def EE(axes, fig : plt.Figure, files):

    SvN, _, _, seffs, *vals = dataloader(['EE', 'effE'], files)

    for i, file in enumerate(files):
        
        S : np.ndarray = SvN[i]
        L = (S.shape[-1] - 9)//2 

        ranges = {
            'S' : [0, L] ,
            'D' : [L + 9, S.shape[-1] + 1],
            'Arr1' : [L, L+3],
            'Arr2' : [L+3, L+6],
            'Arr3' : [L+6, L+9],
        }

        for key in ['S', 'D']:
            astr = 'EE' + file + key
            ax : plt.Axes = axes[astr]
            crange = ranges[key]
            S_cur = S[ :, crange[0] : crange[1]]
            im = ax.imshow( S_cur.transpose(), origin='lower', rasterized=False)
            ax.set_title(key + file, fontsize=20)
            ax.set_xlabel('Time')
            ax.set_ylabel('Site number')
            fig.colorbar(im, orientation="horizontal")

        for key in ['Arr1', 'Arr2', 'Arr3']:

            astr = 'EE' + file + key
            ax : plt.Axes = axes[astr]
            crange = ranges[key]
            S_cur = S[ :, crange[0] : crange[1]]

            for i  in range(S_cur.shape[-1]):
                ax.plot(S_cur[:, i], label=str(crange[0] + i + 1 - L), rasterized=False)
            
            ax.set_title(key, fontsize=20)
            ax.set_xlabel('Time')
            ax.set_ylabel('SvN')
            ax.legend()

        seff = seffs[i]
        ax : plt.Axes = axes[ 'EE' + file + 'Seff']
        ax.plot( seff[i])
        ax.set_xlabel('Time')
        ax.set_ylabel('Effective Entropy')


        axes['EE' + file + "dum"].axis('off')

    return fig





def compare_SvN(groups, verbose =False, filetitle = ''):
    
    groupnames = groups.keys()
    fontsize = 10
    fig, ax = plt.subplots(figsize=(8, 5))
    ax : plt.Axes

    inset_dim=['50%', '40%']
    ax_inset : plt.Axes =   inset_axes(ax, 
                                        width=inset_dim[0], height=inset_dim[1],loc='lower right', 
                                bbox_to_anchor=[0, 0, 1, 1], bbox_transform=ax.transAxes
                                )            

    


    for i, groupname in enumerate(groupnames):


        group = groups[groupname]

        files = list(filter(lambda x: 'mixed' in x and 'biasAinit1.0' not in x, group['files']))


        SvNs, times,  _, _, occs, *vals  = dataloader([], files)

        
        for j, file in enumerate(files):

            time = times[j]
            #current = currents[j]
            #CC = CCs[j]
            SvN = SvNs[j]
            occ = occs[j]


            print(searchkey('biasAinit', file))

            try:
                res = int(searchkey('Ls', file))
            except:
                res = 64

            arr =int(np.sqrt(occ.shape[-1] - res * 2))
            marker, _, color, label = labels_determine(file, -1, arr=None, plottag=r'$\mu_{{init}} = $' +searchkey('biasAinit', file), verbose=verbose, colordetermine='init')

            plotax = ax 
            plotfont = fontsize 
            
            for i, svn in enumerate(SvN[:20]):
                #svn = svn - np.diag(np.diag(svn))

                #offset = np.zeros(svn.shape)

                # offset[res:res + arr ** 2, res:res+ arr**2] = svn[res:res + arr **2, res:res+ arr**2] 

                #svn = svn - offset

                svn[res:res + arr**2] = 0
                t = time[i]
                maxvals = np.sort(np.abs(svn), axis=None)[-1:]
                
                plotax.scatter([t for _ in range(maxvals.shape[0])], maxvals, label=label, marker=marker, color=color)

                plotax.set_ylabel('Max $S_{vN}$ on bipartile cut', fontsize = plotfont)
                plotax.set_xlabel('time step', fontsize = plotfont)
                plotax.set_title('Entanglement increase', fontsize = plotfont)

            plotax = ax_inset
            plotfont = fontsize / 2

            overlap = np.loadtxt(file  + '/productstateoverlap')

            plotax.set_ylabel('Overlap with product state', fontsize = plotfont)
            plotax.set_xlabel('time step', fontsize = plotfont)

            plotax.scatter(range(overlap.shape[0]), overlap, label=label, marker=marker, color=color)
    
    random_utils.no_duplicate_label(ax, fontsize, ncol=1)
    fig.tight_layout()
    fig.savefig('plots/{}.pdf'.format(filetitle))


def work(files=sys.argv[1:]):


    stack = lambda cols, *f : [ ''.join([*f]) + col for col in cols]
    #mosaic_template = [ stack(cols, 'EE', file) for file in files for cols in (["S", "Arr1", "D", "Seff"], ["S", "Arr2", "D", "Seff"], ["S", "Arr3", "D", "Seff"], ["dum", "dum", "dum", "dum"]) ]
    mosaic_template = [ ["S", "Arr1", "D", "Seff"], ["S", "Arr2", "D", "Seff"], ["S", "Arr3", "D", "Seff"]]

    width_ratios = [1, 0.5, 1, 1]
    fontsize = 20

    fig = plt.figure(
                    layout="constrained", 
                     figsize = ( 5*len(mosaic_template[0]), 5 * len(files))
                     )
    subfigs = fig.subfigures(nrows=len(files), ncols=1)

    axes = {}

    SvN, _, _, seffs, *vals = dataloader(['EE', 'effE'], files)


    for i, file in enumerate(files):
        
        S : np.ndarray = SvN[i]
        L = (S.shape[-1] - 9)//2 

        subfig : plt.Figure = subfigs[i]
        mosaic = [ [  'EE' + file + m for m in cols] for cols in mosaic_template ]
        axes |= subfig.subplot_mosaic(mosaic,  width_ratios=width_ratios)
        subfig.suptitle(file.split('/')[-1], fontsize=fontsize * 1.5)

        ranges = {
            'S' : [0, L] ,
            'D' : [L + 9, S.shape[-1] + 1],
            'Arr1' : [L, L+3],
            'Arr2' : [L+3, L+6],
            'Arr3' : [L+6, L+9],
        }

        for key in ['S', 'D']:
            astr = 'EE' + file + key
            ax : plt.Axes = axes[astr]
            crange = ranges[key]
            S_cur = S[ :, crange[0] : crange[1]]
            im = ax.imshow( S_cur.transpose(), origin='lower', rasterized=False)
            ax.set_title(key + ' entanglement', fontsize=20)
            ax.set_xlabel('Time')
            ax.set_ylabel('Site number')
            fig.colorbar(im, orientation="horizontal")

        for key in ['Arr1', 'Arr2', 'Arr3']:

            astr = 'EE' + file + key
            ax : plt.Axes = axes[astr]
            crange = ranges[key]
            S_cur = S[ :, crange[0] : crange[1]]

            for j  in range(S_cur.shape[-1]):
                ax.plot(S_cur[:, j], label=str(crange[0] + j + 1 - L), rasterized=False)
            
            ax.set_title(key, fontsize=20)
            ax.set_xlabel('Time')
            ax.set_ylabel('SvN')
            ax.legend()

        seff = seffs[i]
        ax : plt.Axes = axes[ 'EE' + file + 'Seff']
        ax.plot( seff)
        ax.set_xlabel('Time')
        ax.set_ylabel('Effective Entropy')


        #axes['EE' + file + "dum"].axis('off')

    # fig, axes = plt.subplot_mosaic(mosaic, constrained_layout=True,
    # gridspec_kw={'width_ratios':[1, 0.5, 1, 1],
    # 'height_ratios':[val for _ in range(len(files)) for val in [1, 1, 1, 0.5] ]}, figsize=( 10 * len(mosaic[0]), 1.5 * len(mosaic)))

    
    fig.savefig('plots/EE{}.pdf'.format( '_'.join([f.split('/')[-1] for f in files])))


