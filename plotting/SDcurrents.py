import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from TNplottingutils.dataloader import dataloader
from TNplottingutils import random_utils
from TNplottingutils.searchkey import searchkey
from SDutils import labels_determine

def occ_inset(axes_inset, groupname, insetoption, res, arr, file, plottag, verbose, time, cur_occ):

    ax_inset : plt.Axes = axes_inset[ 'occ' + groupname + 'inset']

    if insetoption == 'initial':

        for k in range(res, res +arr ** 2):

            marker, _, color, label= labels_determine(file, k - res, arr=arr, plottag=plottag, verbose=verbose, colordetermine='sys')
            idx = np.argwhere(time < 2.0).flatten()[-1]
            o = cur_occ[:, k]
            ax_inset.scatter(time[:idx], o[:idx], label=label, color=color, marker=marker)

        ax_inset.set_yscale('log')

    else:

        if '2x2' in searchkey('config', file) :

            diff14 =  cur_occ[:, res] - cur_occ[:, res + arr**2 - 1]

            marker, _, color, label= labels_determine(file, 0, arr=arr, plottag=plottag, verbose=verbose, colordetermine='sys')
            ax_inset.plot(time, diff14, label=r'$\langle n \rangle_L - \langle n \rangle _R$', color=color)

            diff23 =  cur_occ[:, res + 1] - cur_occ[:, res + 2]
            marker, _, color, label= labels_determine(file, 2, arr=arr, plottag=plottag, verbose=verbose, colordetermine='sys')
            ax_inset.plot(time, diff23, label='offdiagonals', color=color)


           
        ax_inset.set_title('Difference between sites')
        ax_inset.legend()


def currents( groups : dict,  mosaic_template=[], width_ratios=[], filetitle='', factor=1.0, verbose=False, insetoption='compare degen'):




    
    ncol = figscale = 1 
    fontsize = figscale * 10

    groupnames = groups.keys()

    fig = plt.figure(
                    layout="constrained", 
                     figsize = ( figscale * 8*len(mosaic_template), figscale* 3 * len(groupnames))
                     )
    subfigs = fig.subfigures(nrows=len(groupnames), ncols=1)

    if len(groupnames) == 1:
        subfigs = [subfigs]

    axes = {}

    for i, groupname in enumerate(groupnames):
        
        group = groups[groupname]

        
        files = group['files']

        print(files)
        
        if 'plottags' in group:
            plottags = group['plottags']
        else:
            plottags = ['' for _ in files]

        grouptitle = group['grouptitle']


        #print("files :", files)
        #print("refs: ", ref)

        subfig : plt.Figure = subfigs[i]
        subfig.suptitle(grouptitle, fontsize=fontsize * 1.5)

        mosaic = [[ m + groupname for m in mosaic_template]]
        axes |= subfig.subplot_mosaic(mosaic,  width_ratios=width_ratios)
        
        SvNs, times, bonds, _, occs, currents, *vals, CCs  = dataloader([], files)
        
        ref = group['ref']

        if ref != []:

            if len(ref) > 1:
                raise ValueError("More than one ref!")

            _, reftime, _, _, refocc, refcurrent, *vals = dataloader([], ref)
            reftime = reftime[0]
            refcurrent = refcurrent[0].transpose()
            refocc = refocc[0]



        # get inset axes
        axes_inset = {}
        for key in [m + groupname for m in mosaic_template]:

            ax = axes[key]
                    #create inset
            
            
            if 'occ' not in key:

                if ref != []:
                    inset_dim=['50%', '40%']
                    ax_inset = inset_axes(ax, 
                                            width=inset_dim[0], height=inset_dim[1],loc='lower left', 
                                    bbox_to_anchor=[0.1, 0.1, 0.9, 0.9], bbox_transform=ax.transAxes
                                    )

                    axes_inset[key + 'inset'] = ax_inset
            
            else:
                inset_dim=['50%', '50%']
                ax_inset = inset_axes(ax, 
                                        width=inset_dim[0], height=inset_dim[1],loc='lower right', 
                                bbox_to_anchor=[0.0, 0.1, 0.9, 0.9], bbox_transform=ax.transAxes
                                )

                axes_inset[key + 'inset'] = ax_inset           


        for j, file in enumerate(files):

            time = times[j]
            current = currents[j]
            CC = CCs[j]
            SvN = SvNs[j]


            try:
                res = int(searchkey('Ls', file))
            except:
                res = 64

            occ = occs[j]
            plottag = plottags[j]

            if 'Electron' in file:
                spins = 2
                spinful = ['Up', 'Down']
                currenttotal = current.shape[-1]
                current = [current[:, :currenttotal//2], current[:, currenttotal//2:]]
            else:
                occ = [occ]
                #CC = [CC]
                SvN = [SvN]
                current = [current]
                spins = 1
                spinful = ['']

            for spin in range(spins):

                cur_occ = occ[spin]
                cur_current = current[spin]
                #cur_cc = CC[spin]
                #cur_svn = SvN[spin]
                
                currentcnt = cur_current.shape[-1]//2

                arr = int(np.sqrt(cur_occ.shape[-1] - res * 2))
                ax : plt.Axes = axes['source' + groupname]
                
                _, linestyle, color, label = labels_determine(file, -1, arr=arr, plottag=plottag + spinful[spin], verbose=verbose, colordetermine='sys')
                
                for k in range(currentcnt):
                    c = cur_current[:, k] * factor
                    ax.plot(time, c, label=label, color=color, linestyle=linestyle)

                    idx = np.argwhere(time < int(searchkey('Ls', file))).flatten()[-1]

                    if ref != []:
                        ax_inset : plt.Axes = axes_inset[ 'source' + groupname + 'inset']
                        ax_inset.plot(time[:idx], np.abs((c[:idx] - refcurrent[:idx, k])/refcurrent[:idx, k]), label=label, color=color, linestyle=linestyle)


                ax.set_title("Current : Source to Array")

                if ref != []:
                    ax_inset : plt.Axes = axes_inset[ 'source' + groupname + 'inset']
                    ax_inset.set_title('Current difference', fontsize=fontsize/2
                            )
                    ax_inset.set_ylabel('$|I/I_{{ref}}|$', fontsize=fontsize/2
                            )
                    ax_inset.set_xlabel('time', fontsize=fontsize/2)
                    ax_inset.tick_params(axis='both',  labelsize=fontsize/2)


                
                ax : plt.Axes = axes['occ' + groupname]

                for k in range(res, res +arr ** 2):

                    marker, linestyle, color, label= labels_determine(file, k - res, arr=arr, plottag=plottag, verbose=verbose, colordetermine='sys')
                    o = cur_occ[:, k]
                    ax.scatter(time, o, label=label + spinful[spin], marker=marker,  c=color, s=5,
                               #linestyle=linestyle,
                               )

                occ_inset(axes_inset, groupname, insetoption, res, arr, file, plottag, verbose, time, cur_occ)
                
                ax.set_yscale('log')
                ax.set_title('Array occupation')



                # ax : plt.Axes = axes['drain' + groupname]

                # marker, linestyle, color, label = labels_determine(file, -1, arr=arr, plottag=plottag + spinful[spin], verbose=verbose)
                # for k in range(currentcnt, currentcnt *2):
                #     c = cur_current[:, k] * factor
                #     ax.plot(time, c, label=plottag  + spinful[spin])

                # ax.set_title('Current: Array to Drain')

        if ref != []:
            ax : plt.Axes = axes['source' + groupname]
                
            for k in range(currentcnt):
                c = refcurrent[:, k]
                ax.plot(reftime, c, c='black', label='ref', linestyle='dotted')
        
            ax : plt.Axes = axes['occ' + groupname]

            arr = int(np.sqrt(refocc.shape[-1]))
            for k in range(arr ** 2):
                o = refocc[:, k]
                marker, linestyle, color, label = labels_determine('ref', k, arr=arr, verbose=verbose, colordetermine='sys')
                ax.plot(reftime, o, label=label, marker=marker, linestyle=linestyle, c=color, markersize=0.1)

        ax : plt.Axes = axes['source' + groupname]
        random_utils.no_duplicate_label(ax, fontsize, ncol=ncol)
        ax : plt.Axes = axes['occ' + groupname]
        random_utils.no_duplicate_label(ax, fontsize, ncol =ncol)

        # ax : plt.Axes = axes['drain' + groupname]
        # for k in range(currentcnt, currentcnt *2):
        #     c = refcurrent[:, k]
        #     ax.plot(reftime, c, c='black', linestyle='dotted', label='ref')
        
        # random_utils.no_duplicate_label(ax, fontsize, ncol=ncol)


        

            
    #fig.tight_layout()
    fig.savefig('/Users/knl20/Desktop/Code/TN/SD/plots/{}.pdf'.format(filetitle))
