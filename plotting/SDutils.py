
import matplotlib.cm as colormap
from TNplottingutils.searchkey import searchkey
import numpy as np

def labels_determine(file, k, arr=1, plottag='', colordetermine= 'invalid', verbose=True):

    systype = searchkey('reservoirtype', file)
    print(systype)
    dim = searchkey('TEdim', file)
    
    valdict = {
        '-0.001': 0.0,
        '-0.01' : 0.4,
        '-0.1' : 0.8,
        '1.0' : 0.0,
        '10.0' : 0.2,
        '100.0' : 0.5,
        '1000.0' : 0.8
    }

    if verbose:
        val = k/(arr ** 2)
        
    elif colordetermine == 'init':
        val = searchkey('biasAinit', file)
        val = valdict[val]
    elif colordetermine == 'coupling':
        val = searchkey('scoupling', file)
        val = valdict[val]

    elif colordetermine == "invalid":
        raise ValueError("unrecognized!")




    linedict = {
        '64': 'dotted',
        '128' : 'solid',
        '256' : 'dashed',
        '512' : 'dotted'
    }

    sysdict = {
        'mixed' : 'blue',
        'spatial' : 'red'
    }

    markerdict = {
        'mixed64' : '+',
        'mixed128' : 'o',
        'mixed256' : 'x',
        'mixed512' : '+',
        'spatial64' : '.',
        'spatial128' : '^',
        'spatial256' : 'v',
        'spatial512' : '.'
    }

    if 'ref' in file:
        linestyle = 'dotted'
        color = 'black'
        label = 'ref'
        marker = ''

    else:
        linestyle = linedict[dim]
        
        if verbose:
            cmap = colormap.hot if systype == 'mixed' else colormap.winter
            color = cmap(val)

            label = plottag + ', ' + systype + 'site ' + str(k + 1)
            marker = markerdict[ systype + dim]

        else:
            if colordetermine !=  'sys':

                cmap = colormap.hot if systype == 'mixed' else colormap.winter
                color = cmap(val)
            
            else:
                color = sysdict[systype]

            label = plottag + ', ' + systype + ', dim =' + dim
            marker = markerdict[ systype + dim]
    
    return marker, linestyle, color, label.strip(' ,')



def get_fermi(file):

    res = int(searchkey('Ls', file))
    arr = int(searchkey('config', file).split('x')[0])
    energy = np.loadtxt( file + '/BIASEDenergiesSORTED')

    LR = np.loadtxt(file + '/LR')
    left = np.argwhere( LR>0).flatten()
    right = np.argwhere( LR<0).flatten()

    leftenergy = energy[left]
    rightenergy = energy[right]


    leftind = np.where(left >= res, left + arr ** 2, left)
    rightind = np.where( right>= res, right + arr ** 2, right)

    print(min(left), max(left), min(right), max(right))

    return leftind, rightind, leftenergy, rightenergy


def gen_energy(k, L, mu):

    return 2 * np.cos( k * np.pi /(L + 1) )  + mu


def ind_determine(which_reservoir, left, right, L):
    """returns the total true index (TN) depending on which part of the system """

    if which_reservoir == 'source':
        inds = left
    
    elif which_reservoir == 'all':
        inds = np.arange( L)

    elif which_reservoir == 'array':
        inds = len(left) + np.arange(L - len(left) - len(right))

    elif which_reservoir =='both':
        inds = np.sort( np.concatenate( (left, right)))
        print(inds)

    else:
        inds = right

    

    return inds