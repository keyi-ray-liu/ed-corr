import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from itertools import product
import glob
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rc
from SDplotter.SDgeneric import *
from SDplotter.SDgeneric import _generic

def phase_range(cat='CS') -> dict:


    if cat == 'CS':
        # biasSDs = np.arange(0, 2.1, 0.1)
        # biasAs = np.arange(0, 2.1, 0.1)

        # Ls = [64]
        # dims = [1, 2, 3, 4]
        # ss = [0.01]
        # ds = [0.001, 0.01, 0.1]      

        biasSDs = [1.0]
        G1 = np.loadtxt('CSScan/G1s')
        G2 = np.loadtxt('CSScan/G2s')

        Ls = [64]
        dims = [3]
        ss = [0.1]
        ds = [0.1]

        paras = {
            'L' : Ls,
            'arr' : dims,
            'biasSD': np.round(biasSDs, 5),
            'G1' : np.round(G1, 5),
            'G2' : np.round(G2, 5),
            's' : ss,
            'd' : ds,
            'biasAinit' : [1000.0],
            'single' : ['true']
        }

    elif cat == 'SD':

        path = 'SDScan/'
        paras = {
            'L' : np.loadtxt( path + 'Ls'),
            'arr' : np.loadtxt( path + 'arrs'),
            'biasSD': np.loadtxt( path + 'biasSDs'),
            'arrhop': np.loadtxt( path + 'arrhops'),
            's' : np.loadtxt(path+'ss'),
            'd' : np.loadtxt(path + 'ds'),
            'biasA' : [0.0],
            'biasAinit' : [1000.0],
            'single' : ['true']
        }


    return paras


def phase_loader(paras : dict, cat='CS'):


    if cat == 'CS':
        path = 'CSScan'
        suffix = ''
        
    elif cat == 'SD':
        path = 'SDScan'
        suffix = ''

    top_path = "/Users/knl20/Desktop/Code/non-interacting/{}/".format(path)
    path = path_gen(paras)
    cur = np.loadtxt( top_path + path + 'currentSD' + suffix)
    times = np.loadtxt( top_path + path + 'times' + suffix)

    return cur, times


# def phase_slider(cat='CS'):


#     def changescale(key, val):

#         print(key, val)
#         if key not in ('arr', 'biasSD'):
#             val = np.log2(val)

#         return val


#     def revertscale(key, val):

#         if key not in ('arr', 'biasSD'):
#             val = 2 ** val

#         return val
    
#     def scalemodifier(key):
        
#         mod = key
#         if key not in ('arr', 'biasSD'):
#             mod += ', $2^x$'

#         return mod
#     # Define the function to plot
#     # Initial parameters for the function

#     # Create the plot
#     fig, axes = plt.subplots(1, 2, figsize=(15, 9))
#     fig.subplots_adjust(left=0.25, bottom=0.4)

#     axes: list[plt.Axes]

#     # Adjust the axes for sliders


#     # Create sliders
#     #params = ["L", "biasSD", "s", "d", "biasA", "arr"]
#     #initial_values = [L_init, bias_init, s_init, d_init, biasA_init, arr_init]
#     paras = phase_range(cat=cat)

#     keys = paras.keys()

#     slider_keys = list(filter(lambda x: len(paras[x]) > 1, keys))

#     step = 0.8/(len(slider_keys) - 1)
#     slider_ax_positions = 0.1 + np.arange(len(slider_keys)) * step
#     sliders : dict[str, Slider] = {}

#     draw_para = {key : paras[key][0] for key in keys}
#     draw_para['N'] = draw_para['L']
#     #ranges = [Ls, biasSDs, ss, ds, biasAs, dims]

#     for i, key in enumerate(slider_keys):

#         val = paras[key]
#         ax_slider = plt.axes([0.25, slider_ax_positions[i] * 0.35, 0.65, 0.03], facecolor="lightgoldenrodyellow")

#         val = changescale(key, val)
#         slider = Slider(ax_slider, scalemodifier(key), min(val), max(val),  valstep=val)
#         sliders[key] = slider

#     # Update function for sliders
#     def update(val):
        
#         for key in slider_keys:
#             draw_para[key] = revertscale(key, sliders[key].val)

#         draw_para['N'] = draw_para['L']

#         cur, times = phase_loader(draw_para, cat=cat)
#         for j in range(2):
#             ax  = axes[j]
#             ax.clear()
#             ax.plot(times, cur[j])

#         #fig.suptitle("L{}_biasSD{}_s{}_d{}_onsite{}_{}x{}".format(L, biasSD, s, d, biasA, dim, dim), #fontsize=20)
#         #fig.canvas.draw_idle()
#         #fig.tight_layout()
        
        

#     # Connect sliders to the update function
#     for slider in sliders.values():
#         slider.on_changed(update)

    
#     #fig.savefig( "/Users/knl20/Desktop/Code/non-interacting/plots/slider.pdf")
#     plt.show()

def path_gen(para: dict):

    def process(key):

        val = para[key]

        if key == 'L' or key == 'N' or key == 'arr':
            val = int(val)
        
        return str(val)
            


    path = '_'.join([ key + process(key) for key in sorted(para.keys())]) + '/'
    return path
    # Ls = [2^i for i in 6:10]
    # biasSDs = [2.0^i for i in -3:3]
    # ss = [trunc(10.0^i, digits=5) for i in -3:1]
    # ds = [trunc(10.0^i, digits=5) for i in -3:1]
    # arrs = [1, 2, 3, 4]
    # biasAs = [0; 2.0.^(-3:4)]


def coarsephase(cat='CS'):

    paras  = phase_range(cat=cat)
    para_names = sorted(paras.keys())
    vals = [ paras[key] for key in para_names]
    combs = product(*vals)

    toppath = "/Users/knl20/Desktop/Code/non-interacting/SDscanphase/" + '_'.join(para_names)

    if not glob.glob(toppath + 'source' + cat):

        source = []
        drain = []
        
        for i, val in enumerate(combs):

            val_dict = { key : val[j] for j, key in enumerate(para_names)}
            val_dict['N'] = val_dict['L']

            print("load : {}".format(i))
            cur, _ = phase_loader(val_dict, cat=cat)

            source.append(cur[0][-1])
            drain.append(cur[1][-1])

        # source = source.reshape((arr, arr))
        # drain = drain.reshape((arr, arr))
        np.savetxt(toppath + 'source' + cat, source)
        np.savetxt(toppath + 'drain' + cat, drain)

    source = np.loadtxt(toppath + 'source' + cat)
    drain = np.loadtxt(toppath + 'drain' + cat)

    return source, drain

def phase(cat):


    def CS():
        ss = paras['s']
        Ls = paras['L']
        dims = paras['arr']
        ds = paras['d']
        biasSDs = paras['biasSD']
        G1s = paras['G1']
        G2s = paras['G2']

        s = ss[0]
        L = Ls[0]
        biasSD = biasSDs[0]
        fig, axes = plt.subplots( len(dims), len(ds), figsize = ( 6 * len(dims), 5 * len(dims)))

        if len(dims) == 1:
            axes = [axes]

        if len(ds) == 1:
            axes = [ [ax] for ax in axes]

        for i, dim in enumerate(dims):
            for j, d in enumerate(ds):

                arr = np.zeros((len(G1s), len(G2s)))
                ax : plt.Axes = axes[i][j]

                for k, G1 in enumerate(G1s):
                    for l, G2 in enumerate(G2s):

                        val_dict = {
                            'L' : L,
                            'arr' : dim,
                            'biasSD': np.round(biasSD, 5),
                            'G1' : G1,
                            'G2' : G2,
                            's' : s,
                            'd' : d,
                            'biasAinit' : 1000.0,
                            'single' : 'true'
                        }
                        
                        val_tuple = tuple( [val_dict[key] for key in para_names])
                        print(val_tuple)

                        ind = para_dict[ val_tuple]
                        arr[k, l] = source[ind]

                im = ax.imshow(arr, cmap='bwr', origin='lower')
                fig.colorbar(im)

                
                ax.set_yticks(range(len(G1s)))
                ax.set_xticks(range(len(G2s)))
                ax.set_yticklabels([str(G1) if i% (len(G1s)//5) == 0 else '' for i, G1 in enumerate(G1s)], fontsize=fontsize)
                ax.set_xticklabels([str(G2) if i% (len(G2s)//5) == 0 else '' for i, G2 in enumerate(G2s)], fontsize=fontsize)

                ax.set_ylabel('G1', fontsize=fontsize)
                ax.set_xlabel('G2', fontsize=fontsize)
                ax.set_title( r'Array: ${}\times{}, s ={}, d ={}$'.format(dim, dim, s, d), fontsize=fontsize)
        
        #fig.suptitle('Weaker couplings: s = {}, L = 64'.format(s))
        fig.tight_layout()
        pdf.savefig(fig)
    
    

    source, drain = coarsephase(cat=cat)
    paras = phase_range(cat=cat)
    para_names = sorted(paras.keys())
    vals = [ paras[key] for key in para_names]

    fontsize = 20

    para_dict = {val : i for i, val in enumerate((product(*vals)))}
    print(para_dict)

    shi = np.amax(source)
    slo = np.amin(source) 

    if cat == 'CS':
        f = 'CS.pdf'

    with PdfPages('plots/Phase' + f) as pdf:

        if cat == 'CS':
            CS()


def slide():


    fig, _ = plt.subplots(figsize=(25, 9))
    fig.clear()
    fig.subplots_adjust(left=0.1, bottom=0.35)

    step = 0.2/5
    constr = '19'
    slider_ax_positions = 0.1 + np.arange(5) * step


    s_ax = plt.axes([0.2, slider_ax_positions[0] , 0.6, 0.03], facecolor="lightgoldenrodyellow")
    s_val = [-100.0, -0.5, -0.1, -0.01]
    #s_slider = Slider(s_ax, 's', min(s_val), max(s_val), valstep=s_val)
    s_slider = Slider(s_ax, 's', 0, len(s_val) - 1, valstep=1)

    d_ax = plt.axes([0.2, slider_ax_positions[1] , 0.6, 0.03], facecolor="lightgoldenrodyellow")
    #d_slider = Slider(d_ax, 'd', min(s_val), max(s_val), valstep=s_val)
    d_slider = Slider(d_ax, 'd', 0, len(s_val) - 1, valstep=1)

    w_ax = plt.axes([0.2, slider_ax_positions[2] , 0.6, 0.03], facecolor="lightgoldenrodyellow")
    #b_val = np.round(np.arange(-2, 2, 0.2), 3) 
    w_val = [-10.0, -1.0, -0.1]
    w_slider = Slider(w_ax, r'$\omega$', 0, len(w_val) - 1, valstep=1)

    L_ax = plt.axes([0.2, slider_ax_positions[3] , 0.6, 0.03], facecolor="lightgoldenrodyellow")
    L_val = [ int(val) for val in 2.0 ** np.arange(9)]
    #L_slider = Slider(L_ax, 'L', min(L), max(L), valstep=L)
    L_slider = Slider(L_ax, 'L', 0, len(L_val), valstep=1)

    t_ax = plt.axes([0.2, slider_ax_positions[4] , 0.6, 0.03], facecolor="lightgoldenrodyellow")
    t_slider = Slider(t_ax, 't', 0, 1000, valstep=1.0)

    # Update function for sliders
    def update(val):

        for subfig in fig.subfigs:
            subfig.clear()

        s = s_val[s_slider.val]
        d = s_val[d_slider.val]
        w = w_val[w_slider.val]
        L = L_val[L_slider.val]

        bias = float(0)
        t = t_slider.val

        f = lambda x: sorted(x, key= lambda y: float(y.split('=')[-1]) )
        cs = lambda ref : current(suppress_label=True, func=f, drawlines=[L], reference=ref, avg=True, xuplim=t)
        cd = lambda ref: current( which_lead=1, func=f, drawlines = [L], reference=ref, avg=True, xuplim=t)
        #socc = lambda x, ref : occ(suppress_label=True, width = 0.5,  func=f, single_site=x, reference=ref, xuplim=t)


        fstr = f'/Users/knl20/Desktop/Code/non-interacting/LengthScale/L{L}_N{L}_arr3_arrhop-1.0_biasA{0.0}_biasAinit1000.0_biasSD0.0_con{constr}_d{d}_modeleft_s{s}_spacing{w}'

        print(fstr)
        #transport = [[cs(None), *[socc(i+j*3, None) for i in range(3)], cd(None)] for j in range(3)]
        transportone = [[cs(None), occ(width =1.0,  func=f, xuplim=t, drawlines = [L]), cd(None)]]
        ed = [
        {
            'paneltitle' : f'Non-int, 3x3, L = {L}, N = {L}, scoupling = {s}, dcoupling = {d}, bias ={bias}, con = {constr}, $\omega $= {w}',
            'panelfuncs' : transportone,
            'panelparameters' : {
                    'dirs': [fstr],
                    'reservoirsizes' : [1],
                    'linestyles': ['solid', 'dashed'],
                    'markers' : ['', 'o'],
                    'labels' : [''],
                    'colors' : ['red', 'blue']
            }
        } 
        ]


        _generic(fig, ed)

            #print("no show")
        

        
    # Connect sliders to the update function
    s_slider.on_changed(update)
    L_slider.on_changed(update)
    d_slider.on_changed(update)
    w_slider.on_changed(update)
    t_slider.on_changed(update)
    

    #fig.savefig( "/Users/knl20/Desktop/Code/non-interacting/plots/slider.pdf")
    plt.show()



def select():

    fig, _ = plt.subplots(figsize=(12, 3.5))
    fig.clear()
    L = 1
    bias = 0.0
    s = -0.1
    d = -0.1
    t = 2.0
    constr = '19'

    fstr = f'/Users/knl20/Desktop/Code/non-interacting/LengthScale/L{L}_N{L}_arr3_arrhop1.0_biasA{bias}_biasAinit1000.0_biasSD0.0_con{constr}_d{d}_s{s}_spacing1.0'


    print(fstr)
    f = lambda x: sorted(x, key= lambda y: float(y.split('=')[-1]) )
    cs = lambda ref : current(suppress_label=True, func=f, drawlines=[L], reference=ref, avg=True, xuplim=t)
    cd = lambda ref: current( which_lead=1, func=f, drawlines = [L], reference=ref, avg=True, xuplim=t)
    #transport = [[cs(None), *[socc(i+j*3, None) for i in range(3)], cd(None)] for j in range(3)]
    transportone = [[cs(None), occ(width =1.0,  func=f, xuplim=t), cd(None)]]
    ed = [
    {
        'paneltitle' : '', #f'Non-int, 3x3, L = {L}, N = {L}, scoupling = {s}, dcoupling = {d}, bias ={bias}',
        'panelfuncs' : transportone,
        'panelparameters' : {
                'dirs': [fstr],
                'reservoirsizes' : [0],
                'linestyles': ['solid', 'dashed'],
                'markers' : ['', 'o'],
                'labels' : [''],
                'colors' : ['red', 'blue']
        }
    } 
    ]


    fig = _generic(fig, ed)
    fig.tight_layout()
    fig.savefig(f'/Users/knl20/Desktop/UMD/March2025/nonintL{L}s{s}d{d}t{t}con{constr}.pdf')


if __name__ == '__main__':
    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True) 

    #phase('CS')
    #phase('all')
    #phase_slider(cat='SD')
    slide()
    #select()
    #coarsephase()