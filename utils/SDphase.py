import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from itertools import product
import glob
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rc


def phase_range(cat='regular'):

    if cat == 'regular':
        biases = np.power(2.0, np.arange(-3, 4))
        mu_tes = np.concatenate(([0], np.power(2.0, np.arange(-3, 5))))

        Ls = [2 ** i  for i in range(6, 10)]
        dims = [1, 2, 3, 4]
        ss = [10.0 ** i for i in range(-3, 2)]
        ds = [10.0 ** i for i in range(-3, 2)]

    elif cat == 'small':
        biases = np.arange(0, 1.1, 0.1)
        mu_tes = np.arange(0, 1.1, 0.1)

        Ls = [64]
        dims = [1, 2, 3, 4]
        ss = [0.01]
        ds = [0.0001, 0.001, 0.01]      

    return Ls, np.round(biases, 5), ss, ds, np.round(mu_tes, 5), dims


def phase_loader(L, bias, s, d, mu_te, dim, cat='regular'):

    if cat == 'regular':
        path = 'SDscanjlnew'
        suffix = '.txt'
    elif cat == 'small':
        path = 'SDscanjlsmall'
        suffix = ''


    top_path = "/Users/knl20/Desktop/Code/non-interacting/{}/".format(path)
    path = path_gen(L, L, bias, s, d, mu_te, dim)
    cur = np.loadtxt( top_path + path + 'currentSD' + suffix)
    times = np.loadtxt( top_path + path + 'times' + suffix)

    return cur, times


def phase_slider():
    # Define the function to plot
    # Initial parameters for the function
    L_init, bias_init, s_init, d_init, mu_te_init, arr_init = 1, 1, 0, 1, 1, 0

    # Create the plot
    fig, axes = plt.subplots(1, 2, figsize=(15, 9))
    fig.subplots_adjust(left=0.25, bottom=0.4)

    axes: list[plt.Axes]

    # Adjust the axes for sliders
    slider_ax_positions = [0.1, 0.25, 0.4, 0.55, 0.7, 0.85]
    sliders : list[Slider] = []

    # Create sliders
    params = ["L", "bias", "s", "d", "mu_te", "arr"]
    initial_values = [L_init, bias_init, s_init, d_init, mu_te_init, arr_init]
    Ls, biases, ss, ds, mu_tes, dims  = phase_range()
    ranges = [Ls, biases, ss, ds, mu_tes, dims]

    for i, (param, init) in enumerate(zip(params, initial_values)):
        ax_slider = plt.axes([0.25, slider_ax_positions[i] * 0.35, 0.65, 0.03], facecolor="lightgoldenrodyellow")
        slider = Slider(ax_slider, param, min(ranges[i]), max(ranges[i]), valinit=init, valstep=ranges[i])
        sliders.append(slider)

    # Update function for sliders
    def update(val):
        L, bias, s, d, mu_te, dim = [slider.val for slider in sliders]

        cur, times = phase_loader(L, bias, s, d, mu_te, dim)
        for j in range(2):
            ax  = axes[j]
            ax.clear()
            ax.plot(times, cur[j])

        fig.suptitle("L{}_bias{}_s{}_d{}_onsite{}_{}x{}".format(L, bias, s, d, mu_te, dim, dim), #fontsize=20
                     )
        #fig.canvas.draw_idle()
        #fig.tight_layout()
        
        

    # Connect sliders to the update function
    for slider in sliders:
        slider.on_changed(update)

    
    #fig.savefig( "/Users/knl20/Desktop/Code/non-interacting/plots/slider.pdf")
    plt.show()

def path_gen(L, N, bias, s, d, mu_te, dim):

    path = "L{}_N{}_{}x{}_bias{}_muinit1000.0_mute{}_S{}_D{}_singletrue_SDbias/".format(
        int(L),
        N,
        int(dim), int(dim),
        bias,
        mu_te,
        s,
        d
    )


    return path
    # Ls = [2^i for i in 6:10]
    # biases = [2.0^i for i in -3:3]
    # ss = [trunc(10.0^i, digits=5) for i in -3:1]
    # ds = [trunc(10.0^i, digits=5) for i in -3:1]
    # arrs = [1, 2, 3, 4]
    # mu_tes = [0; 2.0.^(-3:4)]


def coarsephase(cat='regular'):

    paras  = phase_range(cat=cat)
    paras = list(product(*paras))
    toppath = "/Users/knl20/Desktop/Code/non-interacting/SDscanphase/"

    if not glob.glob(toppath + 'source' + cat):


        arr = int(np.ceil(np.sqrt(len(paras))))
        source = np.zeros(arr ** 2)
        drain = np.zeros(arr ** 2)
        
        for i, (L, bias, s, d, mu_te, dim )in enumerate(paras):

            print("load : {}".format(i))
            cur, _ = phase_loader(L, bias, s, d, mu_te, dim, cat)

            source[i] = cur[0][-1]
            drain[i] = cur[1][-1]

        # source = source.reshape((arr, arr))
        # drain = drain.reshape((arr, arr))
        np.savetxt(toppath + 'source' + cat, source)
        np.savetxt(toppath + 'drain' + cat, drain)

    source = np.loadtxt(toppath + 'source' + cat)
    drain = np.loadtxt(toppath + 'drain' + cat)

    return source, drain

def phase(option):

    def all():
        L = 256
        for dim in dims:

            fig, axes = plt.subplots( len(biases), len(mu_tes), figsize = ( 3 * len(mu_tes), 3 * len(biases)))
            fig.suptitle(r'${}\times{}$'.format(dim, dim), fontsize=fontsize)

            for i, bias in enumerate(biases):
                for j, mu_te in enumerate(mu_tes):

                    arr = np.zeros((len(ss), len(ds)))
                    ax : plt.Axes = axes[i][j]

                    for k, s in enumerate(ss):
                        for l, d in enumerate(ds):

                            ind = para_dict[ (L, bias, s, d, mu_te, dim)]
                            arr[k, l] = source[ind]

                    ax.imshow(arr, cmap='bwr', origin='lower',
                              vmax = shi, vmin=slo)
                    
                    ax.set_yticks(range(len(ss)))
                    ax.set_xticks(range(len(ds)))
                    ax.set_yticklabels([str(s) for s in ss])
                    ax.set_xticklabels([str(d) for d in ds])
                    ax.set_ylabel('source coupling')
                    ax.set_xlabel('drain coupling')
                    ax.set_title( 'bias = {}, $\mu_{{onsite}}$ = {}'.format(bias, mu_te))
            
            
            fig.tight_layout()
            pdf.savefig(fig)
    
    def SDonly():

        s = 0.01
        d = 0.001
        fig, axes = plt.subplots( len(dims), len(Ls), figsize = ( 5 * len(Ls), 3 * len(dims)))

        fig.suptitle('Weaker couplings: s = {}, d = {}'.format(s, d))

        for i, dim in enumerate(dims):
            for j, L in enumerate(Ls):

                arr = np.zeros((len(biases), len(mu_tes)))
                ax : plt.Axes = axes[i][j]

                for k, bias in enumerate(biases):
                    for l, mu_te in enumerate(mu_tes):

                        ind = para_dict[ (L, bias, s, d, mu_te, dim)]
                        arr[k, l] = source[ind]

                im = ax.imshow(arr, cmap='bwr', origin='lower')
                fig.colorbar(im)

                
                ax.set_yticks(range(len(biases)))
                ax.set_xticks(range(len(mu_tes)))
                ax.set_yticklabels([str(bias) for bias in biases])
                ax.set_xticklabels([str(mu_te) for mu_te in mu_tes])
                ax.set_ylabel('bias')
                ax.set_xlabel('onsite potential')
                ax.set_title( r'Array: ${}\times{}, L ={}$'.format(dim, dim, L))
        
        fig.tight_layout()
        pdf.savefig(fig)
    

    def small():

        s = ss[0]
        L = Ls[0]
        fig, axes = plt.subplots( len(dims), len(ds), figsize = ( 5 * len(dims), 3 * len(dims)))

        fig.suptitle('Weaker couplings: s = {}, L = 64'.format(s))

        for i, dim in enumerate(dims):
            for j, d in enumerate(ds):

                arr = np.zeros((len(biases), len(mu_tes)))
                ax : plt.Axes = axes[i][j]

                for k, bias in enumerate(biases):
                    for l, mu_te in enumerate(mu_tes):

                        ind = para_dict[ (L, bias, s, d, mu_te, dim)]
                        arr[k, l] = source[ind]

                im = ax.imshow(arr, cmap='bwr', origin='lower')
                fig.colorbar(im)

                
                ax.set_yticks(range(len(biases)))
                ax.set_xticks(range(len(mu_tes)))
                ax.set_yticklabels([str(bias) for bias in biases])
                ax.set_xticklabels([str(mu_te) for mu_te in mu_tes])
                ax.set_ylabel('bias')
                ax.set_xlabel('onsite potential')
                ax.set_title( r'Array: ${}\times{}, L ={}$'.format(dim, dim, L))
        
        fig.tight_layout()
        pdf.savefig(fig)
    
    if option == 'smallSD':
        cat = 'small'
    else:
        cat = 'regular'

    source, drain = coarsephase(cat=cat)
    Ls, biases, ss, ds, mu_tes, dims = phase_range(cat=cat)
    fontsize = 20

    para_dict = {val : i for i, val in enumerate((product(*phase_range(cat=cat))))}

    shi = np.amax(source)
    slo = np.amin(source) 

    if option == 'all':
        f = 'all.pdf'
    elif option == 'smallSD':
        f = 'SDsmall.pdf'
    else:
        f = 'SDonly.pdf'

    with PdfPages('plots/Phase' + f) as pdf:

        if option == 'all':
            all()

        elif option == 'smallSD':
            small()

        else:
            SDonly()




if __name__ == '__main__':
    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True) 

    phase('smallSD')
    #phase('all')
    #phase_slider()
    #coarsephase()