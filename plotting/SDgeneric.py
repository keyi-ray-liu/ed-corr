import matplotlib.pyplot as plt
import numpy as np
from TNplottingutils.dataloader import dataloader
from TNplottingutils.searchkey import searchkey
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from TNplottingutils.random_utils import *
from matplotlib.widgets import  Slider
from SDplotter.SDutils import get_fermi, gen_energy, ind_determine
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import image
from types import NoneType
from glob import glob
from matplotlib import cm

class plotter:

    def __init__(self, name : str = '', width : float = 1.0, height : float = 1.0):
        self.name = name
        self.width = width
        self.height = height
        pass

    def get_name(self) -> str :
        return self.name 
    
    def get_width(self) -> float :
        return self.width
    
    def get_height(self) -> float:
        return self.height

    def get_data(self) :
        """return SvN, times, bonds, seffs, occs, currents, ee_lo, ee_hi, bd_lo, bd_hi, s_lo, s_hi, tmin, tmax, CCs"""
        raise NotImplementedError
    
    def plot(self) : 
        raise NotImplementedError

class Illustrate(plotter):

    def __init__(self, width=1.3, drive = 'gs'):
        super().__init__(name='illu', width=width)
        self.drive = drive

    def plot(self, ax: plt.Axes, *args, **kwargs):
        
        if self.drive == 'gs':
            im = image.imread('/Users/knl20/Desktop/images/SD/GS.png')
            ax.imshow(im)
            ax.set_axis_off()

        else:
            im = image.imread('/Users/knl20/Desktop/images/SD/Left.png')
            ax.imshow(im)
            ax.set_axis_off()

    


class TotalCharge(plotter):


    def __init__(self, width=1.0, suppress_label=False, func=sorted, height=1.0):
        super().__init__(name='totalocc' , width=width, height=height)
        self.suppress_label = suppress_label
        self.label_func = func

    def get_data(self, parameters : dict):
        return occ().get_data(parameters)

    def plot(self,  ax: plt.Axes,  parameters : dict, **kwargs):

        times, occall, reservoir_sizes = self.get_data(parameters)
        linestyles = ['solid', 'dashed']#parameters['linestyles']
        colors = parameters['colors']
        labels = parameters['labels']
        markers = parameters['markers']

        spindict = {0: 'Up', 1: 'Dn'}
        
        
        for i, occs in enumerate(occall):

            time  = times[i]

            if type(time) == NoneType:
                print("No input occ data!")
                continue
            
            for spin, occ in enumerate(occs):

                
                linestyle = linestyles[spin]
                reservoir_size = reservoir_sizes[i]
                color = colors[i]
                label = spindict[spin] + labels[i]
                marker = markers[i]


                arr =  int(np.sqrt(occ.shape[-1] - reservoir_size * 2))
                o = np.sum(occ[:, reservoir_size :reservoir_size +arr ** 2], axis=1)

                
                ax.plot( time, o, linestyle=linestyle, color=color, label=label, marker=marker)
                    
        ax.set_title(r'$\sum\langle n_{{array}}\rangle $')

        if not self.suppress_label:
            no_duplicate_label(ax, fontsize=10, ncol=1, func=self.label_func)

        ax.set_xlabel(r'time $(1/\bar{t})$')
        ax.set_ylabel(r'$\sum\langle n\rangle$')





class occ(plotter):

    def __init__(self, width=0.7, suppress_label=False, func=sorted, single_site = None, drawlines = [], reference=None, height=1.0, xuplim=None):
        super().__init__(name='occ' + str(single_site), width=width, height=height)
        self.suppress_label = suppress_label
        self.label_func = func
        self.single_site = single_site
        self.reference = reference
        self.xuplim= xuplim
        self.drawlines = drawlines

    def get_data(self, parameters : dict):


        dirs = parameters['dirs']
        _, times, _, _, occs, *_  = dataloader([], dirs, tlim = self.xuplim)

        #print(occs)
        try:
            reservoir_sizes = [ int(searchkey('Ls', f)) for f in dirs]

        except:
            try:
                reservoir_sizes = parameters['reservoirsizes']
            except:
                reservoir_sizes = [ 128 for f in dirs]

        return times, occs, reservoir_sizes
    

    
    def referenceplot(self, ax : plt.Axes, k):
        # L128_N64_arr3_arrhop1.0_biasA0.0_biasAinit1000.0_biasSD0.25_d0.5_s0.5_singletrue_spacing1.0
        reference = self.reference

        times = np.loadtxt(reference + '/times')


        occ = np.loadtxt(reference + '/occ')
        ax.plot(times, occ, color='blue', linestyle='dashed', label='non-interacting')


    def plot(self,  ax: plt.Axes,  parameters : dict, **kwargs):

        reference = self.reference
        times, occall, reservoir_sizes = self.get_data(parameters)
        linestyles = ['solid', 'dashed']#parameters['linestyles']
        colors = parameters['colors']
        labels = parameters['labels']
        markers = parameters['markers']
        markersize = 2
        drawlines = self.drawlines

        spindict = {0: 'Up', 1: 'Dn'}
        single_site = self.single_site
        
        
        for i, occs in enumerate(occall):


            time  = times[i]
            

            if type(time) == NoneType:
                print("No input occ data!")
                continue
            
            for spin, occ in enumerate(occs):

                

                linestyle = linestyles[spin]
                reservoir_size = reservoir_sizes[i]
                color = colors[i]

                if occs.shape[0] == 1:
                    label =  f'site {labels[i]}'

                else:
                    label = spindict[spin] + labels[i]

                marker = markers[spin]

                arr =  int(np.sqrt(occ.shape[-1] - reservoir_size * 2))

                if single_site == None:

                    degen_marker = ['', 'o', 'x', '^']
                    for k in range(reservoir_size, reservoir_size +arr ** 2):

                        o = occ[:, k]

                        k_internal = k - reservoir_size
                        ax.plot( time, o, linestyle=linestyle, label=label + f'{ k_internal+ 1}', color=cm.bwr(k_internal/(arr ** 2)), marker = degen_marker[ k% len(degen_marker)])


                else:
                    k = reservoir_size + single_site
                    o = occ[:, k]

                    if reference != None:
                        self.referenceplot(ax, k)

                    ax.plot( time, o, linestyle=linestyle, color=color, label=label, marker=marker, markersize=markersize)

        ax.set_title(r'$\langle n_{{array}}\rangle $' + '' if single_site == None else f'site {single_site + 1}')

        vmin, vmax = ax.get_ylim()
        ax.vlines(drawlines, vmin, vmax, linestyle='dotted')

        ax.set_ylim( *ax.get_ylim())

        if not self.suppress_label:
            no_duplicate_label(ax, fontsize=10, ncol=1, bbox = (1.05, 0.6), func=self.label_func)

        ax.set_xlabel(r'time $(1/\bar{t})$')
        ax.set_ylabel(r'$\langle n\rangle$')


class Empty(plotter):

    def __init__(self, id=0, width = 0.5):
        super().__init__(name = 'empty' + str(id), width=width)

    def plot(self, ax : plt.Axes, parameters, fig = None):
        ax.set_axis_off()





class EE(plotter):

    def __init__(self, width = 1, id = 0, title='', which_reservoir = 'source'):
        super().__init__(name = 'EE' + str(id) + which_reservoir, width=width)
        # id is used to internally identify panel position
        self.id = id
        self.which_reservoir=which_reservoir
        self.title = title + f'{which_reservoir}: ' + 'Entanglement Entropy on bipartite cuts'

    def get_data(self, parameters : dict):

        dirs = parameters['dirs']
        id = self.id
        SvN, times, *_ = dataloader([], dirs)

        left, right, _, _ = get_fermi(dirs[id])

        return times[id], SvN[id], left, right

    
    def plot(self, ax : plt.Axes, parameters: dict, fig:plt.Figure = None):

        time, SvN, left, right = self.get_data(parameters)

        inds = ind_determine(self.which_reservoir, left, right, SvN.shape[-1])

        extent = [0, time[-1], len(inds), 1]
        
        SvN = SvN[:, inds]
        d1, d2 = SvN.shape
        ylabel = 'time'
        xlabel = f'site number ({self.which_reservoir})'

        if d1 > d2:
            SvN = np.transpose(SvN)
            ylabel, xlabel = xlabel, ylabel

        im = ax.imshow(SvN,  interpolation='none', extent=extent, aspect='auto')
        fig.colorbar(im)
        ax.set_title(self.title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)



class Fermi(plotter):

    def __init__(self, width = 1, id = 0, title='', which_reservoir='source'):
        super().__init__(name = 'Fermi' + str(id) + which_reservoir, width=width)
        # id is used to internally identify panel position
        self.id = id
        self.which_reservoir = which_reservoir
        self.title = title + f'{which_reservoir}: ' + r'$\langle n_{{array}}(t)\rangle  - \langle n_{{array}}(0)\rangle$'

    def get_data(self, parameters : dict):

        dirs = parameters['dirs']
        id = self.id
        _, times, _, _, occs, *_  = dataloader([], dirs)

        left, right, leftenergy, rightenergy = get_fermi(dirs[id])

        return times[id], occs[id], left, right
    
    def plot(self, ax : plt.Axes, parameters: dict, fig:plt.Figure = None):

        time, occ, left, right = self.get_data(parameters)

        occ = occ[0]
        occ = occ - occ[0]

        inds = ind_determine(self.which_reservoir, left, right, occ.shape[-1])

        occ = occ[:, inds]


        d1, d2 = occ.shape
        ylabel = 'time'
        xlabel = f'site number ({self.which_reservoir})'
        #(Equivalent to wavenumber $k$, as coupling is given by $\sim \sin k\pi/L)$'

        if d1 > d2:
            occ = np.transpose(occ)
            ylabel, xlabel = xlabel, ylabel

        extent = [0, time[-1], len(inds), 1]
        im = ax.imshow(occ, interpolation='none', cmap= 'hot', extent=extent, aspect='auto')

        fig.colorbar(im)
        ax.set_title(self.title)
        ax.set_ylabel(ylabel)
        
        ax.set_xlabel(xlabel)

class Occlandscape(plotter):

    def __init__(self, width = 1, id = 0, title='', which_reservoir='source', mu=0.0):
        super().__init__(name = 'Occlandscape' + str(id) + which_reservoir, width=width)
        # id is used to internally identify panel position
        self.id = id
        self.which_reservoir = which_reservoir
        self.title = title + f'{which_reservoir}: ' + r'$\langle n_{{array}}(t_{{last}})\rangle  - \langle n_{{array}}(0)\rangle$'

        # if which_reservoir == 'source' : 
        #     self.mu = mu
        
        # elif which_reservoir == 'drain' :
        #     self.mu = -mu

        # else:
        #     self.mu = 0

    def get_data(self, parameters : dict):

        dirs = parameters['dirs']
        id = self.id
        _, times, _, _, occs, *_  = dataloader([], dirs)

        left, right,  leftenergy, rightenergy  = get_fermi(dirs[id])

        if self.which_reservoir == "source":
            energy =  leftenergy
        
        elif self.which_reservoir == "drain":
            energy = rightenergy

        return times[id], occs[id], left, right,  energy
    
    def plot(self, ax : plt.Axes, parameters: dict, fig:plt.Figure = None):

        time, occ, left, right, energy = self.get_data(parameters)

        occ = occ[0]
        occ = occ[-1] - occ[0]

        inds = ind_determine(self.which_reservoir, left, right, occ.shape[-1])
        occ = occ[ inds]

        threshold = 0.05 * (max(occ) - min(occ))

        peak_ind, _ = find_peaks(np.abs(occ), prominence=threshold)

        ax.scatter(peak_ind, occ[peak_ind])

        for ind in peak_ind:
            ax.annotate(str(ind), (ind - 10, occ[ind]))

        ylabel = r'$\langle n_{{array}}(t_{{last}})\rangle  - \langle n_{{array}}(0)\rangle$'
        xlabel = f'site number ({self.which_reservoir})'

        ax.plot(np.arange(len(occ)), occ, color='black')
        ax.set_title(self.title)
        ax.set_ylabel(ylabel)
        ax.set_xlabel(xlabel)



        ax_inset : plt.Axes = inset_axes(ax, 
                                width=1.0, height=1.0,loc='lower right', 
                                    bbox_to_anchor=[0.1, 0.1, 0.9, 0.9], bbox_transform=ax.transAxes
                                    )
        
        config = searchkey('config', parameters['dirs'][self.id])
        ref = np.loadtxt('/Users/knl20/Desktop/Code/ED/ref/energy' + config)
        
        energies = energy[peak_ind]
        ax_inset.hlines(energies , 0.2, 0.3)
        ax_inset.hlines(ref, 0.7, 0.8)
        ax_inset.set_xlim(0, 1)
        ax_inset.set_ylabel(r'$\bar{t}$')
        #ax_inset.get_xaxis().set_visible(False)

        ax_inset.set_xticks([0.25, 0.75])
        ax_inset.set_xticklabels([ r'$\bar{t}_{{peak}}$', r'$\bar{t}_{{ED}}$'])

class DMRGConvergence(plotter):

    def __init__(self, name = 'dmrgconvergence', width = 1):
        super().__init__(name, width)


    def get_val(self, case : str):

        casestr = searchkey('Na', case)

        

        Nup = int(casestr.split('Up')[0])
        Ndn = int(casestr.split('Dn')[0].split('Up')[-1])
        U = searchkey('U', case)
        print("Nup, Ndn, U:", Nup, Ndn, U)

        return Nup, Ndn, U
    
    def process_lines(self, dir :str):

        dir = glob(dir + '/sl*')[0] 

        with open(dir, "r") as f:
            lines = f.readlines()

        val = []
        for line in lines:
            if 'energy' in line and 'sweep' in line:
                chunks = line.split()

                for chunk in chunks:
                    if 'energy' in chunk:
                        val.append( float(chunk.split('=')[-1]))

        return np.array(val)

    def plot(self, ax : plt.Axes, parameters : dict, **kwargs):

        labels = parameters['labels']
        for i, case in enumerate(parameters['dirs']):

            Nup, Ndn, U = self.get_val(case)

            ed = f'/Users/knl20/Desktop/Code/ED/ref/energy3x3-single{Nup}Up{Ndn}DnU{U}'
            ed = np.loadtxt(ed)[0]

            try:
                vals = self.process_lines(case)
                ax.plot( np.abs(vals - ed), label = labels[i], c=parameters['colors'][i])
            except:
                pass

        ax.set_yscale('log')
        ax.set_xlabel('Iteration number')
        ax.set_ylabel('Difference from ED')
        ax.legend()

        


        


class Occcompare(plotter):


    def __init__(self, width=0.7, suppress_label=False, func=sorted, single_site = None, ):
        super().__init__(name='occcomp' + str(single_site), width=width)
        self.suppress_label = suppress_label
        self.label_func = func
        self.single_site = single_site
        #self.reference = reference

    def get_data(self, parameters : dict):
        return occ().get_data(parameters)
    
    def plot_ref(self, ax : plt.Axes, parameters):

        ref = parameters['ref']
        up = np.loadtxt(ref + 'occup')[0]
        dn = np.loadtxt(ref + 'occdn')[0]


        ax.scatter( np.arange(up.shape[0]), up, color='blue', label = 'Ref Up', s=200)
        ax.scatter( np.arange(dn.shape[0]), dn, color='blue', marker ='x', label = 'Ref Dn', s=200)

    

    def plot(self,  ax: plt.Axes,  parameters : dict, **kwargs):

        _, occall, reservoir_sizes = self.get_data(parameters)
        #linestyles = ['solid', 'dashed']
        markers = ['^', 'v']
        colors = parameters['colors']
        labels = parameters['labels']
        markersize = 2

        spindict = {0: 'Up', 1: 'Dn'}
        self.plot_ref(ax, parameters)
        
        for i, occs in enumerate(occall):
            
            #print(f"occ compare {i}", searchkey('reservoirspacing', parameters['dirs'][i]))
            if type(occs) == NoneType:
                print("No input occ for compare!")
                continue
            
            for spin, occ in enumerate(occs):

                marker = markers[spin]
                reservoir_size = reservoir_sizes[i]
                color = colors[i]
                label = spindict[spin] + labels[i]




                arr =  int(np.sqrt(occ.shape[-1] - reservoir_size * 2))


                o = occ[0, reservoir_size: reservoir_size+ arr ** 2]

                #ax.plot( np.arange(o.shape[0]), o, linestyle=linestyle, color=color, label=label, linewidth = 0.5 * (len(occall)  - i))
                ax.scatter( np.arange(o.shape[0]), o, marker=marker, color=color, label=label, s = markersize * (len(occall)  - i), )


        ax.set_title('Compare array occ at t=0')

        

        if not self.suppress_label:
            no_duplicate_label(ax, fontsize=10, ncol=1, func=self.label_func)

        ax.set_xlabel('site number')
        ax.set_ylabel(r'$\langle n\rangle$')




class EffectiveEntanglement(plotter):

    def __init__(self, name = 'ee', width = 1.0, start = 0, end = 1, func=sorted):
        super().__init__(name=name + f'start{start}end{end}', width = width)
        self.start = start
        self.end = end
        self.func = func


    def get_data(self, parameters : dict):


        dirs = parameters['dirs']
        SvN, times,  *_  = dataloader([], dirs)
        return times, SvN
    

    def plot(self, ax : plt.Axes, parameters : dict, **kwargs):

        times, SvNs = self.get_data(parameters)
        labels = parameters['labels']
        colors = parameters['colors']
        start = self.start 
        end = self.end

        for i, SvN in enumerate(SvNs):

            if type(SvN) == NoneType:
                continue

            ee = SvN[:,  start : end]
            seff =  np.log( np.power( np.sum( np.exp( 3 * ee) , axis=1) / (ee.shape[-1] ), 1/3))

            ax.plot(times[i], seff, label = labels[i], color = colors[i])

        ax.set_title(f'Effective entanglement : start site = {start + 1}, end site = {end }')
        no_duplicate_label(ax, fontsize=10, ncol=1, func=self.func)

    

class current(plotter):

    def __init__(self, width=1.0, xuplim=None, which_lead=0, suppress_label = False, func = sorted, drawlines = [], reference = None, avg = False):
        super().__init__(name='current' + str(which_lead), width=width)
        self.xuplim = xuplim
        self.which_lead = which_lead
        self.lead_title = 'Source to Array' if which_lead == 0 else 'Array to Drain'
        self.suppress_label = suppress_label
        self.label_func = func
        self.drawlines = drawlines
        self.reference = reference
        self.avg = avg

    def get_data(self, parameters : dict):


        dirs = parameters['dirs']
        _, times, _, _, _, currents, *_  = dataloader([], dirs, tlim = self.xuplim)
        return times, currents

    def referenceplot(self, ax : plt.Axes, hmax):
        # L128_N64_arr3_arrhop1.0_biasA0.0_biasAinit1000.0_biasSD0.25_d0.5_s0.5_singletrue_spacing1.0
        reference = self.reference

        times = np.loadtxt(reference + '/times')

        inds = np.argwhere(times < hmax).flatten()
        cur = -np.loadtxt(reference + '/currentSD')[self.which_lead][inds]
        times = times[inds]
        

        ax.plot(times, cur, color='blue', linestyle='dashed', label='non-interacting')




    def plot(self,  ax: plt.Axes,  parameters : dict, **kwargs):

        times, currents = self.get_data(parameters)
        linestyles = parameters['linestyles']
        colors = parameters['colors']
        labels = parameters['labels']
        markers = parameters['markers']
        drawlines = self.drawlines
        reference = self.reference
        avg = self.avg
        markersize = 2
        
        if 'plot_funcs' in parameters:
            plot_funcs = parameters['plot_funcs']

        else:
            plot_funcs = [ 'plot' for _ in range(len(labels))]

        
        for i, current in enumerate(currents):

            time  = times[i]

            if type(time) == NoneType:
                print("No input current!")
                continue

            c = current[:, self.which_lead]

            if avg:
                c = np.concatenate(([0], c))
                c = (c[:-1] + c[1:])/2

            if current.shape[-1] == 4:
                cdn = current[:, self.which_lead + 2]
            else:
                cdn = None

            color = colors[i]
            label = labels[i] 
            plot_func = plot_funcs[i]

            if plot_func == 'plot':
                ax.plot( time, c, linestyle=linestyles[0], color=color, label=label, marker=markers[0], markersize=markersize)

                if type(cdn) != NoneType:
                    ax.plot( time, cdn, linestyle = linestyles[1], color=color, label=label + ' Dn', marker=markers[1], markersize=markersize)


            else:
                ax.scatter(time, c,  color=color, label=label)

        vmin, vmax = ax.get_ylim()
        hmin, hmax = ax.get_xlim()
        ax.vlines(drawlines, vmin, vmax, linestyle='dotted')

        if reference:
            self.referenceplot(ax, hmax)

        ax.set_xlim(hmin, hmax)

        ax.set_title("Current : " + self.lead_title)

        if not self.suppress_label:
            no_duplicate_label(ax, fontsize=10, ncol=1, func=self.label_func)

        ax.set_xlabel(r'time $(1/\bar{t})$')
        ax.set_ylabel(r'Current $(\bar{t})$')

    


def generic(pdf : PdfPages , panels : list[dict], width=5, length=6, save_individual = False, filename = '', figsuptitle = '') -> None :


    fig = plt.figure(
                    layout="constrained", 
                     figsize = ( width * max([len(panel['panelfuncs'][0]) for panel in panels]), length * len(panels))
                     )
    

    fig : plt.Figure = _generic(fig, panels, figsuptitle=figsuptitle)
    
    pdf.savefig(fig, dpi=300)

    if save_individual:
        fig.savefig(f'/Users/knl20/Desktop/Code/TN/SD/plots/SDindividual/{filename}.pdf' )





def _generic(fig :plt.Figure,  panels : list[dict], figsuptitle = ''):

    fontsize = 10
    fig.suptitle(figsuptitle, fontsize=fontsize * 4)
    subfigs = fig.subfigures(nrows=len(panels), ncols=1)


    if len(panels) == 1:
        subfigs = [subfigs]


    for i, panel in enumerate(panels):

        paneltitle : str = panel['paneltitle']
        panelparameters : dict = panel['panelparameters']
        panelfuncs : list[list[plotter]] = panel['panelfuncs']

        axes = {}

        subfig : plt.Figure = subfigs[i]
        subfig.suptitle(paneltitle, fontsize=fontsize * 2)

        mosaic = [[ func.get_name() + paneltitle for func in panelfuncsrows] for panelfuncsrows in panelfuncs]
        width_ratios = [ func.get_width() for func in panelfuncs[0]]
        height_ratios = [ func[0].get_height() for func in panelfuncs]

        print("w ratio, h ratio: ", width_ratios, height_ratios)

        axes |= subfig.subplot_mosaic(mosaic,  width_ratios=width_ratios, height_ratios=height_ratios)


        for _, func in enumerate(set(flatten2d(panelfuncs))):

            ax : plt.Axes = axes[ func.get_name() + paneltitle]
            func.plot(ax, panelparameters, fig=subfig)

    return fig



def slide():

    Nup = 2
    Ndn = 2
    fig, _ = plt.subplots(figsize=(25, 9))
    fig.clear()
    fig.subplots_adjust(left=0.1, bottom=0.35)

    step = 0.2/5
    slider_ax_positions = 0.1 + np.arange(5) * step

    #ranges = [Ls, biasSDs, ss, ds, biasAs, dims]
    U_ax = plt.axes([0.2, slider_ax_positions[0] , 0.6, 0.03], facecolor="lightgoldenrodyellow")
    U_slider =  Slider(U_ax, 'U', 0, 100,  valstep=10.0)

    s_ax = plt.axes([0.2, slider_ax_positions[1] , 0.6, 0.03], facecolor="lightgoldenrodyellow")
    s_val = [-1.0, -0.5, -0.1, -0.01, -0.001]
    s_slider = Slider(s_ax, 's', min(s_val), max(s_val), valstep=s_val)

    d_ax = plt.axes([0.2, slider_ax_positions[2] , 0.6, 0.03], facecolor="lightgoldenrodyellow")
    d_slider = Slider(d_ax, 'd', min(s_val), max(s_val), valstep=s_val)

    b_ax = plt.axes([0.2, slider_ax_positions[3] , 0.6, 0.03], facecolor="lightgoldenrodyellow")
    b_val = np.round(np.arange(-2, 2, 0.2), 3) 

    b_slider = Slider(b_ax, 'b', min(b_val), max(b_val), valstep=b_val)


    t_ax = plt.axes([0.2, slider_ax_positions[4] , 0.6, 0.03], facecolor="lightgoldenrodyellow")
    t_slider = Slider(t_ax, 't', 0, 1000, valstep=1.0)

    # Update function for sliders
    def update(val):

        for subfig in fig.subfigs:
            subfig.clear()

        U = U_slider.val
        s = s_slider.val
        d = d_slider.val
        bias = b_slider.val

        if abs(bias) == 0:
            bias = float(0)
        t = t_slider.val

        f = lambda x: sorted(x, key= lambda y: float(y.split('=')[-1]) )
        cs = lambda ref : current(suppress_label=True, func=f, drawlines=[], reference=ref, avg=True, xuplim=t)
        cd = lambda ref: current( which_lead=1, func=f, reference=ref, avg=True, xuplim=t)
        socc = lambda x, ref : occ(suppress_label=True, width = 0.5,  func=f, single_site=x, reference=ref, xuplim=t)


        fstr = f'/Users/knl20/Desktop/Code/ED/scan/SD3x3-single_{Nup}Up{Ndn}DnU{float(U)}_1emultiplier{float(bias)}_sc{s}_dc{d}'
        print(fstr)
        #transport = [[cs(None), *[socc(i+j*3, None) for i in range(3)], cd(None)] for j in range(3)]
        transportone = [[cs(None), occ(width =1.0,  func=f, xuplim=t), cd(None)]]
        ed = [
        {
            'paneltitle' : f'3x3, L = {1}, N = {1}, scoupling = {s}, dcoupling = {d}, U = {U}',
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
    d_slider.on_changed(update)
    U_slider.on_changed(update)
    b_slider.on_changed(update)
    t_slider.on_changed(update)
    

    #fig.savefig( "/Users/knl20/Desktop/Code/non-interacting/plots/slider.pdf")
    plt.show()

# def test():

#     fig = plt.figure(
#                     layout="constrained", 
#                      figsize = (20, 20)
#                      )

#     row = 3
#     subcol = 3
#     subrow = 3
#     col = 1
#     subfigs = fig.subfigures(nrows=row, ncols=col)
#     axes = {}

#     overall = ['A', 'B', 'C']
#     names = [ ['a', *[str(i + j * subcol) for i in range(subcol)], 'b', 'c'] for j in range(subrow)]

#     for i in range(row):

#         subfig : plt.Figure = subfigs[i]

#         mosaic = [[ overall[i] + name for name in namerow] for namerow in names]
#         width_ratios = [ 1.0, 0.5, 0.5, 0.5, 1.0, 1.0]
#         height_ratios = [ 1.0, 1.0, 1.0]

#         axes |= subfig.subplot_mosaic(mosaic,  width_ratios=width_ratios, height_ratios=height_ratios)

#         for ax_name in axes:

#             axes[ax_name].set_title(ax_name)


#     fig.savefig('test.png')
    
