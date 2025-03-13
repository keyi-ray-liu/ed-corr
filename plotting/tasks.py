from SDcurrents import currents
from glob import glob
from SDfermi import compare_fermi
from SDgeneric import *
import matplotlib.pyplot as plt

from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm



def get_biases(U):

    biases = {
        2.0 : 0.5761793031432951,
        50.0 : 1.8257251433266721,
        100.0 : 1.7916179881169065
    }

    return biases[U]

def biase_fraction(cur_bias, U):

    ref_bias = get_biases(U)
    return str(np.round( cur_bias /ref_bias, decimals=5)) 
        



def plotspacing():

    save_individual = False

    sysstr = '/Users/knl20/Desktop/Code/TN/SD/Feb_18_spacing/'
    file = 'spacingNF.pdf'
    exampledim = 128
    exampleconfig = '3x3-single'
    examplecoupling = 1/32
    spacings = np.power(2.0, -np.arange( 9))
    Us = [2.0, 100.0]
    Nstrings = ['3Up4Dn', '4Up4Dn']
    exampleL = 128
    coul = 0.0

    mode_str = { 'bias' : 'Drive mode: non-zero bias across source-drain',
                'left': 'Drive mode: source at half-filling (both up and down electron), then released'}
    

    fermifuncrelease = [EE(which_reservoir='source'), EE(which_reservoir='drain'), Empty(), Fermi(which_reservoir='source'), Fermi(which_reservoir='drain'), Empty(id=1), Occlandscape(which_reservoir='source'), Occlandscape(which_reservoir='drain')]


    f = lambda x: sorted(x, key= lambda y: float(y.split('=')[-1]) )

    illu = Illustrate(drive = 'left')
    cs = lambda ref : current(suppress_label=True, func=f, drawlines=[], reference=ref, avg=True)
    cd = lambda ref: current(which_lead=1, func=f, reference=ref, avg=True)
    socc = lambda x, ref : occ(suppress_label=True, width = 0.5, func=f, single_site=x, reference=ref)

    #ref = lambda L, N, arr, s, spacing: f'{refstr}/L{L}_N{N}_arr{arr}_arrhop1.0_biasA0.0_biasAinit1000.0_biasSD0.25_d{s}_s{s}_singletrue_spacing{spacing}'



    paneltransition = [
        [
            
            {
                'paneltitle' : f'3x3, {mode_str['left']}, coupling = {convert_to_frac(examplecoupling)}, U = {convert_to_frac(U)}',
                'panelfuncs' : [[cs(None), *[socc(i+j*3, None) for i in range(3)], cd(None)] for j in range(3)],
                'panelparameters' : {
                        'dirs': [
                            sysstr + f'Ls{exampleL}_Ns{exampleL//2}_Na{Nstring}_Nd0_scoupling-{examplecoupling}_systypeElectron_coul{coul}_TEdim{exampledim}_reservoirtypemixed_U{U}_config{exampleconfig}_biasAinit0.0_biasA0.0_biasS-3.0_biasD3.0_modeBiasRelease_reservoirspacing{spacing}' for spacing in spacings
                        ],
                        'linestyles': ['solid' for _ in range(len(spacings))],
                        'labels' : [f'spacing = {spacing}' for spacing in spacings],
                        'colors' : [ cm.hot(i/len(spacings)) for i in range(len(spacings))]
                }
                
            } for U in Us
        ] for Nstring in Nstrings
    ]
    
    panelcompare = [
        [
            
            {
                'paneltitle' : f'3x3, {mode_str['left']}, coupling = {convert_to_frac(examplecoupling)}, U = {U}',
                'panelfuncs' : [[ Occcompare()]],
                'panelparameters' : {
                        'dirs': [
                            sysstr + f'Ls{exampleL}_Ns{exampleL//2}_Na{Nstring}_Nd0_scoupling-{examplecoupling}_systypeElectron_coul{coul}_TEdim{exampledim}_reservoirtypemixed_U{U}_config{exampleconfig}_biasAinit0.0_biasA0.0_biasS-3.0_biasD3.0_modeBiasRelease_reservoirspacing{spacing}' for spacing in spacings
                        ],
                        'ref' : f'/Users/knl20/Desktop/Code/TN/SD/TNreference/NF{Nstring}3x3U{U}',
                        'linestyles': ['solid' for _ in range(len(spacings))],
                        'labels' : [f'spacing = {spacing}' for spacing in spacings],
                        'colors' : [ cm.hot(i/len(spacings)) for i in range(len(spacings))]
                }
                
            } for U in Us
        ] for Nstring in Nstrings
    ]

    exampleL = 128
    paneltransitionfermi = [
        [
            [
                {
                    'paneltitle' : f'3x3, L =128, {mode_str['left']}, coupling = {convert_to_frac(examplecoupling)}, spacing = {convert_to_frac(spacing)}',
                    'panelfuncs' : [fermifuncrelease],
                    'panelparameters' : {
                            'dirs': [
                                sysstr + f'Ls{exampleL}_Ns{exampleL//2}_Na{Nstring}_Nd0_scoupling-{examplecoupling}_systypeElectron_coul{coul}_TEdim{exampledim}_reservoirtypemixed_U{U}_config{exampleconfig}_biasAinit0.0_biasA0.0_biasS-3.0_biasD3.0_modeBiasRelease_reservoirspacing{spacing}'
                            ]
                    }
                    
                } for spacing in spacings
            ] for U in Us ]
        for Nstring in Nstrings
    ]

    paneldmrgconv = [
            
                {
                    'paneltitle' : 'DMRG convergence',
                    'panelfuncs' : [[DMRGConvergence()]],
                    'panelparameters' : {
                            'cases': [
                                f'{Nstring}_U{U}' for Nstring in Nstrings for U in Us
                            ],
                            'colors' : [ cm.hot(i/4) for i in range(4)]
                    }
                    
                } 
            
    ]

    with PdfPages(f'/Users/knl20/Desktop/Code/TN/SD/plots/{file}') as pdf:
        
        #generic(pdf, paneltest)
        generic(pdf, paneldmrgconv, width = 10)
        for i, Nstring in enumerate(Nstrings):
            generic(pdf, paneltransition[i], save_individual=False, filename='GS_noadjust', figsuptitle=Nstring)
            generic(pdf, panelcompare[i], width =15, save_individual=False, filename='GS_noadjust', figsuptitle=Nstring)

        
        # for i, Nstring in enumerate(Nstrings):
        #     for j, U in enumerate(Us):
        #         generic(pdf, paneltransitionfermi[i][j], save_individual=False, filename='GS_noadjust', figsuptitle=Nstring + ',  ' + str(U))

        

     

def plottransport1():

    save_individual = False

    sysstr = '/Users/knl20/Desktop/Code/TN/SD/Feb_26_transport2_largerpotential/'
    
    exampledim = 128
    exampleconfig = '3x3-single'

    examplespacing = 1/16
    Us = [2.0, 50.0, 100.0]
    Nstring = '3Up4Dn'
    exampleL = 128
    exampleNs = 64
    coul = 0.0

    

    mode_str = { 'bias' : 'Drive mode: non-zero bias across source-drain',
                'left': 'Drive mode: source at half-filling (both up and down electron), then released'}
    

    fermifuncrelease = [EE(which_reservoir='source'), EE(which_reservoir='drain'), Empty(), Fermi(which_reservoir='source'), Fermi(which_reservoir='drain'), Empty(id=1), Occlandscape(which_reservoir='source'), Occlandscape(which_reservoir='drain')]


    f = lambda x: sorted(x, key= lambda y: float(y.split('=')[-1]) )
    #f = sorted

    illu = Illustrate(drive = 'left')
    cs = lambda ref : current(suppress_label=True, func=f, drawlines=[], reference=ref, avg=True)
    cd = lambda ref: current( which_lead=1, func=f, reference=ref, avg=True)
    socc = lambda x, ref : occ(suppress_label=True, width = 0.5,  func=f, single_site=x, reference=ref)
    effe = EffectiveEntanglement(start = 0, end = exampleL * 2 + 3 * 3, func=f)

    transport = [[cs(None), *[socc(i+j*3, None) for i in range(3)], cd(None)] for j in range(3)]

    totalmeasure = [[effe, TotalCharge(func=f) ]]
    #ref = lambda L, N, arr, s, spacing: f'{refstr}/L{L}_N{N}_arr{arr}_arrhop1.0_biasA0.0_biasAinit1000.0_biasSD0.25_d{s}_s{s}_singletrue_spacing{spacing}'

    biases = { U : sorted(set([float(searchkey('biasA', f)) for f in glob( sysstr + f'*Ns{exampleNs}*_U{U}_*' )])) for U in Us}

    paneltransport = lambda examplecoupling : [
        
            
            {
                'paneltitle' : f'3x3, L = {exampleL}, N = {exampleNs}, coupling = {convert_to_frac(examplecoupling)}, U = {convert_to_frac(U)}',
                'panelfuncs' : transport,
                'panelparameters' : {
                        'dirs': [
                            sysstr + f'Ls{exampleL}_Ns{exampleNs}_Na{Nstring}_Nd0_scoupling-{examplecoupling}_systypeElectron_coul{coul}_TEdim{exampledim}_reservoirtypemixed_U{U}_config{exampleconfig}_biasAinit0.0_biasA{bias}_biasS-3.0_biasD3.0_modeBiasRelease_reservoirspacing{examplespacing}' for bias in biases[U]
                        ],
                        'linestyles': ['solid' for _ in range(len(biases[U]))],
                        'labels' : [r'bias ($\times \Delta_{{e^-}})$ =' + f'{biase_fraction(bias, U)}' for bias in biases[U]],
                        'colors' : [ cm.hot(i/len(biases[U])) for i in range(len(biases[U]))]
                }
                
            } for U in Us
        
    ]

    paneltotal = lambda examplecoupling : [
        
            
            {
                'paneltitle' :  f'3x3, L = {exampleL}, N = {exampleNs}, coupling = {convert_to_frac(examplecoupling)}, U = {convert_to_frac(U)}',
                'panelfuncs' : totalmeasure,
                'panelparameters' : {
                        'dirs': [
                            sysstr + f'Ls{exampleL}_Ns{exampleNs}_Na{Nstring}_Nd0_scoupling-{examplecoupling}_systypeElectron_coul{coul}_TEdim{exampledim}_reservoirtypemixed_U{U}_config{exampleconfig}_biasAinit0.0_biasA{bias}_biasS-3.0_biasD3.0_modeBiasRelease_reservoirspacing{examplespacing}' for bias in biases[U]
                        ],
                        'linestyles': ['solid' for _ in range(len(biases[U]))],
                        'labels' : [r'bias ($\times \Delta_{{e^-}})$ =' + f'{biase_fraction(bias, U)}' for bias in biases[U]],
                        'colors' : [ cm.hot(i/len(biases[U])) for i in range(len(biases[U]))]
                }
                
            } for U in Us
        
    ]




    
    panelcompare = lambda examplecoupling : [
        
            
        {
            'paneltitle' : f'3x3, {mode_str['left']}, coupling = {convert_to_frac(examplecoupling)}, U = {U}',
            'panelfuncs' : [[ Occcompare()]],
            'panelparameters' : {
                    'dirs': [
                        sysstr + f'Ls{exampleL}_Ns{exampleNs}_Na{Nstring}_Nd0_scoupling-{examplecoupling}_systypeElectron_coul{coul}_TEdim{exampledim}_reservoirtypemixed_U{U}_config{exampleconfig}_biasAinit0.0_biasA{bias}_biasS-3.0_biasD3.0_modeBiasRelease_reservoirspacing{examplespacing}' for bias in biases[U]
                    ],
                    'ref' : f'/Users/knl20/Desktop/Code/ED/ref/3x3-single{Nstring}U{U}', #f'/Users/knl20/Desktop/Code/TN/SD/TNreference/NF{Nstring}3x3U{U}',
                    'linestyles': ['solid' for _ in range(len(biases[U]))],
                    'labels' : [r'bias ($\times \Delta_{{e^-}})$ =' +  f'{biase_fraction(bias, U)}' for bias in biases[U]],
                    'colors' : [ cm.hot(i/len(biases[U])) for i in range(len(biases[U]))]
            }
            
        } for U in Us

    ]



    paneldmrgconv = lambda examplecoupling : [
            
                {
                    'paneltitle' : f'DMRG convergence, U = {U}',
                    'panelfuncs' : [[DMRGConvergence()]],
                    'panelparameters' : {
                            'dirs': [
                            sysstr + f'Ls{exampleL}_Ns{exampleNs}_Na{Nstring}_Nd0_scoupling-{examplecoupling}_systypeElectron_coul{coul}_TEdim{exampledim}_reservoirtypemixed_U{U}_config{exampleconfig}_biasAinit0.0_biasA{bias}_biasS-3.0_biasD3.0_modeBiasRelease_reservoirspacing{examplespacing}' for bias in biases[U]
                        ],
                            'labels' : [  r'bias ($\times \Delta_{{e^-}})$ = ' + biase_fraction(bias, U) for bias in biases[U]],
                            'colors' : [ cm.hot(i/len(biases[U])) for i in range(len(biases[U]))]
                    }
                    
                } for U in Us
            
    ]
    # paneltransitionfermi = [
    #     [
    #         [
    #             {
    #                 'paneltitle' : f'3x3, L =128, {mode_str['left']}, coupling = {convert_to_frac(examplecoupling)}, spacing = {convert_to_frac(spacing)}',
    #                 'panelfuncs' : [fermifuncrelease],
    #                 'panelparameters' : {
    #                         'dirs': [
    #                             sysstr + f'Ls{exampleL}_Ns{exampleL//2}_Na{Nstring}_Nd0_scoupling-{examplecoupling}_systypeElectron_coul{coul}_TEdim{exampledim}_reservoirtypemixed_U{U}_config{exampleconfig}_biasAinit0.0_biasA0.0_biasS-3.0_biasD3.0_modeBiasRelease_reservoirspacing{spacing}'
    #                         ]
    #                 }
                    
    #             } for spacing in spacings
    #         ] for U in Us ]
    #     for Nstring in Nstrings
    # ]

    
    file = f'transport1Ls{exampleL}Ns{exampleNs}largepotential.pdf'

    with PdfPages(f'/Users/knl20/Desktop/Code/TN/SD/plots/{file}') as pdf:
        
        generic(pdf, paneldmrgconv(1/8), width=15, figsuptitle=Nstring)
        generic(pdf, panelcompare(1/8), width = 15, save_individual=False, filename='GS_noadjust', figsuptitle=Nstring)
        #generic(pdf, paneltransport(1/2), length=10, save_individual=False, filename='GS_noadjust', figsuptitle=Nstring)
        #generic(pdf, paneltotal(1/2), width = 10, save_individual=False, filename='GS_noadjust', figsuptitle=Nstring)
        generic(pdf, paneltransport(1/8), length=10, save_individual=False, filename='GS_noadjust', figsuptitle=Nstring)
        generic(pdf, paneltotal(1/8), width = 10, save_individual=False, filename='GS_noadjust', figsuptitle=Nstring)




        
        # for i, Nstring in enumerate(Nstrings):
        #     for j, U in enumerate(Us):
        #         generic(pdf, paneltransitionfermi[i][j], save_individual=False, filename='GS_noadjust', figsuptitle=Nstring + ',  ' + str(U))

        

def plottest():

    f = lambda x: sorted(x, key= lambda y: float(y.split('=')[-1]) )
    cs = lambda ref : current(suppress_label=True, func=f, drawlines=[], reference=ref, avg=True)
    cd = lambda ref: current( which_lead=1, func=f, reference=ref, avg=True)
    socc = lambda x, ref : occ(suppress_label=True, width = 0.5,  func=f, single_site=x, reference=ref)

    transport = [[cs(None), *[socc(i+j*3, None) for i in range(3)], cd(None)] for j in range(3)]


    refdir = '/Users/knl20/Desktop/Code/TN/SD/TNreference/reversetest/'
    forward_backward = [
        {
            'paneltitle' : f'3x3, L = {128}, N = {64}, coupling = {convert_to_frac(0.125)}, U = {2.0}',
            'panelfuncs' : transport,
            'panelparameters' : {
                    'dirs': [refdir + 'ref_lohigh', refdir + 'ref_highlo'],
                    'linestyles': ['solid', 'dashed'],
                    'labels' : ['sort: low high', 'sort: high low'],
                    'markers' : ['', 'o'],
                    'colors' : ['red', 'blue']
            }
        } 
    ]


    docc = lambda x, ref : occ(suppress_label=False, width = 0.5,  func=f, single_site=x, reference=ref)

    edempty = [
        {
            'paneltitle' : f'3x3, L = {0}, N = {0}, coupling = {convert_to_frac(0.125)}, U = {4.0}',
            'panelfuncs' : transport,
            'panelparameters' : {
                    'dirs': ['/Users/knl20/Desktop/Code/ED/ref/SD3x3-single1Up1DnU4.0bias0.0', '/Users/knl20/Desktop/Code/TN/SD/TNreference/NF1Up1Dn3x3U4.0/work' ],
                    'reservoirsizes' : [1, 1],
                    'linestyles': ['solid', 'dashed'],
                    'markers' : ['', 'o'],
                    'labels' : ['ED', 'MPS'],
                    'colors' : ['red', 'blue']
            }
        } 
    ]



    ed = [
        {
            'paneltitle' : f'3x3, L = {1}, N = {1}, coupling = {convert_to_frac(0.125)}, U = {4.0}',
            'panelfuncs' : transport,
            'panelparameters' : {
                    'dirs': ['/Users/knl20/Desktop/Code/ED/ref/SD3x3-single2Up2DnU4.0bias0.0', '/Users/knl20/Desktop/Code/TN/SD/TNreference/NF2Up2DnU4.0ED/work' ],
                    'reservoirsizes' : [1, 1],
                    'linestyles': ['solid', 'dashed'],
                    'markers' : ['', 'o'],
                    'labels' : ['ED', 'MPS'],
                    'colors' : ['red', 'blue']
            }
        } 
    ]

    ed_transport = [
        {
            'paneltitle' : f'3x3, L = {1}, N = {1}, coupling = {convert_to_frac(0.125)}, U = {100.0}, bias = 1e',
            'panelfuncs' : [[ *[docc(i+j*3, None) for i in range(3)]] for j in range(3)],
            'panelparameters' : {
                    'dirs': ['/Users/knl20/Desktop/Code/ED/ref/SD3x3-single4Up5DnU100.0bias-1.79161798', '/Users/knl20/Desktop/Code/TN/SD/TNreference/NF3Up4Dn3x3U100.0bias/work' ],
                    'reservoirsizes' : [1, 1],
                    'linestyles': ['solid', 'dashed'],
                    'markers' : ['', 'o'],
                    'labels' : ['ED', 'MPS'],
                    'colors' : ['red', 'blue']
            }
        } 
    ]

    portera = [
        {
            'paneltitle' : f'ERA, 3x3, L = {1}, N = {1}, coupling = {convert_to_frac(0.125)}, U = {0.0}, $\gamma$= {0.1}',
            'panelfuncs' : transport,
            'panelparameters' : {
                    'dirs': ['/Users/knl20/Desktop/github/portera/work/SDtest' ],
                    'reservoirsizes' : [1],
                    'linestyles': ['solid', 'dashed'],
                    'markers' : ['', 'o'],
                    'labels' : ['ED', 'MPS'],
                    'colors' : ['red', 'blue']
            }
        } 
    ]

    with PdfPages(f'/Users/knl20/Desktop/Code/TN/SD/plots/testMar3.pdf') as pdf:
        
        generic(pdf, forward_backward, width=5, figsuptitle='3Up4Dn')
        generic(pdf, edempty, width=5, figsuptitle='0Up0Dn')
        generic(pdf, ed, width=5, figsuptitle='1Up1Dn')
        generic(pdf, ed_transport, width=5, figsuptitle='3Up4Dn')
        generic(pdf, portera, width=5, figsuptitle='2Up2Dn')