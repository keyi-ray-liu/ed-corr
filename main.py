import numpy as np
from matplotlib import rc
from solver import *
from scipy.linalg import expm
from systems import *
from solver import *
from utils import *
from itertools import product
import multiprocessing
from joblib import Parallel, delayed
import os
import time

# def set_U(**kwargs):

#     N = kwargs["N"]

#     U = np.zeros((N, N))

#     for j in range(N):
#         for k in range(N):
#             U[j ,k] = u(j +1, k + 1, N)

#     return U


def LSR():
    key = 'LSR'

    timestep = 0.25
    fin = 200
    left=0.25
    right=-0.25

    inter = 1/np.sqrt(2)
    on_site = 1.0
    factor = 4*np.pi * inter
    N = 128
    solver = BiasSolver()
    
    system = LSR_SIAM(N, left, right, inter, on_site)
    INIT = LSR_SIAM(N, 0.0, 0.0, inter, on_site)

    fig, ax = plt.subplots( figsize=(6, 5))

    h0 = INIT.set_hij()
    h = system.set_hij()

    current, _ = solver.solve(h0, h, timestep, fin, N, factor=factor)
    ax : plt.Axes = ax
    ax.plot( np.arange(current.shape[0]) * timestep, current)

    title = "LSR Current with correlation matrix,  N={}, timestep={}".format(N, timestep)
    #desc = "N{}t{}.pdf".format(N, timestep)

    ax.set_xlabel('Time')
    ax.set_ylabel(' $ I$')
    ax.set_title(title)

    fig.tight_layout()
    fig.savefig('plots/' + 'LSR.pdf')

def cal_DPT():

    key = 'DPT'

    timestep = 0.25
    fin = 128
    factor = 1.0
    Us = np.arange(0, 4.125, 0.125)

    solver = BiasSolver()
    #Ns = [34, 66]
    Ns = [130]

    fig, axes = plt.subplots(1, len(Us), figsize=(6 * len(Us), 5))
    for i, U in enumerate(Us):

        for j, N in enumerate(Ns):

        # sys = LSR_SIAM(N, left, right, inter, on_site)
        # INIT = LSR_SIAM(N, 0.0, 0.0, inter, on_site)
            for mu  in [#0.0,
                         0.5, 2.0
                        ]:

                left = mu/2
                right = -mu/2

                system = DPT(N, left, right, U)
                INIT = DPT(N, 0.0, 0.0, U)

                h0 = INIT.set_hij()
                h = system.set_hij()

                current, time = solver.solve(h0, h, timestep, fin, N, factor=factor)

                np.savetxt( str(N) + 'currentCC' + str(U) + str(mu), current)
                np.savetxt( str(N) + 'timeCC' + str(U) + str(mu), time)

                # ax : plt.Axes = axes[i]
                # ax.plot( np.arange(current.shape[0]) * timestep, current)
            
                # title = "DPT Current with correlation matrix, U={}, N={}, timestep={}".format(U, N, timestep)
                # #desc = "N{}t{}.pdf".format(N, timestep)

                # ax.set_xlabel('Time')
                # ax.set_ylabel(' $ I$')
                # ax.set_title(title)

    fig.tight_layout()
    fig.savefig('plots/' + 'DPT.pdf')

def SDprod():
    key = 'SDprod'

    timestep = 0.2
    fin = 200
    L = 64
    N = L//2
    sd = SD(L, N, 0, 0, -0.1)
    #N = L //2


    init = sd.set_init()
    solver = SDSolver()
    #solver = ExponentiateSolver()
    currents, occs = solver.solve(sd, init, timestep, fin)

    sd.save_result(currents, occs, [], timestep, fin)
    #sd.plot_current(currents, timestep)

    comp1 = "/Users/knl20/Desktop/Code/TN/SD/64_32_3x3_prod_spin_dim64_nocoul_tun0.1"
    comp2 = "/Users/knl20/Desktop/Code/TN/SD/64_32_3x3_prod_spin_dim64_nocoul_tun0.1_highdim"
    sd.plot_occupation( comp1, comp2)


def SDGSnoint():
    key = 'SDGSnoint'   
    L, N = 1, 1


    

    combs = [
        (1, 1),
        [1, 0.01],
        [0.1, 0.1],
        [0.1, 0.001],
        [0.01, 0.01],
        [0.01, 0.0001]
    ]

    for s, d in combs:

        factor = s
        timestep = 0.01/factor
        fin = 10/factor

        name = 'L{}_N{}_S{}_D{}'.format(L, N, abs(s), abs(d))
        path = name + '_' + key

        system = SD(L, N, 0, 0, 0, s, d)
        INIT = SD(L, N, -1000.0, 0, 0, s, d)

        h0 = INIT.set_hij(path, 'init')
        h = system.set_hij(path, 'full')


        solver = BiasSolver()
        #solver = ExponentiateSolver()
        CCs, _ = solver.solve(h0, h, timestep, fin, N, fullCC=True)


        currents = system.current(CCs, fullCC=True)
        system.save_result(currents, np.array([np.real(np.diag((CC))) for CC in CCs]), CCs, timestep, fin, path)
    #sd.plot_current(currents, timestep)

    # comp1 = "/Users/knl20/Desktop/Code/TN/SD/6432/6432_LGS_tonly_0.1_mixed"
    # comp2 = "/Users/knl20/Desktop/Code/TN/SD/6432/6432_LGS_tonly_0.1_spatial"
    # system.plot_occupation('GS', comp1, comp2)

def SDbias(L=128, N=128, fin=128.0, bias=0.25, s=0.1, d=0.1, arr=1, single=True, mu_init=1000, mu_te=0, prefix = '', cpu=1):

    begin = time.time()
    key = 'SDbias'
    
    timestep = 0.125
    name = 'L{}_N{}_{}x{}_muinit{}_mute{}_S{}_D{}_single{}'.format(L, N, arr, arr, mu_init, mu_te, abs(s), abs(d), single)
    path = prefix + name + '_' + key

    system = SD(L, N, arr, bias, mu_te, -bias, s, d, single=single)
    INIT = SD(L, N, arr, 0.0, mu_init, 0, s, d, single=single)

    h0 = INIT.set_hij(path, 'init')
    h = system.set_hij(path, 'full')


    solver = BiasSolver()
    #solver = ExponentiateSolver()
    CCs, _ = solver.solve(h0, h, timestep, fin, N, fullCC=True)


    currents = system.current(CCs, fullCC=True)
    # if arr == 1:
    #     system.plot_current(currents, path)

    occ = np.array([np.abs(np.diag((CC))[L:L + arr ** 2]) for CC in CCs])
    system.save_result(currents, occ, CCs, timestep, fin, path)

    avgtime = (time.time() - begin)/cpu

    print("avg time : {}".format(avgtime))

# comp1 = '/Users/knl20/Desktop/Code/TN/SD/LRbiasNoCoul/L64_N64_bias0.5_Spinless_0.1_0.1_mixed'
# comp2 = '/Users/knl20/Desktop/Code/TN/SD/LRbiasNoCoul/L64_N64_bias0.5_Spinless_0.1_0.1_spatial'
# system.plot_occupation(path, comp1, comp2)



def phase():

    
    N = 1000
    W = np.zeros((N, 11))
    grid = np.linspace(0, 5, N)
    for i, s in enumerate(grid):

        d = s
        system = SD(1, 1, 0, 0, 0, s, d)
        h = system.set_hij('temp', '')

        #hstr = to_mathematica(h, {-1: 't', -s : 's', 2 * -s : '2s'})
        #print(hstr)
        w, _ = eigh(h)
        W[i] = w

    fig, ax = plt.subplots(figsize=(8, 5))
    ax : plt.Axes

    for j in range(11):
        ax.scatter( grid, W[:, j], s=0.1)

    ax.set_xlabel('coupling (abs)')
    ax.set_ylabel('Energies')
    fig.savefig('plots/phase.pdf')



def parallelrun(kwargs):
    return SDbias(**kwargs)

def SDscan(cpu=1):


    # Get the number of CPUs
    num_cpus = os.cpu_count()

    print("Number of CPUs available:", num_cpus)

    
    Ls = np.power(2, np.arange(6, 11), dtype=int)
    biases = np.power(2.0, np.arange(-3, 4))
    ss = np.power(10.0, np.arange(-3, 2))
    ds = np.power(10.0, np.arange(-3, 2))
    arrs = [1, 2, 3, 4]
    mu_tes = np.concatenate( ([0], np.power(2.0, np.arange(-3, 5))) )

    combs = list(product(Ls, biases, ss, ds, arrs, mu_tes))
    single = True

    kwargs = [{
            "L" : L,
            "N" : L,
            "fin" : L * 0.9,
            "bias" : bias,
            "s" : s,
            "d" : d,
            "arr" : arr,
            "single" : single,
            "mu_init" : 1000,
            "mu_te" : mu_te,
            "prefix" : 'SDscan/',
            "cpu" : cpu
    } for L, bias, s, d, arr, mu_te in combs]

    # pool = multiprocessing.Pool(1)
    # pool.map(parallelrun, kwargs)

    Parallel(n_jobs=cpu)(delayed(parallelrun)(kwarg) for kwarg in kwargs)


def Res():
    key = 'Res'
    L = 64
    N = 32

    res = Reservoir(L, N)
    res.test_transform()

    
#reference(np.diag(hmn), energy)
#sol = solve(hmn, cmn, timestep=timestep, fin=fin)
#LSR.current(sol, momentum)


def test():
    Ls = [64]
    biases = [0.25]
    ss = [0.5]
    ds = [0.5]
    arrs = [1, 2, 3, 4]
    mu_tes = [0, 1]

    combs = list(product(Ls, biases, ss, ds, arrs, mu_tes))
    single = True

    kwargs = [{
            "L" : L,
            "N" : L,
            "fin" : L * 0.9,
            "bias" : bias,
            "s" : s,
            "d" : d,
            "arr" : arr,
            "single" : single,
            "mu_init" : 1000,
            "mu_te" : mu_te,
            "prefix" : 'SDscan/',
            "cpu" : 1
    } for L, bias, s, d, arr, mu_te in combs]

    # pool = multiprocessing.Pool(1)
    # pool.map(parallelrun, kwargs)

    Parallel(n_jobs=1)(delayed(parallelrun)(kwarg) for kwarg in kwargs)

if __name__ == '__main__':

    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True) 

    #SDbias()
    #SDscan(cpu=int(sys.argv[1]))
    test()
    #LSR()
    #SDGSnoint()
    #phase()
    #cal_DPT()

