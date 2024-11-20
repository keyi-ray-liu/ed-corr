import numpy as np
from itertools import product
from joblib import Parallel, delayed
import time 
import os
from numpy.linalg import eigh



class BiasSolver:

    def __init__(self) -> None:
        pass


    def solve(self, h0, h, timestep, fin, N, fullCC =False, factor=1.0):

        _, U0 = eigh(h0)
        w, U = eigh(h)


        total = w.shape[0]
        E = np.zeros(( total, total), dtype=complex)

        for i in range(total):
            for j in range(total):
                E[i, j] = np.exp( -1j * (w[i] - w[j]) * timestep)

        #E = np.transpose(E)

        cur_E = E
        U0 = U0[:, :N]
        Ud = np.conj(np.transpose(U))
        U0d = np.conj(np.transpose(U0))



        #print("Ud, U, U0d, U0 shapes:", Ud.shape, U.shape, U0d.shape, U0.shape)
        steps = int(fin/timestep)

        if fullCC:
            CCs = np.zeros( [ steps + 1] + list(U.shape), dtype=complex)
            CCs[0] = np.transpose(Ud) @ ( np.multiply(Ud @ U0 @ U0d @ U ,np.ones(cur_E.shape))) @ np.transpose(U)

        else:
            current = np.zeros(steps + 1)

        for i in range(steps):

            #print(i)

            if fullCC:
                CC = np.transpose(Ud) @ ( np.multiply(Ud @ U0 @ U0d @ U ,cur_E)) @ np.transpose(U)
                CCs[i + 1] = CC

            else:
                cc = Ud[:, N - 1] @ ( np.multiply(Ud @ U0 @ U0d @ U ,cur_E)) @ U[N, :]
                current[i + 1] =  2 * np.imag(cc) * factor

            cur_E = np.multiply(cur_E, E)

            #print(cur_E)

        if fullCC:
            return CCs, np.arange(CCs.shape[0])*timestep
        else:

            return current, np.arange(current.shape[0])*timestep



class SD():

    def __init__(self, L, N, arr, left_bias, arr_bias, right_bias, source_coupling, drain_coupling, single = False, t=-1.0) -> None:
        super().__init__()

        self.L = L
        self.name = 'L{}_N{}_S{}_D{}'.format(L, N, abs(source_coupling), abs(drain_coupling))
        self.source_coupling = source_coupling
        self.drain_coupling = drain_coupling
        self.t = t
        self.N = N
        self.left_bias = left_bias
        self.arr_bias = arr_bias
        self.right_bias = right_bias
        self.arr = arr
        self.contact_source = L - 1
        self.contact_drain = L + arr ** 2 

        if not single:
            self.arr_source = self.contact_source + np.arange( 1, arr **2 + 1, arr, dtype=int)
            self.source_scaling = np.ones(arr)
            self.arr_drain = self.contact_source + np.arange( arr, arr **2 + 1, arr, dtype=int)
            self.drain_scaling = np.ones(arr)

            if arr == 3:
                self.source_scaling[1] = 2
                self.drain_scaling[1] = 2

        else:
            self.arr_source = np.array([self.contact_source + 1 + arr * ((arr - 1)//2)])
            self.source_scaling = np.ones(1)
            self.arr_drain = np.array([self.contact_drain - 1 - arr * ((arr - 1)//2)])
            self.drain_scaling = np.ones(1)

        #print(self.arr_source, self.arr_drain)

    # 0 to L - 1 : left
    # L to L + 8: arr

    # 0 1 2
    # 3 4 5
    # 6 7 8

    def res_size(self):
        return self.L

    # L + 9 to end: right
    def set_hij(self, path, mod, savehij = False):

        t = self.t
        L = self.L
        arr = self.arr
        arr_size = arr ** 2
        left = self.left_bias
        cent = self.arr_bias
        right = self.right_bias
        source_coupling = self.source_coupling
        drain_coupling = self.drain_coupling
        contact_source = self.contact_source
        contact_drain = self.contact_drain
        arr_source = self.arr_source
        arr_drain = self.arr_drain
        source_scaling = self.source_scaling
        drain_scaling = self.drain_scaling

        total = L + L + arr_size
        M = np.zeros((total, total))

        # Left
        for i in range(L - 1):

            M[i + 1, i] = M[i, i + 1] = t
            M[i, i] = left
        
        M[L -1, L -1] = left

        # coupling

        for i, site in enumerate(arr_source):
            M[contact_source, site] = M[site, contact_source] = source_coupling * t * source_scaling[i]

        for i, site in enumerate(arr_drain):

            M[contact_drain, site] = M[site, contact_drain] = drain_coupling * t * drain_scaling[i]  

        # on-site, right
        for i in range(L + arr_size, total - 1):

            M[i + 1, i] = M[i, i+ 1] = t
            M[i, i] = right

        M[total - 1, total -1] = right

        # arr
        for i in range(arr):
            for j in range(arr):
                

                idx = i * arr + j + L 
                M[idx, idx] += cent
                if i < arr - 1:
                    M[idx, idx + arr] = t
                    M[idx + arr, idx] = t

                if j < arr - 1:
                    M[idx, idx + 1] = t
                    M[idx + 1, idx] = t

        if not np.allclose(M, np.conj(np.transpose(M))):
            raise(ValueError("M not hermitian!"))
        
        try:
            os.mkdir(path)
        except:
            pass
        
        if savehij:
            np.savetxt(path + '/M' + mod + '.txt', M, fmt='%6.4f')

        return M
    
    def set_init(self):

        N = self.N
        total = self.L * 2 + self.arr ** 2
        init = np.zeros(total, dtype=complex)
        init[:N] = 1.0#/np.sqrt(N)

        return init

    def current(self, states, fullCC=False):

        contact_source = self.contact_source
        contact_drain = self.contact_drain
        arrsize = self.arr
        
        arr_source = self.arr_source
        arr_drain = self.arr_drain

        #print(arr_source, arr_drain, contact_source, contact_drain)

        #print(contact_source, contact_drain, arr_source, arr_drain)
        currents = np.zeros((arrsize * 2, states.shape[0]))

        for i in range(1, states.shape[0]):

            for j, arr in enumerate(arr_source):
                
                if not fullCC:
                    cc = np.conj( states[i, contact_source]) * states[i, arr]
                    currents[j, i] = -2 * np.imag(cc) #* factor
                else:
                    currents[j, i] = -2 * np.imag(states[i, contact_source, arr]) 
                

            for k, arr in enumerate(arr_drain):
                
                if not fullCC:
                    cc = np.conj( states[i, arr]) * states[i, contact_drain]
                    currents[k + len(arr_source), i] = -2 * np.imag(cc) #* factor

                else:
                    currents[k + len(arr_source), i ] = -2 * np.imag(states[i, arr, contact_drain]) 


        return currents
                
        
    def save_result(self, currents, occs, CCs, timestep, fin, path, saveCC=False):

        #print('saving')
    
        try:
            os.mkdir(path)
        except:
            pass
        
        if saveCC:
            np.save(path + '/CC', CCs)
        np.savetxt( path + '/times', np.arange(0, fin + timestep, timestep))
        np.savetxt( path + '/currentSD', currents)
        np.savetxt( path + '/occ', occs)
    




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

def parallelrun(kwargs):
    return SDbias(**kwargs)

def SDscan(cpu=1):


    # Get the number of CPUs
    num_cpus = os.cpu_count()

    print("Number of CPUs available:", num_cpus)

    
    Ls = np.power(2, np.arange(7, 11), dtype=int)
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
