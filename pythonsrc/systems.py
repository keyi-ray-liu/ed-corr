import numpy as np
from numpy.linalg import eigh
import glob
import os
import sys
import matplotlib.pyplot as plt

sys.path.append(os.path.abspath("/Users/knl20/Desktop/Code/TN"))
from occdirect import *


class SYS:

    def __init__(self) -> None:

        self.name = None
        self.L = None
    
    def set_hij(self):
        raise NotImplementedError
    
    def set_init(self):
        raise NotImplementedError


    def get_name(self):
        return self.name
    
    
    
    def u(self, k, m):

        L = self.L
        return np.sqrt(2/ (L + 1)) * np.sin(  k * m * np.pi/ (L + 1))
    
    
    
class LSR_SIAM(SYS):

    def __init__(self, L, left_bias, right_bias, v, sys_onsite) -> None:

        self.L = L
        self.left_bias = left_bias
        self.right_bias = right_bias
        self.v = v
        self.sys_onsite = sys_onsite
        self.total = L + 1 + L
        self.name = 'LSR_SIAM'




    def set_init(self):

        # cal energy
        w, v = self.cal_energy()

        # select basis vec
        gs = v[:, 0].astype(complex)

        return np.diag(gs)


    
    # def set_cmn(self):
    #     '''set initial Cmn. For LSR-SIAM, the initial state was prepared s.t. it's in the unbiased GS. Since H is setup simply according to k's, for the initial state of the Cmn matrix, the smaller k's are occupied. '''
        

    #     L = self.L
    #     total = self.total

    #     diag = np.zeros(total, dtype=complex)

    #     diag[:L] = 1.0

    #     #print(np.count_nonzero(diag))

    #     init_mat = np.diag(diag)
    #     return init_mat

    def set_hij(self):


        L = self.L
        left = self.left_bias
        right = self.right_bias
        inter = self.v
        on_site = self.sys_onsite


        total = L + 1 + L
        M = np.zeros((total, total))

        # left
        for i in range(L - 1):

            M[i + 1, i] = 1.0
            M[i, i + 1] = 1.0

        # left to SYS

        M[ L - 1, L] = inter
        M[ L, L - 1] = inter

        # SYS to right
        M[ L, L + 1] = inter
        M[ L + 1, L] = inter

        # right
        for i in range(L + 1, L + L):
            M[i + 1, i] = 1.0
            M[i, i + 1] = 1.0

        # on-site, left
        for i in range(L):
            M[i, i] = left

        # on-site, SYS
        M[ L, L] = on_site

        # on-site, left
        for i in range(L + 1, L + 1 + L):
            M[i, i] = right


        return M
    


class DPT(SYS):

    def __init__(self, L, left, right, U, t=1.0) -> None:
        super().__init__()

        self.L = L
        self.U = U
        self.t = t
        self.left_bias = left
        self.right_bias = right

    def set_hij(self):

        U = self.U
        L = self.L
        left = self.left_bias
        right = self.right_bias

        total = L + L
        M = np.zeros((total, total))

        # all
        for i in range(total - 1):

            M[i + 1, i] = 1.0
            M[i, i + 1] = 1.0


        # on-site, left
        for i in range(L):
            M[i, i] = left

        # on-site, left
        for i in range(L, total):
            M[i, i] = right


        for i in range(L - 2, L):
            M[i, i] += U/2

        for i in range(L, L + 2):
            M[i, i] +=  U/2

        return M
    

class SD(SYS):

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
                
        

    def plot_current(self, currents, path):
        
        arr_source = self.arr_source
        arr_drain = self.arr_drain
        arrsize = self.arr
        fig, axes = plt.subplots(1, 2, figsize=(8, 2))

        ax : plt.Axes = axes[0]
        for i, arr in enumerate(arr_source):
            ax.plot(  currents[i], label=str(arr- self.L + 1))
        ax.legend()

        ax : plt.Axes = axes[1]
        for i, arr in enumerate(arr_drain):
            ax.plot( currents[i + arrsize], label=str(arrsize - self.L + 1))

        ax.legend()
        fig.tight_layout()
        fig.savefig('plots/{}.pdf'.format(path))


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
    

    def plot_occupation(self, path, *comp):

        
        occ_direct_general(files = [path, *comp])


class Reservoir(SYS):

    def __init__(self, L, N) -> None:
        super().__init__()
        self.L = L
        self.N = N
        self.t = 1

    def set_init(self):

        N = self.N
        total = self.L
        init = np.zeros(total, dtype=complex)
        init[:N] = 1.0#/np.sqrt(N)

        return init
    
    def set_hij(self):

        M = np.zeros((self.L, self.L))

        for i in range(self.L - 1):
            M[i, i+1] = M[i+1, i] = self.t


        return M

    def test_transform(self):

        M = self.set_hij()
        w, v = eigh(M)

        U = np.array([[self.u(k + 1, m + 1) for k in range(self.L)] for m in range(self.L)])

        # print(np.allclose(v, U))
        # print(np.allclose(*map(np.abs, (v, U))))
        # init = self.set_init()

        # print(np.sum(np.absolute(v@init) ** 2))
        # plt.plot( v@ init)
        # plt.show()
