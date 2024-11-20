import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.linalg import expm
from numpy.linalg import eigh
from systems import *
# def solve( hmn : np.ndarray, cmn: np.ndarray, **kwargs):


#     def gen_RHS(t, y: np.ndarray):

#         y = y.reshape(hmn.shape)
#         RHS_mat = - 1j * (hmn @ y - y @ hmn)
#         return RHS_mat.flatten()

#     timestep = kwargs["timestep"]

#     fin = kwargs["fin"]
#     sol = solve_ivp( gen_RHS, y0 =cmn, t_span=[0, fin], t_eval=np.arange(0, fin, timestep))

#     print(sol.y.shape)
#     return sol





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


class ExponentiateSolver:

    def __init__(self) -> None:
        pass

    def solve(self, system : sys, init, timestep, fin, N, factor=1.0):

        h = system.set_hij() 
        steps = int(fin/timestep)
        M = expm( -1j * h * timestep)

        state = init
        states = np.zeros((steps + 1, init.shape[0]), dtype=complex)
        states[0] = state

        for i in range(steps):
            
            print(i)
            state = M @ state
            states[i + 1] = state

        return states




class SDSolver:

    def __init__(self) -> None:
        pass


    def solve(self, system : SD, init, timestep, fin, factor=1.0):

        h = system.set_hij()
        # for this solver we only have one U
        w, U = eigh(h)

        total = init.shape[0]
        E = np.zeros(( total, total), dtype=complex)

        for i in range(total):
            for j in range(total):
                E[i, j] = np.exp( -1j * (w[i] - w[j]) * timestep)

        #Uc = U[:, :N]
        Ud = np.conj(np.transpose(U))
        #Ucd = np.conj(np.transpose(Uc))

        #print("Ud, Uc, Ucd, U shapes:", Ud.shape, Uc.shape, Ucd.shape, U.shape)
        steps = int(fin/timestep)

        # we calculate 6 currents: 1, 4, 7 in; 3, 6, 9 out
        current = np.zeros(( system.arr * 2, steps + 1))
        occs = np.zeros((steps + 1, total))
        contact_source = system.contact_source
        contact_drain = system.contact_drain

        arr_source = system.arr_source
        arr_drain = system.arr_drain
        occs[0] = np.abs(init) ** 2
        

        #CC_init = np.kron(np.conj(init), init).reshape((total, total))
        CC_init = np.diag(init)
        AA_init = Ud @ CC_init @ U
        AA = AA_init

        # plt.imshow(np.real(AA_init), origin='lower')
        # plt.show()

        for i in range(steps):

            print(i)
            AA = np.multiply(E, AA)
            CC = U @ AA @ Ud

            for j, arr in enumerate(arr_source):
                
                cc = CC[contact_source, arr]
                current[j, i + 1] = -2 * np.imag(cc) * factor

            for k, arr in enumerate(arr_drain):

                cc = CC[arr, contact_drain]
                current[k + len(arr_source), i + 1] = -2 * np.imag(cc) * factor


            occs[i + 1] = np.diag(CC)

        return current, occs

