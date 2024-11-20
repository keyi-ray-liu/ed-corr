import numpy as np
import matplotlib.pyplot as plt





def g(t, E, ELR):
    return np.sqrt( 4 * (t ** 2) - (E - ELR)**2)

def transmission(E, mu, t):


    EL = mu/2
    ER = -mu/2

    gg = t ** 2

    gL = g(t, E, EL)
    gR = g(t, E, ER)

    mb = ( E - EL - gg * (E- ER))
    val = 4 * gg * gL * gR / (mb ** 2 + (gL + gg * gR) ** 2)
    return val


def plot():

    T = np.vectorize(transmission)
    mu =0.5
    
    for t in (0.2, 0.5, 1):

        Es = np.linspace(-2 * t, 2 * t, 1000)
        plt.plot( Es/t, T(Es, mu, t))

    plt.show()

if __name__ == '__main__':

    plot()