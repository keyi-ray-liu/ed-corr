import matplotlib.pyplot as plt
import numpy as np


def reference_energy():


    emu = np.loadtxt('ref/finiteMUenergy')
    lrmu = np.loadtxt('ref/finiteMULR')

    eleft = np.loadtxt('ref/Leftenergy')
    lrleft = np.loadtxt('ref/LeftLR')

    refs =[(emu, lrmu), (eleft, lrleft)]
    

    titles = ['mu', 'left']
    mus = [0.25, -3.0]
    fig, axes = plt.subplots( 1, len(refs), figsize = (7 * len(refs), 5))

    for i, (e, lr) in enumerate(refs):

        mu = mus[i]
        ax : plt.Axes = axes[i]
        

        leftind = np.argwhere( lr > 0).flatten()
        rightind = np.argwhere( lr < 0).flatten()
        ax.plot( leftind, e[leftind], c ='red', label='Source')
        ax.plot( rightind, e[rightind], c ='blue', label='Drain')

        ax.plot( leftind, e[leftind] - mu, c ='red', label='Source zero bias', linestyle = 'dashed')
        ax.plot( rightind, e[rightind] + mu, c ='blue', label='Drain zero bias', linestyle = 'dashed')

        ax.set_title(titles[i])
        ax.legend()

        ax.hlines( [-5.04, 0.625], 0, ax.get_xlim()[-1], linestyle='dotted')

    fig.savefig('/Users/knl20/Desktop/Code/TN/SD/plots/reference.pdf')