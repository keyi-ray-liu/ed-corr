
import matplotlib.pyplot as plt

def test():

    fig = plt.figure(
                    layout="constrained", 
                     figsize = (20, 20)
                     )

    row = 3
    subcol = 3
    subrow = 3
    col = 1
    subfigs = fig.subfigures(nrows=row, ncols=col)
    axes = {}

    overall = ['A', 'B', 'C']
    names = [ ['a', *[str(i + j * subcol) for i in range(subcol)], 'b', 'c'] for j in range(subrow)]

    for i in range(row):

        subfig : plt.Figure = subfigs[i]

        mosaic = [[ overall[i] + name for name in namerow] for namerow in names]
        width_ratios = [ 1.0, 0.5, 0.5, 0.5, 1.0, 1.0]
        height_ratios = [ 1.0, 1.0, 1.0]

        axes |= subfig.subplot_mosaic(mosaic,  width_ratios=width_ratios, height_ratios=height_ratios)

        for ax_name in axes:

            axes[ax_name].set_title(ax_name)


    fig.savefig('test.png')

test()