import numpy as np

def to_mathematica(M, val_dict):


    M = str(np.array2string(M, separator=', ')).replace('[', '{').replace(']','}')
    for val in val_dict:

        val_str = val_dict[val]
        if float(val) == int(val):
            val = str(int(val)) + '.'

        else:
            val = str(val)

        M = M.replace( val, val_str)

    return M

