import numpy as np
import pyfits

def interpolate(filename):


    d,h = pyfits.getdata(filename,header=True)
    shape = d.shape
    #print(shape)
    # -- setup --
    #shape = (10,10)
    dim = len(shape)
    #data = np.random.random(shape)
    flag = np.zeros(shape, dtype=bool)
    #t_ct = int(data.size/5)

    flag = np.where(np.isnan(d),False,True)
    #flag.flat[np.random.randint(0, flag.size, t_ct)] = True
    # True flags the data
    # -- end setup --
    data = d
    slcs = [slice(None)]*dim

    while np.any(~flag): # as long as there are any False's in flag
        for i in range(dim): # do each axis
            # make slices to shift view one element along the axis
            slcs1 = slcs[:]
            slcs2 = slcs[:]
            slcs1[i] = slice(0, -1)
            slcs2[i] = slice(1, None)

            # replace from the right
            repmask = np.logical_and(~flag[slcs1], flag[slcs2])
            data[slcs1][repmask] = data[slcs2][repmask]
            flag[slcs1][repmask] = True

            # replace from the left
            repmask = np.logical_and(~flag[slcs2], flag[slcs1])
            data[slcs2][repmask] = data[slcs1][repmask]
            flag[slcs2][repmask] = True

    pyfits.writeto("TempFile.fits",data,h,clobber=True)
    return("TempFile.fits")