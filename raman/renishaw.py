## Small library for reading images and spectra exported as txt from Renishaw software

import numpy as np
import pylab as pl
import operator as op
from scipy import stats,ndimage,signal

bclose, bopen = ndimage.binary_closing, ndimage.binary_opening

def load_spectrum(fname):
    "Read spectrum exported as text from Renishaw software"
    spec = np.loadtxt(fname)
    knu = np.argsort(spec[:,0])
    return spec[knu,:]
              
def load_image(fname):
    "Loads image exported as text from Renishaw software"
    d = np.loadtxt(fname)
    nx = len(np.unique(d[:,0]))
    ny = len(np.unique(d[:,1]))
    nnu = len(np.unique(d[:,2]))

    d2 = d.reshape((ny,nx,nnu,-1), order = 'C')
    y,x = d2[0,:,0,0], d2[:,0,0,1]
    nu = d2[0,0,:,2]

    knu = np.argsort(nu)

    pre_out = d2[:,:,:,3]
    out = np.zeros((nx,ny,nnu))
    print pre_out.shape, out.shape
    for j,k in enumerate(knu):
        out[:,:,j] = pre_out[:,:,k].T
    return out, nu[knu], x, y




