## Small library for reading images and spectra exported as txt from Renishaw software

import numpy as np

def load_image(fname):
    "Loads image exported as text from Renishaw software"
    d = np.loadtxt(fname)
    nx = len(np.unique1d(d[:,0]))
    ny = len(np.unique1d(d[:,1]))
    nnu = len(np.unique1d(d[:,2]))

    d2 = d.reshape(ny,nx,nnu,-1)
    x,y = d2[0,:,0,0], d2[:,0,0,1]
    nu = d2[0,0,:,2]

    knu = np.argsort(nu)

    out = d2[:,:,:,3]
    return out[:,:,knu], nu[knu], x, y

def process_image(arr, fn):
    "Apply a function to each spectrum in an image"
    figsh = arr.shape[:2]
    new_shape= figsh + fn(arr[0,0]).shape
    out = zeros(new_shape)
    nrow,ncol = figsh
    for r in xrange(nrow):
        for c in xrange(ncol):
            out[r,c] = fn(arr[r,c])
    return out

def inspect_image(arr, nu):
    "Inspect an image: show spectra for a pixel mouse hovers on"
    import pylab as pl
    fh  = pl.figure();
    ax2 = pl.subplot(122);
    plh = pl.plot(nu, arr[0,0], 'k-')[0]
    axis((nu[0], nu[-1], arr.min(), arr.max()))

    def _on_hover(event):
        if event.inaxes == ax1:
            x = int(event.xdata)
            y = int(event.ydata)
            plh.set_data((nu, arr[y,x]))
            fh.canvas.draw()
            pass
    ax1 = pl.subplot(121);
    pl.imshow(range_sum(arr, nu, (nu.min(), nu.max())),
              aspect='equal', interpolation='nearest')
    fh.canvas.mpl_connect('motion_notify_event', _on_hover)

def in_range(vec, range):
    return (vec > range[0]) * (vec < range[1])

def range_sum(arr, nu, nurange):
    "Sum intensities in the given nu range"
    return process_image(arr, lambda x: np.sum(x[in_range(nu,nurange)]))
    
def load_spectrum(fname):
    "Read spectrum exported as text from Renishaw software"
    spec = np.loadtxt(fname)
    knu = np.argsort(spec[:,0])
    return spec[knu,:]
              

