## Small library for reading images and spectra exported as txt from Renishaw software

import numpy as np
import pylab as pl
from scipy import stats

def load_image(fname):
    "Loads image exported as text from Renishaw software"
    d = np.loadtxt(fname)
    nx = len(np.unique1d(d[:,0]))
    ny = len(np.unique1d(d[:,1]))
    nnu = len(np.unique1d(d[:,2]))

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

def process_image(arr, fn):
    "Apply a function to each spectrum in an image"
    figsh = arr.shape[:2]
    new_shape= figsh + fn(arr[0,0]).shape
    out = np.zeros(new_shape)
    nrow,ncol = figsh
    for r in xrange(nrow):
        for c in xrange(ncol):
            out[r,c] = fn(arr[r,c])
    return out

def range_sum(arr, nu, nurange):
    "Sum intensities in the given nu range"
    return process_image(arr, lambda x: np.sum(x[in_range(nu,nurange)]))

def range_rms(arr, nu, nurange):
    "RMS intensities in the given nu range"
    return process_image(arr, lambda x: pl.rms_flat(x[in_range(nu,nurange)]))
    
def in_range(vec, range):
    return (vec > range[0]) * (vec < range[1])

def range_sum(arr, nu, nurange):
    "Sum intensities in the given nu range"
    return process_image(arr, lambda x: np.sum(x[in_range(nu,nurange)]))

def range_rms(arr, nu, nurange):
    "RMS intensities in the given nu range"
    return process_image(arr, lambda x: pl.rms_flat(x[in_range(nu,nurange)]))
    



def toggle(val):
    if val : return False
    else : return True
    
def inspect_image(arr, nu, 
                  fn = range_rms, nuint=None,
                  **kwargs):
    "Inspect an image: show spectra for a pixel mouse hovers on"
    fh  = pl.figure();
    ax2 = pl.subplot(122);
    plh = pl.plot(nu, arr[0,0], 'k-')[0]
    mn = stats.scoreatpercentile(arr.flatten(), 0.1)
    mx = stats.scoreatpercentile(arr.flatten(), 99.9)
    print mn, mx
    pl.axis((nu[0], nu[-1], mn, mx))

    if nuint is None:
        nuint = (nu.min(), nu.max())

    hoverp = [True, None]

    def _on_hover(event):
        if (event.inaxes == ax1) and hoverp[0]:
            x = int(event.xdata)
            y = int(event.ydata)
            plh.set_data((nu, arr[y,x]))
            fh.canvas.draw()
            pass
    def _on_press(event):
        if (event.inaxes != ax1):
            return 
        hoverp[0] = toggle(hoverp[0])
        if not hoverp[1]:
            xh = ax1.axhline(event.ydata)
            yh = ax1.axvline(event.xdata)
            hoverp[1] = (xh, yh)
            fh.canvas.draw()
        else:
            xh,yh = hoverp[1]
            xh.remove(); yh.remove()
            hoverp[1] = None
            fh.canvas.draw()
        
        
    ax1 = pl.subplot(121);
    pl.imshow(fn(arr, nu, nuint), aspect='equal', **kwargs)
    fh.canvas.mpl_connect('motion_notify_event', _on_hover)
    fh.canvas.mpl_connect('button_press_event', _on_press)



def load_spectrum(fname):
    "Read spectrum exported as text from Renishaw software"
    spec = np.loadtxt(fname)
    knu = np.argsort(spec[:,0])
    return spec[knu,:]
              

