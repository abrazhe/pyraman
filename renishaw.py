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


def mirrorpd(k, L):
    if 0 <= k < L : return k
    else: return -(k+1)%L


def bspline_denoise(sig, phi = np.array([1./16, 1./4, 3./8, 1./4, 1./16])):
    L = len(sig) 
    padlen = len(phi)
    assert L > padlen
    indices = map(lambda i: mirrorpd(i, L),
                  range(-padlen, 0) + range(0,L) + range(L, L+padlen))
    padded_sig = sig[indices]
    apprx = np.convolve(padded_sig, phi, mode='same')[padlen:padlen+L]
    return apprx


def locextr(v, x = None, mode = 1, **kwargs):
   from scipy.interpolate import splrep,splev
   "Finds local maxima when mode = 1, local minima when mode = -1"
   if x is None:
       x = np.arange(len(v))
       xfit = np.linspace(1,len(v), len(v)*10)
   else:
       xfit = np.linspace(x[0], x[-1], len(v)*10)
   tck = splrep(x, v, **kwargs)
   vals = splev(xfit, tck)
   dersign = mode*np.sign(splev(xfit, tck, der=1))
   extrema = dersign[:-1] - dersign[1:] > 1.5
   return xfit[extrema], vals[extrema]


def find_peak(x, band,nhood=6):
    """
    given nu values and band to search around, return a function to map
    spectrum to location of the highest peak in the nhood of band
    """
    def _(y):
	yn = y/np.std(y)
	peaks = [lm for lm in zip(*locextr(bspline_denoise(y),x))
		 if np.abs(lm[0]-band) <=nhood]
	if len(peaks):
	    loc, amp = peaks[np.argmax([p[1] for p in peaks])]
	else:
	    loc, amp = -1,-1
	return np.array((loc,amp))
    return _
	

def imagemap(fn, arr):
    "Map a function on an array, wrapper for process_image"
    return process_image(arr, fn)


def in_range(vec, range):
    return (vec > range[0]) * (vec < range[1])

def range_imagemap(fn, arr, nu, nurange):
    return imagemap(lambda x: fn(x[in_range(nu,nurange)]))

#def range_rms(arr, nu, nurange):
#    return range_imagemap(pl.rms_flat, arr, nu, nurange)

def range_rms(arr, nu, nurange):
    "RMS intensities in the given nu range"
    def _fn(x) : return pl.rms_flat(x[in_range(nu,nurange)])
    return imagemap(_fn, arr)
    

def range_sum(arr, nu, nurange):
    "Sum intensities in the given nu range"
    def _fn(x) : return np.sum(x[in_range(nu,nurange)])
    return imagemap(_fn, arr)

def range_mean(arr, nu, nurange):
    "Sum intensities in the given nu range"
    def _fn(x) : return np.mean(x[in_range(nu,nurange)])
    return imagemap(_fn, arr)


def valid_loc(loc,shape):
    "location not outside bounds"
    return reduce(op.__and__, [(0 <= x < s) for x,s in zip(loc,shape)])


def neighbours_x(loc):
    n = len(loc)
    d = np.diag(np.ones(n))
    return map(tuple, np.concatenate((d,-d)) + loc)


def adaptive_median_filter(m):
    out = m.copy()
    mask = in_range(m, (percentile(m,1.0), percentile(m,99.0)))
    for loc in zip(*np.where(mask < 1)):
        out[loc] = np.median([m[x] for x in neighbours_x(loc) if valid_loc(x, m.shape)])
    return out


def simple_peak_ratio(arr, nu, nuint_n, nuint_d, fn=range_mean, filt_output=False):
    nominator = fn(arr, nu, nuint_n)
    denominator = fn(arr, nu, nuint_d)
    out = np.ma.masked_where(denominator < np.finfo(float).eps, nominator/denominator)
    if filt_output:
        return adaptive_median_filter(out)
    else:
        return out

def peak_ratio(arr, nu, nuint_n, nuint_d, fn = range_mean, filt_output=True):
    arr2 = arr
    nominator = fn(arr2, nu, nuint_n)
    denominator = fn(arr2, nu, nuint_d)
    mask_n = (nominator > np.median(nominator))*(nominator>0)
    mask_d = (denominator > np.median(denominator))*(denominator > 0)
    mask = bopen(bclose(mask_n * mask_d))
    out = mask*(nominator/denominator)
    if filt_output:
	out = adaptive_median_filter(out)
    return out

def peak_ratio2(arr, nu, nuint_n, nuint_d, fn = range_mean,
		mask = None,
                filt_output = True,
		as_masked_array = False,
                baseline_region = (1750,2000)):
    nominator = fn(arr, nu, nuint_n)
    denominator = fn(arr, nu, nuint_d)
    x = range_rms(arr, nu, baseline_region)
    th_rms = x.mean() + 5.0*x.std()
    x = fn(arr, nu, baseline_region)
    th_fn = x.mean() + 5.0*x.std()
    if mask is None:
	mask_n = (range_rms(arr,nu, nuint_n) > th_rms)*(nominator > th_fn)
	mask_d = (range_rms(arr,nu, nuint_d) > th_rms)*(denominator >th_fn)
	mask = mask_n*mask_d
    mask = mask*(denominator > th_fn)*(nominator > th_fn)
    out = mask*(nominator/denominator)
    if filt_output:
        out = adaptive_median_filter(out)
    if as_masked_array:
	out = np.ma.masked_where(out <=0, out)
    return out


    
def toggle(val):
    if val : return False
    else : return True

def percentile(arr, p):
    return stats.scoreatpercentile(arr.flatten(), p)
    
def inspect_image(arr, nu, 
                  fn = range_rms, nuint=None,
                  filt_output = False,
                  **kwargs):
    "Inspect an image: show spectra for a pixel mouse hovers on"
    fh  = pl.figure();
    ax2 = pl.subplot(122);
    plh = pl.plot(nu, arr[0,0], 'k-')[0]
    mn, mx = [percentile(arr,p) for p in (0.5, 99.5)]
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
    if filt_output:
        arrshow = ndimage.median_filter(fn(arr, nu, nuint), 3) # looks nicer
    else:
        arrshow = fn(arr, nu, nuint)
    pl.imshow(arrshow, aspect='equal', **kwargs)
    fh.canvas.mpl_connect('motion_notify_event', _on_hover)
    fh.canvas.mpl_connect('button_press_event', _on_press)




