import operator as op
import numpy as np
import pylab as pl

from scipy import stats,ndimage,signal
bclose, bopen = ndimage.binary_closing, ndimage.binary_opening

#### Raman spectral images are orgaized such as arr[i,j,:] gives a spectrum at
#### location i,j (i is rows or Y, j is columns or X)


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

def imagemap(fn, arr):
    "Map a function on an array, wrapper for process_image"
    return process_image(arr, fn)

def range_imagemap(fn, arr, nu, nurange):
    return imagemap(lambda x: fn(x[_in_range(nu,nurange)]))


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

def adaptive_median_filter(m):
    out = m.copy()
    mask = _in_range(m, (np.percentile(m,1.0), np.percentile(m,99.0)))
    for loc in zip(*np.where(mask < 1)):
        out[loc] = np.median([m[x] for x in _neighbours_x(loc) if _valid_loc(x, m.shape)])
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

    
def inspect_image(arr, nu, 
                  fn = range_rms, nuint=None,
                  filt_output = False,
                  **kwargs):
    "Inspect an image: show spectra for a pixel mouse hovers on"
    fh  = pl.figure();
    ax2 = pl.subplot(122);
    plh = pl.plot(nu, arr[0,0], 'k-')[0]
    mn, mx = [np.percentile(arr,p) for p in (0.5, 99.5)]
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


def _in_range(vec, range):
    return (vec > range[0]) * (vec < range[1])

def _neighbours_x(loc):
    n = len(loc)
    d = np.diag(np.ones(n))
    return map(tuple, np.concatenate((d,-d)) + loc)


def _valid_loc(loc,shape):
    "location not outside bounds"
    return reduce(op.__and__, [(0 <= x < s) for x,s in zip(loc,shape)])
