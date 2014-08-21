# License: GPL
# Copyright: Alexey Brazhe, 2010-2013

import numpy as np
import pylab as pl
import renishaw as rshaw

from scipy.interpolate import splrep, splev
#from scipy.optimize import leastsq

def unique_tag(tags, max_tries = 1e4):
    n = 0
    while n < max_tries:
        tag = np.random.rand()
        n += 1
        if not tag in tags:
            return tag
    return "Err"

def rezip(coll):
    import itertools as itt
    return itt.izip(*coll)

def in_range(n, region):
    return (n > region[0]) * (n < region[1])

def ind_to_val(ind, vrange, nsteps, dv=None):
    if dv is None: dv = (vrange[1] - vrange[0])/nsteps
    return vrange[0] + ind*dv

def val_to_ind(val, vrange, nsteps, dv=None):
    if dv is None: dv = (vrange[1] - vrange[0])/nsteps
    return (val-vrange[0])/dv


def locextr(v, x = None, mode = 1, **kwargs):
   "Finds local maxima when mode = 1, local minima when mode = -1"
   res = 0.05
   if x is None:
       x = np.arange(len(v))
       xfit = np.arange(0,len(v), res)
   else:
       xfit = np.arange(x[0], x[-1], res)
   tck = splrep(x, v, **kwargs)
   vals = splev(xfit, tck)
   dersign = mode*np.sign(splev(xfit, tck, der=1))
   extrema = dersign[:-1] - dersign[1:] > 1.5
   return xfit[extrema], vals[extrema]

def any_obj_contains(objects, event):
    "Checks if event is contained by any ROI"
    if len(objects) < 1 : return False
    return reduce(lambda x,y: x or y,
                  [o.obj.contains(event)[0]
                   for o in objects.values()])




def bspline_denoise(sig, phi = np.array([1./16, 1./4, 3./8, 1./4, 1./16])):
    def _mirrorpd(k, L):
	if 0 <= k < L : return k
	else: return -(k+1)%L
    L = len(sig) 
    padlen = len(phi)
    assert L > padlen
    indices = map(lambda i: _mirrorpd(i, L),
                  range(-padlen, 0) + range(0,L) + range(L, L+padlen))
    padded_sig = sig[indices]
    apprx = np.convolve(padded_sig, phi, mode='same')[padlen:padlen+L]
    return apprx


def onpress_peaknotifier(event, ax, x,y, coll):
    tb = pl.get_current_fig_manager().toolbar
    if event.inaxes != ax or tb.mode != '': return
    if event.button != 1 : return
    if any_obj_contains(coll, event) : return
    ex,ey = event.xdata, event.ydata
    max_pos, max_vals = locextr(y, x)
    j = np.argmin(abs(max_pos - ex))
    peak_x, peak_y = max_pos[j], max_vals[j]
    #print ' '.join(['%3.3f'%v for v in (x, y, peak_x, peak_y)])
    label = unique_tag(coll.keys())
    lp =ax.plot(peak_x, peak_y, 'ro', alpha=0.5, label=label)[0]
    coll[label] = LabeledPoint(lp, coll)
    #ax.text(peak_x, peak_y, '%3.3f, %3.3f'%(peak_x, peak_y), size='x-small')
    pl.draw()
        

def plot_with_peaks(x, y, **kwargs):
    peak_points = {}
    def print_coll(event):
        if event.key == 'e':
            for xy in sorted([lp.get_xy() for lp in peak_points.values()],
                             key=lambda x:x[0]):
                print '%3.3f, %3.3f' % xy
        return
    newax = pl.figure().add_subplot(111)
    newax.plot(x,y)
    newax.set_title("Click on peaks to select, press 'e' to export selected...")
    canvas = newax.figure.canvas
    canvas.mpl_connect('button_press_event',
                       lambda e: onpress_peaknotifier(e,newax,x,y,peak_points))
    canvas.mpl_connect('key_press_event', print_coll)
    return peak_points

def loc_max_pos(v):
    print len(v)
    return [i for i in xrange(1,len(v)-1)
            if (v[i] >= v[i-1]) and (v[i] > v[i+1])]


import pickle

class RamanCooker():
    linestyle = 'k-'
    def cook(self, nu, spectrum, mode='spans', **kwargs):
        if mode =='spans':
            self.cook_spans(nu, spectrum, **kwargs)
        elif mode == 'knots':
            self.cook_knots(nu, spectrum, **kwargs)
        else:
            print "Unknown mode: ", mode
    def shared_setup(self, nu, spectrum):
        self.connected = False
        L = min(len(nu), len(spectrum))
        self.nu = nu[:L]
        self.sp = spectrum[:L]
        self.axspl = pl.figure().add_subplot(111)
        self.figspl = self.axspl.figure
        self.axspl.plot(self.nu, self.sp, self.linestyle)
        self.pressed = None
        self.plfit = self.axspl.plot([],[],'m--')[0]


    def cook_spans(self, nu, spectrum, s=200, 
                   xspans=[], bands=None):
        """
        UI for removing baseline with splines,
        
        Inputs:
        -------
        - nu : values for wavelength shift
        - spectrum: intensity values
        - s : spline smoothing parameter(200 by default)
        - xspans: Spans (list of pairs of left,right nu values )
                  to exclude from fitting (empty list by default)
        """
        self.shared_setup(nu, spectrum)
        self.spans = {}
        self.bands = bands
	self.mode = 'spans'
        if bands:
            _ = [pl.axvline(band, color='r') for band in bands]
        self.spl_k= 3
        self.curr_span = None
        self.spankw = {'facecolor':(0.9,0.9,0.9), 'alpha':0.9}
        self.load_recipe({'s':s, 'xspans':xspans})
        if not self.connected:
            canvas = self.figspl.canvas
            canvas.mpl_connect('button_press_event',self.onpress_spans)
            canvas.mpl_connect('motion_notify_event',self.onmotion_spans)
            canvas.mpl_connect('button_release_event',self.onrelease_spans)
            canvas.mpl_connect('key_press_event',self.onkeypress_spans)
            canvas.mpl_connect('scroll_event',self.onpress_spans)
            self.connected = True
        return self.axspl

    def cook_knots(self, nu, spectrum, bands = None, nuspan=20, xlocs=None):
        self.shared_setup(nu, spectrum)
        self.bands = bands
        self.nuspan = nuspan
        self.knots = {}
	self.mode = 'knots'
	if xlocs is not None:
	    self.load_knots(xlocs)
        if bands:
            _ = [pl.axvline(band, color='r') for band in bands]
        if not self.connected:
            canvas = self.figspl.canvas
            canvas.mpl_connect('button_press_event',self.onpress_knots)
            canvas.mpl_connect('key_press_event',self.onkeypress_knots)
            self.connected = True
        return self.axspl

    def cook_als(self, nu, spectrum, lam = None, p = None):
	if lam is None:
	    lam = len(spectrum)**2
	if p is None:
	    p = 0.1
	self.shared_setup(nu, spectrum)
	baseline = als(spectrum, lam, p)
	self.plfit.set_data(nu, baseline)
	return self.axspl
	
		
 
    def update_smooth_hint(self,renewp=False):
        if (not hasattr(self, 'smooth_hint')) or renewp:
            self.smooth_hint = pl.text(1.0,1.1,
                                       "Smoothing: %3.2e"% self.spl_smooth,
                                       horizontalalignment = 'right',
                                       verticalalignment = 'top',
                                       transform = self.axspl.transAxes)
        else:
            self.smooth_hint.set_text("Smoothing: %3.2e"%self.spl_smooth)            
            
    def redraw(self, smooth_hint = True):
        if smooth_hint:
            self.update_smooth_hint()
        self.figspl.canvas.draw()
    def onkeypress_spans(self,event):
        if event.key == 'a':
            self.apply_recipe(show=True)
    def onkeypress_knots(self,event):
        if event.key == 'a':
            self.apply_recipe(show=True)
    def onpress_spans(self, event):
        tb = pl.get_current_fig_manager().toolbar
        if event.inaxes != self.axspl or tb.mode != '': return
        x,y = event.xdata, event.ydata
        axrange = self.axspl.get_xbound() + self.axspl.get_ybound()
        if event.button is 1 :
            pass
        elif event.button is 3:
            self.pressed = event.xdata,event.ydata
            x0 = event.xdata
            self.curr_span = self.addspan((x0,x0))
            self.axspl.axis(axrange)
            pass
        elif event.button is 'up':
            self.spl_smooth *= 1.1
        elif event.button is 'down':
            self.spl_smooth /= 1.1
        self.apply_recipe()
        self.axspl.axis(axrange)
        self.redraw()
    def onpress_knots(self, event):
        tb = pl.get_current_fig_manager().toolbar
        if event.inaxes != self.axspl or tb.mode != '': return
        if any_obj_contains(self.knots, event) : return
        x,y = event.xdata, event.ydata
        axrange = self.axspl.get_xbound() + self.axspl.get_ybound()
        if event.button is 1:
            kp = self.addknot(event.xdata)
            self.knots[kp.tag] = kp
        self.apply_recipe()
        self.axspl.axis(axrange)
        self.redraw(smooth_hint = False)
        
    def onmotion_spans(self,event):
        if (self.pressed is None) or (event.inaxes != self.axspl):
            return
        x0 = self.pressed[0]
        x = event.xdata
        if self.curr_span:
            self.curr_span.set_xy([[x0,0], [x0,1], [x,1], [x,0], [x0,0]])
        self.redraw()

    def onrelease_spans(self,event):
        if (not self.curr_span) or (event.inaxes != self.axspl):
            return
        vert = self.curr_span.get_xy()
        w = abs(vert[0][0] - vert[2][0])
        if w>5:
            label = self.curr_span.get_label()
            self.spans[label]  = Span(self.curr_span, self.spans)
        else:
            self.curr_span.remove()
        self.curr_span = None
        self.apply_recipe()

    def addspan(self, (x1,x2)):
        label = unique_tag(self.spans.keys())
        spanh = pl.axvspan(x1,x2,**self.spankw)
        spanh.set_label(label)
        return spanh
    def addknot(self, x):
        label = unique_tag(self.knots.keys())
        nux = in_range(self.nu, (x-self.nuspan, x+self.nuspan))
        if np.any(nux):
            y = np.mean(self.sp[nux]) # todo: use median as an option
            lp =self.axspl.plot(x, y, #xerr=self.nuspan,
                                color='red', marker='o',
                                alpha=0.5, label=label)[0]
            kp = KnotPoint(lp, self.knots,  self.nu, self.sp, self.nuspan)
            return kp

    def get_weights(self, nu):
        weights = np.ones(len(nu))
        mask_regs = [s.xspan() for s in self.spans.values()]
        if len(mask_regs) > 0:
            masked = reduce(lambda a,b:a+b,
                            map(lambda x: in_range(nu,x), mask_regs))
            weights -= masked*(1-1e-6)
        return weights
   
    def knots2pars(self, nu = None, sp = None):
        xlocs = sorted([k.get_xloc() for k in self.knots.values()])
        if nu is None : nu = self.nu
        if sp is None : sp = self.sp
        nuxs = [in_range(nu, (x-self.nuspan, x + self.nuspan)) for x in xlocs]
        y = [np.mean(sp[nx]) for nx in nuxs]
        return xlocs,y

    def apply_recipe(self, show=False):
        if self.mode == 'spans':
            smooth_hint = True
            process = self.process_spans
        elif self.mode=='knots':
            smooth_hint = False
            process = self.process_knots
        else:
            print "[RamanCooker] applying recipe: unknown mode"
            return
        out = process(self.nu, self.sp, 'full')
        if out:
            self.plfit.set_data(self.nu, out[0])
            self.redraw(smooth_hint=smooth_hint)
            if show:
                plot_with_peaks(self.nu, out[1])
                if self.bands:
                    _ = [pl.axvline(b, color='r') for b in self.bands] 
    def process_spans(self, nu, sp, ret = None):
        w = self.get_weights(nu)
        nux,spx = nu[w>0.5], sp[w>0.5]
        spsum = spx.sum()
        tck = splrep(nux,spx,s=self.spl_smooth*spsum,k=self.spl_k)
        sp_fit = splev(nu, tck)
        if ret == 'full':
            return sp_fit, sp - sp_fit
        else:
            return sp-sp_fit
        
    def process_knots(self, nu=None, sp=None, ret = None):
        if nu is None: nu = self.nu
        if sp is None: sp = self.sp
        xlocs, y = self.knots2pars(nu, sp)
        if len(xlocs) > 5:
            tck = splrep(xlocs, y)
            sp_fit = splev(nu, tck)
            diff = self.sp - sp_fit
            if ret == 'full':
                return sp_fit, sp - sp_fit
            else:
                return sp-sp_fit
        else:
            return None

    def export_recipe(self, out=None):
        if self.mode == 'spans':
            rec =  {'s': self.spl_smooth,
                    'xspans':self.export_xspans()}
        elif self.mode == 'knots':
            rec = {'nuspan': self.nuspan,
                   'xlocs': self.export_knots()}
        else:
            print "[RamanCooker] unknown mode: ", self.mode
        if isinstance(out, str):
            fo = file(out, 'w')
            pickle.dump(rec, fo)
            fo.close()
        return rec
    def load_recipe(self, recipe, mode=None):
        if isinstance(recipe, dict):
            rec = recipe
        elif isinstance(recipe, str):
            fo = file(recipe,'r')
            rec = pickle.load(fo)
            fo.close()
        else:
            print "[RamanCooker] Can't recognize type of recipe!"
            return
        if mode is None:
            if not hasattr(self, 'mode'):
                print "[RamanCooker] Can't determine mode"
                return
            mode = self.mode
        else:
            self.mode = mode
        if mode == 'spans':
            self.spl_smooth = rec['s']
            self.load_spans(rec['xspans'])
        elif mode == 'knots':
            self.nuspan = rec['nuspan']
            self.load_knots(rec['xlocs'])
            
        
        
    def export_xspans(self):
        return [s.xspan() for s in self.spans.values()]
    def export_knots(self):
        return sorted([k.get_xloc() for k in self.knots.values()])
    def load_spans(self, xspans):
        for xsp in xspans:
            newspan = self.addspan(xsp)
            label = newspan.get_label()
            self.spans[label] = Span(newspan, self.spans)
        self.apply_recipe()
    def load_knots(self, xlocs):
        for x in xlocs:
            kp = self.addknot(x)
            self.knots[kp.tag] = kp
        self.apply_recipe()
        

# --- ROIS

class DraggableObj(object):
    verbose = True
    def __init__(self, obj, coll):
        self.coll = coll
        self.obj = obj
        self.connect()
        self.pressed = None
        self.tag = obj.get_label()
    def redraw(self):
        self.obj.axes.figure.canvas.draw()
    def event_ok(self, event, should_contain=False):
        containsp = True
        if should_contain:
            containsp, _ = self.obj.contains(event)
        return event.inaxes == self.obj.axes and \
               containsp and \
               pl.get_current_fig_manager().toolbar.mode ==''
    def connect(self):
        "connect all the needed events"
        cf = self.obj.axes.figure.canvas.mpl_connect
        self.cid = {
            'press': cf('button_press_event', self.on_press),
            'release': cf('button_release_event', self.on_release),
            'motion': cf('motion_notify_event', self.on_motion),
            'scroll': cf('scroll_event', self.on_scroll),
            'type': cf('key_press_event', self.on_type)
            }
    def on_motion(self, event):
        if not (self.event_ok(event, False) and self.pressed):
            return
        p = event.xdata,event.ydata
        self.move(p)
        self.redraw()
    def on_release(self, event):
        if not self.event_ok(event):
            return
        self.pressed = None
        self.redraw()
    def on_scroll(self, event): pass
    def on_type(self, event): pass
    def on_press(self, event): pass
    def disconnect(self):
        map(self.obj.axes.figure.canvas.mpl_disconnect,
            self.cid.values())
    def get_color(self):
        return self.obj.get_facecolor()

class Span(DraggableObj):
    def on_press(self, event):
        if not self.event_ok(event, True): return
        #x0,y0 = self.obj.get_xdata(), self.obj.get_ydata()
        if event.button is 1:
            pass
        #    self.pressed = event.xdata, event.ydata, x0, y0
        elif event.button is 2:
            pass
          
        elif event.button is 3:
            p = self.coll.pop(self.tag)
            self.obj.remove()
            self.disconnect()
            self.redraw()
    def xspan(self):
        vert = self.obj.get_xy()
        l,r = vert[0][0], vert[2][0]
        return min(l,r), max(l,r)

class DraggablePoint(DraggableObj):
    "Labeled Point"
                              
    def on_press(self, event):
        if not self.event_ok(event, True): return
        x0,y0 = self.obj.get_xdata(), self.obj.get_ydata()
        if event.button is 1:
            self.pressed = event.xdata, event.ydata, x0, y0
        elif event.button is 2:
            pass
        elif event.button is 3:
            p = self.coll.pop(self.tag)
            self.obj.remove()
            self.disconnect()
            self.redraw()

    def move(self, p):
        "Move the ROI when the mouse is pressed over it"
        xp,yp, x0, y0 = self.pressed
        dx = p[0] - xp
        dy = p[1] - yp
        #x0, y0 = self.obj.get_xdata(), self.obj.get_ydata()
        self.obj.set_data(x0+dx, y0+dy)
        self.redraw()

    def set_xy(self,(x,y)):
        self.obj.set_data([x],[y])
        self.redraw()
        
    def get_xy(self):
        return self.obj.get_xdata()[0], self.obj.get_ydata()[0]



class KnotPoint(DraggablePoint):
    "Knot point"
    def __init__(self, obj, coll, nu, sp, nuspan):
        super(KnotPoint, self).__init__(obj, coll)
        self.nu = nu
        self.sp = sp
        self.nuspan = nuspan
    def move(self, p):
        "Move the KnotPoint when the mouse is pressed on it"
        xp,yp, x0, y0 = self.pressed
        dx = p[0] - xp
        x = x0+dx
        nx = in_range(self.nu, (x- self.nuspan, x+self.nuspan))
        y = np.mean(self.sp[nx])
        self.obj.set_data(x,y)
        self.redraw()
    def get_xloc(self):
        return self.obj.get_xdata()[0]

        

class LabeledPoint(DraggablePoint):
    "Labeled Point"
    def __init__(self, obj, coll):
        super(LabeledPoint, self).__init__(obj, coll)
        x,y = self.get_xy()
        self.textlabel = obj.axes.text(x,y,"%3.3f, %3.3f"%(x,y))
    def on_press(self,event):
        DraggablePoint.on_press(self,event)
        if self.event_ok(event, True) and event.button is 3:
            self.textlabel.remove()

    def redraw(self):
        xy = self.get_xy()
        self.textlabel.set_position(xy)
        self.textlabel.set_text("%3.3f, %3.3f"%xy)
        DraggableObj.redraw(self)

                              
### --- Non-interactive (global) baseline algorithms --
from scipy import sparse
from scipy.sparse.linalg import spsolve

def als(y, lam, p, niter=20, tol=1e-5):
    """Implements an Asymmetric Least Squares Smoothing
    baseline correction algorithm
    (P. Eilers, H. Boelens 2005)
    """
    L = len(y)
    D = sparse.csc_matrix(np.diff(np.eye(L),2))
    w = np.ones(L)
    zprev = None
    for i in xrange(niter):
	W = sparse.spdiags(w, 0, L, L)
	Z = W + lam*D.dot(D.T)
	z = spsolve(Z,w*y)
	w = p*(y>z) + (1-p)*(y<=z)
	if zprev is not None:
	    err = np.sum((z-zprev)**2)
	    if err < tol:
		return z
	zprev = z
    return z

def fillpeaks(y, lam, hwi, niter, nint):
    """Implements and iterative baseline correction algorithm based on
    mean suppression (originally written for R by Kristian Hovde Liland) """
    L = len(y)
    D = sparse.csc_matrix(np.diff(np.eye(L),2))
    ww = np.ones(L)
    # make decreasing windows:
    if niter > 1:
	d1 = np.log10(hwi)
	d2 = 0.0
	_x = np.arange(0,niter-1)*(d2-d1)/(np.floor(niter)-1)
	w = np.ceil(10**(np.concatenate((d1+_x, [d2]))))
    else:
	w = [hwi]
    # Primary smoothing
    W = sparse.spdiags(ww, 0, L, L)
    Z = W + lam*D.dot(D.T)
    z = spsolve(Z,y)
    #Center points
    lims = np.linspace(0,L-1, nint+1)
    lefts = np.ceil(lims[:-1])
    rights = np.floor(lims[1:])
    minip = np.round((lefts+rights)/2)
    xx = np.array([np.min(z[l:r+1]) for l,r in zip(lefts,rights)])
    for k in range(niter):
	w0 = w[k]
	# to the right
	for i in range(1,nint-1):
	    l,r = max(i-w0,0), min(i+w0+1, nint)
	    a = np.mean(xx[l:r])
	    xx[i] = min(a,xx[i])
	# to the left
	for i in range(1,nint-1):
	    j = nint-i
	    l,r = max(j-w0,0), min(j+w0+1, nint)
	    a = np.mean(xx[l:r])
	    xx[j] = min(a, xx[j])
    tck = splrep(minip,xx)
    return splev(np.arange(L),tck)
	
    
    
def aws(v, **kwargs):
    import wavelets
    z = wavelets.asymmetric_smooth(v,verbose=True,**kwargs)
    return z


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
