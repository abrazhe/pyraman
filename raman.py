# License: GPL
# Copyright: Alexey Brazhe, 2010

import numpy as np
import pylab as pl
import renishaw as rshaw

from scipy.interpolate import splrep, splev
from scipy.optimize import leastsq


def unique_tag(tags, max_tries = 1e4):
    n = 0
    while n < max_tries:
        tag = np.random.rand()
        n += 1
        if not tag in tags:
            return tag
    return "Err"

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
    def cook(self, nu, spectrum, s = 200, xspans = []):
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
        #if not hasattr(self, 'connected'):
        #    self.connected = False
        self.connected = False
        L = min(len(nu), len(spectrum))
        self.nu = nu[:L]
        self.sp = spectrum[:L]
        self.axspl = pl.figure().add_subplot(111)
        self.figspl = self.axspl.figure
        self.spans = {}
        self.spl_k= 3
        self.curr_span = None
        self.axspl.plot(self.nu, self.sp)
        self.pressed = None
        self.plfit = self.axspl.plot([],[],'m--')[0]
        self.spankw = {'facecolor':(0.9,0.9,0.9), 'alpha':0.9}
        self.load_recipe({'s':s, 'xspans':xspans})
        #self.update_smooth_hint(True)
        if not self.connected:
            canvas = self.figspl.canvas
            canvas.mpl_connect('button_press_event',self.onpress_spl)
            canvas.mpl_connect('motion_notify_event',self.onmotion_spl)
            canvas.mpl_connect('button_release_event',self.onrelease_spl)
            canvas.mpl_connect('key_press_event',self.onkeypress_spl)
            canvas.mpl_connect('scroll_event',self.onpress_spl)
            self.connected = True
        return self.axspl
    def cook_knots(self, nu, spectrum):
        self.connected = False
        L = min(len(nu), len(spectrum))
        self.nu = nu[:L]
        self.sp = spectrum[:L]
        self.axspl = pl.figure().add_subplot(111)
        self.figspl = self.axspl.figure
        self.nuspan = 5.0
        self.knots = {}
        self.axspl.plot(self.nu, self.sp)
        self.pressed = None
        self.plfit = self.axspl.plot([],[],'m--')[0]
        if not self.connected:
            canvas = self.figspl.canvas
            canvas.mpl_connect('button_press_event',self.onpress_spl_knots)
            #canvas.mpl_connect('motion_notify_event',self.onmotion_spl)
            #canvas.mpl_connect('button_release_event',self.onrelease_spl)
            canvas.mpl_connect('key_press_event',self.onkeypress_spl_knots)
            #canvas.mpl_connect('scroll_event',self.onpress_spl)
            self.connected = True
        return self.axspl
 
    def update_smooth_hint(self,renewp=False):
        #print "Renewp: ", renewp
        #print "Hasattr: ", hasattr(self, 'smooth_hint')
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
    def onkeypress_spl(self,event):
        if event.key == 'a':
            self.apply_spl2(show=True)
        pass
    def onkeypress_spl_knots(self,event):
        if event.key == 'a':
            self.apply_spl_knots(show=True)
        pass
    def onpress_spl(self, event):
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
        self.apply_spl2()
        self.axspl.axis(axrange)
        self.redraw()
    def onpress_spl_knots(self, event):
        tb = pl.get_current_fig_manager().toolbar
        if event.inaxes != self.axspl or tb.mode != '': return
        if any_obj_contains(self.knots, event) : return
        x,y = event.xdata, event.ydata
        axrange = self.axspl.get_xbound() + self.axspl.get_ybound()
        if event.button is 1:
            self.pressed = x,y
            x0 = event.xdata
            label = unique_tag(self.knots.keys())
            nux = in_range(self.nu, (x0-self.nuspan, x0+self.nuspan))
            if np.any(nux):
                y = np.mean(self.sp[nux])
                lp =self.axspl.plot(x, y, 'ro', alpha=0.5, label=label)[0]
                self.knots[label] = DraggablePoint(lp, self.knots)
                self.axspl.axis(axrange)
                pl.draw()
            pass
        self.apply_spl_knots()
        self.axspl.axis(axrange)
        self.redraw(smooth_hint = False)
        
    def onmotion_spl(self,event):
        if (self.pressed is None) or (event.inaxes != self.axspl):
            return
        x0 = self.pressed[0]
        x = event.xdata
        if self.curr_span:
            self.curr_span.set_xy([[x0,0], [x0,1], [x,1], [x,0], [x0,0]])
        self.redraw()

    def onrelease_spl(self,event):
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
        self.apply_spl2()
    def addspan(self, (x1,x2)):
        label = unique_tag(self.spans.keys())
        spanh = pl.axvspan(x1,x2,**self.spankw)
        spanh.set_label(label)
        return spanh


    def get_weights(self, nu):
        weights = np.ones(len(nu))
        mask_regs = [s.xspan() for s in self.spans.values()]
        if len(mask_regs) > 0:
            masked = reduce(lambda a,b:a+b,
                            map(lambda x: in_range(nu,x), mask_regs))
            weights -= masked*(1-1e-6)
        return weights
   
    def apply_spl2(self,show=False):
        nu,sp = self.nu, self.sp
        sp_fit, diff = self.process(nu, sp, 'full')
        self.plfit.set_data(nu, sp_fit)
        self.redraw()
        if show:
            plot_with_peaks(nu, diff)
        return

    def knots2pars(self, nu = None, sp = None):
        xlocs = sorted([k.get_xy()[0] for k in self.knots.values()])
        if nu is None : nu = self.nu
        if sp is None : sp = self.sp
        self.nuspan = 5.0
        nuxs = [in_range(nu, (x-self.nuspan, x + self.nuspan)) for x in xlocs]
        y = [np.mean(sp[nx]) for nx in nuxs]
        return xlocs,y

    def update_knots(self):
        for knot in self.knots.values():
            x,y = knot.get_xy()
            nux = in_range(self.nu, (x- self.nuspan, x+self.nuspan))
            if np.any(nux):
                y = np.mean(self.sp[nux])
            knot.set_xy((x,y))
            self.redraw(smooth_hint=False)
    
    def apply_spl_knots(self, show=False):
        nu,sp = self.nu, self.sp
        self.update_knots()
        if len(self.knots.values()) > 5:
            sp_fit, diff = self.process_knots(nu, sp, 'full')
            self.plfit.set_data(nu, sp_fit)
            self.redraw(smooth_hint = False)
            if show:
                plot_with_peaks(nu, diff)


    def export_recipe(self, out=None):
        rec =  {'s': self.spl_smooth,
                'xspans':self.export_xspans()}
        if isinstance(out, str):
            fo = file(out, 'w')
            pickle.dump(rec, fo)
            fo.close()
        return rec
    def load_recipe(self, recipe):
        if isinstance(recipe, dict):
            rec = recipe
        elif isinstance(recipe, str):
            fo = file(recipe,'r')
            rec = pickle.load(fo)
            fo.close()
        else:
            print "Don't recognize type of recipe!"
            print "Doing nothing"
            return
        self.spl_smooth = rec['s']
        self.load_spans(rec['xspans'])
        
    def process(self, nu, sp, ret = None):
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
        
    def export_xspans(self):
        return [s.xspan() for s in self.spans.values()]
    def load_spans(self, xspans):
        for xsp in xspans:
            newspan = self.addspan(xsp)
            label = newspan.get_label()
            self.spans[label] = Span(newspan, self.spans)
        self.apply_spl2()
        

# --- ROIS

class DraggableObj:
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

class LabeledPoint(DraggablePoint):
    "Labeled Point"
    def __init__(self, obj, coll):
        DraggableObj.__init__(self, obj, coll)
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

                              

