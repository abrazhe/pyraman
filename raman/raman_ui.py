import wx
import matplotlib
matplotlib.use('WxAgg')

import pylab as pl

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg  as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.backends.backend_wx import NavigationToolbar2Wx as NToolbar

from traits.api import *
from traitsui.api import *


from traitsui.menu \
     import Action, CloseAction, Menu, MenuBar, OKCancelButtons, Separator

from traitsui.wx.editor import Editor
from traitsui.basic_editor_factory import BasicEditorFactory

class _MPLFigureEditor(Editor):
    scrollable  = True
    def init(self, parent):
        self.control = self._create_canvas(parent)
        self.set_tooltip()
    def update_editor(self):
        pass
    def _create_canvas(self, parent):
        """ Create the MPL canvas. """
        # The panel lets us add additional controls.
        panel = wx.Panel(parent, -1, style=wx.CLIP_CHILDREN)
        sizer = wx.BoxSizer(wx.VERTICAL)
        panel.SetSizer(sizer)
        # matplotlib commands to create a canvas
        mpl_control = FigureCanvas(panel, -1, self.value)
        sizer.Add(mpl_control, 1, wx.LEFT | wx.TOP | wx.GROW)
        toolbar = NToolbar(mpl_control)
        sizer.Add(toolbar, 0, wx.EXPAND)
        self.value.canvas.SetMinSize((10,10))
        return panel


class MPLFigureEditor(BasicEditorFactory):
    klass = _MPLFigureEditor

from .baselines import als

class ALS_traits(HasTraits):
    figure = Instance(Figure, ())
    smoother = Float()
    p = Float()
    view = View(VGroup(Item('figure', editor=MPLFigureEditor(),
                            show_label=False),
                       HGroup(Item('smoother'), Item('p'))),
                width=800,
                height=500,
                title = 'Baseline with ALS',
                resizable=True)
    def __init__(self, nu, sp,smooth=2.0, p=0.01):
        self.nu = nu
        self.sp = sp
        self.smoother = smooth
        print (smooth, self.smoother)
        self.p = p
    def _figure_default(self):
        figure = Figure()
        self.axes = figure.add_subplot(111)
        self.axes.set_title('ALS baseline coorection')
        self._update_baseline()
        self.pl1 = self.axes.plot(self.nu, self.sp)
        self.pl2 = self.axes.plot(self.nu,self.baseline, 'm--')[0]
        return figure
    def _update_baseline(self):
        self.baseline = als(self.sp, 10**self.smoother, self.p)
        if hasattr(self, 'pl2'):
            self.pl2.set_data(self.nu,self.baseline)
            self.figure.canvas.draw()
    def _p_changed(self):
        self._update_baseline()
    def _smoother_changed(self):
        self._update_baseline()

        
def main(nu, sp):
    Raman_als(nu,sp).configure_traits()
