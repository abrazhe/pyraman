from pylab import *

verbose = 1

def defined(obj):
    return obj is not None

def provided(argument):
    return defined(argument)

def best (scoref, lst):
    n,winner = 0, lst[0]
    for i, item in enumerate(lst):
        if  scoref(item, winner): n, winner = i, item
    return n,winner

def min1(scoref, lst):
    return best(less, map(scoref, lst))

def in_range(n, region):
    return (n > region[0]) and (n < region[1])

def nearest_item_ind(items, x, fn = lambda a: a):
    return min1(lambda p: abs(fn(p) - x), items)[0]

def cpoints_cmp(one, other):
    diff = one.x - other.x
    if diff < 0: return -1
    elif diff > 0: return 1
    else: return 0

def regions_from_cpoints(cpoints):
    """Takes ordered list of coloured points and returnes
    a list of regions
    """
    regions = []
    for n, p in enumerate(cpoints[:-1]):
        regions.append((p.x, cpoints[n+1].x))
    return regions

def paired_regions_from_cpoints(cpoints):
    return  [map(lambda p: p.x, pair)
             for pair in group(cpoints,2)
             if len(pair) > 1]

def group(lst, n):
    if n == 0: return lst
    acc = []
    n = int(n)
    if lst:
        for obj in arange(len(lst))[::n]:
            acc.append(lst[obj : obj+n])
    return acc


def slots_from_cpoints2(cpoints):
    """makes range slots from sorted cpoints"""
    state = 'zero'
    slots = []
    curr_slot_dict = {'blue':[], 'start':0, 'stop':0}
    for point in cpoints:
        if point.color is 'r': 
            if state is not 'zero':
                curr_slot_dict['stop'] = point
                slots.append(Slot(curr_slot_dict))
                curr_slot_dict['blue'] = []
                pass
            curr_slot_dict['start'] = point
            state = 'appending'
            pass
        elif point.color is 'b': 
            if state is not 'zero':
                curr_slot_dict['blue'].append(point)
    return slots


def print_point(ev):
    print "You clicked: %f, %f" % (ev.xdata, ev.ydata )






class ColoredPoint:
    def __init__(self, color, x, plh):
        self.color = color
        self.x = x
        self.plh = plh
    def clear(self):
        setp(self.plh, 'Visible', False, 'data', (None, None))
    def __str__(self):
        return str( self.color + ':' + str(self.x) )
    def __cmp__(self, other):
        return cpoints_cmp(self, other)
        



class Slot:
    def __init__(self, slot_dict, order=3):
        self.order = order
        self.start = slot_dict['start'].x
        self.stop = slot_dict['stop'].x
        self.blue_regions = []
        if slot_dict.has_key('blue'):
            self.blue_regions = paired_regions_from_cpoints(slot_dict['blue'])
            print self.blue_regions
        self.vector = self.make_vector()
    def make_vector(self):
        vector = range(self.start, self.stop)
        if len(self.blue_regions)>0:
            for br in self.blue_regions:
                vector = [i for i in vector if i < br[0] or i > br[1]]
                pass
        return vector

    def fit(self, nu, data):
        x2fit = nu[self.vector]
        data2fit = data[self.vector]
        pc = polyfit(x2fit, data2fit, self.order)
        self.yfitted = polyval(pc, nu[self.start:self.stop])
        return self.yfitted

    def in_regionp(self,x):
        return (x > self.start) and (x < self.stop)

    def set_order(self, new_order):
        if verbose:
            print "Changing polynomial order from ", self.order, 'to', new_order
        self.order = new_order

    def printv(self):
        print self.vector

    def __str__(self):
        str(self.vector)
                


class SlotsCollection:
    def __init__(self):
        self.red_points = []
        self.blue_points = []
        self.points = []
        pass
    def push_cpoint(self, cpoint):
        self.points.append(cpoint)
        if len(self.points) > 1:
            self.points.sort(cpoints_cmp)
            self.slots = slots_from_cpoints2(self.points)
    def first(self):
        if len(self.slots)> 0:
            return self.slots[0].start
        else:
            return 0
    def last(self):
        if len(self.slots) > 0:
            return self.slots[-1].stop
        else: return -1

    def apply(self, nu, data):
        y_fits = []
        for slot in self.slots:
            y_fits.append(slot.fit(nu,data))
        nufit = nu[self.first():self.last()]
        total_y_fit = [i for i in flatten(y_fits)]
        return [nufit, total_y_fit]
        pass
    def pop_cpoint(self, x):
        #n = nearest_cpoint_ind(self.points,x)
        n = nearest_item_ind(self.points, x, lambda a: a.x)
        cpoint = self.points.pop(n)
        cpoint.clear()
        self.slots = slots_from_cpoints2(self.points)
        pass

    def set_order(self, order, x):
        "Sets order of polynomial for fitting"
        for slot in self.slots:
            if slot.in_regionp(x):
                slot.set_order(order)
                break

class BasicRaman:
    def __init__(self, data = None, data2d=None):
        if provided(data):
            if verbose: print "Provided with data"
            self.take_data(data[0], data[1])
        if provided(data2d):
            if verbose: print "Provided with data2d"
            self.take_data_2d(data2d)

        self.macophilf = False
        self.pointcolors = {1:'r', 2:'b', 3:'b'}

        pass
    
    def take_data(self,nu,spectrum):
        N = min(len(nu), len(spectrum))
        self.pre_spectr = spectrum[:N]
        self.nu = nu[:N]

    def take_from_files(self, nu_file, spectrum_file):
        nu = load(nu_file)
        spectrum = load(spectrum_file)
        self.take_data(nu, spectrum)

    def give_results(self):
        return self.new_nu,self.new_spectrum
        
    def take_data_2d(self,data):
        self.nu = data[:,0]
        self.pre_spectr = data[:,1]


class RamanCooker(BasicRaman):
    def start_cooking(self):
        self.slots_collection = SlotsCollection()
        self.plots = {}
        self.events_handlers = {}
        self.read_mouse_eventsp = True
        self.new_nu, self.new_spectrum = None,None

        figure(1000); hold(True)
        self.main_axes = gca()
        main_plot = plot(self.nu, self.pre_spectr, 'k-')
        self.events_handlers['mouse'] = connect('button_press_event', self.click)
        self.events_handlers['keys'] = connect('key_press_event', self.type)
        #self.plots['main'] = main_plot
        
    def type(self, event):
        print "you typed:", event.key

        if event.key == 't':
            self.read_mouse_eventsp = not self.read_mouse_eventsp

        elif event.key == '-':
            #p = nearest_point_ind(self.nu, event.xdata)
            p = nearest_item_ind(self.nu, event.xdata)
            self.slots_collection.pop_cpoint(p)
            print "Xdata:", event.xdata

        elif event.key in '1234567890':
            order = int(event.key)
            p = nearest_item_ind(self.nu, event.xdata)
            self.slots_collection.set_order(order,p)

        elif event.key is 'm':
            self.macophilf = not self.macophilf

        elif event.key == 'a':
            xfit,yfit = self.slots_collection.apply(self.nu, self.pre_spectr)
            if self.plots.has_key('fits'):
                setp(self.plots['fits'], 'data', (None, None), 'visible', False)
            if self.macophilf:
                self.plots['fits'] = plot(xfit, yfit-0.5*mean(self.pre_spectr), 'g--')
            else:
                self.plots['fits'] = plot(xfit,yfit,'m--')
                axis((self.nu[0], self.nu[-1],
                      min(self.pre_spectr), max(self.pre_spectr)))

            self.new_nu = xfit
            self.new_spectrum = self.pre_spectr[self.slots_collection.first():
                                                self.slots_collection.last()] - yfit

            print "New lengths:", len(self.new_nu), len(yfit)

            figure();
            plot(self.new_nu, self.new_spectrum,'b-')
            connect("button_press_event", print_point)
            figure(1000)
            
    def click(self, event):
        if defined(event.inaxes) and self.read_mouse_eventsp:
            print "In the main axes?", event.inaxes == self.main_axes
            p = nearest_item_ind(self.nu, event.xdata)
            if verbose:
                print "you clicked:", event.xdata, event.ydata
                print "Most close nu:", p, self.nu[p]
            color = self.pointcolors[event.button]
            plh = plot([self.nu[p]],[self.pre_spectr[p]], 'o' + color)
            self.slots_collection.push_cpoint(ColoredPoint(color, p, plh))
                    

import pywt


def ith_details(coeffs, i, w):
    N = len(coeffs)
    j = N - i
    return pywt.waverec([None, coeffs[j]] + [None]*j, w)

def ith_averages(coeffs, i, w):
    N = len(coeffs)
    print "N = ", N
    j = N - i
    return pywt.waverec([coeffs[j], None] + [None]*j, w)

class RamanWavelets(BasicRaman):
    def new_spectr(self, coeffs):
        return pywt.waverec(coeffs, self.w)

    def start(self, w = 'db1'):

        self.w = pywt.Wavelet(w)

        N = pywt.dwt_max_level(data_len = len(self.pre_spectr), 
                          filter_len = self.w.dec_len)

        print N

        coeffs = pywt.wavedec(self.pre_spectr, w, level=N)
        rec_d = []
        

        #figure(2000); #hold(True)
        self.main_fig = figure(); 
        ax_main = subplot(N+2,1,1)
        ax_main.plot(self.nu, self.pre_spectr,'k-')
        xlim(self.nu[0], self.nu[-1])

        print len(rec_d)
        
        self.ax_details = {}

        for i in xrange(1,N+1):
            ax = subplot(N+2, 1, i+1)
            y = ith_details(coeffs, i, w)
            plh = ax.plot(y, 'g-')
            connect('scroll_event', lambda event : self.click(event, coeffs) )
            self.ax_details[ax] = [i, plh]
            xlim(0, len(y)-1)
            ylabel('D %d' % i)
        
        ax_aver = subplot(N+2,1,N+2)
        plh = ax_aver.plot(ith_averages(coeffs,N+1, w), 'b-')
        self.ax_details[ax_aver] = [N + 1, plh]
        ylabel('A %d' % (N+1))
        ax_aver.draw()

        figure(2000)
        hold(False)
        self.reconstr_axes = gca()
        #hold(True)
        #plot(self.nu,self.pre_spectr,'k-')
        self.pl_spec_rec = plot(self.nu, self.new_spectr(coeffs), 'b-')
        
    
    def click(self, event, coeffs):
        self.up_factor = 0.7**-0.25
        self.down_factor = 0.7**0.25
        scale_factor = 1
        if defined(event.inaxes):
            curr_axes = event.inaxes
            N = len(coeffs)
            i = self.ax_details[curr_axes][0]
            plh = self.ax_details[curr_axes][1]
            if event.button == 'up':
                scale_factor = self.up_factor
            elif event.button == 'down':
                scale_factor = self.down_factor
            #print N,i, scale_factor
            
            coef_prev = coeffs[N - i]
            coeffs[N - i] = [k * scale_factor for k in coef_prev]
            #print array(coef_prev)/array(coeffs[N-i])
            #print coef_prev
            #print '----------------'
            y = []
            if i <= N:
                y = ith_details(coeffs, i, self.w)
            else:
                y = ith_averages(coeffs, i, self.w)
                #print y
            setp(plh, 'ydata', y, 'color', 'm')
            setp(curr_axes, 'ylim', [min(y), max(y)])
            curr_axes.draw()
            figure(2000)
            #plot(self.nu, self.new_spectr(coeffs))
            setp(self.pl_spec_rec, 'ydata', self.new_spectr(coeffs))

            


import random
def make_evolve(genotype,
                 mutate_gene_fn,
                 score_genotype_fn,
                 population_size=512,
                 mutation_rate = 0.1):
    """Basic evolutionary optimization.
    genotype: a list of genes
    mutate_gene_fn: a function to mutate a gene
    score_genotype_fn: a function to score a genotype"""

    new_population = range(population_size)
    def mutate(gene):
        if random.random() < mutation_rate:
            gene = mutate_gene_fn(gene)
        return gene

    def new_generation(genotype):
        for i in xrange(population_size):  # iteration here only for the sake of speed
            new_population[i] = (map(mutate, genotype))
        return new_population

    def fittest(population):
        x = best(less, map(score_genotype_fn, population))
        return population[x[0]], x[1]

    def update(genotype):
        return fittest(new_generation(genotype))
    return update

def mutate_number(number):
    return abs(random.gauss(number, 0.5*number))

def mutate_number_uniform(n):
    if n == 0: n = 1
    return random.uniform(0, 2*n)

def lorentzian(x, (intensity, mean, gamma)):
    return intensity/(1. + ((x - mean)/gamma)**2)

def gaussian(x, (a, b, c)):
    return a * exp(-(x-b)**2/(2*c**2))

def lorentz2str(p):
       return "{'I': %e, 'x0': %3.3f, 'g': %3.3f}" % tuple(p)

def gauss2str(p):
    return "{'a' %e, 'b': %3.3f, 'c': %3.3f}" %tuple(p)

def lorentzians2str(params):
    strs = map(lorentz2str,group(params,3))
    return '[' + str.join(', ', strs)

def curves(xdata, func, params, n=3):
    return [func(xdata,p) for p in group(params,n)]

def lorentzians(xdata, params_vect):
    return [lorentzian(xdata,p) for p in group(params_vect,3)]

def squared_deviations(ydata, yfit):
    return sum((ydata-yfit)**2)

def score_curve(xdata, ydata, func, params_vect):
    return squared_deviations(ydata,
                              reduce(add,
                                     curves(xdata, func, params_vect),
                                     0))

def score_lorentzians(xdata,ydata,params_vect):
    return squared_deviations(ydata,
                              reduce(add,
                                     lorentzians(xdata, params_vect),
                                     0))

import os,sys
def curve_fitter(xdata,
                 ydata,
                 ranges = None,
                 func = lorentzian,
                 nfunc = 6,
                 tolerance = 1e-7,
                 max_iter = 1e6):

    offrange_penalty = 1000
    def score_genotype(g):
        score = score_curve(xdata,ydata,func,g)
        if defined(ranges):
            for i,p in enumerate(group(g,3)):
                if not in_range(p[1], ranges[i]):
                    score += offrange_penalty
        return score*1000
        
    L = len(xdata)
    if ranges is None:
        means = linspace(xdata[int(L/6)], xdata[-int(L/6)],nfunc)
    else:
        means = []
        nfunc = len(ranges)
        for irange in ranges:
            means.append(mean(irange))
    xfit = linspace(xdata[0], xdata[-1], len(xdata)*10)
    start_i = mean(ydata)
    start_gamma = 4.0
    genotype = []
    for i in xrange(nfunc):
        genotype.append(start_i)
        genotype.append(means[i])
        genotype.append(start_gamma)


    single_fits = curves(xfit, func, genotype)

    yfit = reduce(add,single_fits)

    figure(); grid(); hold(True); axh = gca()

    plh_data = plot(xdata,ydata,'b.')
    plh_fit = plot(xfit,yfit,'r-')

    plh_single_fits = [plot(xfit,yfit,'g-') for yfit in single_fits]
    draw()

    my_evolve = make_evolve(genotype,
                            mutate_number_uniform,
                            score_genotype,
                            population_size=512)
    old_score  = 1e5
    diff  = tolerance*1000

    n_draw = 10
    n_genotype_print = 100
    for n in xrange(int(max_iter)):
        genotype, score = my_evolve(genotype)

        if not (n+1)%n_draw:
            single_fits = curves(xfit, func, genotype)
            yfit = reduce(add,single_fits)

            setp(plh_fit, 'ydata', yfit)
            for i, plh_i in enumerate(plh_single_fits):
                setp(plh_i, 'ydata', single_fits[i])
            draw()

            sys.stderr.write("\rN: %d, Score: %e, Score diff: %e"
                             %(n+1, score, score-old_score))

        if not (n+1)%(n_genotype_print):
            print "\nGenotype:", lorentzians2str(genotype)
                
        old_score = score
        if score < tolerance:
            print "\nNice score:", score
            break

    return genotype
        


def test_fit():
    x = load("nu_part.txt")[10:]
    data =load("rbc_ph7_mean_norm_part.txt")[10:]

    if min(data) < 0:
        data -= min(data)

    """
    Variant 1:
    1545-1552, 1560-1565, 1582-1588, 1620-1626, 1635-1640
    """
    ranges1 = ((1545, 1552),
               (1560, 1565),
               (1582, 1588),
               (1620, 1626),
               (1635, 1640))
    ranges2 = ((1540, 1552), (1560, 1565), (1582, 1588), (1620, 1626), (1635, 1640))

    """
    Variant 2:
    1545-1647, 1550-1553, 1582-1583, 1586-1588,  1620-1626, 1638-1640
    """
    ranges3 = ((1544, 1547),
               (1550, 1555),
               (1582, 1583),
               (1586, 1588),
               (1620, 1626),
               (1638, 1640))

    return lorentz_fitter(x,data,
                          ranges = ranges2,
                          max_iter=200)
