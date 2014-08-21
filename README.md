<div id="table-of-contents">
<h2>Table of Contents</h2>
<div id="text-table-of-contents">
<ul>
<li><a href="#sec-1">1. Tutorial</a>
<ul>
<li><a href="#sec-1-1">1.1. Getting started</a>
<ul>
<li><a href="#sec-1-1-1">1.1.1. Load PyRaman</a></li>
</ul>
</li>
<li><a href="#sec-1-2">1.2. Removing baseline from a spectrum</a>
<ul>
<li><a href="#sec-1-2-1">1.2.1. Loading a spectrum</a></li>
<li><a href="#sec-1-2-2">1.2.2. Using <code>RamanCooker</code></a></li>
<li><a href="#sec-1-2-3">1.2.3. Exporting a `recipe'</a></li>
<li><a href="#sec-1-2-4">1.2.4. Loading a `recipe'</a></li>
</ul>
</li>
</ul>
</li>
</ul>
</div>
</div>


# Tutorial<a id="sec-1" name="sec-1"></a>

## Getting started<a id="sec-1-1" name="sec-1-1"></a>

First, we need to start a python shell with loaded `pylab` and `numpy`
libraries. A good way of doing that is to use `ipython` shell. If you're on
Linux, use the following command in your shell:

    ipython -pylab

If you're on Windows, use `PythonXY` distribution.

In the rest of the document it's assumed that `pylab` and `numpy` are
imported and usable.

### Load PyRaman<a id="sec-1-1-1" name="sec-1-1-1"></a>

    1  import raman

## Removing baseline from a spectrum<a id="sec-1-2" name="sec-1-2"></a>

### Loading a spectrum<a id="sec-1-2-1" name="sec-1-2-1"></a>

A spectrum is a (&nu;, *I*) pair, where &nu; is wavelength shift and *I* is intensity
of a signal. Values for &nu; and *I* may be stored separately, e.g. in text files:

    2  nu = np.loadtxt('nu-filename.txt')
    3  sp = np.loadtxt('specturm-values.txt')

As another option, the &nu; and *I* values may be stored together, e.g. when
exported from Renishaw software:

    4  spectrum = raman.rshaw.load_spectrum("filename.txt")

In the latter case the spectrum will be stored in an `Nx2` variable
`spectrum`, first column of which contains &nu; values and the second one
&#x2014; *I* values.

### Using `RamanCooker`<a id="sec-1-2-2" name="sec-1-2-2"></a>

To use `RamanCooker`, first create an instance of its class:

    5  rc1 = raman.RamanCooker()

The basic user interface is realised through the `cook` method of
`RamanCooker` class:

    6  rc1.cook(nu, sp)

To use the keyword arguments:

    7  rc1.cook(nu, sp, s = 200, xspans = [])

Using the variable loaded with `rshaw.load_spectrum()`:

    8  rc1.cook(*spectrum.T)

Note that asterisk `*` is a standard way of arguments expansion in Python
and transposing the spectrum array `spectrum.T` is needed to get the
correct order of arguments.

To shade spectrum regions (generally, containing the peaks) from spline
fitting, right click-drag on the axes in the desired locations. To change
the smoothing parameter, use scrolling. Scrolling up increases smoothness,
scrolling down allows for a more detailed fit.

### Exporting a \`recipe'<a id="sec-1-2-3" name="sec-1-2-3"></a>

1.  To a variable

        9  recipe1 = rc1.export_recipe()

2.  To a file

        10  rc1.export_recipe('Recipe-filename.txt')

### Loading a \`recipe'<a id="sec-1-2-4" name="sec-1-2-4"></a>

Can do that for the same instance, or for a new one:

    11  rc2 = raman.RamanCooker()

1.  From start

        12  rc2.cook(nu,sp,**recipe1)
    
    Note that double asterisk `**` expands the keyword arguments in Python.

2.  From a variable

        13  rc2.load_recipe(recipe1)

3.  From a pre-saved file

        14  rc2.load_recipe('Recipe-filename.txt')
