import pyfits
import pylab

def plot(filename):

  f = open(filename)

  for l in f:
    g = pyfits.open(l)
    d = g[0].data
    s = d.sum(0)

    pylab.figure()
    pylab.imshow(d)
    pylab.title(l)

    pylab.figure()
    pylab.plot(s)
    pylab.title(l)

    