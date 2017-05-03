import pyfits
import pylab

def plot(filename)

  f = open(filename)

  for l in f:
    g = pyfits.open(l)
    d = g[0].data
    s = d.sum(0)

    pylab.figure(1)
    pylab.imshow(d)
    title(l)
    
    pylab.figure(2)
    pylab.plot(s)
    title(l)