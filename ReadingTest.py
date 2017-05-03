import pyfits
import pylab

def plot(filename):
  i = 1
  f = open(filename)
  for l in f:
    g = pyfits.open(l)
    d = g[0].data
    pylab.subplot(5,7,i)
    pylab.imshow(d)
    i = i + 1