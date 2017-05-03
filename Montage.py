import pyfits
import pylab

def plot(filename,x,y):

  f = open(filename)
  i = 1

  for l in f:
    g = pyfits.open(l)
    d = g[0].data
    s = d.sum(0)
    

    pylab.figure(1)
    pylab.subplot(x,y,i)
    pylab.imshow(d)
    

    pylab.figure(2)
    pylab.subplot(x,y,i)
    pylab.plot(s)
    i = i + 1