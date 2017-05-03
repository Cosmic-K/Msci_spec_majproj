import pyfits
import pylab
import numpy as np
import minuit

def plot(filename):
  
  data = np.load(filename)
  s=data.sum(0)
  x = pylab.arange(0,len(s),1)
  
  
  
  yerr = pylab.sqrt(s)

  pylab.errorbar(x,s,yerr)

  chi2 = lambda bg, a, mu, sigma: ((s-a*pylab.exp(-(x-mu)**2/sigma**2)-bg)**2/yerr**2).sum()

  print chi2(16500,50000,25,10)

  m = minuit.Minuit(chi2, sigma = 10, bg = 16500, mu = 25, a = 50000)
  m.printMode = 1
  m.migrad()
  m.hesse()
  m.values

  bg = m.values['bg']
  a  = m.values['a']
  mu = m.values['mu']
  sigma = m.values['sigma']

  yfit = a*pylab.exp(-((x-mu)**2)/(sigma**2))+bg

  pylab.plot(x,yfit)
  
  print m.values
  print m.errors
  print m.covariance