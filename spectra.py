import pyfits
import pylab

def plot(fileName,fileName1) :
    #double argument input, file and dark-frame 
    f = pyfits.open(fileName)
    g = pyfits.open(fileName1)

    print f[0].data
    print g[0].data
    #print contents to screen
    pylab.figure(1)
    pylab.imshow(f[0].data)
    #show image 
    
    h = f[0].header
    print h['EXPTIME']
    print h['CCD-TEMP']
    #print ccd temp ad exposure t to screen
    
    lf = f[0].data
    df = g[0].data
    lf.sum(0)
    df.sum(0)
    pf=(2.3)*((lf.sum(0))-(df.sum(0)))
    #subtract background
    pylab.figure(2)
    pylab.plot(pf)
