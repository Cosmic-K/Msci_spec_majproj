import pylab
import pyfits
import numpy
import minuit
import matplotlib
import math
import itertools
import scipy
from scipy import integrate

# load files and background subtract
# find lines
# fit lines


def determineCuts(spec2d, dCutFactor = 4.0, dAddPix = 10, bPlot = True) :
    
    vproj = spec2d.sum(1)
    dvproj = vproj-pylab.roll(vproj,1)#derivitive
    
    index  = pylab.arange(0,len(vproj),1)
    index1 = index[dvproj > dvproj.max()/dCutFactor]
    start  = index1[0]
    end    = index1[-1]
    
    startWide = start-dAddPix
    endWide   = end+dAddPix 
    
    if bPlot :
        pylab.figure(12)
        pylab.clf()

        pylab.subplot(2,2,1)
        pylab.plot(vproj)

        pylab.subplot(2,2,2)
        pylab.plot(pylab.absolute(dvproj))
        pylab.axhline(dvproj.max()/dCutFactor)

        pylab.subplot(2,2,3)
        pylab.plot(vproj)
        pylab.axvline(start)
        pylab.axvline(end)

        pylab.axvline(startWide)
        pylab.axvline(endWide)

    return [startWide, endWide]

def cutBlankFrame(frame, limits) :
    return frame[limits[0]:limits[1],:]


def backgroundSub(fileSignal = '../spectradata/2013_02_06/sirius/sirius_1_3000_60s.FIT',fileBackground = '../spectradata/2013_02_06/sirius/dark_3000_30s.FIT', bPlot = True) :

    ''' Takes two files and returns a background subtracted profile
        in horizontal axis'''

    fSignal     = pyfits.open(fileSignal)
    fBackground = pyfits.open(fileBackground)
    dSignal     = fSignal[0].data
    dBackground = fBackground[0].data

    # test spectrum cut
    cuts = determineCuts(dSignal,bPlot=False)
    dSignalCut     = cutBlankFrame(dSignal,cuts)
    dBackgroundCut = cutBlankFrame(dBackground,cuts)
    
    # Get this from the header
    h  = fSignal[0].header
    dh = fBackground[0].header
    gain     = h['EGAIN']
    expTime  = h['EXPTIME']
    expTimed = dh['EXPTIME']
    
    #if exposures time do not match change to equate(dark is smaller than data in this case)
    if h['EXPTIME']!= dh['EXPTIME']:
        dSignalCut = (dSignalCut/h['EXPTIME'])*dh['EXPTIME']
    
    # Apply Gain and project
    dBgPE     = gain*dBackgroundCut.sum(0)
    dSignalPE = gain*dSignalCut.sum(0)
    
    # Calculate error
    dBgPEErr    = pylab.sqrt(dBgPE)
    dSignalPEErr= pylab.sqrt(dSignalPE)
    
    
    # Background subtract
    dBgSubCut   = dSignalPE-dBgPE
    dBgSubCutErr= pylab.sqrt(dBgPEErr**2+dSignalPEErr**2) #not sure if i believe this
    x=pylab.arange(0,len(dBgSubCut),1)
    
    if bPlot :
        pylab.figure(10)
        pylab.imshow(dSignal)
        
        pylab.figure(12)
        pylab.clf()
        pylab.errorbar(x,dBgSubCut,dBgSubCutErr)

        pylab.title(fileSignal)
        pylab.xlabel('pixel position')
        pylab.ylabel('no.photoelectrons')

    return [dBgSubCut, dBgSubCutErr, gain]


def findspectralLines(data, cutLevel =0.866, bPlot=True) :#0.02 is good up to cad 4000->0.01
    # moving average 3, smooth data
    data = 1/3.0*(data + pylab.roll(data,1) + pylab.roll(data,2))
    #normalise
    ndata   = data/data.max()
    #where values are abover cut assign value to -1
    #same basic idea as findcallibration except we define the area of interest as below rather than abover cut value
    lines=numpy.where(ndata>cutLevel,-1,ndata)
    coords = []
    
        #find start and end values of each line
    for index, item in enumerate(lines):
        if lines[index-1]==-1 and item !=-1:
            coords.append(index)
        elif  index==764:
            break
        elif lines[index+1]==-1 and item !=-1:
            coords.append(index)

    coordPairs = pylab.reshape(coords,(len(coords)/2,2))

    
    if bPlot:
        pylab.figure(10)
        pylab.clf()
        pylab.plot(ndata,label='3 point moving average, 1/3.0*(data + pylab.roll(data,1) + pylab.roll(data,2))')
        pylab.xlabel('Pixels')
        pylab.ylabel('Normalised "Smoothed" yValues')
        pylab.axhline(cutLevel,color='r',label='Threshold=%s'%(cutLevel))
        pylab.title('Spectral Lines Found')
        pylab.legend(loc='best',prop={'size':8})

        for lineRange in coordPairs :
            pylab.axvline(lineRange[0],color='r',ls='--')
            pylab.axvline(lineRange[1],color='r',ls='--')

        print 'findCalibrationLines> line coordinates :',coords
        print 'findCalibratoinLines> length of list   :',len(coords)
        print 'findCalibrationLines> pairs of coordinates :\n',coordPairs

    return[coordPairs]

def findBasicLineData(data, coordPairs, bPlot = True) :
    xdata =pylab.arange(0,len(data),1)
    
    basicData = []
    for lineRange in coordPairs :
        #print 'findBasicLineData> lineRange', lineRange
        xmean  = (data[lineRange[0]:lineRange[1]]*xdata[lineRange[0]:lineRange[1]]).sum()/(data[lineRange[0]:lineRange[1]]).sum()
        ymax   = data[lineRange[0]:lineRange[1]].min()
        #amplitude set as the minimum value between the defined coords
        xwidth = lineRange[1]-lineRange[0]
        
        basicData.append([xmean,ymax,xwidth])

   
    if bPlot:
        pylab.figure(20)
        pylab.clf()
        pylab.plot(data)
        i=1
        pylab.xlabel('Pixels')
        pylab.ylabel('no. Photoelectrons')
        pylab.title('Guess Parameters')
        for d in basicData :
            a=pylab.axvline(d[0],color='r',ls='--')
            pylab.text(d[0],data.max(),i,fontsize=10)
            b=pylab.axhline(d[1],color='k',ls='--')
            i+=1
        pylab.figlegend((a,b),('$\mu$','Amplitude'),'best',prop={'size':10})
    return pylab.array(basicData)


def CutLines(data, dataErr, coordPairs):
    #setting line widths
    #cutting data and errors

    CutLinesValues = []
    pixelad =89
    #57
    for index, lines in enumerate(coordPairs):
    
        if lines[0]<100:
            start  = lines[0]
        else:
            start  = lines[0]-pixelad
    
        end = lines[1]+pixelad
    
        dataCrop     = data[start:end]
        dataErrCrop  = dataErr[start:end]
        cutdata_x    = pylab.arange(start,end,1)
        CutLinesValues.append([cutdata_x,dataCrop,dataErrCrop])

    #print 'CutLines> CutLinesVaues',CutLinesValues

    return[CutLinesValues]

def circ(xp,R,x):
    if R<math.fabs(x-xp):
        return 0
    else:
        return (1/(math.pi*R**2))*2*pylab.sqrt(R**2-(x-xp)**2)
#circle func

def lor(xp,gamma):
    if gamma <= 0:
        return 0
    else:
        return (1/math.pi)*(gamma)/((xp)**2+gamma**2)


def circLor(xp, R, gamma, x):
        return circ(xp,R,x)*lor(xp,gamma)

def convolutionfuncL(a, mu, gamma, x, R):
    C = []
    for n in x:
        result, error = integrate.quad(circLor, -R*4.5, R*4.5, args = (R, gamma, n-mu))
        C.append(result)
    
    CC = pylab.array(C)
    CCN = CC/CC.max()
    CC = CCN*a
    return pylab.array(CC)
#lorentz

def voiSb(xp, sigma, gamma) :
    z = (xp+1j*gamma)/(sigma*pylab.sqrt(2))
    ff = scipy.special.wofz(z)
    vf = pylab.real(ff)

    #CC = pylab.array(vf)
    #CCNormalised = CC/CC.max()
    #CC = CCNormalised*a

    return pylab.array(vf)

def circvoi(xp, R, sigma, gamma, x):
    return circ(xp,R,x)*voiSb(xp, sigma, gamma)

def convolutionfuncvoi(a, x, mu, sigma, gamma, R):
    Convolve = []
    for n in x:
        result, error = integrate.quad(circvoi, -R*3, R*3, args = (R, sigma, gamma, n-mu))
        Convolve.append(result)

    CC = pylab.array(Convolve)
    CCNormalised = CC/CC.max()
    CC = CCNormalised*a
    return pylab.array(CC)
#voigt
def gauss(xp,sigma):
    if sigma <= 0:
        return 0
    else:
        return (1/(sigma*pylab.sqrt(2*math.pi)))*pylab.exp(-(xp)**2/(2*sigma**2))

def circGauss(xp, R, sigma, x):
    return circ(xp,R,x)*gauss(xp,sigma)

def convolutionfuncG(a, mu, sigma, x, R):
    Convolve = []
    for n in x:
        result, error = integrate.quad(circGauss, -R*5,R*5, args = (R, sigma, n-mu))
        Convolve.append(result)
    
    CC = pylab.array(Convolve)
    CCNormalised = CC/CC.max()
    CC = CCNormalised*a
    return pylab.array(CC)
    #gaussian

def fitLines(coordPairs,lineBasicData, CutLinesValues, bPlot = False) :
    
    fitData = []
    j = 1
    for i in range(0,len(coordPairs),1) :
        print i
        print coordPairs[i]
        print lineBasicData[i]
       
    for Values, basicData in itertools.izip(CutLinesValues[1:],lineBasicData[1:]):
        #chnage slicer values to skip lines which cannot be minimised i.e edge effects or detected false lines
        x   = Values[0]
        f   = Values[1]
        Err = Values[2]

        #chi2 = lambda bg, a, mu, sigma, R: (((f-convolutionfuncG(a, mu, sigma, x, R)+bg)/Err)**2).sum()
        #chi2 = lambda bg,a,mu,sigma: ((f -(a/sigma*pylab.sqrt(2*math.pi))*pylab.exp(-(x-mu)**2/sigma**2)+bg)**2/Err**2).sum() #just gauss
        #chi2 = lambda bg, a, mu, gamma,R: ((f-convolutionfuncL(a, mu, gamma, x, R)+bg)**2/Err**2).sum()	#circleLorentz
        chi2  = lambda bg, a, mu, gamma, sigma, R:(((f-convolutionfuncvoi(a, x, mu, sigma, gamma, R)+bg)/Err)**2).sum()	#circlevoigt
        #chi2  = lambda bg, a, mu, gamma, sigma:(((f-voiSb((x-mu), sigma, gamma,a)+bg)/Err)**2).sum()	#voigt

        #m=minuit.Minuit(chi2, bg=-505864.835856, a = -408608, mu= 666.961236213, gamma=11.1507399367, sigma=0.750795227728, R=4.87)
        #m=minuit.Minuit(chi2, bg=-510099.995856, a = -890000, mu= 668.961236213, gamma=11.1507399367, sigma=0.750795227728)#voigt beta extended
        #m=minuit.Minuit(chi2, bg=-505864.835856, a = -408608, mu= 666.961236213, gamma=11.1507399367, sigma=0.750795227728)#voigt shallow
        #m=minuit.Minuit(chi2, bg = -500000, a = -basicData[1], mu = basicData[0], gamma = basicData[2], sigma = basicData[2], R=4.87)#tester
        #m=minuit.Minuit(chi2, bg = -603807.71721, a = -801908.440558, mu = 389.77133122, gamma = 10.133591096, sigma = 10.161872412, R=4.255771982)
        m=minuit.Minuit(chi2, bg = -494532.214627, a = -275192.356194, mu = 389.755898794, gamma = 9.06574437224, sigma = 2, R=11.0398160531)

        m.tol = 10000000000
        m.printMode = 1
        m.migrad()
        m.hesse()
        m.values
        bg    = m.values['bg']
        a     = m.values['a']
        mu    = m.values['mu']
        gamma = m.values['gamma']
        sigma = m.values['sigma']
        R     = m.values['R']
        #fit =(a/sigma*pylab.sqrt(2*math.pi))*pylab.exp(-((x-mu)**2)/sigma**2)-bg
        #fit =convolutionfuncL(a, mu, gamma, x, R)-bg
        #fit = convolutionfuncG(a, mu, sigma, x, R)-bg
        fit = convolutionfuncvoi(a, x, mu, sigma, gamma,R)-bg
        #fit = voiSb((x-mu), sigma, gamma,a)-bg
        #fitData.append([a,mu,bg,sigma,m.errors])
        fitData.append([a,mu,bg,gamma,sigma,R,m.errors])
        print'fitLines> Covariance Matrix:\n', m.covariance
        

        if bPlot:
            pylab.figure()
            pylab.errorbar(x,f,Err,label='Data')
            pylab.xlabel('Pixels')
            pylab.ylabel('no. Photo electrons')
            pylab.plot(x,fit,label='Fitted Voigt X Circle')
            pylab.legend(loc='best')
            pylab.figure(108)

            pylab.subplot(5,2,j)
            pylab.tight_layout()
            ax = pylab.gca()
            for tick in ax.xaxis.get_major_ticks():
                tick.label1.set_fontsize(9)
            for tick in ax.yaxis.get_major_ticks():
                tick.label1.set_fontsize(9)
            pylab.errorbar(x,f,Err)
            pylab.plot(x,fit,label='line:%s'%(j))
            pylab.legend(loc='best',prop={'size':8})
            j+=1
                        


    fitData = pylab.reshape(fitData,(len(fitData),7))
    #out_txt(fitData,'')
    print 'fitLines> fitData:\n',fitData


    return[fitData]

def Plotlines(data, dataErr, coordPairs, CutLinesValues, montage=True, fastplot=False):
    
    #fast plot lines Turn callibrationlinedata plot to True if want full spectra plot
    #plot subplot of lines
    i = 1

    if fastplot:
        for Values in CutLinesValues:
            pylab.figure()
            pylab.errorbar(Values[0], Values[1], Values[2])
            pylab.xlabel('Pixels')
            pylab.ylabel('no.photoelectrons')
            

    if montage:
        for Values in CutLinesValues:
        
            pylab.figure(71)
            pylab.subplot(5,2,i)
            pylab.plot(Values[0],Values[1])
            pylab.title('subplot,lines')
            i+=1


def out_csv(mydata, filename):
    with open(filename, 'w') as out_handle:
        for line in mydata:
            out_handle.write(','.join(str(line)))

def out_txt(mylist,filename):
    with open(filename,'w') as file:
        for item in mylist:
            print>>file, item


def setting_x_to_nanometer(data,coordPairs,Lambda,dLambda,Lambdax):
    
    startpixel   = (Lambda-(Lambdax*dlambda))
    x_wavelength = pylab.arange(startpixel,dlambda*len(data),1)
    
    return[x_wavelength]


def RUN():
    [data, dataErr, gain]                               = backgroundSub(bPlot=False)
    #subtracts background and returns data and error
    [coordPairs]                                        = findspectralLines(data, bPlot=True)
    #Find spectral lines and returns the coords, works the same a callibration except region of interest has been redefined
    lineBasicData                                       = findBasicLineData(data, coordPairs, bPlot=True)
    #establishes and returns guess parameters 
    [CutLinesValues]                                    = CutLines(data, dataErr, coordPairs)
    #cuts raw data around defined coords and returns these cut lines to be minimised
    plot                                                = Plotlines(data, dataErr, coordPairs, CutLinesValues, montage=False, fastplot=True)
    lineFits                                            = fitLines(coordPairs,lineBasicData, CutLinesValues, bPlot = True)
    #loops through list of cut lines and minimises values to gaussian, can skip lines using slicers.


