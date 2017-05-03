import scipy.integrate as si
import pylab as pl 

def arrayConvolution(func1 = 'gaus', func2 = 'circle', 
                     f1params = [], f2params = [], 
                     x = 0) :
    pass

def integralConvolution(func1 = 'gauss', func2 = 'circle', 
                        f1params= [], f2params = [] ,
                        x = 0) :     
    if func1 == 'gauss' : 
        sigma = f1params[0]
        f1 = lambda xp, sigma : 1/(pl.sqrt(2*pl.pi)*sigma)*pl.exp(-xp**2/sigma**2)
    elif func1 == 'lorentz' : 
        pass
    if func2 == 'circle' : 
        R = f2params[0]
        f2 = lambda xp, R : 2*pl.sqrt(R**2-xp**2)/(pl.pi*R**2) if abs(xp)<=R else 0.0    
    elif func2 == 'gauss' : 
        pass

    fc = lambda xp : f1(x-xp,sigma)*f2(xp,R)

    if func2 == 'circle' : 
        return si.quad(fc,-R ,R,limit=200)[0]
    elif func2 == 'gauss' : 
        pass

def plotIntegralConvolution(s = 1.0, r = 1.0) :
    x = pl.linspace(-10,10,1000)
    y = []

    for xx in x :
        yy = integralConvolution('gauss','circle',[s],[r], xx)
        y.append(yy)
#        print xx,yy
        
    y = pl.array(y)

    pl.plot(x,y)

def linearRes(p,x,y,yerr = pl.array([])) : 
    if len(yerr) == 0 : 
        res = y-linear(p,x)
    else :
        res = (y-linear(p,x))/yerr
    return res

def linear(p,x) :
    y   = p[0]*x+p[1]
    return y


