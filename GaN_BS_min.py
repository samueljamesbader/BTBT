# -*- coding: utf-8 -*-
"""
Author: Sam Bader
Module: 
"""
import numpy as np
import scipy.linalg as la
from Constants import eps0,q,hbar,m0,kT,eVToJoules
import matplotlib.pyplot as mpl
import scipy.interpolate as sciint


def GaN_imagBS():
    mepar=0.2*m0
    print "REDEFINING ELECTRON MASS"
    mepar=0.189*m0
    meperp=0.18*m0
    A1=-6.98
    A2=-0.56
    A3=6.42
    A4=-3.21
    A5=-2.9
    A6=-3.66
    A7=0.0 /1e8 #eV cm
    Eppar=17.292 #eV
    Epperp=16.265
    P1=0;P2=0
    P1=np.sqrt(hbar**2 *q*Eppar  /(2*m0)) * 100 # eV cm
    P2=np.sqrt(hbar**2 *q*Epperp /(2*m0)) * 100 # eV cm
    Eg=3.44-.034
    DCR=.034
    
    Ap1=q*hbar**2*10000/2*(1/mepar -1/m0) - P1**2/Eg # eV cm^2
#    print Ap1;stop
    Ap2=q*hbar**2*10000/2*(1/meperp-1/m0) - P2**2/Eg # eV cm^2
    Lp1=q*hbar**2*10000/(2*m0)*(A2+A4+A5-1) + P1**2/Eg # eV cm^2
    Lp2=q*hbar**2*10000/(2*m0)*(A1-1) + P2**2/Eg # eV cm^2
    M1 =q*hbar**2*10000/(2*m0)*(A2+A4-A5-1) # eV cm^2
    M2 =q*hbar**2*10000/(2*m0)*(A1+A3-1) # eV cm^2
    M3 =q*hbar**2*10000/(2*m0)*(A2-1) # eV cm^2
    Np1=q*hbar**2*10000/(2*m0)*(2*A5) + P1**2/Eg# eV cm^2
    Np2=q*hbar**2*10000/(2*m0)*(np.sqrt(2)*A6) + P1*P2/Eg# eV cm^2
    Np3=1j*np.sqrt(2)*A7
    
    
    Egp=Eg+DCR
    me=1/(1/m0+2/hbar**2*Ap1/(10000*q))
    mh=1/(1/m0+2/hbar**2*Lp2/(10000*q))
    mr=1/(1/me+1/mh)
    
    a=lambda kp: 1+4.6397e-16*kp**2+1.5231e-9*kp
    b=lambda kp: 1.9663e-8*kp**2+.009221*kp
    kvl=lambda l: np.sqrt(((Egp*hbar**2*q*10000/(2*mh)-P1**2-l*hbar**2*q*10000/(2*mr))\
        +np.sqrt((Egp*hbar**2*q*10000/(2*mh)-P1**2-l*hbar**2*q*10000/(2*mr))**2+hbar**4*q**2*10000**2*l/(me*mh)*(Egp-l)))/
        (hbar**4*q**2*10000**2/(2*me*mh)))
    
    return lambda l,kp: np.nan_to_num(kvl((l)/a(kp))+b(kp))
    
kvl=GaN_imagBS()

def GaN_imagBS2():
    mepar=0.2*m0
    print "REDEFINING ELECTRON MASS"
    mepar=0.189*m0
    meperp=0.18*m0
    A1=-6.98
    A2=-0.56
    A3=6.42
    A4=-3.21
    A5=-2.9
    A6=-3.66
    A7=0.0 /1e8 #eV cm
    Eppar=17.292 #eV
    Epperp=16.265
    P1=np.sqrt(hbar**2 *q*Eppar  /(2*m0)) * 100 # eV cm
    P2=np.sqrt(hbar**2 *q*Epperp /(2*m0)) * 100 # eV cm
    Eg=3.44-.034
    DCR=.034
    Egp=Eg+DCR
    
    Ap1=q*hbar**2*10000/2*(1/mepar -1/m0) - P1**2/Eg # eV cm^2
    Ap2=q*hbar**2*10000/2*(1/meperp-1/m0) - P2**2/Eg # eV cm^2
    Lp1=q*hbar**2*10000/(2*m0)*(A2+A4+A5-1) + P1**2/Eg # eV cm^2
    Lp2=q*hbar**2*10000/(2*m0)*(A1-1) + P2**2/Eg # eV cm^2
    M1 =q*hbar**2*10000/(2*m0)*(A2+A4-A5-1) # eV cm^2
    M2 =q*hbar**2*10000/(2*m0)*(A1+A3-1) # eV cm^2
    M3 =q*hbar**2*10000/(2*m0)*(A2-1) # eV cm^2
    Np1=q*hbar**2*10000/(2*m0)*(2*A5) + P1**2/Eg# eV cm^2
    Np2=q*hbar**2*10000/(2*m0)*(np.sqrt(2)*A6) + P1*P2/Eg# eV cm^2
    Np3=1j*np.sqrt(2)*A7
    
    
    m0z=1/(1/m0+2/hbar**2*Ap1/(10000*q))
    m3z=1/(1/m0+2/hbar**2*Lp2/(10000*q))
    m0x=1/(1/m0+2/hbar**2*Ap2/(10000*q))
    m3x=1/(1/m0+2/hbar**2*M3/(10000*q))
    print "Ap1", Ap1
    print "m0z", m0z/m0
    
    mez=1/((2*P1**2/(q*10000*hbar**2*Egp))+1/m0z)
    mSOz=1/((2*P1**2/(q*10000*hbar**2*Egp))-1/m3z)
    mex=1/((2*P2**2/(q*10000*hbar**2*Eg))+1/m0x)
    mSOx=-m3x
    print "mez", mez/m0
    print "mSOz", mSOz/m0
    mu2=m0z*(-m3z)
    
    a=lambda kp:1
    b=lambda kp:0
    akz2=lambda l: (mu2/(q*10000*hbar**2))*(((Egp-l)/mSOz+l/mez)-np.sqrt(((Egp-l)/mSOz+l/mez)**2-(4/mu2)*(Egp-l)*l))
#    akz2=lambda l: (((Egp-l)/mSOz+l/mez))
    kvl=lambda l: np.sqrt(akz2(l))
#    kvl=lambda l: np.sqrt(((Egp*hbar**2*q*10000/(2*mh)-P1**2-l*hbar**2*q*10000/(2*mr))\
#        +np.sqrt((Egp*hbar**2*q*10000/(2*mh)-P1**2-l*hbar**2*q*10000/(2*mr))**2+hbar**4*q**2*10000**2*l/(me*mh)*(Egp-l)))/
#        (hbar**4*q**2*10000**2/(2*me*mh)))
    alpha=lambda l: (l/Egp)*(mez/mex-mSOz/mSOx)+mSOz/mSOx
    return lambda l,kp: np.nan_to_num(np.sqrt(akz2(l)+alpha(l)*(kp)**2))
kvl_2=GaN_imagBS2()



if __name__ == "__main__":
    mpl.close('all')
    import matplotlib.pyplot as mpl
    E=np.linspace(-.5,4,100)
    mpl.plot(E,kvl(E,0)/(np.pi/5.189e-8),'b.')
    mpl.plot(E,kvl_2(E,0)/(np.pi/5.189e-8),'g-')
    mpl.plot(E,kvl(E,(.1*np.pi/5.189e-8))/(np.pi/5.189e-8),'b.')
    mpl.plot(E,kvl_2(E,(.1*np.pi/5.189e-8))/(np.pi/5.189e-8),'g-')
    mpl.xlim(-.5,4)
#    mpl.ylim(0,1)
