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
mpl.close('all')


a=3.19e-8 #cm
c=5.189e-8 #cm


#######
###### RINKE
#######
mepar=0.186*m0
meperp=0.209*m0
A1=-5.947
A2=-0.528
A3=5.414
A4=-2.512
A5=-2.510
A6=-3.202
A7=.046 /1e8 #eV cm
Eppar=17.292 #eV
Epperp=16.265
P1=np.sqrt(hbar**2 *q*Eppar  /(2*m0)) * 100 # eV cm
P2=np.sqrt(hbar**2 *q*Epperp /(2*m0)) * 100 # eV cm
Eg=3.24
DCR=.034

#######
###### ME
#######
if 1:
    mepar=0.2*m0
#    print "REDEFINING ELECTRON MASS"
#    mepar=0.189*m0
    
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
print Ap1;stop
Ap2=q*hbar**2*10000/2*(1/meperp-1/m0) - P2**2/Eg # eV cm^2
Lp1=q*hbar**2*10000/(2*m0)*(A2+A4+A5-1) + P1**2/Eg # eV cm^2
Lp2=q*hbar**2*10000/(2*m0)*(A1-1) + P2**2/Eg # eV cm^2
M1 =q*hbar**2*10000/(2*m0)*(A2+A4-A5-1) # eV cm^2
M2 =q*hbar**2*10000/(2*m0)*(A1+A3-1) # eV cm^2
M3 =q*hbar**2*10000/(2*m0)*(A2-1) # eV cm^2
Np1=q*hbar**2*10000/(2*m0)*(2*A5) + P1**2/Eg# eV cm^2
Np2=q*hbar**2*10000/(2*m0)*(np.sqrt(2)*A6) + P1*P2/Eg# eV cm^2
Np3=1j*np.sqrt(2)*A7



def bs(k,bad=False):
    #print k
    kx,ky,kz=k
    k=np.sqrt(kx**2+ky**2+kz**2)
    H1= np.array([ [Eg+DCR+q*hbar**2*10000*k**2/(2*m0), 1j*P2*kx, 1j*P2*ky, 1j*P1*kz],
               [-1j*P2*kx, DCR+q*hbar**2*10000*k**2/(2*m0), 0, 0],
               [-1j*P2*ky,0, DCR+q*hbar**2*10000*k**2/(2*m0),0],
               [-1j*P1*kz, 0, 0, q*hbar**2*10000*k**2/(2*m0)]])
    H2= np.array([[Ap2*(kx**2+ky**2)+Ap1*kz**2,0,0,0],
                  [0,Lp1*kx**2+M1*ky**2+M2*kz**2,Np1*kx*ky,Np2*kx*kz-Np3*kx],
                  [0,Np1*ky*kx,M1*kx**2+Lp1*ky**2+M2*kz**2,Np2*ky*kz-Np3*ky],
                  [0,Np2*kz*kx+Np3*kx,Np2*kz*ky+Np3*ky,M3*(kx**2+ky**2)+Lp2*kz**2]])
    
    if bad:
        Np2p=Np2*0
        Np3p=Np3*0
        P2p=P2
        
        H1= np.array([ [Eg+DCR+q*hbar**2*10000*k**2/(2*m0), 1j*P2p*kx, 0, 1j*P1*kz],
                   [-1j*P2p*kx, DCR+q*hbar**2*10000*k**2/(2*m0), 0, 0],
                   [0,0, DCR+q*hbar**2*10000*k**2/(2*m0),0],
                   [-1j*P1*kz, 0, 0, q*hbar**2*10000*k**2/(2*m0)]])
        H2= np.array([[Ap2*(kx**2)+Ap1*kz**2,0,0,0],
                      [0,Lp1*kx**2+M2*kz**2,0,Np2p*kx*kz-Np3p*kx],
                      [0,0,M1*kx**2+M2*kz**2,0],
                      [0,Np2p*kz*kx+Np3p*kx,0,M3*kx**2+Lp2*kz**2]])
    H=H1+H2
    w=la.eigvals(H)
    w=np.sort(np.choose(abs(w.imag)<1e-8,[w.real*np.NaN,w.real]))
    
    return w

def bs2bm(k):
    kx,ky,kz=k
    k=np.sqrt(kx**2+ky**2+kz**2)
    mr=.0799*m0
    P=(hbar/m0)*m0/2*np.sqrt(q*(Eg+DCR)/mr)*100
    H1= np.array([ [Eg+DCR+q*hbar**2*10000*k**2/(2*m0), 0, 0, 1j*P*kz],
               [0, DCR+q*hbar**2*10000*k**2/(2*m0), 0, 0],
               [0,0, DCR+q*hbar**2*10000*k**2/(2*m0),0],
               [-1j*P*kz, 0, 0, q*hbar**2*10000*k**2/(2*m0)]])
    H2= np.array([[0,0,0,0],
                  [0,M2*kz**2,0,0],
                  [0,0,M2*kz**2,0],
                  [0,0,0,0]])
    H=H1+H2
    w=la.eigvals(H)
    w=np.sort(np.choose(abs(w.imag)<1e-8,[w.real*np.NaN,w.real]))
    
    return w


n=200
def path(pts):
    return np.concatenate(
        [s+(e-s)*np.reshape(np.linspace(0,1,n-1,endpoint=False),[n-1,1])
            for s,e in zip(pts[:-1],pts[1:])]
        +[[pts[-1]]])

GammaPoint=np.array([0,0,0])
KPoint=np.array([4*np.pi/(3*a),0,0])
APoint=np.array([0,0,np.pi/c])



if 1:
    figure()
    karr=path([KPoint/8,GammaPoint,APoint/4])
    b=np.array([bs(k) for k in karr])
    b2=np.array([bs2bm(k) for k in karr])
    subplot(211)
    plot(b[:,3]-np.min(b[:,3]),'k.')
    plot(b2[:,3]-np.min(b2[:,3]),'r')
    xlim([0,2*n-2])
    xticks([0,n-1,2*n-2],['$K/8$','$\Gamma$','$A/4$'])
    ylim([0,.45])
    yticks([0,.1,.2,.3,.4],['0','','.2','','.4'])
    mpl.axvline(n-1,color='k')
    ylabel('Energy above $E_C$ [eV]')
    subplot(212)
    plot(np.arange(2*n-1),b[:,0:3]-np.max(b[:,0:3]),'k.')
    plot(np.arange(2*n-1),b2[:,0:3]-np.max(b2[:,0:3]),'r')
    xlim([0,2*n-2])
    xticks([0,n-1,2*n-2],['$K/8$','$\Gamma$','$A/4$'])
    ylim([-.6,0])
    yticks([0,-.1,-.2,-.3,-.4,-.5,-.6],['0','','-.2','','-.4','','-.6'])
    xlabel('Wavevector $\\vec{k}$')
    ylabel('Energy below $E_V$ [eV]')
    mpl.axvline(n-1,color='k')
if 1:
    n=20
#    figure()
    karr=path([-APoint/80,GammaPoint,APoint/80])
    b=np.array([bs(k) for k in karr])
    dk=np.max(abs(karr[1]))-np.max(abs(karr[0]))
    m=1/(np.gradient(np.gradient(b[:,0]*q*10000,dk),dk)/(hbar*q*10000)**2)
    print "M: "
    print m[20]/m0
    n=10000


figure()
karr=path([1j*APoint/2,GammaPoint,APoint/2])+.01*KPoint
b=np.array([bs(k) for k in karr])
plot(b[:,0],'k',linewidth=2)
xlim([0,n-1])
xticks([0,n-1,2*n-2],['$iA/2$','$\Gamma$','$A/2$'])
xlabel('Wavevector $k_z$')
ylabel('Energy [eV]')




#b2=np.array([bs2bm(k) for k in karr])
#plot(b2[:10000,2:4],'g--',linewidth=2)

#bbad=np.array([bs(k,bad=True) for k in karr])
#plot(bbad[:10000,2:4],'ro',linewidth=2)


plot(b[:,1:4],'k',linewidth=2)
mpl.axvline(n-1,color='k')
stop

#karr=path([1j*APoint,GammaPoint,APoint])+.1*KPoint
#b=np.array([bs(k,False) for k in karr])
#plot(b[:,:],'--')

karr=path([1j*APoint,GammaPoint,APoint])+.2*KPoint
b=np.array([bs(k,True) for k in karr])
plot(b[:,:],'.')
theb=b


def gen_im_EK(k,bad=False):
    b=np.array([bs(ki,bad) for ki in k])
    k=k[np.where(np.logical_not(np.any(np.isnan(b),1)))[0],:]
    b=np.array([bs(ki,bad) for ki in k])
    
    ki=np.concatenate([k[:,2].imag,np.flipud(k[:,2].imag)])
    Ei=np.concatenate([b[:,2],np.flipud(b[:,3])])
    
    Ebounds=[np.min(Ei),np.max(Ei)]
    return sciint.interp1d(Ei,ki,kind='linear'),Ebounds
figure()
karr=path([GammaPoint,1j*APoint])+0*KPoint
EK,Ebounds=gen_im_EK(karr)
Earr=np.linspace(Ebounds[0],Ebounds[1],1000)
plot(Earr,EK(Earr)/(np.pi/c),'b')
EK,Ebounds=gen_im_EK(karr,True)
Earr=np.linspace(Ebounds[0],Ebounds[1],1000)
plot(Earr,EK(Earr)/(np.pi/c),'b--',linewidth=2)

#karr=path([GammaPoint,1j*APoint])+.2*KPoint
#EK,Ebounds=gen_im_EK(karr)
#Earr=np.linspace(Ebounds[0],Ebounds[1],100)
#plot(Earr,EK(Earr)/(np.pi/c))

karr=path([GammaPoint,1j*APoint])+.5*KPoint
EK,Ebounds=gen_im_EK(karr)
Earr=np.linspace(Ebounds[0],Ebounds[1],100)
plot(Earr,EK(Earr)/(np.pi/c),'g')
EK,Ebounds=gen_im_EK(karr,True)
Earr=np.linspace(Ebounds[0],Ebounds[1],100)
plot(Earr,EK(Earr)/(np.pi/c),'g--',linewidth=2)

karr=path([GammaPoint,1j*APoint])+.025*KPoint
EK,Ebounds=gen_im_EK(karr)
Earr=np.linspace(Ebounds[0],Ebounds[1],100)
plot(Earr,EK(Earr)/(np.pi/c),'r')
EK,Ebounds=gen_im_EK(karr,True)
Earr=np.linspace(Ebounds[0],Ebounds[1],100)
plot(Earr,EK(Earr)/(np.pi/c),'r--',linewidth=2)

b=np.array([bs(k,True) for k in karr])
#plot(b,'.')

karr=path([1j*APoint,GammaPoint,APoint])+0*KPoint
b=np.array([bs(k) for k in karr])
#plot(b,'--')


if 0:
    #figure()
    Egp=Eg+DCR
    me=1/(1/m0+2/hbar**2*Ap1/(10000*q))
    mh=1/(1/m0+2/hbar**2*Lp2/(10000*q))
    mr=1/(1/me+1/mh)
    
    kvl=lambda l: np.sqrt(((Egp*hbar**2*q*10000/(2*mh)-P1**2-l*hbar**2*q*10000/(2*mr))\
        +np.sqrt((Egp*hbar**2*q*10000/(2*mh)-P1**2-l*hbar**2*q*10000/(2*mr))**2+hbar**4*q**2*10000**2*l/(me*mh)*(Egp-l)))/
        (hbar**4*q**2*10000**2/(2*me*mh)))

    def a(kp):
        return 1+4.6397e-16*kp**2+1.5231e-9*kp
    def b(kp):
        return 3.2478e-16*kp**2+1.5231e-10*kp
    
    print "DUMMY"
    print ([a(.025*KPoint[0]),b(.025*KPoint[0])])    
    print ([a(.05*KPoint[0]),b(.05*KPoint[0])])    
    
    plot(Earr,kvl(Earr)/(np.pi/c),'x')
    #plot(Earr*1.08,kvl(Earr)/(np.pi/c)+.05,'o')
    #plot(Earr,kvl(Earr/a(.2*KPoint[0]))/(np.pi/c)+b(.2*KPoint[0]),'o')
    plot(Earr,kvl(Earr/a(.025*KPoint[0]))/(np.pi/c)+b(.025*KPoint[0]),'o')
    plot(Earr,kvl(Earr/a(.05*KPoint[0]))/(np.pi/c)+b(.05*KPoint[0]),'o')
    
    plot(Earr,sqrt(2*(2*m0/(1/.2+1/1.78))*(Earr*Eg-Earr**2)/Eg)/(hbar*sqrt(q*10000))/(np.pi/c),'s')
    
    #stop
    
    l=lambda kz: .5*((Egp+hbar**2*q*10000*kz**2/(2*mr))*np.array([[1],[1]])\
        +np.array([[-1],[1]])*np.sqrt((Egp+hbar**2*q*10000*kz**2/(2*mr))**2-4*(-P1**2*kz**2+(Egp+hbar**2*q*10000*kz**2/(2*me))*hbar**2*q*10000*kz**2/(2*mh)))).T
    
    l3d=lambda kx,ky,kz: .5*((Egp+hbar**2*q*10000*kz**2/(2*mr)+(Ap2+M3)*(kx**2+ky**2))*np.array([[1],[1]])\
        +np.array([[-1],[1]])*np.sqrt((Egp+hbar**2*q*10000*kz**2/(2*mr)+(Ap2+M3)*(kx**2+ky**2))**2-4*(-P1**2*kz**2+(Egp+hbar**2*q*10000*kz**2/(2*me)+Ap2*(kx**2+ky**2))*(hbar**2*q*10000*kz**2/(2*mh)+M3*(kx**2+ky**2))))).T
        
    lbad=lambda kz: .5*((Egp+hbar**2*q*10000*kz**2/(2*mr))*np.array([[1],[1]])\
        +np.array([[-1],[1]])*np.sqrt((Egp+hbar**2*q*10000*kz**2/(2*mr))**2-4*(-P1**2*abs(kz)**2+(Egp+hbar**2*q*10000*kz**2/(2*me))*hbar**2*q*10000*kz**2/(2*mh)))).T
    
    #l=lambda kz: .5*((Egp+hbar**2*q*10000*kz**2/(2*mr))*np.array([[1],[1]])\
    #    +np.array([[-1],[1]])*np.sqrt((Egp+hbar**2*q*10000*kz**2/(2*mr))**2)).T
    
    #lvals=l(karr[:,2])
    lvals=l3d(karr[:,0],karr[:,1],karr[:,2])
    lvals=np.choose(abs(lvals.imag)<1e-10,[np.NaN,lvals.real])
    plot(lvals,'.')
    #plot(lbad(karr[:,2]))
    mpl.xlim([0,n-1])
    mpl.xticks([0,n-1,2*n-2],['$iA$','$\Gamma$','$A$'])
    mpl.xlabel('Wavevector $k_z$')
    mpl.ylabel('Energy [eV]')
    mpl.axvline(n-1,color='k')
    
if 1:
    figure()
    k=[np.array([2*x*np.pi/a,2*y*np.pi/a,0]) for x,y in zip(*[np.ravel(x) for x in np.meshgrid(np.linspace(-1,1),np.linspace(-1,1))])]
    kx=np.reshape([ki[0] for ki in k],[50,50])
    ky=np.reshape([ki[1] for ki in k],[50,50])
    b=np.reshape([bs(ki)[2] for ki in k],[50,50])
    mpl.contour(kx,ky,b)

if 0:
    figure()
    title('kx,kz')
    k=[np.array([2*x*np.pi/a,0,2*z*np.pi/a]) for x,z in zip(*[np.ravel(x) for x in np.meshgrid(np.linspace(0,1),np.linspace(0,1))])]
    kx=np.reshape([ki[0] for ki in k],[50,50])
    kz=np.reshape([ki[1] for ki in k],[50,50])
    b=np.reshape([bs(ki)[2] for ki in k],[50,50])
    #mpl.den(kx,ky,b)

if 1:
    figure()
    isim=1
    
    title('vs kx')
    karr=path([-.3*KPoint,GammaPoint,.3*KPoint])+0*isim*APoint
    b=np.array([bs(k) for k in karr])
    plot(b[:,1]-np.max(b[:,1]),'b')
    
    karr=path([-.3*KPoint,GammaPoint,.3*KPoint])+.05*isim*APoint
    b=np.array([bs(k) for k in karr])
    plot(b[:,1]-np.max(b[:,1]),'g')
    
    karr=path([-.3*KPoint,GammaPoint,.3*KPoint])+.25*isim*APoint
    b=np.array([bs(k) for k in karr])
    plot(b[:,1]-np.max(b[:,1]),'r')
    
    karr=path([-.3*KPoint,GammaPoint,.3*KPoint])+.5*isim*APoint
    b=np.array([bs(k) for k in karr])
    plot(b[:,1]-np.max(b[:,1]),'k')
    
    karr=path([-.3*KPoint,GammaPoint,.3*KPoint])
    plot(-hbar**2*(q*10000)*karr[:,0]**2/(2*.15*m0),'--k')
    
    xlim([0,n-1])
    ylim([.5,-1])
    xticks([0,n-1,2*n-2],['$-K/3$','$\Gamma$','$K/3$'])
    xlabel('Wavevector $k_x$')
    ylabel('Energy [eV]')
    mpl.axvline(n-1,color='k')