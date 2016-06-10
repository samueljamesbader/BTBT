# -*- coding: utf-8 -*-
from __future__ import division
import numpy as np
import matplotlib.pyplot as mpl
from Constants import eps0,q,hbar,m0,kT,eVToJoules
from SpikedPINJunction import SpikedPINJunction
from GaN_BS_min import GaN_imagBS
from scipy.interpolate import interp1d
kvl=GaN_imagBS2()

def tunnel_current(junc,justD=False):
        
    vwindowmin=max(junc.Ec[-1],junc.mu[-1])
    vwindowmax=min(junc.Ev[0],junc.mu[0])
    
    print vwindowmax
    print vwindowmin
    
    Fbound=1e4
    izero=np.argmax(junc.x==0)
    iwindowstart=np.argmax(np.logical_and(junc.F[:izero]<-Fbound,junc.Ev[:izero]<vwindowmax))
    iwindowend=np.argmax(np.logical_or(junc.F[iwindowstart:]>-Fbound,junc.Ec[iwindowstart:]<vwindowmin))+iwindowstart

    x=junc.x[iwindowstart:iwindowend]
    Ec=junc.Ec[iwindowstart:iwindowend]
    Ev=junc.Ev[iwindowstart:iwindowend]
    
    
        
#    print Ez
#    print Ev[Ev>vwindowmin][-1]
    
    
#    close('all')
#    figure()
#    plot(junc.x*1e7,junc.Ec)
#    plot(junc.x*1e7,junc.Ev)
#    plot(junc.x*1e7,junc.mu,'r')
    
    
    mpl.figure()
    mpl.plot(x*1e7,Ec,'.')
    mpl.plot(x*1e7,Ev,'.')

    if justD:
        Evx=interp1d(np.flipud(Ev),np.flipud(x))
        Ecx=interp1d(np.flipud(Ec),np.flipud(x))
        E=np.linspace(np.min(Ec),np.max(Ev))
        ds=np.min(Ecx(E)-Evx(E))
        return ds
        
        
        
    E=vwindowmax-Ev[Ev>vwindowmin]
#    print "Ex goes from to ", np.min(Ex),np.max(Ex)
#    figure()
    
    Ep=E[:100];
    if len(E): Ep=Ep-E[0];
    #Ep=[0,.04,.08]
    expts=[
        [ -2*np.trapz(kvl(-Ei-(Ev-junc.Ev[0]),np.sqrt(2*m0*Epi/(hbar**2*q*10000))),x)
                #if Exi+2*Epi<(vwindowmax-vwindowmin) else -np.inf
            for Ei in E]
        for Epi in Ep]
    print 'hi'
    print expts
    print np.shape(expts)
#    if not(np.prod(np.shape(expts))):
#        hello
    iint=np.trapz(np.trapz(np.exp(expts),Ep,axis=0),E,axis=0)\
        if np.prod(np.shape(expts)) else 0
    figure();plot(E,np.transpose(np.array(expts)));title('V: -'+str(junc.Vr))
    
    pref=q*(m0)/(2*np.pi**2*hbar**3*eVToJoules) #A/m^2eV^2
    print "IMAG BS IS WONKY NEAR 0"
    print pref*iint / 1e4, " A/cm^2" # A/cm^2
    
    
    if iint>0:
        print "max epxt: ", np.max(expts)

    return pref*iint    
    
    figure()
    for Exi in Ex:
        kxvx=kvl(-Exi-(Ev-junc.Ev[0]),0)
        expt=-2*trapz(kxvx,x)
#        print expt
        plot(x*1e7,kxvx/1e7)
    
    
def barrier_tunnel_current(junc):
        
    vwindowmin=np.min(junc.Ev)
    vwindowmax=junc.Ev[-1]
    
    print vwindowmax
    print vwindowmin
    
    Fbound=1e4
    izero=np.argmax(junc.x==0)
    iwindowstart=np.argmax(np.logical_and(junc.F[:izero]<-Fbound,junc.Ev[:izero]<vwindowmax))
    iwindowend=np.argmax(abs(junc.F[iwindowstart:])<Fbound)+iwindowstart

    x=junc.x[iwindowstart:iwindowend]
    Ec=junc.Ec[iwindowstart:iwindowend]
    Ev=junc.Ev[iwindowstart:iwindowend]
    

#    print Ez
#    print Ev[Ev>vwindowmin][-1]
    
    if 0:
        close('all')
        figure()
        plot(junc.x*1e7,junc.Ec)
        plot(junc.x*1e7,junc.Ev)
        plot(junc.x*1e7,junc.mu,'r')
        
        
        figure()
        plot(x*1e7,Ec,'.')
        plot(x*1e7,Ev,'.')

    
    
    E=junc.Ev[0]-Ev[0:np.argmin(Ev)]
#    print "Ex goes from to ", np.min(Ex),np.max(Ex)
#    figure()
    
    Ep=E[:100];
    if len(E): Ep=Ep-E[0];
    #Ep=[0,.04,.08]
    
    barkvl=lambda l,Epi: np.nan_to_num(np.sqrt(2*.15*m0*(l+Epi)/(hbar**2*q*10000)))
    print E
    print -E[0]-(Ev-junc.Ev[0])
#    expts=[
#        [ -2*np.trapz(barkvl(-Ei-(Ev-junc.Ev[0]),Epi),x)
#                #if Exi+2*Epi<(vwindowmax-vwindowmin) else -np.inf
#            for Ei in E]
#        for Epi in Ep]
#    iint=np.trapz(np.trapz(np.exp(expts),Ep,axis=0),E,axis=0)\
#        if np.prod(np.shape(expts)) else 0
#    figure();plot(E,np.transpose(np.array(expts)));title('V: -'+str(junc.Vr))
    
    
    wexpts=[
        [ -2*np.trapz(barkvl(-Ei-(Ev-junc.Ev[0]),Epi),x)+(-(Ei+Epi)/kT)
                #if Exi+2*Epi<(vwindowmax-vwindowmin) else -np.inf
            for Ei in E]
        for Epi in Ep]
    iint=np.trapz(np.trapz(np.exp(wexpts),Ep,axis=0),E,axis=0)\
        if np.prod(np.shape(wexpts)) else 0
    
        
    
    #pref=q*(m0)/(2*np.pi**2*hbar**3*eVToJoules) #A/m^2eV^2
    pref=np.exp((junc.Ev[0]-junc.Ev[-1])/kT)/(kT)**2
    print pref*iint , " (weighted T)"
    
#    
#    if iint>0:
#        print "max epxt: ", np.max(expts)

    return pref*iint    


if __name__=="__main__":    
    #    figure()
    #    for Ei in E:
    #        kxvx=kvl(-Ei-(Ev-junc.Ev[0]),0)
    #        expt=-2*trapz(kxvx,x)
    ##        print expt
    #        plot(x*1e7,kxvx/1e7)
        
    
    #SpikedPINJunction("GaN",
    #           "Mg",8e19,.1e-4,
    #            0*.005e-4,0e13,
    #            "Si",2e19,.2e-4,
    #            -.2).plot_band_diagram()
    #stop
    ###### FOR EXAMINING THE INTERBAND TUNNELING CURRENT
    V=np.concatenate([np.linspace(-5,-.5,51),np.linspace(-.5,0,20)])
    #V=[-5]
    curs=np.array(
        [tunnel_current(
            SpikedPINJunction("GaN",
               "Mg",8e19,.1e-4,
                0*.005e-4,0*8e13,
                "Si",2e19,.2e-4,
                Vi))
        for Vi in V])
    figure()
    plot(V,curs/1e4)
    yscale('log')
        
        
    
    #### FOR EXAMINING THE BARRIER TUNNELING CURRENT    
    #V=np.linspace(0,2,11)
    ##V=[0]
    #juncs=[SpikedPINJunction("GaN",
    #   "Mg",4e19,.1e-4,
    #    0*.005e-4,5e13,
    #    "Si",1e19,.2e-4,
    #    Vi) for Vi in V]
    #curs=np.array([barrier_tunnel_current(j) for j in juncs])
    #Vns=np.array([j.Ev[-1]-np.min(j.Ev) for j in juncs])
    #figure()
    #plot(Vns/kT,curs)
    #plot(Vns/kT,exp(-Vns/(2.8*kT)),'.')
    #title('vs vns')
    #yscale('log')
    #
    #print (np.log(curs[-1]/curs[0])/((Vns[-1]-Vns[0])/kT))**-1
    
    
    
    
    #figure()
    #plot(V,curs)
    #title('vs va')
    #yscale('log')