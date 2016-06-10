# -*- coding: utf-8 -*-
"""
Author: Sam Bader
Module: 
"""
from __future__ import division
from Constants import kT,q
from SpikedPINJunction import SpikedPINJunction
import numpy as np
import matplotlib.pyplot as mpl
from Material import Material
from FermiDirac import int_fermi_dirac_one_half as G12
from scipy.stats import linregress
mpl.close('all')

def Vn_vs_VA_curve(NA,ND,sigma,Va):
    pside=Material("GaN",{"Mg":NA})
    Nv=pside.Nv
    EA=pside.dopants.values()[0]['E']
    EVp=pside.Ev-pside.mu_bulk
    NAp=NA*np.log(1+4*np.exp((EVp+EA)/kT))
    Gp=G12(EVp/kT)
    
    
    nside=Material("GaN",{"Si":ND})
    Nc=nside.Nc
    ED=nside.dopants.values()[0]['E']
    ECn=nside.Ec-nside.mu_bulk
    NDp=ND*np.log(1+.5*np.exp((ECn-ED)/kT))
    Gn=G12(-ECn/kT)
    
    eps=nside.eps
    Vbi=nside.mu_bulk-pside.mu_bulk
    
    Sigma=np.sqrt(q*sigma**2/(2*eps*kT*NA))
    
    Vsigma=Vbi-Sigma**2*kT+(Nc/NA*kT*Gn-NAp/NA*kT-Nv/NA*Gp*kT)
    print "Vsigma: ",Vsigma
    
    def Fn_of_Vn(Vn):
        return np.sqrt(2*q*kT*Nc/eps)*\
            np.sqrt(G12((Vn-ECn)/kT)-Gn-NDp/Nc)
    def Fp_of_Vn(Vn,Va):
        Vp=Vn+Vbi-Va
        return -np.sqrt(2*q*kT*Nc/eps)*\
            np.sqrt(NA/Nc*Vp/kT+G12((Vn-ECn)/kT)-NAp/Nc-Nc/Nc*Gp)
    
    
    def Vn_of_Va(Vai):
        print ""
        print "VA: ",Vai
        g=SpikedPINJunction("GaN","Mg",NA,.1e-4,0*.005e-4,sigma,"Si",ND,.2e-4,Vai)
    #    g.plot_all()
        Vni=g.Ec[-1]-np.min(g.Ec)
        print "Predicted Fn: "+str(Fn_of_Vn(Vni)/1e6)
        print "Predicted Fp: "+str(Fp_of_Vn(Vni,Vai)/1e6)

        if Vni<.3*kT:
            return np.nan
        return Vni
        
    
    Vn=np.vectorize(Vn_of_Va,otypes=['float64'])(Va)
    
    
    ind=np.argmin(np.abs(Va))
    alpha_exp,Vn0=linregress(Va[ind-2:ind+3],Vn[ind-2:ind+3])[0:2]
    alpha= 1/((2/(np.pi)**.25*Sigma*np.sqrt(Nc/NA)*((Vn0-ECn)/kT)**.25)-1)
    
    
    mpl.figure()
    mpl.plot(Va,Vn,'.')
    
    LHS=4*Nc/NA*Sigma**2*(G12((Vn-ECn)/kT)-Gn-NDp/Nc)
    RHS=((Vn+Vsigma-Va)/kT+NDp/NA)**2
    #print "LHS ", LHS
    #print "RHS ", RHS
    print "LHS RHS Error: ", (LHS-RHS)/RHS
    
    S=(32/(15*np.sqrt(np.pi))*(Nc/NA)*Sigma**2)**.4
    DV=Va-Vsigma-NDp/NA*kT
    
    
    Vnpre=kT*((1/S)*(DV/kT)**.8+ECn/kT)/(1+(4/(5*S)*(kT/DV)**.2))
    
    return Vn,Vnpre,Vsigma
    
    #nu=(Va-Vsigma)/kT
    #A=2/(25*(nu)**(6/5))
    #B=((32/(15*np.sqrt(np.pi))*(Nc/NA)*Sigma**2)**.4+(4/(5*nu**.2)))
    #C=(-nu**.8-(32/(15*np.sqrt(np.pi))*(Nc/NA)*Sigma**2)**.4*ECn/kT)
    #Vnpre=kT*(-B+np.sqrt(B**2-4*A*C))/(2*A)
    
    
    
    #LHSapp=S*(Vn-ECn)/kT
    #RHSapp=((nu*kT)**.8-.8*Vn/(kT*nu)**.2-.08*Vn**2/(kT*nu)**1.2)/kT**.8
    
    
    #print "FirstTerm"
    #print np.sqrt(q/(2*eps*kT*Nc))*sigma
    #print "SecondTerm"
    #print np.sqrt(G12((Vn-ECn)/kT)-Gn-NDp/Nc)
    #print "ThirdTerm"
    #print np.sqrt(NA/Nc*(Vn+Vbi-Va)/kT+G12((Vn-ECn)/kT)-NAp/Nc-Nv/Nc*Gp)
    #
    #
    #print "LHS"
    #print G12((Vn-ECn)/kT)
    #print "RHS"
    #print A*((Vn-Va)/kT+B)**2+C

if 1:
    Va=np.linspace(-1,1,9)
    
    NA=1e19
    ND=1e19
    sigma=5e13
    Vn15,Vnpre15,Vsigma=Vn_vs_VA_curve(NA,ND,sigma,Va)
    title("NA="+str(NA)+" sig="+str(sigma)+" Vsig="+str(Vsigma))
    
    NA=4e19
    ND=1e19
    sigma=5e13
    Vn45,Vnpre45,Vsigma=Vn_vs_VA_curve(NA,ND,sigma,Va)
    title("NA="+str(NA)+" sig="+str(sigma)+" Vsig="+str(Vsigma))
    
    NA=7e19
    ND=1e19
    sigma=5e13
    Vn75,Vnpre75,Vsigma=Vn_vs_VA_curve(NA,ND,sigma,Va)
    title("NA="+str(NA)+" sig="+str(sigma)+" Vsig="+str(Vsigma))
    
    NA=1e19
    ND=1e19
    sigma=6e13
    Vn16,Vnpre16,Vsigma=Vn_vs_VA_curve(NA,ND,sigma,Va)
    title("NA="+str(NA)+" sig="+str(sigma)+" Vsig="+str(Vsigma))
    
    NA=4e19
    ND=1e19
    sigma=6e13
    Vn46,Vnpre46,Vsigma=Vn_vs_VA_curve(NA,ND,sigma,Va)
    title("NA="+str(NA)+" sig="+str(sigma)+" Vsig="+str(Vsigma))
    
    NA=7e19
    ND=1e19
    sigma=6e13
    Vn76,Vnpre76,Vsigma=Vn_vs_VA_curve(NA,ND,sigma,Va)
    title("NA="+str(NA)+" sig="+str(sigma)+" Vsig="+str(Vsigma))

figure()
plot(Va,Vn15,'bo-',label='$N_A=1e19$, $\sigma=5e13$')
plot(Va,Vn45,'go-',label='$N_A=4e19$, $\sigma=5e13$')
plot(Va,Vn75,'ro-',label='$N_A=7e19$, $\sigma=5e13$')

plot(Va,Vn16,'b',label='$N_A=1e19$, $\sigma=6e13$')
plot(Va,Vn46,'g',label='$N_A=4e19$, $\sigma=6e13$')
plot(Va,Vn76,'r',label='$N_A=7e19$, $\sigma=6e13$')

plot(Va,Vnpre15,'b.',label='$N_A=1e19$, $\sigma=5e13$')
plot(Va,Vnpre45,'g.',label='$N_A=4e19$, $\sigma=5e13$')
plot(Va,Vnpre75,'r.',label='$N_A=7e19$, $\sigma=5e13$')

plot(Va,Vnpre16,'b.',label='$N_A=1e19$, $\sigma=6e13$')
plot(Va,Vnpre46,'g.',label='$N_A=4e19$, $\sigma=6e13$')
plot(Va,Vnpre76,'r.',label='$N_A=7e19$, $\sigma=6e13$')

#plot(Va,[Vn15,Vn45,Vn75,Vn16,Vn46,Vn76])
ylim(0,1)
ylabel('$V_n$')
xlabel('$V_A$')
#legend(loc='upper center',prop={'size':16},ncol=2)