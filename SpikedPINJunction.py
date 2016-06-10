# -*- coding: utf-8 -*-
"""
Author: Sam Bader
Module: 
"""
from __future__ import division
from __future__ import print_function
import numpy as np
import scipy.optimize as sciopt
from Constants import eps0,q,hbar,m0,kT,eVToJoules
from Material import Material,materialparameters
import matplotlib.pyplot as mpl
import warnings
#from OutputRedirection import stdout_redirected, merged_stderr_stdout

class SpikedPINJunction:
    def __init__(self, mat_name,
                acceptor_name, Na, p_len,
                i_len, sigma,
                donor_name, Nd, n_len,
                Va, sigma_name="Si"):
        
        # The Material class will find the bulk energy levels
        mp=Material(mat_name,{acceptor_name: Na})
        mn=Material(mat_name,{donor_name: Nd})
        
        # Band edges far on the p and n side
        Ecp=mp.Ec-mp.mu_bulk-Va/2
        Evp=mp.Ev-mp.mu_bulk-Va/2
        Ecn=mn.Ec-mn.mu_bulk+Va/2
        Evn=mn.Ev-mn.mu_bulk+Va/2    
        
        Vr=-Va
        
        Edsigma=materialparameters[mat_name]["dopants"][sigma_name]["E"]
        
        
        #minfield=000 #V/cm
        def Fp(dEFIp):
            intrhodEFI=-mp.intrho(EFn=mn.mu_bulk-(Ecp-Ecn)-dEFIp,EFp=mp.mu_bulk-dEFIp)
#            print (intrhodEFI,mn.mu_bulk-(Ecp-Ecn)-dEFIp,mp.mu_bulk-dEFIp)
            intrhodEFI=intrhodEFI*(intrhodEFI>0)
            
            #if np.all(np.abs(np.sqrt((2/mn.eps)*intrhodEFI))<minfield): return 0*dEFIp
            #return np.sign(mp.rho(mp.mu_bulk-dEFIp))*np.sqrt((2/mp.eps)*intrhodEFI)
            return np.sign(dEFIp)*np.sqrt((2/mp.eps)*intrhodEFI)
        def Fn(dEFIn,x=0):
            intrhodEFI=-mn.intrho(EFn=mn.mu_bulk-dEFIn,EFp=mp.mu_bulk+(Ecp-Ecn)-dEFIn)
            intrhodEFI=intrhodEFI*(intrhodEFI>0)
            
            #if np.all(np.abs(np.sqrt((2/mn.eps)*intrhodEFI))<minfield): return 0*dEFIn
            #return -np.sign(mn.rho(mn.mu_bulk-dEFIn))*np.sqrt((2/mn.eps)*intrhodEFI)
            #print(x,"    ",dEFIn,"    ",-np.sign(dEFIn)*np.sqrt((2/mn.eps)*intrhodEFI))
            return -np.sign(dEFIn)*np.sqrt((2/mn.eps)*intrhodEFI)
        
        def FDiff(dEFIp):
            fp=Fp(dEFIp)
            dEFIn=dEFIp+fp*i_len+Ecp-Ecn
            fn=Fn(dEFIn)
#            sigma_eff=sigma*(1/(1+2*np.exp((mn.mu_bulk-dEFIn-mn.Ec+Edsigma)/kT)))
#            print ("srat",1/(1+2*np.exp((mn.mu_bulk-dEFIn-mn.Ec+Edsigma)/kT)))
            dF=q*sigma/mn.eps
#            dF=q*sigma_eff/mn.eps
#            print("dF: ",dF/1e6)
            return np.abs(fn-fp-dF)
            
        dEFIp_bracket=[-(Ecp-Ecn),-Nd/(Nd+Na)*(Ecp-Ecn)]
        res=sciopt.minimize_scalar(FDiff,bracket=dEFIp_bracket)
        dEFIp=res.x
        fp=Fp(dEFIp)
        dEFIn=dEFIp+fp*i_len+Ecp-Ecn
        fn=Fn(dEFIn)
        
        print("dEFIp: ",dEFIp)
        print("dEFIn: ",dEFIn)
        print ("Fn: ",fn/1e6," MV/cm")
        print ("Fp: ",fp/1e6," MV/cm")
        
        
        from scipy.integrate import odeint
        
        # x-grid for calculation
        nr=1000
        ni=500
        nl=1000
        xr=np.concatenate([[0],np.logspace(-10,np.log10(n_len),nr)])
        xi=np.linspace(0,-i_len,ni) if i_len else np.array([0])
        xl=np.concatenate([xi,-np.logspace(-10,np.log10(p_len),nl)-i_len])
        x=np.concatenate([xl[::-1],xr[1:]])        
        
        Ecr=Ecn+np.ravel(odeint(\
            lambda E,x: Fn(E,x),\
            dEFIn,xr,rtol=1e-10,atol=1e-8))
        print('int\'d r')
        
        
#        something,fo=odeint(\
#            lambda E,x: Fn(E,x),\
#            dEFIn,xr,rtol=1e-10,atol=1e-8,hmax=1e-7,full_output=1)#,printmessg=1)
#        Ecr=Ecn+np.ravel(something)
#        print('int r')
#        print(fo)
        
        
        FpdEFIp=Fp(dEFIp)
        Ecl=Ecp+np.ravel(odeint(\
            lambda E,x: np.choose(x<-i_len,[FpdEFIp+0*x,Fp(E)]),\
            dEFIp+i_len*fp,xl,rtol=1e-10,atol=1e-8,hmax=1e-7))
        print('int\'d l')
        
        rhor=mn.rho(EFn=-(Ecr-Ecn-mn.mu_bulk),EFp=-(Ecr-Ecn-mp.mu_bulk-(Ecp-Ecn))) ;  print("cheating rhor EFp near junction")
        rhol=np.choose(xl<-i_len,[[0]*len(xl),mp.rho(EFn=-(Ecl-Ecp-mn.mu_bulk+(Ecp-Ecn)),EFp=-(Ecl-Ecp-mp.mu_bulk))]); print("cheating rhol EFn far from junction")
        #rhol=np.choose(xl<-i_len,[[0]*len(xl),mp.rho(EFn=-(Ecl-Ecp-mp.mu_bulk),EFp=-(Ecl-Ecp-mp.mu_bulk))]); print("cheating BAD rhol EFn")
        Fr=Fn(Ecr-Ecn)
        Fl=np.choose(xl<-i_len,[(Fp(dEFIp))+0*xl,Fp(Ecl-Ecp)])
        Ec=np.concatenate([Ecl[:0:-1],Ecr])
        mu=(x<0)*(Vr/2)+(x>0)*(-Vr/2)+np.choose(x==0,[0,np.nan])
        rho=np.concatenate([rhol[:0:-1],[np.nan],rhor[1:]])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore",RuntimeWarning)
            rho=rho*(abs(rho)/q>1e11)
        F=np.concatenate([Fl[:0:-1],[(Fl[0]+Fr[0])/2],Fr[1:]])
        Ev=Ec-mn.Eg
        
        self.x=x+i_len/2
        self.Ec=np.ravel(Ec)
        self.Ev=np.ravel(Ev)
        self.mu=mu
        self.F=F
        self.rho=rho
        self.Vr=Vr
        self.mp=mp
        self.mn=mn
    
    def plot_band_diagram(self,fignum=None,legend=True, dot=False, **kwargs):
        fig=mpl.figure(fignum)
        mpl.plot(self.x*1e7,self.Ev,'g',markersize=3,label='$E_v$', **kwargs)
        mpl.plot(self.x*1e7,self.Ec,'b',markersize=3,label='$E_c$', **kwargs)
        mpl.plot(self.x*1e7,self.mu,'r',markersize=1,label='$\\mu$', **kwargs)
        if dot:
            mpl.plot(0,self.Ec[np.argmin(np.abs(self.x))]-self.mn.dopants['Si']["E"],'o')
        mpl.xlabel("$x$ [nm]")
        mpl.ylabel("$E$ [eV]")
        if legend: mpl.legend(loc='best')
        return fig
    
    def plot_electric_field(self,fignum=None, **kwargs):
        fig=mpl.figure(fignum)
        mpl.plot(self.x*1e7,self.F/1e6,'g',markersize=3, **kwargs)
        mpl.xlabel("$x$ [nm]")
        mpl.ylabel("$F$ [MV/cm]")
        return fig
    
    def plot_charge(self,fignum=None,**kwargs):
        fig=mpl.figure(fignum)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore",RuntimeWarning)
            mpl.plot(self.x*1e7,np.choose(self.rho>0,[np.nan,self.rho/q]),'b',markersize=5,label='+',**kwargs)
            mpl.plot(self.x*1e7,np.choose(self.rho<0,[np.nan,-self.rho/q]),'g',markersize=5,label='-',**kwargs)
        mpl.xlabel("$x$ [nm]")
        mpl.ylabel("$|\\rho|$ [$q$/cm$^3$]")
        mpl.legend(loc='best')
        mpl.yscale('log')
        mpl.ylim(1e17,1e21)
        return fig
    
    def plot_all(self):
        #fig=mpl.figure(figsize=(6,9))
        fig,ax=mpl.subplots(3,sharex=True,figsize=(6,9))
        mpl.axes(ax[0])
        self.plot_band_diagram(fignum=fig.number)
#        mpl.xticks([])
        mpl.xlabel("")
        mpl.axes(ax[1])
        self.plot_electric_field(fignum=fig.number)
#        mpl.xticks([])
        mpl.xlabel("")
        mpl.axes(ax[2])
        self.plot_charge(fignum=fig.number)
        mpl.xlim([-15,15])
        mpl.tight_layout()

if __name__ == "__main__":
    mpl.rc('font',**{'size':20})
    mpl.rc('legend',**{'fontsize':24})
    mpl.close('all')

#    g=SpikedPINJunction("GaN","Mg",8e19,.1e-4,0,8e13,"Si",2e19,.2e-4,-4)
#    g.plot_all()
#    g=SpikedPINJunction("GaN","Mg",8e19,.1e-4,0,8e13,"Si",2e19,.2e-4,-4.5)
#    g.plot_all()
#    #g=SpikedPINJunction("GaN","Mg",8e19,.1e-4,0,8e13,"Si",2e19,.2e-4,-6)
#    #g.plot_all()
#    stop
    

    
    # GOOD FOR VISUALIZATION PURPOSES
#    g=SpikedPINJunction("GaN","Mg",4e19,.1e-4,0*.005e-4,5e13,"Si",1e19,20e-7,-5)
#    g.plot_all()
#    sigma=5.202e13;
    sigma=4.02e13
    
    NA=8e19
    ND=1e19
    Va=0
    
    
    g=SpikedPINJunction("GaAs","C",NA,100e-7,0*.005e-4,sigma,"S",ND,100e-7,Va,sigma_name="S")    
#    g=SpikedPINJunction("GaN","Mg",NA,100e-7,0*.005e-4,sigma,"Si",ND,100e-7,Va)
    
    
#    g.plot_charge()
    f=g.plot_band_diagram(legend=False,linewidth=3)
    
#    sigma=0
#    g=SpikedPINJunction("GaN","Mg",NA,100e-7,0*.005e-4,sigma,"Si",ND,100e-7,Va)    
#    g.plot_band_diagram(fignum=f.number,linestyle='--',legend=False)
    mpl.xlim(-20,20)
#    mpl.axes([.25, .2, .25, .25])
#    g.plot_band_diagram(f.number,legend=False, dot=True)
#    xlabel('')
#    ylabel('')
#    xlim(-5,5)
#    ylim(-3*kT+g.Ec[np.argmin(np.abs(g.x))],3*kT+g.Ec[np.argmin(np.abs(g.x))])
    
#    stop
#    figure()
#    mpl.plot(g.x*1e7,g.Ev,'g',markersize=3,label='$E_v, \, \\mathrm{numeric}$')
#    mpl.xlabel("$x$ [nm]")
#    mpl.ylabel("$E$ [eV]")
#    stop
    
####
    
#    
#    from TunnelingCurrent import tunnel_current as tc
#    print ("D: ",tc(g,justD=False)*1e7)
#    stop
#    g.Ec-
    
####    
    
    from FermiDirac import int_fermi_dirac_one_half as G12
    Nc=g.mn.Nc
    eps=g.mp.eps
    ED=g.mn.dopants.values()[0]['E']
    
    ECn=g.mn.Ec-g.mn.mu_bulk
    Nv=g.mn.Nv
    NDp=ND*np.log(1+.5*np.exp((ECn-ED)/kT))
    Gn=G12(-ECn/kT)
    Vbi=g.mn.mu_bulk-g.mp.mu_bulk
    Sigma=np.sqrt(q*sigma**2/(2*eps*kT*NA))
    
    EA=g.mp.dopants.values()[0]['E']
    EVp=g.mp.Ev-g.mp.mu_bulk
    NAp=NA*np.log(1+4*np.exp((EVp+EA)/kT))
    Gp=G12(EVp/kT)
    
    #Vn=g.Ev[-1]-np.min(g.Ev)
    S=(32/(15*np.sqrt(np.pi))*(Nc/NA)*Sigma**2)**.4
    Vsigma=Vbi-Sigma**2*kT+(Nc/NA*kT*Gn-NAp/NA*kT-Nv/NA*Gp*kT)
    DV=Va-Vsigma-NDp/NA*kT
    Vnpre=kT*((1/S)*(DV/kT)**.8+ECn/kT)/(1+(4/(5*S)*(kT/DV)**.2))
    Vn=Vnpre
#    Vn=g.Ev[-1]-np.min(g.Ev)
    
    stop
    print ("FrG", sqrt(2*q*kT*Nc/eps*(G12((Vn-ECn)/kT))) / 1e6)
    Fr=np.sqrt(16*q*kT*Nc/(15*eps*np.pi**.5)*((Vn-ECn)/kT)**2.5-Gn-NDp/Nc)
    Fl=(q*sigma/eps-Fr)
    print("Fr",Fr/1e6)
    print("Fl",Fl/1e6)
    app=(np.choose(g.x>0,[
            (Fl*(-g.x)+(g.Ev[-1]-Vn)),
            (g.Ev[-1]-ECn-kT*(np.sqrt(q*Nc/(15*eps*kT*np.pi**.5))*(g.x)+((Vn-ECn)/kT)**-.25)**-4)]))
    plot(g.x*1e7,app,'b-',label='$E_v,\, \\mathrm{analytic}$',linewidth=3)
    xlim(-2,20)
    ylim(-3.2,-2.4)
#    ylim(-4,0)
    legend(['$E_V,\, \\mathrm{ numeric}$','$E_V,\, \\mathrm{ analytic}$'],loc='lower right',prop={'size':24})
    
    mpl.axes([.55, .4, .35, .4])
    mpl.plot(g.x*1e7,app,linewidth=3)
    mpl.plot(g.x*1e7,g.Ev,'g',linewidth=3)
    xlim(-.07,.07)
    xticks([-.05,.05])
    ylim(-3.167,-3.156)
    yticks([-3.167,-3.162,-3.157])
    
    
#    g.plot_all()
#    g=SpikedPINJunction("GaN","Mg",4e19,.1e-4,0*.005e-4,5e13,"Si",1e19,20e-7,2.5)
#    g.plot_all()
#    bdplot=g.plot_band_diagram(linewidth=2)
#    mpl.title('Spiked PN Junction')
#    mpl.xlim(-50,50)
#        
    #g=SpikedPINJunction("GaN","Mg",8e19,.1e-4,0*.005e-4,0e14,"Si",4e20,.2e-4,-1)
#    g.plot_all()
#    mpl.title('PN Junction')
    #mpl.xlim(-50,50)
    

    
#    g.plot_electric_field(marker='.')
#    g.plot_charge()    
    
    stop
    
    
#    g=SpikedPINJunction("GaN","Mg",8e19,.1e-4,0*.005e-4,7e13,"Si",2e19,.2e-4,-.6)
    
    g.plot_all()
    