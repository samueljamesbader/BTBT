from __future__ import division
from __future__ import print_function
import numpy as np
import scipy.optimize as sciopt
from Constants import eps0,q,hbar,m0,kT,eVToJoules
from Material import Material
from OutputRedirection import stdout_redirected, merged_stderr_stdout

# Generates a PIN band diagram of straight lines
#
# mat: the material name (eg "Si") of this homojunction
# Lp, Aname, Na: the length, dopant name, and concentration (cm^-3) of the p-region
# Li: the length of the intrinsic region
# Ln, Dname, Nd: the length, dopant name, and concentration (cm^-3) of the n-region
# Va: the applied voltage
#
# Returns a dictionary with x (the position grid), Ec and Ev (the band edges),
# EB (the Fermi level in the doped regions), and Vr (reverse bias) and mat
# (the Material, from the p side, stored for convenience).
def generate_flatpin_bandstructure(mat, Lp,Aname,Na, Li, Ln,Dname,Nd, Va):
    
    # The Material class will find the bulk energy levels
    mp=Material(mat,{Aname:Na})
    mn=Material(mat,{Dname:Nd})
    
    # Band edges far on the p and n side
    Ecp=mp.Ec-mp.EBulk-Va/2
    Evp=mp.Ev-mp.EBulk-Va/2
    Ecn=mn.Ec-mn.EBulk+Va/2
    Evn=mn.Ev-mn.EBulk+Va/2
    
    # x-grid for the calculation
    x=np.linspace(-Lp-Li/2,Ln+Li/2,50000)
    return {"x":x,\
            
        # np.piecewise linear functions for the conduction band and valence band
        "Ec": np.piecewise(x,
            [x<-Li/2, (x>=-Li/2)&(x<Li/2), (x>=Li/2)],
            [lambda x: Ecp, lambda x: Ecp+(Ecn-Ecp)*(x+Li/2)/Li, lambda x: Ecn]),\
        "Ev": np.piecewise(x,
            [x<-Li/2, (x>=-Li/2)&(x<Li/2), (x>=Li/2)],
            [lambda x: Evp, lambda x: Evp+(Evn-Evp)*(x+Li/2)/Li, lambda x: Evn]),\
        
        # The Fermi levels, only given in the p and n regions
        "EB": np.piecewise(x,
            [x<-Li/2, (x>=-Li/2)&(x<Li/2), (x>=Li/2)],
            [lambda x: -Va/2, lambda x: np.NaN, lambda x: Va/2]),\
        
        # The reverse bias
        "Vr": -Va,
        
        # The material
        "material": mp}

# np.vectorize the above function so it actually expects an array for the applied
# voltage, and returns an array of the above return values
generate_flatpin_bandstructure=np.vectorize(generate_flatpin_bandstructure)


# Generates a PIN band diagram using the depletion approximation
#
# mat: the material name (eg "Si") of this homojunction
# Lp, Aname, Na: the length, dopant name, and concentration (cm^-3) of the p-region
# d: the length of the intrinsic region
# Ln, Dname, Nd: the length, dopant name, and concentration (cm^-3) of the n-region
# Va: the applied voltage
#
# Returns a dictionary with x (the position grid), Ec and Ev (the band edges),
# EB (the Fermi level in the bulk regions), and Vr (reverse bias) and mat
# (the Material, from the p side, stored for convenience).
def generate_pin_bandstructure(mat, Lp,Aname,Na,d, Ln,Dname,Nd, Va):
    
    # The Material class will find the bulk energy levels
    mp=Material(mat,{Aname:Na})
    mn=Material(mat,{Dname:Nd})
    
    # Band edges far on the p and n side
    Ecp=mp.Ec-mp.EBulk-Va/2
    Evp=mp.Ev-mp.EBulk-Va/2
    Ecn=mn.Ec-mn.EBulk+Va/2
    Evn=mn.Ev-mn.EBulk+Va/2
    
    # Depletion edges, assuming depletion fully within device.
    # Note, these are locations in x (with appropriate signs), where depletion ends
    # ...not the absolute widths of the respective regions!
    xn= (-Na*d+np.sqrt(Na**2*d**2+2*mn.eps*(Ecp-Ecn)/q*Na/Nd*(Na+Nd)))/(Na+Nd) + d/2
    xp=-(-Nd*d+np.sqrt(Nd**2*d**2+2*mn.eps*(Ecp-Ecn)/q*Nd/Na*(Na+Nd)))/(Na+Nd) - d/2
    
    # P-side depletion should stay within the device
    assert (xp > -Lp-d/2), "p-side depletion region extends outside device!"
    
    # x-grid for calculation
    x=np.linspace(-Lp-d/2,Ln+d/2,50000)
    
    # If the n-side depletion region is fully within the device
    if (xn<Ln+d/2):
        
        # Return np.piecewise band diagram from depletion approximation
        return {"x":x,\
            "Ec":np.piecewise(x,
                [x<xp, (x>=xp)&(x<-d/2), (x>=-d/2)&(x<d/2), (x>=d/2)&(x<xn), x>=xn],
                [lambda x: Ecp, lambda x: Ecp-.5*q*Na/mp.eps*(x-xp)**2,
                 lambda x: Ecp-.5*q*Na/mp.eps*(-d/2-xp)**2-q*Na/mp.eps*(-d/2-xp)*(x+d/2),
                 lambda x: Ecn+.5*q*Nd/mn.eps*(xn-x)**2, lambda x: Ecn]),\
            "Ev":np.piecewise(x,
                [x<xp, (x>=xp)&(x<-d/2), (x>=-d/2)&(x<d/2), (x>=d/2)&(x<xn), x>=xn],
                [lambda x: Evp, lambda x: Evp-.5*q*Na/mp.eps*(x-xp)**2,
                 lambda x: Evp-.5*q*Na/mp.eps*(-d/2-xp)**2-q*Na/mp.eps*(-d/2-xp)*(x+d/2),
                 lambda x: Evn+.5*q*Nd/mn.eps*(xn-x)**2, lambda x: Evn]),\
            "EB": np.piecewise(x,
                [x<xp, (x>=xp)&(x<xn), (x>=xn)],
                [lambda x: -Va/2, lambda x: np.NaN, lambda x: Va/2]),\
            "Vr": -Va,\
            "material": mp}
    
    # If the n-side depletion region hits the device edge
    else:
        # Assume a sheet charge after the n-side to rebalance charge
        # Recalculate xp, and the sheet charge.
        xp=-(-(d+Ln)+np.sqrt((d+Ln)**2+(Nd/Na*Ln**2+2*mp.eps*(Ecp-Ecn)/(q*Na))))-d/2
        sigma_n=Na*(-d/2-xp)-Nd*Ln
        
        # Again, make sure the p-side depletion is still inside the device
        assert (xp > -Lp), "p-side depletion region extends outside device!"
        
        # Return a np.piecewise bandstructure given our maxed out depletion
        return {"x":x,\
            "Ec":np.piecewise(x,
                [x<xp, (x>=xp)&(x<-d/2), (x>=-d/2)&(x<d/2), (x>=d/2)&(x<Ln+d/2), x>=Ln+d/2],
                [lambda x: Ecp, lambda x: Ecp-.5*q*Na/mp.eps*(x-xp)**2,
                 lambda x: Ecp-.5*q*Na/mp.eps*(-d/2-xp)**2-q*Na/mp.eps*(-d/2-xp)*(x+d/2),
                 lambda x: Ecn+.5*q*Nd/mn.eps*(Ln+d/2-x)**2+q*sigma_n*(Ln+d/2-x)/mn.eps, lambda x: Ecn]),\
            "Ev":np.piecewise(x,
                [x<xp, (x>=xp)&(x<-d/2), (x>=-d/2)&(x<d/2), (x>=d/2)&(x<Ln+d/2), x>=Ln+d/2],
                [lambda x: Evp, lambda x: Evp-.5*q*Na/mp.eps*(x-xp)**2,
                 lambda x: Evp-.5*q*Na/mp.eps*(-d/2-xp)**2-q*Na/mp.eps*(-d/2-xp)*(x+d/2),
                 lambda x: Evn+.5*q*Nd/mn.eps*(Ln+d/2-x)**2+q*sigma_n*(Ln+d/2-x)/mn.eps, lambda x: Evn]),\
            "EB": np.piecewise(x,
                [x<xp, (x>=xp)&(x<xn), (x>=Ln+d/2)],
                [lambda x: -Va/2, lambda x: np.NaN, lambda x: Va/2]),\
            "Vr": -Va,\
            "material": mn}
    
# np.vectorize the above function so it actually expects an array for the applied
# voltage, and returns an array of the above return values
generate_pin_bandstructure=np.vectorize(generate_pin_bandstructure)

# Generates a PN band diagram using the depletion approximation
#
# mat: the material name (eg "Si") of this homojunction
# Lp, Aname, Na: the length, dopant name, and concentration (cm^-3) of the p-region
# Ln, Dname, Nd: the length, dopant name, and concentration (cm^-3) of the n-region
# Va: the applied voltage
#
# Returns a dictionary with x (the position grid), Ec and Ev (the band edges),
# EB (the Fermi level in the bulk regions), and Vr (reverse bias) and mat
# (the Material, from the p side, stored for convenience).
def generate_pn_bandstructure(mat, Lp,Aname,Na, Ln,Dname,Nd, Va):
    
    # The Material class will find the bulk energy levels
    mp=Material(mat,{Aname:Na})
    mn=Material(mat,{Dname:Nd})
    
    # Band edges far on the p and n side
    Ecp=mp.Ec-mp.EBulk-Va/2
    Evp=mp.Ev-mp.EBulk-Va/2
    Ecn=mn.Ec-mn.EBulk+Va/2
    Evn=mn.Ev-mn.EBulk+Va/2    
    
    # Depletion edges, assuming depletion fully within device.
    # Note, these are locations in x (with appropriate signs), where depletion ends
    # ...not the absolute widths of the respective regions!
    xn=np.sqrt(2*mn.eps*(Ecp-Ecn)/q*Na/(Nd*(Na+Nd)))
    xp=-np.sqrt(2*mp.eps*(Ecp-Ecn)/q*Nd/(Na*(Na+Nd)))
    
    # P-side depletion should stay within the device
    assert (xp > -Lp), "p-side depletion region extends outside device!"
    
    # x-grid for calculation
    x=np.linspace(-Lp,Ln,50000)
    
    # If the n-side depletion region is fully within the device
    if (xn<Ln):
        
        # Return np.piecewise band diagram from depletion approximation
        return {"x":x,\
            "Ec":np.piecewise(x,
                [x<xp, (x>=xp)&(x<0), (x>=0)&(x<xn), x>=xn],
                [lambda x: Ecp, lambda x: Ecp-.5*q*Na/mp.eps*(x-xp)**2,
                 lambda x:Ecn+.5*q*Nd/mn.eps*(xn-x)**2, Ecn]),\
            "Ev":np.piecewise(x,
                [x<xp, (x>=xp)&(x<0), (x>=0)&(x<xn), x>=xn],
                [lambda x: Evp, lambda x: Evp-.5*q*Na/mp.eps*(x-xp)**2,
                 lambda x:Evn+.5*q*Nd/mn.eps*(xn-x)**2, Evn]),\
            "EB": np.piecewise(x,
                [x<xp, (x>=xp)&(x<xn), (x>=xn)],
                [lambda x: -Va/2, lambda x: np.NaN, lambda x: Va/2]),\
            "Vr": -Va,\
            "material": mp}
    
    # If the n-side depletion region hits the device edge
    else:
        # Assume a sheet charge after the n-side to rebalance charge
        # Recalculate xp, and the sheet charge.        
        xp=-Ln*(np.sqrt(1+Nd/Na+2*mn.eps*(Ecp-Ecn)/(q*Na*Ln**2))-1)
        sigma_n=-Na*xp-Nd*Ln
        
        # Again, make sure the p-side depletion is still inside the device
        assert (xp > -Lp), "p-side depletion region extends outside device!"
        
        # Return a np.piecewise bandstructure given our maxed out depletion
        return {"x":x,\
            "Ec":np.piecewise(x,
                [x<xp, (x>=xp)&(x<0), (x>=0)&(x<xn), x>=xn],
                [lambda x: Ecp, lambda x: Ecp-.5*q*Na/mp.eps*(x-xp)**2,
                 lambda x:Ecn+.5*q*Nd/mn.eps*(Ln-x)**2+q*sigma_n*(Ln+-x)/mn.eps, Ecn]),\
            "Ev":np.piecewise(x,
                [x<xp, (x>=xp)&(x<0), (x>=0)&(x<xn), x>=xn],
                [lambda x: Evp, lambda x: Evp-.5*q*Na/mp.eps*(x-xp)**2,
                 lambda x:Evn+.5*q*Nd/mn.eps*(Ln-x)**2+q*sigma_n*(Ln-x)/mn.eps, Evn]),\
            "EB": np.piecewise(x,
                [x<xp, (x>=xp)&(x<xn), (x>=Ln)],
                [lambda x: -Va/2, lambda x: np.NaN, lambda x: Va/2]),\
            "Vr": -Va,\
            "material": mp}
    
# np.vectorize the above function so it actually expects an array for the applied
# voltage, and returns an array of the above return values
generate_pn_bandstructure=np.vectorize(generate_pn_bandstructure)

# Generates a PN band diagram with a donor delta spike using the depletion approximation
#
# mat: the material name (eg "Si") of this homojunction
# Lp, Aname, Na: the length, dopant name, and concentration (cm^-3) of the p-region
# sigmaD: sheet doping at the interface (cm^-2)
# Ln, Dname, Nd: the length, dopant name, and concentration (cm^-3) of the n-region
# Va: the applied voltage
#
# Returns a dictionary with x (the position grid), Ec and Ev (the band edges),
# EB (the Fermi level in the bulk regions), and Vr (reverse bias) and mat
# (the Material, from the p side, stored for convenience).
def generate_spikedpn_bandstructure(mat, Lp,Aname,Na, sigmaD, Ln,Dname,Nd, Va):
    
    # The Material class will find the bulk energy levels
    mp=Material(mat,{Aname:Na})
    mn=Material(mat,{Dname:Nd})
    
    # Band edges far on the p and n side
    Ecp=mp.Ec-mp.EBulk-Va/2
    Evp=mp.Ev-mp.EBulk-Va/2
    Ecn=mn.Ec-mn.EBulk+Va/2
    Evn=mn.Ev-mn.EBulk+Va/2    
    
    # x-grid for calculation
    x=np.linspace(-Lp,Ln,50000)
    
    # Depletion edges, assuming spike is undepleted
    sigma=np.sqrt(2*mp.eps*Na*(Ecp-Ecn)/q)
    print ("Va",Va," sigma", sigma)
    xp=-sigma/Na
    
    # If the spike is not fully depleted
    if(sigma<sigmaD):
        
        print ("not dep")
    
        # P-side depletion should stay within the device
        assert (xp > -Lp), "p-side depletion region extends outside device!"
           
        
        # Depletion edges, assuming depletion fully within device.
        # Note, these are locations in x (with appropriate signs), where depletion ends
        # ...not the absolute widths of the respective regions!
        return {"x":x,\
                "Ec":np.piecewise(x,
                    [x<xp,(x>=xp)&(x<0),(x>=0)],
                    [lambda x: Ecp, lambda x: Ecp-.5*q*Na/mp.eps*(x-xp)**2,lambda x: Ecn]),\
                "Ev":np.piecewise(x,
                    [x<xp,(x>=xp)&(x<0),(x>=0)],
                    [lambda x: Evp, lambda x: Evp-.5*q*Na/mp.eps*(x-xp)**2,lambda x: Evn]),\
                "EB": np.piecewise(x,
                    [x<xp, (x>=xp)&(x<0), (x>=0)],
                    [lambda x: -Va/2, lambda x: np.NaN, lambda x: Va/2]),\
                "Vr": -Va,\
                "material": mp}
    
    # If the spike is fully depleted, we actually have to solve for xn
    else:
        xn=-sigmaD/(Nd+Na)+np.sqrt(sigmaD**2/(Nd+Na)**2-2/(Nd+Na)*(sigmaD**2/(2*Nd)-Na*mn.eps*(Ecp-Ecn)/(q*Nd)))
        xp=-(Nd/Na*xn+sigmaD/Na)
        
        # N-side depletion should stay within the device, because I haven't considered full depletion
        assert (xn < Ln), "n-side depletion region extends outside device, and I haven't coded that"
                
        # Return np.piecewise band diagram from depletion approximation
        return {"x":x,\
            "Ec":np.piecewise(x,
                [x<xp, (x>=xp)&(x<0), (x>=0)&(x<xn), x>=xn],
                [lambda x: Ecp, lambda x: Ecp-.5*q*Na/mp.eps*(x-xp)**2,
                 lambda x:Ecn+.5*q*Nd/mn.eps*(xn-x)**2, Ecn]),\
            "Ev":np.piecewise(x,
                [x<xp, (x>=xp)&(x<0), (x>=0)&(x<xn), x>=xn],
                [lambda x: Evp, lambda x: Evp-.5*q*Na/mp.eps*(x-xp)**2,
                 lambda x:Evn+.5*q*Nd/mn.eps*(xn-x)**2, Evn]),\
            "EB": np.piecewise(x,
                [x<xp, (x>=xp)&(x<xn), (x>=xn)],
                [lambda x: -Va/2, lambda x: np.NaN, lambda x: Va/2]),\
            "Vr": -Va,\
            "material": mp}
    
# np.vectorize the above function so it actually expects an array for the applied
# voltage, and returns an array of the above return values
generate_spikedpn_bandstructure=np.vectorize(generate_spikedpn_bandstructure,otypes=[dict])



def generate_exactspikedpn_bandstructure(mat, Lp,Aname,Na, sigmaD, Ln,Dname,Nd, Va):
    import time
    t0=time.time()
    
    # The Material class will find the bulk energy levels
    mp=Material(mat,{Aname:Na})
    mn=Material(mat,{Dname:Nd})
    
    # Band edges far on the p and n side
    Ecp=mp.Ec-mp.EBulk-Va/2
    Evp=mp.Ev-mp.EBulk-Va/2
    Ecn=mn.Ec-mn.EBulk+Va/2
    Evn=mn.Ev-mn.EBulk+Va/2    
    Vr=-Va
    
    print ("EBulk N: ", mn.EBulk    )
    def build_interpolated_rho_ints(bounds,n=1000000):
        if n<10000000: print ("Low precision")
        from scipy.interpolate import interp1d
        bd=bounds[1]-bounds[0]
        l=bounds[0]-bd/100
        r=bounds[1]+bd/100
        Envals=np.linspace(l,r,n)
        dEn=Envals[1]-Envals[0]
        nrhos=mn.rho(Envals,ignoreMinority=True)
        intnrho=interp1d(Envals,-np.cumsum(nrhos)*dEn)
        
        nrhobulk=intnrho(mn.EBulk)
        
        Epvals=Envals+Vr
        dEp=dEn
        prhos=mp.rho(Epvals,ignoreMinority=True)
        intprho=interp1d(Epvals,-np.cumsum(prhos)*dEp)
        prhobulk=intprho(mp.EBulk) 
        
        return intnrho,intprho,nrhobulk,prhobulk
        
    bounds=[mn.EBulk-(Ecp-Ecn),mn.Ec+10*kT]
    intnrho,intprho,nrhobulk,prhobulk\
        =build_interpolated_rho_ints(bounds)
        
    print ("Bounds:")
    print (bounds)
    t1=time.time()
    print ("Building time: ", t1-t0)
    t0=t1
    
    def Fn(E0n,x=None):
        if x is not None:
            pass#print (E0n[0],E0n[-1]),"  ",x
        try:
            insqrt=2/mn.eps*(intnrho(E0n)-nrhobulk)
            #if np.where(np.abs(insqrt)<
            #assert np.all(insqrt>=-1), "Wrong sign in Fn"

            insqrt=np.clip(insqrt,0,None)            
            
            #fd=np.sign(E0n-mn.EBulk)*np.sqrt(insqrt)
            #print (fd[0],fd[-1])
        except Exception as e:
            print (e, E0n)
            return 0
        #    exit()
            
        #if x is not None:
        #    fd=np.sign(E0n-mn.EBulk)*np.sqrt(insqrt)
        #    print "    ",(fd[0],fd[-1])
            
            
            
        return np.sign(E0n-mn.EBulk)*np.sqrt(insqrt)
    def Fp(E0p,x=None):
        insqrt=2/mn.eps*(intprho(E0p)-prhobulk)
        insqrt=np.clip(insqrt,0,None)            
        #assert np.all(insqrt>=0), "Wrong sign in Fp"
        return np.sign(mp.EBulk-E0p)*np.sqrt(insqrt)
    
    
    
    dF=q*sigmaD/mn.eps
    def fdiff(E0n):
        #print "fdiff at E0n=",E0n
        try:
            fd= abs(Fn(E0n)-Fp(E0n+Vr)-dF)
            #print "      fdiff: ",fd, " (Fn: ",Fn(E0n), " )"
            return fd
        except:
            return np.inf

    import matplotlib.pyplot as mpl
    if 0:
        E0n=np.linspace(bounds[0],bounds[1],50)
        #E0n=np.linspace(1,2,50)
        fdiffs=[Fn(e)-Fp(e+Vr)-dF for e in E0n]
        
        #mpl.close('all')
        mpl.figure()
        mpl.plot(E0n,np.abs(fdiffs))
        mpl.yscale('log')

    if 1:
        #print "bounds: ",[mn.EBulk+mn.Eg,mn.EBulk-(Ecp-Ecn)]
        bd=bounds[1]-bounds[0]
        while (bd>1e-5):
            bd=bounds[1]-bounds[0]
            mincheck=min(.001,bd/10)
            res=sciopt.minimize_scalar(fdiff,bounds=bounds,method="bounded")
            E0n=res.x
            fd=fdiff(E0n)        
            E0p=E0n+Vr
            print ("min: ", fd)
            print ("Exact... E0n: ",E0n," E0p: ",E0p)
            print ("below: ",fdiff(E0n-mincheck))
            print ("above: ",fdiff(E0n+mincheck))
            assert ((fdiff(E0n-mincheck)>fd) and (fdiff(E0n+mincheck)>fd)), "Bad field minimization!"
            bounds=[E0n-mincheck,E0n+mincheck]
            bd=0
            #intnrho,intprho,nrhobulk,prhobulk\
                #=build_interpolated_rho_ints(bounds)
                
    t1=time.time()
    print ("Charge control Time: ", t1-t0)
    t0=t1
    
    print ("Fn: ",Fn(E0n)/1e6," MV/cm")
    print ("Fp: ",Fp(E0p)/1e6," MV/cm")
    #time.sleep(4)
    
    
    #Fn(E0n)
    
    #print "ED: ",Ecp-Ecn


    if 0:
        print ("Dep Approx")
        bs=generate_pn_bandstructure(mat,Lp,Aname,Na,Ln,Dname,Nd,[Va])[0]
        i0=np.argmin(np.abs(bs["x"]))
        E0n=Ecn-bs["Ec"][i0]+mn.EBulk
        E0p=Ecp-bs["Ec"][i0]+mp.EBulk
        Fmax=np.min(np.diff(bs["Ec"])/np.mean(np.diff(bs["x"])))
        
        print ("EBp", mp.EBulk, " EBn ", mn.EBulk)
        print ("Fmax ", Fmax)
        print ("E0n ", E0n, " Fn ", Fn(E0n))
        print ("E0p ", E0p, " Fp ", Fp(E0p))
        
    from scipy.integrate import odeint
    #def fn(y,x):
        #En=y[0]
        #F=y[1]
        #return [y[1],mn.rho(y[0],ignoreMinority=(y[0]<0))]
   #     return Fn(y)
    #odeint(fn,[E0n,Fn(E0n)])
    
    # x-grid for calculation
    x=np.linspace(-Lp,Ln,500000)
    zeroi=np.argmin(np.abs(x))
    zeroadded=False
    if(x[zeroi]):
        x=np.insert(x,zeroi if x[zeroi]>0 else zeroi+1,0)
        zeroadded=True
    zeroi=np.argmin(np.abs(x))
    
    xr=x[zeroi:]
    xl=x[zeroi::-1]
    x=np.linspace(-Lp,Ln,500000)
    
    
    with stdout_redirected(), merged_stderr_stdout():
        Ecr=Ecn+mn.EBulk-odeint(lambda E,x: -Fn(E,x),E0n,xr,rtol=1e-12,atol=1e-12)
        Ecl=Ecp+mp.EBulk-odeint(lambda E,x: -Fp(E,x),E0p,xl,rtol=1e-12,atol=1e-12)
    Ec=np.concatenate([Ecl[:0:-1],Ecr[zeroadded:]])
    Ef=(x<0)*(Vr/2)+(x>0)*(-Vr/2)
    
    Ev=Ec-mn.Eg
    
    
    #Ep=x[zeroi::-1]*0
    
    #Ep=odeint(Fp,E0p,x[zeroi::-1])[::-1]
    
    return {
        "x":x,
        "Ec": np.ravel(Ec),
        "Ev": np.ravel(Ev),
        "EB": Ef,
        "Vr": Vr,
        "material": mp
        }
    
    
# np.vectorize the above function so it actually expects an array for the applied
# voltage, and returns an array of the above return values
generate_exactspikedpn_bandstructure=np.vectorize(generate_exactspikedpn_bandstructure,otypes=[dict])




######################
######################
######################
def generate_exactspikedpin_bandstructure(mat, Lp,Aname,Na, Li, sigmaD, Ln,Dname,Nd, Va):
    import time
    t0=time.time()
    
    # The Material class will find the bulk energy levels
    mp=Material(mat,{Aname:Na})
    mn=Material(mat,{Dname:Nd})
    
    # Band edges far on the p and n side
    Ecp=mp.Ec-mp.EBulk-Va/2
    Evp=mp.Ev-mp.EBulk-Va/2
    Ecn=mn.Ec-mn.EBulk+Va/2
    Evn=mn.Ev-mn.EBulk+Va/2    
    Vr=-Va
    
    print ("EBulk N: ", mn.EBulk    )
    def build_interpolated_rho_ints(bounds,n=1000000):
        from scipy.interpolate import interp1d
        bd=bounds[1]-bounds[0]
        l=bounds[0]-bd/100
        r=bounds[1]+bd/100
        Envals=np.linspace(l,r,n)
        dEn=Envals[1]-Envals[0]
        nrhos=mn.rho(Envals,ignoreMinority=True)
        intnrho=interp1d(Envals,-np.cumsum(nrhos)*dEn)
        
        nrhobulk=intnrho(mn.EBulk)
        
        Epvals=Envals+Vr
        dEp=dEn
        prhos=mp.rho(Epvals,ignoreMinority=True)
        intprho=interp1d(Epvals,-np.cumsum(prhos)*dEp)
        prhobulk=intprho(mp.EBulk) 
        
        return intnrho,intprho,nrhobulk,prhobulk
        
    bounds=[mn.EBulk-(Ecp-Ecn),mn.Ec+15*kT]
    intnrho,intprho,nrhobulk,prhobulk\
        =build_interpolated_rho_ints(bounds)
        
    print ("Bounds:")
    print (bounds)
    t1=time.time()
    print ("Building time: ", t1-t0)
    t0=t1
    
    def Fn(E0n,x=None):
        if x is not None:
            pass#print (E0n[0],E0n[-1]),"  ",x
        try:
            insqrt=2/mn.eps*(intnrho(E0n)-nrhobulk)
            #if np.where(np.abs(insqrt)<
            #assert np.all(insqrt>=-1), "Wrong sign in Fn"

            insqrt=np.clip(insqrt,0,None)            
            
            #fd=np.sign(E0n-mn.EBulk)*np.sqrt(insqrt)
            #print (fd[0],fd[-1])
        except Exception as e:
            print (e, E0n)
            return 0
        #    exit()
            
        #if x is not None:
        #    fd=np.sign(E0n-mn.EBulk)*np.sqrt(insqrt)
        #    print "    ",(fd[0],fd[-1])
            
            
            
        return np.sign(E0n-mn.EBulk)*np.sqrt(insqrt)
    def Fp(E0p,x=None):
        insqrt=2/mn.eps*(intprho(E0p)-prhobulk)
        insqrt=np.clip(insqrt,0,None)            
        #assert np.all(insqrt>=0), "Wrong sign in Fp"
        return np.sign(mp.EBulk-E0p)*np.sqrt(insqrt)
    
    
    
    dF=q*sigmaD/mn.eps
    def fdiff(E0n):
        #print "fdiff at E0n=",E0n
        try:
            fn=Fn(E0n)
            fp=Fp(E0n+Vr+(fn-dF)*Li)
            #fd= abs(Fn(E0n)-Fp(E0n+Vr)-dF)
            fd= abs(fn-fp-dF)
            #print "      fdiff: ",fd, " (Fn: ",Fn(E0n), " )"
            return fd
        except:
            return np.inf

    import matplotlib.pyplot as mpl
    if 1:
        E0n=np.linspace(bounds[0]+.01,bounds[1]-.01,50)
        #E0n=np.linspace(1,2,50)

        #fdiffs=[Fn(e)-Fp(e+Vr)-dF for e in E0n]
        fdiffs=fdiff(E0n[5:-5])
        
        #mpl.close('all')
        mpl.figure()
        mpl.plot(E0n[5:-5],np.abs(fdiffs))
        mpl.yscale('log')

    if 1:
        #print "bounds: ",[mn.EBulk+mn.Eg,mn.EBulk-(Ecp-Ecn)]
        bd=bounds[1]-bounds[0]
        while (bd>1e-5):
            bd=bounds[1]-bounds[0]
            mincheck=min(.001,bd/10)
            res=sciopt.minimize_scalar(fdiff,bounds=bounds,method="bounded")
            E0n=res.x
            fd=fdiff(E0n)        
            E0p=E0n+Vr
            print ("min: ", fd)
            print ("Exact... E0n: ",E0n," E0p: ",E0p)
            print ("below: ",fdiff(E0n-mincheck))
            print ("above: ",fdiff(E0n+mincheck))
            assert ((fdiff(E0n-mincheck)>fd) and (fdiff(E0n+mincheck)>fd)), "Bad field minimization!"
            bounds=[E0n-mincheck,E0n+mincheck]
            bd=0
            #intnrho,intprho,nrhobulk,prhobulk\
                #=build_interpolated_rho_ints(bounds)
                
    t1=time.time()
    print ("Charge control Time: ", t1-t0)
    t0=t1
    
    print ("Fn: ",Fn(E0n)/1e6," MV/cm")
    print ("Fp: ",Fp(E0p)/1e6," MV/cm")
    #time.sleep(4)
    
    
    #Fn(E0n)
    
    #print "ED: ",Ecp-Ecn


    if 0:
        print ("Dep Approx")
        bs=generate_pn_bandstructure(mat,Lp,Aname,Na,Ln,Dname,Nd,[Va])[0]
        i0=np.argmin(np.abs(bs["x"]))
        E0n=Ecn-bs["Ec"][i0]+mn.EBulk
        E0p=Ecp-bs["Ec"][i0]+mp.EBulk
        Fmax=np.min(np.diff(bs["Ec"])/np.mean(np.diff(bs["x"])))
        
        print ("EBp", mp.EBulk, " EBn ", mn.EBulk)
        print ("Fmax ", Fmax)
        print ("E0n ", E0n, " Fn ", Fn(E0n))
        print ("E0p ", E0p, " Fp ", Fp(E0p))
        
    from scipy.integrate import odeint
    #def fn(y,x):
        #En=y[0]
        #F=y[1]
        #return [y[1],mn.rho(y[0],ignoreMinority=(y[0]<0))]
   #     return Fn(y)
    #odeint(fn,[E0n,Fn(E0n)])
    
    # x-grid for calculation
    x=np.linspace(-Lp-Li,Ln,500000)
    zeroi=np.argmin(np.abs(x))
    zeroadded=False
    if(x[zeroi]):
        x=np.insert(x,zeroi if x[zeroi]>0 else zeroi+1,0)
        zeroadded=True
    zeroi=np.argmin(np.abs(x))
    
    xr=x[zeroi:]
    xl=x[zeroi::-1]
    x=np.linspace(-Lp-Li,Ln,500000)
    
    print("care")
    print(E0p)
    print((lambda E,x: -Fp(E,x))(E0p,xl))
    print((lambda E,x: np.choose(x<-Li,[[-Fn(E0n,0)+dF]*len(x),-Fp(E,x)]))(E0p,xl))
    with stdout_redirected(), merged_stderr_stdout():
        Ecr=Ecn+mn.EBulk-odeint(lambda E,x: -Fn(E,x),E0n,xr,rtol=1e-12,atol=1e-12)
            
        Ecl=Ecp+mp.EBulk-odeint(lambda E,x: np.choose(x<-Li,[(-Fn(E0n,0)+dF)*(1+0*x),-Fp(E,x)]),E0p,xl,rtol=1e-12,atol=1e-12)
        #Ecl=Ecp+mp.EBulk-odeint(lambda E,x: -Fp(E,x),E0p,xl,rtol=1e-12,atol=1e-12)
    Ec=np.concatenate([Ecl[:0:-1],Ecr[zeroadded:]])
    Ef=(x<0)*(Vr/2)+(x>0)*(-Vr/2)
    
    Ev=Ec-mn.Eg
    
    
    #Ep=x[zeroi::-1]*0
    
    #Ep=odeint(Fp,E0p,x[zeroi::-1])[::-1]
    
    return {
        "x":x+Li/2,
        "Ec": np.ravel(Ec),
        "Ev": np.ravel(Ev),
        "EB": Ef,
        "Vr": Vr,
        "material": mp
        }
    
    
# np.vectorize the above function so it actually expects an array for the applied
# voltage, and returns an array of the above return values
generate_exactspikedpin_bandstructure=np.vectorize(generate_exactspikedpin_bandstructure,otypes=[dict])

##########################
##########################
##########################




















import re
def import_bd(filename,xmin=-1e10,xmax=1e10):
    with open(filename) as f:
        data=[]
        for l in f:
            mo=re.match("\s*([\d\.e+-]+)\s+([\d\.e+-]+)\s+([\d\.e+-]+)\s+([\d\.e+-]+)\s+"\
                        "([\d\.e+-]+)\s+([\d\.e+-]+)\s+([\d\.e+-]+)\s+([\d\.e+-]+)\s+", l)
            if mo:
                x=float(mo.groups()[0])
                if(x>=xmin and x<=xmax):
                    data+=[[float(a) for a in mo.groups()]]
                    data[-1][0]=(x-(xmin+xmax)/2)/1e8
        data=np.array(data)
        #print data
        bd={}
        bd["x"]=data[:,0]
        bd["Ec"]=data[:,1]
        bd["Ev"]=data[:,2]
        bd["E"]=data[:,3]
        bd["Ef"]=data[:,4]
        bd["n"]=data[:,5]
        bd["p"]=data[:,6]
        bd["material"]=Material("GaN")
    return bd