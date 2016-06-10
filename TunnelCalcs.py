from __future__ import division
from Constants import eps0,q,hbar,m0,kT,eVToJoules
import numpy as np
from GaN_BS_min import kvl


# Returns the tunneling current from a band diagram bs
def J(bs):
    
    # Returns the integrand (np.meant to be integrated over Ep and Ex) for the tunneling
    # current, including both the exponential and the prefactor.
    #
    # Ex: (hbar k_x)^2/2m
    # Ep: (hbar k_perpendicular)^2/2m
    # bs: the band diagram
    def integrand(Ex,Ep,bs):
        mr=bs["material"].mrtunnel        
        
        # Find the indices of the classical turning points
        E=Etop-Ex-Ep-Ebot
        try:
            il=np.where(Etop-Ex>bs["Ev"])[0][0]
            ir=np.where(Ebot+E-Ep>bs["Ec"])[0][0]
        except:
            print 'no int bounds'
            return 0
        i=np.arange(il,ir)
        
        # Log of WKB tunneling probability
        lt= -2*np.trapz(
        np.sqrt(2*mr*(bs["Ec"][i]-Ebot+Ep-E)*eVToJoules)\
            /(hbar*eVToJoules) / 100,\
            bs["x"][i])
        
        # Log of prefactor
        lpref=np.log(q*mr/(2*np.pi**2*hbar**3*eVToJoules)) # log ( A m^-2 eV^-2)
        
        # Combine the prefactor and exponential and return
        return np.exp(lt+lpref)
    
    # Returns the inner integral (ie over Ep) for the tunneling current
    #
    # Ex: (hbar k_x)^2/2m
    # bs: the band diagram
    def perp_int(Ex,bs):
        
        # Ep gets integrated up to 1/2 * (Emax-Ex)
        Ep=np.linspace(0,.5*(Etop-Ebot-Ex),100,endpoint=False)
        return np.trapz([integrand(Ex,ep,bs) for ep in Ep],Ep)
    
    # Find the top of top and bottom of the energy range
    v=bs["Ev"]
    c=bs["Ec"]
    Etop=v[np.where(v<max(v))[0][0]]
    Ebot=c[np.where(c>min(c))[0][-1]]
    Vr=bs["Vr"]
    
    if(Etop>Vr/2 or Ebot<-Vr/2):
        print "boundaries of energy integration are invalid because of degeneracy"
    
    # Ex can run over the entire range
    Ex=np.linspace(0,Etop-Ebot,100,endpoint=False)
    
    # Return the outer integral
    return np.trapz([perp_int(ex,bs) for ex in Ex],Ex)

# Returns the tunneling current from a band diagram bs
def Jv2(bs):
    
    def integrand(E,Ep):
        mr=bs["material"].mrtunnel        
        
        # Find the indices of the classical turning points
        try:
            il=np.where(E+Ep>bs["Ev"])[0][0]
            ir=np.where(E-Ep>bs["Ec"])[0][0]
        except:
            print 'no int bounds for E=', E, ' Ep= ', Ep
            return 0
        i=np.arange(il,ir)
        
        # Log of WKB tunneling probability
        lt= -2*np.trapz(
        np.sqrt(2*mr*(bs["Ec"][i]+Ep-E)*eVToJoules)\
            /(hbar*eVToJoules) / 100,\
            bs["x"][i])
        print( "E: {:6.3f} Ep: {:6.3f} logT: {:6.3f}".format(E,Ep,lt/np.log(10)))
        
        # Log of prefactor
        lpref=np.log(q*mr/(2*np.pi**2*hbar**3*eVToJoules)) # log ( A m^-2 eV^-2)
        
        # Combine the prefactor and exponential and return
        return np.exp(lt+lpref)
    

    def int_at_fixed_Ep(Ep):
        
        E=np.linspace(Emin+Ep,Emax-Ep,101,endpoint=False)[1:]
        #print('about to do E for Ep= '+str(Ep))
        #print(E)
        return np.trapz(np.vectorize(lambda e: integrand(e,Ep),otypes=['float64'])(E),E)
    
    # Find the top of top and bottom of the energy range
    Emax=min(bs["Ev"][0],bs["Vr"]/2)
    Emin=max(bs["Ec"][-1],-bs["Vr"]/2)
    print Emax-Emin
    
    Ep=np.linspace(0,.5*(Emax-Emin),100,endpoint=False)
    return np.trapz(np.vectorize(int_at_fixed_Ep,otypes=['float64'])(Ep),Ep)   
    
    '''
    Etop=v[np.where(v<max(v))[0][0]]
    Ebot=c[np.where(c>min(c))[0][-1]]
    
    if(Etop>Vr/2 or Ebot<-Vr/2):
        print "boundaries of energy integration are invalid because of degeneracy"
    
    # Ex can run over the entire range
    Ex=np.linspace(0,Etop-Ebot,100,endpoint=False)
    
    # Return the outer integral
    return np.trapz([perp_int(ex,bs) for ex in Ex],Ex)
    '''


# Returns the maximum field in a band diagram bs
def max_field(bs):
    return np.max(-np.diff(bs["Ec"])/np.diff(bs["x"]))

# Returns the tunneling current applying our analytical np.expression to the maximum field
def J_tri(bs):
    xi=max_field(bs)
    mr=bs["material"].mrtunnel
    Eg=bs["material"].Eg
    Emax=min(bs["Ev"][0],bs["Vr"]/2)
    Emin=max(bs["Ec"][-1],-bs["Vr"]/2)
    return np.sqrt(2*mr)*q*xi*100*(Emax-Emin)/(16*np.pi**2*hbar**2*(Eg*eVToJoules)**.5)\
        *np.exp(-4*np.sqrt(2*mr)*Eg**1.5/(3*xi*100*hbar*eVToJoules**.5))

# Returns the tunneling current applying our analytical np.expression to the maximum field
def J_para(bs):
    xi=max_field(bs)
    mr=bs["material"].mrtunnel
    Eg=bs["material"].Eg
    Emax=min(bs["Ev"][0],bs["Vr"]/2)
    Emin=max(bs["Ec"][-1],-bs["Vr"]/2)
    return np.sqrt(2*mr)*q*xi*100*(Emax-Emin)/(3*np.pi**3*hbar**2*(Eg*eVToJoules)**.5)\
        *np.exp(-np.pi*np.sqrt(mr)*Eg**1.5/(2*np.sqrt(2)*xi*100*hbar*eVToJoules**.5))

# Returns the tunneling current applying our analytical np.expression to the np.mean field
# within the junction region
def J_trimean(bs):
    Efield=np.diff(bs["Ec"])/np.diff(bs["x"])
    xi=np.mean(-Efield[np.where(Efield!=0)])
    mr=bs["material"].mrtunnel
    Eg=bs["material"].Eg
    Emax=min(bs["Ev"][0],bs["Vr"]/2)
    Emin=max(bs["Ec"][-1],-bs["Vr"]/2)
    return np.sqrt(2*mr)*q*xi*100*(Emax-Emin)/(16*np.pi**2*hbar**2*(Eg*eVToJoules)**.5)\
        *np.exp(-4*np.sqrt(2*mr)*Eg**1.5/(3*xi*100*hbar*eVToJoules**.5))

def J_exp(DE,xi,bs):
    mr=bs["material"].mrtunnel
    Eg=bs["material"].Eg
    if not(DE):
        print DE, xi, bs
        Emax=min(bs["Ev"][0],bs["Vr"]/2)
        Emin=max(bs["Ec"][-1],-bs["Vr"]/2)
        DE=Emax-Emin
    return np.sqrt(2*mr)*q*xi*100*DE/(16*np.pi**2*hbar**2*(Eg*eVToJoules)**.5)\
        *np.exp(-4*np.sqrt(2*mr)*Eg**1.5/(3*xi*100*hbar*eVToJoules**.5))
        
def J_d(bs):
    mr=bs["material"].mrtunnel
    Eg=bs["material"].Eg
    Emax=min(bs["Ev"][0],bs["Vr"]/2)
    Emin=max(bs["Ec"][-1],-bs["Vr"]/2)
    E=np.linspace(Emin,Emax,10001,endpoint=False)[1:]
    d=np.min(
        [    bs["x"][np.argmax(e>bs["Ec"])]
            -bs["x"][np.argmax(e>bs["Ev"])]
        for e in E])
    xi=Eg/d
    return np.sqrt(2*mr)*q*xi*100*(Emax-Emin)/(16*np.pi**2*hbar**2*(Eg*eVToJoules)**.5)\
        *np.exp(-4*np.sqrt(2*mr)*Eg**1.5/(3*xi*100*hbar*eVToJoules**.5))
    
def Jb(bs):
    mr=bs["material"].mrtunnel
    Eg=bs["material"].Eg
    Emax=min(bs["Ev"][0],bs["Vr"]/2)
    Emin=max(bs["Ec"][-1],-bs["Vr"]/2)
    

# Convenience function, calls the appropriate function of the above three when supplied
# a band diagram and a method: "integral", "np.meanfield", or "maxfield"
def current(bs,method="integral"):
    method={"orig": J, "v2": Jv2,"para":J_para, "intpara":J_intpara, "imag":Jrb, "d": J_d, "integral": J, "meanfield": J_trimean, "maxfield": J_tri}[method]
    return method(bs)
    
    
def dOfE(bs):
    Emax=min(bs["Ev"][0],bs["Vr"]/2)
    Emin=max(bs["Ec"][-1],-bs["Vr"]/2)
    E=np.linspace(Emin,Emax,10001,endpoint=False)[1:]
    d=np.array(
        [    bs["x"][np.argmax(e>bs["Ec"])]
            -bs["x"][np.argmax(e>bs["Ev"])]
        for e in E])
    return E,d
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
# Returns the tunneling current from a band diagram bs
def Jrb(bs):
    
    def integrand(E,Ep):
        mr=bs["material"].mrtunnel        
        
        # Find the indices of the classical turning points
        try:
            il=np.where(E+Ep>bs["Ev"])[0][0]
            ir=np.where(E-Ep>bs["Ec"])[0][0]
        except:
            print 'no int bounds for E=', E, ' Ep= ', Ep
            return 0
        i=np.arange(il,ir)
        
        # Log of WKB tunneling probability
        lt= -2*np.trapz(
        kvl(E-bs["Ev"][i],Ep),
            bs["x"][i])
        #print( "E: {:6.3f} Ep: {:6.3f} logT: {:6.3f}".format(E,Ep,lt/np.log(10)))
        
        # Log of prefactor
        lpref=np.log(q*mr/(2*np.pi**2*hbar**3*eVToJoules)) # log ( A m^-2 eV^-2)
        
        # Combine the prefactor and exponential and return
        return np.exp(lt+lpref)
    

    def int_at_fixed_Ep(Ep):
        
        E=np.linspace(Emin+Ep,Emax-Ep,101,endpoint=False)[1:]
        #print('about to do E for Ep= '+str(Ep))
        #print(E)
        return np.trapz(np.vectorize(lambda e: integrand(e,Ep),otypes=['float64'])(E),E)
    
    # Find the top of top and bottom of the energy range
    Emax=min(bs["Ev"][0],bs["Vr"]/2)
    Emin=max(bs["Ec"][-1],-bs["Vr"]/2)
    print Emax-Emin
    
    Ep=np.linspace(0,.5*(Emax-Emin),100,endpoint=False)
    return np.trapz(np.vectorize(int_at_fixed_Ep,otypes=['float64'])(Ep),Ep)   
    
    '''
    Etop=v[np.where(v<max(v))[0][0]]
    Ebot=c[np.where(c>min(c))[0][-1]]
    
    if(Etop>Vr/2 or Ebot<-Vr/2):
        print "boundaries of energy integration are invalid because of degeneracy"
    
    # Ex can run over the entire range
    Ex=np.linspace(0,Etop-Ebot,100,endpoint=False)
    
    # Return the outer integral
    return np.trapz([perp_int(ex,bs) for ex in Ex],Ex)
    '''



# Returns the tunneling current from a band diagram bs
def J_intpara(bs):
    
    def integrand(E,Ep):
        mr=bs["material"].mrtunnel
        Eg=bs["material"].Eg
        
        # Find the indices of the classical turning points
        try:
            il=np.where(E+Ep>bs["Ev"])[0][0]
            ir=np.where(E-Ep>bs["Ec"])[0][0]
        except:
            print 'no int bounds for E=', E, ' Ep= ', Ep
            return 0
        i=np.arange(il,ir)
        
        print "here goes"
        print Eg
        print Ep
        print bs["Ev"][il]
        print bs["Ec"][ir]
        print  np.sqrt(2*mr*(((Eg+2*Ep)**2/4-(E-bs["Ev"][i]+Ep-(Eg+2*Ep)/2)**2)/(Eg+2*Ep))*eVToJoules)\
            /(hbar*eVToJoules) / 100
        # Log of WKB tunneling probability
        lt= -2*np.trapz(
        np.sqrt(2*mr*(((Eg+2*Ep)**2/4-(E-bs["Ev"][i]+Ep-(Eg+2*Ep)/2)**2)/(Eg+2*Ep))*eVToJoules)\
            /(hbar*eVToJoules) / 100,\
            bs["x"][i])
        print( "E: {:6.3f} Ep: {:6.3f} logT: {:6.3f}".format(E,Ep,lt/np.log(10)))
        if(np.isnan(lt)):
            stop
        # Log of prefactor
        lpref=np.log(q*mr/(2*np.pi**2*hbar**3*eVToJoules)) # log ( A m^-2 eV^-2)
        
        # Combine the prefactor and exponential and return
        return np.exp(lt+lpref)
    

    def int_at_fixed_Ep(Ep):
        
        E=np.linspace(Emin+Ep,Emax-Ep,101,endpoint=False)[1:]
        #print('about to do E for Ep= '+str(Ep))
        #print(E)
        return np.trapz(np.vectorize(lambda e: integrand(e,Ep),otypes=['float64'])(E),E)
    
    # Find the top of top and bottom of the energy range
    Emax=min(bs["Ev"][0],bs["Vr"]/2)
    Emin=max(bs["Ec"][-1],-bs["Vr"]/2)
    print Emax-Emin
    
    Ep=np.linspace(0,.5*(Emax-Emin),100,endpoint=False)
    return np.trapz(np.vectorize(int_at_fixed_Ep,otypes=['float64'])(Ep),Ep)   
    
    '''
    Etop=v[np.where(v<max(v))[0][0]]
    Ebot=c[np.where(c>min(c))[0][-1]]
    
    if(Etop>Vr/2 or Ebot<-Vr/2):
        print "boundaries of energy integration are invalid because of degeneracy"
    
    # Ex can run over the entire range
    Ex=np.linspace(0,Etop-Ebot,100,endpoint=False)
    
    # Return the outer integral
    return np.trapz([perp_int(ex,bs) for ex in Ex],Ex)
    '''

