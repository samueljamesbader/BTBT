from __future__ import division
import scipy.optimize as sciopt
from copy import deepcopy
import numpy as np
from Constants import eps0,q,hbar,m0,kT,eVToJoules

# Dictionary of parameters for materials used in these calculations
materialparameters={
                    
    # Gallium Nitride
    "GaN": {
            
        # Permittivity
        "eps": 9.7*eps0,
        
        # Bandgap
        "Eg": 3.4, # eV
        
        # Effective masses for DOS calculations
        "medos": .2*m0, # kg
        "mhdos": 1.5*m0, # kg
        
        # Effective mass for tunneling (Sze convention)
        "mrtunnel": 2/(1/.2+1/.259)*m0, # kg
        
        # Dopant types and depths
        "dopants": {
            "Si": {
                "type": "donor",
                "E":.015 # eV
            },
            "Mg": {
                "type": "acceptor",
                "E": .23 # eV
            }
        }
    },
    
    # Silicon
    "Si": {
           
        # Permittivity
        "eps": 11.7*eps0,
        
        # Bandgap
        "Eg": 1.1, # eV
        
        # Effective masses for DOS calculations
        "medos": 1.08*m0, # kg
        "mhdos": .81*m0, # kg
        
        # Effective mass for tunneling (Sze convention)
        "mrtunnel": .16*m0, # kg
        
        # Dopant types and depths
        "dopants": {
            "P": {
                "type": "donor",
                "E": .045 # eV
            },
            "B": {
                "type": "acceptor",
                "E": .045 #eV
            }
        }
    }
}

# Convenience class that represents a chunk of a specific material with uniform doping
class Material:
    
    # Initialize the material with a name and doping, and oes the bulk energy level calculations
    # 
    # name: name of the material (eg. "Si")
    #     uses this to pull data from materialparamters dictionary
    # doping: dictionary of the form {'dopantname': concentration_per_cm3}
    #     'dopantname' should match a listed dopant for this material in materialparams
    def __init__(self,name,doping={}):
        
        # Take each key from the dictionary and make it an attribute of self
        for k,v in materialparameters[name].iteritems():
            setattr(self,k,deepcopy(v))
            
        # Self should have a dictionary of non-zero dopants with their parameters from 
        # materialparamters as well as their concentrations given at initialization.
        for k,v in self.dopants.items():
            if (k in doping and doping[k]!=0):
                v["conc"]=doping[k]
            else:
                self.dopants.pop(k)
        
        # Calculate intrinsic carrier concentration
        self.Nc=2*(self.medos*kT/(2*np.pi*hbar**2 * eVToJoules))**(3/2)/100**3 # cm^-3
        self.Nv=2*(self.mhdos*kT/(2*np.pi*hbar**2 * eVToJoules))**(3/2)/100**3 # cm^-3
        self.ni=np.sqrt(self.Nc*self.Nv)*np.exp(-self.Eg/(2*kT))
        
        # With zero of energy defined by the *intrinsic* Fermi level...
        
        # Calculate the energy of mid-gap 
        self.Emg=-3/4*kT*np.log(self.mhdos/self.medos)
        
        # Calculate the conduction and valence bands
        self.Ec= self.Eg/2+self.Emg
        self.Ev=-self.Eg/2+self.Emg
        
        # Calculate where the (extrinsic) Fermi level should be in the bulk
        self._recomputeEBulk()
    
    # Enforces charge neutrality to find the location of the (extrinsic) Fermi level
    # referenced to intrinsic Fermi level, sets it as self.EBulk
    def _recomputeEBulk(self):
        
        # If undoped, EBulk=0 (ie intrinsic)
        if not len(self.dopants.values()):
            self.EBulk=0
        else:
            # Find the net doping
            netDope=sum([
                v["conc"]*(-1 if v["type"]=="donor" else 1)
                    for v in self.dopants.values()])
            # Use the net doping with complete ionization to get a starting guess (E0) for EBulk
            if(not netDope): E0=0
            else: E0=np.sign(-netDope)*kT*np.log(abs(netDope)/self.ni)
                
            # Find EBulk to minimize the net charge in the material
            res=sciopt.minimize_scalar(lambda E:abs(self.rho(E)),sorted([0,E0]),method='Golden')
            self.EBulk=res.x
    
    # Returns the electron density when the Fermi level is E (referenced to intrinsic)
    def n(self,E):
        return self.ni*np.exp(E/kT)
        return np.choose((self.Ec-E)<5*kT, [
            self.ni*np.exp(E/kT),
            self.Nc*FOneHalf(np.array(E)-self.Ec)
        ])
        
        
        #return self.Nc*FOneHalf(E,self.Ec)
        #return self.Nc*FOneHalf(E,self.Ec)
    
    # Returns the hole density when the Fermi level is E (referenced to intrinsic)
    def p(self,E):
        return self.ni*np.exp(-E/kT)
        return np.choose((E-self.Ev)<5*kT, [
            self.ni*np.exp(-E/kT),
            self.Nv*FOneHalf(self.Ev-np.array(E))
        ])
        
        
        
        
        #return self.Nv*FOneHalf(self.Ev,E)
    
    # Returns the ionized donor concentration when the Fermi level is E (referenced to intrinsic)
    def Nd_ionized(self,E):
        #return sum([ d["conc"]*(1-1/(1+.5*np.exp((self.Ec-d["E"]-E)/kT))) for d in self.dopants.values()
        #            if d["type"]=="donor"])
        d=self.dopants.values()[0]
        if d["type"]=="donor":
            return d["conc"]*(1-1/(1+.5*np.exp((self.Ec-d["E"]-E)/kT)))
        else: return 0;
    
    # Returns the ionized acceptor concentration when the Fermi level is E (referenced to intrinsic)
    def Na_ionized(self,E):
        #return sum([ d["conc"]*(1/(1+4*np.exp((self.Ev+d["E"]-E)/kT))) for d in self.dopants.values()
                    #if d["type"]=="acceptor"])
        d=self.dopants.values()[0]
        if d["type"]=="acceptor":
            return d["conc"]*(1/(1+4*np.exp((self.Ev+d["E"]-E)/kT)))
        else: return 0;
    
    # Returns the net charge density when the Fermi level is E (referenced to intrinsic)
    def rho(self,E,ignoreMinority=False):
        #assert self.Ec-E>3*kT and E-self.Ev>3*kT, "DEGENERATE!!!"
        #if (np.any(self.Ec-E<3*kT) or np.any(E-self.Ev<3*kT)):
            #print "DEGENERATE!!!"
    
        
            
        if(ignoreMinority):
            if(self.EBulk>0):
                return q*( self.Nd_ionized(E)-self.n(E))
            else:
                return q*(-self.Na_ionized(E)+self.p(E))
        else:
            return q*(self.p(E)+self.Nd_ionized(E)-self.n(E)-self.Na_ionized(E))

def FOneHalf(E):
    #print np.array((E-Eb)/kT)
    try:
        len(E)
    except:
        return fermihalf(E/kT,1)
    return np.vectorize(lambda eta: fermihalf(eta,1))(np.array(E/kT))
    

# http://www.scientificpython.net/pyblog/approximate-fermi-dirac-integrals
def fermihalf(x,sgn):
    """ Series approximation to the F_{1/2}(x) or F_{-1/2}(x) 
        Fermi-Dirac integral """

    f = lambda k: np.sqrt(x**2+np.pi**2*(2*k-1)**2)


    if sgn>0: # F_{1/2}(x)
        a = np.array((1.0/770751818298,-1.0/3574503105,-13.0/184757992,
              85.0/3603084,3923.0/220484,74141.0/8289,-5990294.0/7995))
        g = lambda k:np.sqrt(f(k)-x)

    else:  # F_{-1/2}(x)
        a = np.array((-1.0/128458636383,-1.0/714900621,-1.0/3553038,
                      27.0/381503,3923.0/110242,8220.0/919))
        g = lambda k:-0.5*np.sqrt(f(k)-x)/f(k)

    F = np.polyval(a,x) + 2*np.sqrt(2*np.pi)*sum(map(g,range(1,21)))
    return F
    