from __future__ import division
import scipy.optimize as sciopt
from copy import deepcopy
import numpy as np
from Constants import eps0,q,hbar,m0,kT,eVToJoules
from FermiDirac import fermi_dirac_one_half, int_fermi_dirac_one_half, fermi_dirac_zero
# Dictionary of parameters for materials used in these calculations

#kT=77/300 * kT

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
        "mrtunnel": 2/(1/.2+1/1.78)*m0, # kg
        
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
                    
    # Gallium Arsenide
    "GaAs": {
            
        # Permittivity
        "eps": 12.9*eps0,
        
        # Bandgap
        "Eg": 1.42, # eV
        
        # Effective masses for DOS calculations
        "medos": .063*m0, # kg
        "mhdos": .53*m0, # kg
        
        # Effective mass for tunneling (Sze convention)
        #"mrtunnel": 2/(1/+1/)*m0, # kg
        
        # Dopant types and depths
        "dopants": {
            "S": {
                "type": "donor",
                "E":.006 # eV
            },
            "C": {
                "type": "acceptor",
                "E": .02 # eV
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
    
    # Initialize the material with a name and doping, and do the bulk chemical potential calculations
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
        self._recomputeMuBulk()
    
    # Enforces charge neutrality to find the location of the (extrinsic) chemical potential
    # referenced to intrinsic Fermi level, sets it as self.EBulk
    def _recomputeMuBulk(self):
        
        # If undoped, EBulk=0 (ie intrinsic)
        if not len(self.dopants.values()):
            self.mu_bulk=0
        else:
            # Find the net doping
            netDope=sum([
                v["conc"]*(-1 if v["type"]=="donor" else 1)
                    for v in self.dopants.values()])
            # Use the net doping with complete ionization to get a starting guess (E0) for EBulk
            if(not netDope): mu_0=0
            else: mu_0=np.sign(-netDope)*kT*np.log(abs(netDope)/self.ni)
                
            # Find EBulk to minimize the net charge in the material
            res=sciopt.minimize_scalar(lambda mu:abs(self.rho(mu,mu)),sorted([0,mu_0]),method='Golden')
            self.mu_bulk=res.x
        


            eterm=-q*kT*self.Nc*int_fermi_dirac_one_half((self.mu_bulk-self.Ec)/kT)
            hterm=-q*kT*self.Nv*int_fermi_dirac_one_half((self.Ev-self.mu_bulk)/kT)
            if(len(self.dopants)!=1): return
            d=self.dopants.values()[0]
            if d["type"]=="donor":
                dterm=-q*kT*d["conc"]*fermi_dirac_zero(-(self.mu_bulk-self.Ec+d["E"]+np.log(2)*kT)/kT)
            elif d["type"]=="acceptor":
                dterm=-q*kT*d["conc"]*fermi_dirac_zero(-(self.Ev+d["E"]-self.mu_bulk+np.log(4)*kT)/kT)
            self.intrho_ref=eterm+hterm+dterm
        
        #assert(len(self.dopants)==1),"Should have one dopant"
            
#        if(len(self.dopants)!=1): return
#        d=self.dopants.values()[0]
#        if d["type"]=="donor":
#            self._intrho_ref=\
#                -q*kT*(self.Nc*int_fermi_dirac_one_half((self.mu_bulk-self.Ec)/kT)
#                    +d["conc"]*fermi_dirac_zero(-(self.mu_bulk-self.Ec+d["E"]+np.log(2)*kT)/kT))
#            #self._minorrho=q*(self.p(self.mu_bulk)-self.Na_ionized(self.mu_bulk))
#            
#        elif d["type"]=="acceptor":
#            self._intrho_ref=\
#                -q*kT*(self.Nv*int_fermi_dirac_one_half((self.Ev-self.mu_bulk)/kT)
#                    +d["conc"]*fermi_dirac_zero(-(self.Ev+d["E"]-self.mu_bulk+np.log(4)*kT)/kT))
#            self._minorrho=-q*(self.n(self.mu_bulk)-self.Nd_ionized(self.mu_bulk))


            
#        if d["type"]=="acceptor":
#            self._intrho_ref=\
#                q*(self.Nv*int_fermi_dirac_one_half((self.Ev-self.mu_bulk)/kT)
#                    -d["conc"]*fermi_dirac_zero((self.Ec-d["Ed"]-self.mu_bulk)/kT))
    
    # Returns the electron density when the Fermi level is E (referenced to intrinsic)
    def n(self,mu):
        return self.Nc*fermi_dirac_one_half((mu-self.Ec)/kT)
    
    # Returns the hole density when the Fermi level is E (referenced to intrinsic)
    def p(self,mu):
        return self.Nv*fermi_dirac_one_half((self.Ev-mu)/kT)
    
    # Returns the ionized donor concentration when the Fermi level is E (referenced to intrinsic)
    def Nd_ionized(self,mu):
        #return sum([ d["conc"]*(1-1/(1+.5*np.exp((self.Ec-d["E"]-E)/kT))) for d in self.dopants.values()
        #            if d["type"]=="donor"])
        assert(len(self.dopants)==1),"Should have one dopant"
        d=self.dopants.values()[0]
        if d["type"]=="donor":
            return d["conc"]*(1/(1+2*np.exp((mu-self.Ec+d["E"])/kT)))
        else: return 0;
    
    # Returns the ionized acceptor concentration when the Fermi level is E (referenced to intrinsic)
    def Na_ionized(self,mu):
        #return sum([ d["conc"]*(1/(1+4*np.exp((self.Ev+d["E"]-E)/kT))) for d in self.dopants.values()
                    #if d["type"]=="acceptor"])
        assert(len(self.dopants)==1),"Should have one dopant"
        d=self.dopants.values()[0]
        if d["type"]=="acceptor":
            return d["conc"]*(1/(1+4*np.exp((self.Ev+d["E"]-mu)/kT)))
        else: return 0;
    
    # Returns the net charge density when the Fermi level is E (referenced to intrinsic)
    def rho(self,EFn,EFp):            
        #print "ignore minority should probs include minority but eval'd at mu_bulk"
#        if(ignoreMinority):
#            if(self.mu_bulk>0):
#                return q*( self.Nd_ionized(mu)-self.n(mu)
#                    -self.Na_ionized(self.mu_bulk)+self.p(self.mu_bulk))
#            else:
#                return q*(-self.Na_ionized(mu)+self.p(mu)
#                    +self.Nd_ionized(self.mu_bulk)-self.n(self.mu_bulk))
#        else:
#            return q*(self.p(mu)+self.Nd_ionized(mu)-self.n(mu)-self.Na_ionized(mu))
        return q*(self.p(EFp)+self.Nd_ionized(EFn)-self.n(EFn)-self.Na_ionized(EFp))

    # int dmu from majority,minority=mu_bulk to majority,minority=given values
    def intrho(self,EFn,EFp):   
        assert(len(self.dopants)==1),"Should have one dopant"
        #d=self.dopants.values()[0]
#        if d["type"]=="donor":
#            intrho_maj=\
#                -q*kT*(self.Nc*int_fermi_dirac_one_half((mu-self.Ec)/kT)
#                    +d["conc"]*fermi_dirac_zero(-(mu-self.Ec+d["E"]+np.log(2)*kT)/kT))
#        elif d["type"]=="acceptor":
#            intrho_maj=\
#                -q*kT*(self.Nv*int_fermi_dirac_one_half((self.Ev-mu)/kT)
#                    +d["conc"]*fermi_dirac_zero(-(self.Ev+d["E"]-mu+np.log(4)*kT)/kT))
            #print 'sup: ',intrho_maj[-1], self._intrho_ref, self.Nc*int_fermi_dirac_one_half((mu-self.Ec)/kT)[-1]



#        intrho=intrho_maj-self._intrho_ref+self._minorrho*(mu-self.mu_bulk)

        eterm=-q*kT*self.Nc*int_fermi_dirac_one_half((EFn-self.Ec)/kT)
        hterm=-q*kT*self.Nv*int_fermi_dirac_one_half((self.Ev-EFp)/kT)
        if(len(self.dopants)!=1): return
        d=self.dopants.values()[0]
        if d["type"]=="donor":
            dterm=-q*kT*d["conc"]*fermi_dirac_zero(-(EFn-self.Ec+d["E"]+np.log(2)*kT)/kT)
        elif d["type"]=="acceptor":
            dterm=-q*kT*d["conc"]*fermi_dirac_zero(-(self.Ev+d["E"]-EFp+np.log(4)*kT)/kT)
        intrho=eterm+hterm+dterm - self.intrho_ref
        

        return intrho
#        if d["type"]=="acceptor":
#            self._intrho=\
#                q*(self.Nv*int_fermi_dirac_one_half((self.Ev-self.mu_bulk)/kT)
#                    -d["conc"]*fermi_dirac_zero(-(self.Ec-d["Ed"]-self.mu_bulk+log(2))/kT))
    
#if __name__ == "__main__":
#    
#    uno=Material("GaN",{"Mg":8e19})
#    print uno.mu_bulk
#    dos=Material("GaN",{"Si":2e19})
#    print dos.mu_bulk
#    print "phib: ",uno.mu_bulk-dos.mu_bulk
#    
#    print "N"
#    tmat=Material("GaN",{"Si":2e19})
#    #print tmat.Ec-tmat.mu_bulk
#    defi=np.linspace(0,1.5,1000000)
#    rhos=tmat.rho(tmat.mu_bulk-defi,tmat.mu_bulk-defi)
#    #print rhos[-1]/q
#    print (np.cumsum(rhos)*(defi[1]-defi[0]))[-1]/q
#    print -tmat.intrho(tmat.mu_bulk-defi,tmat.mu_bulk-defi)[-1]/q
#    
#    print "P"
#    tmat=Material("GaN",{"Mg":8e19})
#    print tmat.Ev-tmat.mu_bulk
#    defi=-np.linspace(0,1.5,1000000)
#    rhos=tmat.rho(tmat.mu_bulk-defi,tmat.mu_bulk-defi)
#    print rhos[-1]/q
#    print (np.cumsum(rhos)*(defi[1]-defi[0]))[-1]/q
#    print -tmat.intrho(tmat.mu_bulk-defi,tmat.mu_bulk-defi)[-1]/q    
#    
#    print "SPACING"
#    nds=tmat.Nd_ionized(tmat.mu_bulk-defi)
#    print (np.cumsum(nds)*(defi[1]-defi[0]))[-1]
#    print 2e19*kT*fermi_dirac_zero(-(tmat.mu_bulk-defi[-1]-tmat.Ec+.015+np.log(2)*kT)/kT)-2e19*kT*fermi_dirac_zero(-(tmat.mu_bulk-tmat.Ec+.015+np.log(2)*kT)/kT)
#    
#    print "MORE SPACING"
#    ns=tmat.n(tmat.mu_bulk-defi)
#    ns2=tmat.Nc*fermi_dirac_one_half((tmat.mu_bulk-defi-tmat.Ec)/kT)
#    print (np.cumsum(ns)*(defi[1]-defi[0]))[-1]
#    print tmat.Nc*kT*(int_fermi_dirac_one_half((tmat.mu_bulk-defi[-1]-tmat.Ec)/kT)-int_fermi_dirac_one_half((tmat.mu_bulk-tmat.Ec)/kT))
#    