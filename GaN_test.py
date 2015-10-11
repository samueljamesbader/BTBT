# -*- coding: utf-8 -*-
"""
Author: Sam Bader
Module: 
"""
from __future__ import division
import matplotlib.pyplot as mpl
mpl.close('all')

from Constants import eps0,q,hbar,m0,kT,eVToJoules
from Material import Material
from BandDiagrams import import_bd, generate_exactspikedpn_bandstructure,generate_exactspikedpin_bandstructure, generate_pn_bandstructure, generate_flatpin_bandstructure
from TunnelCalcs import current, max_field

# vs 1dp
if 1:
    print "Mine ", Material("GaN",{"Si":2e17}).Ec-Material("GaN",{"Si":2e17}).EBulk
    
    #myev=Material("GaN",{"Mg":2e20}).Ev-Material("GaN",{"Mg":2e20}).EBulk
    #myec=Material("GaN",{"Si":2e20}).Ec-Material("GaN",{"Si":2e20}).EBulk
    
    g=import_bd("/home/sam/1DPoisson/1D Poisson Beta 8g Linux Distribution/Input_Files_Examples/GaNDelta_Out.txt")#,750,1250)
    g["x"]=g["x"]-.1e-4
    print "1DP ", g["Ec"]
    
    
    # THESE ARE GOOD
    #mpl.plot(g["x"]*1e4,g["Ev"],'b--')
    #mpl.plot(g["x"]*1e4,g["Ec"],'g--')
    #mpl.plot(g["x"]*1e4,np.gradient(np.gradient(g["Ec"]))/np.gradient(g["x"])**2,'g--')
    
    
    
    #plot(g["x"]*1e4,g["x"]*0+myev,'r')
    #mpl.plot(g["x"]*1e4,g["x"]*0+myec,'r')
    
    #ylim([40,50])
    #xlim([-.05,.05])
    #figure()
    #dx=mean(diff(g["x"]))
    #plot(g["x"][:-1]*1e4,diff(g["Ev"])/dx/1e6)
    #ylabel("MV/cm")





# Generate a depletion approximation band diagram for a GaN PIN junction
#gan_pin_u=generate_exactspikedpn_bandstructure("GaN",.1e-4,"Mg",1e19,0e14,.2e-4,"Si",1e18,[-.25,-.5,-1,-5,-10])
#gan_pin_u=generate_exactspikedpn_bandstructure("GaN",.1e-4,"Mg",5e19,.5e14,.2e-4,"Si",2e19,[-5])

gan_pin=generate_exactspikedpin_bandstructure("GaN",.1e-4,"Mg",2e19,.005e-4,1.3e14,.2e-4,"Si",2e19,[-5])
gan_pn=generate_exactspikedpn_bandstructure  ("GaN",.1e-4,"Mg",2e19,        1.3e14,.2e-4,"Si",2e19,[-5])

#gan_pin_upn=generate_pn_bandstructure("GaN",.1e-4,"Mg",1e19,.2e-4,"Si",1e18,[-0])

#gan_pin=generate_exactspikedpn_bandstructure("GaN",.4e-4,"Mg",2e20,1e14,.4e-4,"Si",2e20,[-1])
#gan_pin_hs=generate_exactspikedpn_bandstructure("GaN",.4e-4,"Mg",2e20,1e14,.4e-4,"Si",2e20,[-.3,-.5,-1])
#gan_pin_ls=generate_exactspikedpn_bandstructure("GaN",.4e-4,"Mg",2e19,1e14,.4e-4,"Si",2e17,[-.3,-.5,-1])
#gan_pin_hu=generate_exactspikedpn_bandstructure("GaN",.4e-4,"Mg",2e20,0e14,.4e-4,"Si",2e20,[-.3,-.5,-1])
#gan_pin_lu=generate_exactspikedpn_bandstructure("GaN",.4e-4,"Mg",2e19,0e14,.4e-4,"Si",2e17,[-.3,-.5,-1])

#gan_pin=generate_flatpin_bandstructure("GaN",.01e-4,"Mg",2e19,.0035e-4,.01e-4,"Si",2e17,[-1])

#gan_pin=generate_pn_bandstructure("GaN",.4e-4,"Mg",2e19,.4e-4,"Si",2e17,[-80])
jv2={}
jd={}
jm={}
#for name,gan_pin in [("hs",gan_pin_hs),("ls",gan_pin_ls),("hu",gan_pin_hu),("lu",gan_pin_lu)]:
for name,g in [("pn",gan_pn),("pin",gan_pin)]:
    if 1:
        #figure()
        gp=g
        g=g[0]
        print name
        print g
        #figure(name)
        #for g in gs:
        plot(g["x"]*1e4,g["Ev"],label=name)
        plot(g["x"]*1e4,g["Ec"],label=name)
        plot(g["x"]*1e4,g["EB"],'r.',markersize=1,label=name)
        #figure()
        #for g in gan_pin:
        #    plot(g["x"]*1e7,abs(-np.gradient(np.gradient(g["Ec"]))/np.gradient(g["x"])**2*g["material"].eps/q),label=("Vr="+str(g["Vr"])+"V"))
        #g=gan_pin_upn[0]
        #plot(g["x"]*1e4,np.gradient(np.gradient(g["Ec"]))/np.gradient(g["x"])**2,'b')
        print "Max field: ", max_field(g)/1e6
        #legend(loc='upper right')
        #xlabel('$x$ [nm]')
        #ylabel('net charge density [e/cm$^{3}$]')
        #yscale('log')
        g=gp
    if 0:
        figure()
        #jo=array([current(b,"orig") for b in gan_pin])
        jv2[name]=array([current(b,"v2") for b in g])
        #jd[name]=array([current(b,"d") for b in gan_pin])
        jm[name]=array([current(b,"maxfield") for b in g])
        #print array([jo,jv2,jm])/1e4
        #print array([jv2,jm])/1e4
        #print array([jv2[name],jd[name],jm[name]])/1e4
        print array([jv2[name],jm[name]])/1e4
        
        
        #plot([max_field(b)/1e6 for b in gan_pin],j/1e4,'b',label='v2')

