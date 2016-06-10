# -*- coding: utf-8 -*-
"""
Author: Sam Bader
Module: 
"""


close('all')


for name,g in [("unspiked",gan_upin)]:#,("pin",gan_pin)]:
    if 1:
        #figure()
        gp=g
        g=g[0]
        print name
        print g
        #figure(name)
        #for g in gs:
        mpl.plot(g["x"]*1e7,g["Ev"],'g--',markersize=3,label='$E_v \\mathrm{ pn}$')#label=name+" Ev")
        mpl.plot(g["x"]*1e7,g["Ec"],'b--',markersize=3,label='$E_c \\mathrm{ pn}$')#,label=name+" Ec")
        #mpl.plot(g["x"]*1e7,g["EB"],'r.',markersize=1,label='$\\mu$')#label=name+" EB")
        mpl.xlabel("$z$ [nm]")
        mpl.ylabel("$E$ [eV]")


for name,g in [("spiked",gan_pin)]:#,("pin",gan_pin)]:
    if 1:
        #figure()
        gp=g
        g=g[0]
        print name
        print g
        #figure(name)
        #for g in gs:
        mpl.plot(g["x"]*1e7,g["Ev"],'g',markersize=3,label='$E_v \\mathrm{ spiked}$')#label=name+" Ev")
        mpl.plot(g["x"]*1e7,g["Ec"],'b',markersize=3,label='$E_c \\mathrm{ spiked}$')#,label=name+" Ec")
        mpl.plot(g["x"]*1e7,g["EB"],'r.',markersize=1)#,label='$\\mu$')#label=name+" EB")
        
        mpl.legend(loc='best')
        

xlim(-40,90)
ylim(-4.5,4)