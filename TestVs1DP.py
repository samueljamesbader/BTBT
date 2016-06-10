# -*- coding: utf-8 -*-
"""
Author: Sam Bader
Module: 
"""


from __future__ import division
import numpy as np
from Constants import eps0,q,hbar,m0,kT,eVToJoules
from Material import Material
import matplotlib.pyplot as mpl
from SpikedPINJunction import SpikedPINJunction

mpl.close('all')

import re
def import_bd(filename,xmin=-1e10,xmax=1e10,center=0):
    with open(filename) as f:
        data=[]
        for l in f:
            mo=re.match("\s*([\d\.e+-]+)\s+([\d\.e+-]+)\s+([\d\.e+-]+)\s+([\d\.e+-]+)\s+"\
                        "([\d\.e+-]+)\s+([\d\.e+-]+)\s+([\d\.e+-]+)\s+([\d\.e+-]+)\s+", l)
            if mo:
                x=float(mo.groups()[0])
                if(x>=xmin and x<=xmax):
                    data+=[[float(a) for a in mo.groups()]]
                    data[-1][0]=(x-center)/1e8
        data=np.array(data)
        #print data
        bd={}
        bd["x"]=data[:,0]
        bd["Ec"]=data[:,1]
        bd["Ev"]=data[:,2]
        bd["F"]=data[:,3]
        bd["Ef"]=data[:,4]
        bd["n"]=data[:,5]
        bd["p"]=data[:,6]
        bd["material"]=Material("GaN")
    return bd

#gP=import_bd("/home/sam/1DPoisson/1D Poisson Beta 8g Linux Distribution/"\
#        "Input_Files_Examples/gd_withoutSch.txt",0,center=1000)
#
#mpl.plot(gP['x']*1e7,gP['Ev'],'g',markersize=3,label='$E_v p$')
#mpl.plot(gP['x']*1e7,gP['Ec'],'b',markersize=3,label='$E_c p$')
#mpl.plot(gP['x']*1e7,gP['Ef'],'r',markersize=1,label='$\\mu p$')
#
#gP=import_bd("/home/sam/1DPoisson/1D Poisson Beta 8g Linux Distribution/"\
#        "Input_Files_Examples/gd_withSch2.txt",0,center=1000)
#        
#mpl.plot(gP['x']*1e7,gP['Ev'],'g--',markersize=3,label='$E_v p$')
#mpl.plot(gP['x']*1e7,gP['Ec'],'b--',markersize=3,label='$E_c p$')
#mpl.plot(gP['x']*1e7,gP['Ef'],'r--',markersize=1,label='$\\mu p$')
#        
#stop

gM=SpikedPINJunction("GaN","Mg",8e16,.1e-4,0*.005e-4,1.2e13,"Si",2e16,100e-7,0)
#gM=SpikedPINJunction("GaAs","C",8e19,.1e-4,0*.005e-4,0e13,"S",1e19,100e-7,0,sigma_name="S")
#gM.plot_electric_field()
gM.plot_band_diagram()


gP=import_bd("/home/sam/1DPoisson/1D Poisson Beta 8g Linux Distribution/"\
        "Input_Files_Examples/GaNDelta_Out.txt",center=1000)

#fig,ax=mpl.subplots(3,sharex=True,figsize=(6,9))
#mpl.axes(ax[0])
#mpl.plot(gP['x']*1e7,gP['Ev'],'g',markersize=3,label='$E_v p$')
#mpl.plot(gP['x']*1e7,gP['Ec'],'b',markersize=3,label='$E_c p$')
#mpl.plot(gP['x']*1e7,gP['Ef'],'r',markersize=1,label='$\\mu p$')
#mpl.axes(ax[1])
#mpl.plot(gP['x']*1e7,gP['F']/1e6,'g',markersize=3,label='$F$')
#mpl.axes(ax[2])
#mpl.plot(gP['x']*1e7,gP['n'],'g',markersize=3,label='$n$')
#mpl.plot(gP['x']*1e7,gP['p'],'g',markersize=3,label='$p$')
#mpl.yscale('log')
#mpl.ylim(1e16,1e21)
#mpl.xlim(70,130)

#gP=import_bd("/home/sam/1DPoisson/1D Poisson Beta 8g Linux Distribution/"\
#        "Input_Files_Examples/GIG_Out_Class.txt",0,center=000)
#
#mpl.axes(ax[0])
#mpl.plot(gP['x']*1e7,gP['Ev'],'g--',markersize=3,label='$E_v p$')
#mpl.plot(gP['x']*1e7,gP['Ec'],'b--',markersize=3,label='$E_c p$')
#mpl.plot(gP['x']*1e7,gP['Ef'],'r--',markersize=1,label='$\\mu p$')
#mpl.axes(ax[1])
#mpl.plot(gP['x']*1e7,gP['F']/1e6,'g--',markersize=3,label='$F$')
#mpl.axes(ax[2])
#mpl.plot(gP['x']*1e7,gP['n'],'g--',markersize=3,label='$n$')
#mpl.plot(gP['x']*1e7,gP['p'],'g--',markersize=3,label='$p$')
#mpl.yscale('log')
#mpl.ylim(1e16,1e21)
#mpl.xlim(70,130)


#
mpl.plot(gP['x']*1e7,gP['Ev'],'g--',markersize=3,label='$E_v p$')
mpl.plot(gP['x']*1e7,gP['Ec'],'b--',markersize=3,label='$E_c p$')
mpl.plot(gP['x']*1e7,gP['Ef'],'r--',markersize=1,label='$\\mu p$')
