# -*- coding: utf-8 -*-
"""
Author: Sam Bader
Module: 
"""
import numpy as np
import matplotlib.pyplot as mpl
mpl.close('all')

def import_iv(filename):
    with open(filename) as f:
        vals=np.array([map(float,l.split(',')[1:3]) for i,l in enumerate(f) if i>0])
        return vals[:,1],vals[:,0]
            

hdhs_v,hdhs_i=import_iv('/home/sam/My_All/Cornell/Jena/Delta Tunneling/Data/1003b_110um_3.csv')
mpl.plot(hdhs_v,abs(hdhs_i))

hdns_v,hdns_i=import_iv('/home/sam/My_All/Cornell/Jena/Delta Tunneling/Data/1003a_110_1.csv')
mpl.plot(hdns_v,abs(hdns_i))

ldhs_v,ldhs_i=import_iv('/home/sam/My_All/Cornell/Jena/Delta Tunneling/Data/1003c_110_2.csv')
mpl.plot(ldhs_v,abs(ldhs_i))

yscale('log')

figure()

mpl.plot(hdhs_v,1/gradient(log10(abs(hdhs_i)),mean(diff(hdhs_v)))/.06)

mpl.plot(hdns_v,1/gradient(log10(abs(hdns_i)),mean(diff(hdns_v)))/.06)

mpl.plot(ldhs_v,1/gradient(log10(abs(ldhs_i)),mean(diff(ldhs_v)))/.06)

ylim(0,50)