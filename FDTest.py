# -*- coding: utf-8 -*-
"""
Author: Sam Bader
Module: 
"""
from __future__ import division
import numpy as np

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


if 1:#__name__ == '__main__':

    from pylab import *
    n = 350
    x = np.linspace(-6,10,n)
    f = lambda s:fermihalf(s,1)

    y = map(f,x)

    dex = np.where(x>0)

    low = np.exp(x)
    high = np.tile(None,n)
    high[dex] = 4*x[dex]**1.5/np.sqrt(9*np.pi)

    p1, = semilogy(x,y,'r-',lw=2)
    p2, = semilogy(x,low,'g--',lw=2)
    p3, = semilogy(x,high,'b--',lw=2)
    legend([p1,p2,p3],[r"$F_{1/2}(\eta)$",r"$e^{\eta}$",
            r"$\frac{4\eta^{3/2}}{3\sqrt{\pi}} $"],loc=2)
    title('Approximate Fermi-Dirac 1/2 integral')
    xlabel(r'$\eta$',fontsize=16)
    grid()
    show()
