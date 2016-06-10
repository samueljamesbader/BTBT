# -*- coding: utf-8 -*-
"""
Author: Sam Bader
Module: 
"""
from __future__ import division
import numpy as np
import warnings

# https://en.wikipedia.org/wiki/Complete_Fermi%E2%80%93Dirac_integral
def fermi_dirac_zero(x):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore",RuntimeWarning)
        return np.log(1+np.exp(x))

# Wong et al Solid-State Electronics Vol. 37, No. I, pp. 61~64, 1994
a=np.array([1,.353568,.192439,.122973,.077134,.036228,.008346])
b1=np.array([.76514793,.60488667,.19003355,2.00193968e-2,-4.12643816e-3,-4.70958992e-4,1.50071469e-4])
b2=np.array([.78095732,.57254453,.21419339,1.38432741e-2,-5.54949386e-3,6.48814900e-4,-2.84050520e-5])
c=np.array([.752253,.928195,.680839,25.7829,-553.636,3531.43,-3254.65])

def fermi_dirac_one_half(x):
#    print "ABOUT TO FD ",x
    with warnings.catch_warnings():
        warnings.simplefilter("ignore",RuntimeWarning)
        return np.choose((0+(x>0)+(x>2)+(x>5)), [
            sum(map(lambda i: (-1)**i*a[i]*np.exp((i+1)*x),range(7)),0),
            np.polyval(b1[::-1],x),
            np.polyval(b2[::-1],x),
            sum(map(lambda i: c[i]*abs(x+(x==0))**(1.5-2*i),range(7)),0)
        ])

def bad_fermi_dirac_one_half(x):
    return (4*x**(3/2)/(3*np.sqrt(np.pi))+.75/(1+.1*x))
def bad_int_fermi_dirac_one_half(x):
    return (8*x**(5/2)/(15*np.sqrt(np.pi)))
    #return (np.sqrt(8/(15*np.sqrt(np.pi)))*x**(5/4))**2+int_fermi_dirac_one_half(100)-(np.sqrt(8/(15*np.sqrt(np.pi)))*100**(5/4))**2


# Direct integration of the above expressions
int_0=sum(map(lambda i: (-1)**i/(i+1)*a[i],range(7)))
int_2=int_0+sum(map(lambda i: 1/(i+1)*b1[i]*2**(i+1),range(7)))
int_5=int_2+sum(map(lambda i: 1/(i+1)*b2[i]*(5**(i+1)-2**(i+1)),range(7)))
def int_fermi_dirac_one_half(x):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore",RuntimeWarning)
        return np.choose((0+(x>0)+(x>2)+(x>5)), [
            sum(map(lambda i: (-1)**i/(i+1)*a[i]*np.exp((i+1)*x),range(7)),0),
            int_0+sum(map(lambda i: 1/(i+1)*b1[i]*x**(i+1),range(7)),0),
            int_2+sum(map(lambda i: 1/(i+1)*b2[i]*(x**(i+1)-2**(i+1)),range(7)),0),
            int_5+sum(map(lambda i: 1/(2.5-2*i)*c[i]*(abs(x+(x==0))**(2.5-2*i)-5**(2.5-2*i)),range(7)),0)
        ])


if __name__ == "__main__":
    close('all')
    x=np.linspace(-5,10,10000)
    fd=fermi_dirac_one_half(x)

    fd2=bad_fermi_dirac_one_half(x)        
    title('FD')
    plot(x,fd)
    plot(x,fd2,'.')
    #yscale('log')

    figure()
    title('FD Int')
    fdsum=np.cumsum(fd)*(x[1]-x[0])
    fdint=int_fermi_dirac_one_half(x)-int_fermi_dirac_one_half(-5)
    fdint2=bad_int_fermi_dirac_one_half(x)
    plot(x,fdsum,'.')
    plot(x,fdint)
    plot(x,fdint2,'r.')
    