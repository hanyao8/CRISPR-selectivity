# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 23:01:37 2018

@author: Choon
"""

import numpy as np
import scipy as sp

avgs=np.array([1.0,1.0])
devs=np.array([1.0,1.0])

def Gaussian(x,mu,sigma):
    return((np.exp(-((x-mu)**2)/2/(sigma**2)))/sigma/np.sqrt(2*np.pi) )
    
def f(x,z,mu,sigma):
    return(abs(x)*Gaussian(x,mu[0],sigma[0])*Gaussian(x/z,mu[0],sigma[0]) )
    
def rho_z(z,mu,sigma):
    lo_lim=min(mu[0]-5*sigma[0],z*(mu[1]-5*sigma[1]))
    hi_lim=max(mu[0]+5*sigma[0],z*(mu[1]+5*sigma[1]))
    ans=1/(z**2)
    ans*=sp.integrate.quad(lambda x: f(x,z,mu,sigma),lo_lim,hi_lim)
    return(ans)
    
print(rho_z(1.0,avgs,devs))