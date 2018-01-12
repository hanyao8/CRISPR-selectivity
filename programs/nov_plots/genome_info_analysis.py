# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 23:09:24 2017

@author: Choon
"""

#genome info analysis

import numpy as np
import matplotlib.pyplot as plt

txt_data=np.genfromtxt(r'C:\Users\Choon\Desktop\CRISPR\Code\results\genome_entropy3.txt',delimiter=',')
info=txt_data[:,1]
Z=txt_data[:,2]
S_real=txt_data[:,3]
Boltz_prob=txt_data[:,4]

N_datapts=len(info)

offbind_prob=np.ones(N_datapts)-Boltz_prob

f1=plt.figure()
f2=plt.figure()
f3=plt.figure()
#f4=plt.figure()

ax1=f1.add_subplot(121)
ax12=f1.add_subplot(122)
ax2=f2.add_subplot(121)
ax22=f2.add_subplot(122)
ax3=f3.add_subplot(121)
ax32=f3.add_subplot(122)
#ax4=f4.add_subplot(111)

ax1.scatter(info,Z,s=1)
ax12.scatter(info,np.log(Z)/np.log(10),s=1)
f1.suptitle("Genome information vs Z")
ax1.set_ylabel("Z")
ax1.set_xlabel("info_gnm")
ax12.set_ylabel("log_10(Z)")
ax12.set_xlabel("info_gnm")

ax2.scatter(info,S_real,s=1)
ax22.scatter(info,np.log(S_real)/np.log(10),s=1)
f2.suptitle("Genome information vs S_real")
ax2.set_ylabel("S_real")
ax2.set_xlabel("info_gnm")
ax22.set_ylabel("log_10(S_real)")
ax22.set_xlabel("info_gnm")

ax3.scatter(info,offbind_prob,s=1)
ax32.scatter(info,np.log(offbind_prob)/np.log(10),s=1)
f3.suptitle("Genome information vs 1-Boltz_prob")
ax3.set_ylabel("1-Boltz_prob")
ax3.set_xlabel("info_gnm")
ax32.set_ylabel("log_10(1-Boltz_prob)")
ax32.set_xlabel("info_gnm")

plt.show()