# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 13:52:59 2017

@author: Choon
"""

import matplotlib.pyplot as plt
import numpy as np


k_B=1.38e-23
N_avo=6.022e23
"""
G_dpx=[]
for i in range(0,4):
    G_dpx.append([])
    for j in range(0,4):
        G_dpx[i].append([])
        for k in range(0,4):
            G_dpx[i][j].append([0,0,0,0])
"""

G_dpx=np.zeros((4,4,4,4))
#Energies displayed in kcal/mol, converted to joules
G_dpx[0][0][1][1]=G_dpx[1][1][0][0]=-1.00*4184/N_avo #AA/TT      
G_dpx[0][1][1][0]=-0.88*4184/N_avo #AT/TA
G_dpx[1][0][0][1]=-0.58*4184/N_avo #TA/AT
G_dpx[2][0][3][1]=G_dpx[1][3][0][2]=-1.45*4184/N_avo #CA/GT
G_dpx[3][1][2][0]=G_dpx[0][2][1][3]=-1.44*4184/N_avo #GT/CA
G_dpx[2][1][3][0]=G_dpx[0][3][1][2]=-1.28*4184/N_avo #CT/GA
G_dpx[3][0][2][1]=G_dpx[1][2][0][3]=-1.30*4184/N_avo #GA/CT
G_dpx[2][3][3][2]=-2.17*4184/N_avo #CG/GC
G_dpx[3][2][2][3]=-2.24*4184/N_avo #GC/CG
G_dpx[3][3][2][2]=G_dpx[2][2][3][3]=-1.84*4184/N_avo #GG/CC

#single defects
G_dpx[0][0][1][0]=G_dpx[0][1][0][0]=0.61*4184/N_avo
G_dpx[0][0][1][2]=G_dpx[2][1][0][0]=0.88*4184/N_avo
G_dpx[0][0][1][3]=G_dpx[3][1][0][0]=0.14*4184/N_avo
G_dpx[0][1][1][1]=G_dpx[1][1][1][0]=0.69*4184/N_avo
G_dpx[0][1][1][2]=G_dpx[2][1][1][0]=0.73*4184/N_avo
G_dpx[0][1][1][3]=G_dpx[3][1][1][0]=0.07*4184/N_avo
G_dpx[0][2][1][0]=G_dpx[0][1][2][0]=0.77*4184/N_avo
G_dpx[0][2][1][1]=G_dpx[1][1][2][0]=0.64*4184/N_avo
G_dpx[0][2][1][2]=G_dpx[2][1][2][0]=1.33*4184/N_avo
G_dpx[0][3][1][0]=G_dpx[0][1][3][0]=0.02*4184/N_avo
G_dpx[0][3][1][1]=G_dpx[1][1][3][0]=0.71*4184/N_avo
G_dpx[0][3][1][3]=G_dpx[3][1][3][0]=-0.13*4184/N_avo

G_dpx[1][0][0][0]=G_dpx[0][0][0][1]=0.69*4184/N_avo
G_dpx[1][0][0][2]=G_dpx[2][0][0][1]=0.92*4184/N_avo
G_dpx[1][0][0][3]=G_dpx[3][0][0][1]=0.42*4184/N_avo
G_dpx[1][1][0][1]=G_dpx[1][0][1][1]=0.68*4184/N_avo
G_dpx[1][1][0][2]=G_dpx[2][0][1][1]=0.75*4184/N_avo
G_dpx[1][1][0][3]=G_dpx[3][0][1][1]=0.34*4184/N_avo
G_dpx[1][2][0][0]=G_dpx[0][0][2][1]=1.33*4184/N_avo
G_dpx[1][2][0][1]=G_dpx[1][0][2][1]=0.97*4184/N_avo
G_dpx[1][2][0][2]=G_dpx[2][0][2][1]=1.05*4184/N_avo
G_dpx[1][3][0][0]=G_dpx[0][0][3][1]=0.74*4184/N_avo
G_dpx[1][3][0][1]=G_dpx[1][0][3][1]=0.43*4184/N_avo
G_dpx[1][3][0][3]=G_dpx[3][0][3][1]=0.44*4184/N_avo

G_dpx[2][0][3][0]=G_dpx[0][3][0][2]=0.43*4184/N_avo
G_dpx[2][0][3][2]=G_dpx[2][3][0][2]=0.75*4184/N_avo
G_dpx[2][0][3][3]=G_dpx[3][3][0][2]=0.03*4184/N_avo
G_dpx[2][1][3][1]=G_dpx[1][3][1][2]=-0.12*4184/N_avo
G_dpx[2][1][3][2]=G_dpx[2][3][1][2]=0.40*4184/N_avo
G_dpx[2][1][3][3]=G_dpx[3][3][1][2]=-0.32*4184/N_avo
G_dpx[2][2][3][0]=G_dpx[0][3][2][2]=0.79*4184/N_avo
G_dpx[2][2][3][1]=G_dpx[1][3][2][2]=0.62*4184/N_avo
G_dpx[2][2][3][2]=G_dpx[2][3][2][2]=0.70*4184/N_avo
G_dpx[2][3][3][0]=G_dpx[0][3][3][2]=0.11*4184/N_avo
G_dpx[2][3][3][1]=G_dpx[1][3][3][2]=-0.47*4184/N_avo
G_dpx[2][3][3][3]=G_dpx[3][3][3][2]=-0.11*4184/N_avo

G_dpx[3][0][2][0]=G_dpx[0][2][0][3]=0.17*4184/N_avo
G_dpx[3][0][2][2]=G_dpx[2][2][0][3]=0.81*4184/N_avo
G_dpx[3][0][2][3]=G_dpx[3][2][0][3]=-0.25*4184/N_avo
G_dpx[3][1][2][1]=G_dpx[1][2][1][3]=0.45*4184/N_avo
G_dpx[3][1][2][2]=G_dpx[2][2][1][3]=0.98*4184/N_avo
G_dpx[3][1][2][3]=G_dpx[3][2][1][3]=-0.59*4184/N_avo
G_dpx[3][2][2][0]=G_dpx[0][2][2][3]=0.47*4184/N_avo
G_dpx[3][2][2][1]=G_dpx[1][2][2][3]=0.62*4184/N_avo
G_dpx[3][2][2][2]=G_dpx[2][2][2][3]=0.79*4184/N_avo
G_dpx[3][3][2][0]=G_dpx[0][2][3][3]=-0.52*4184/N_avo
G_dpx[3][3][2][1]=G_dpx[1][2][3][3]=0.08*4184/N_avo
G_dpx[3][3][2][3]=G_dpx[3][2][3][3]=-1.11*4184/N_avo

plot_array=(np.ndarray.flatten(G_dpx))
plt.hist(plot_array,bins=100)
plt.show()