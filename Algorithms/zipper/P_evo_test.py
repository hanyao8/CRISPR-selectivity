# -*- coding: utf-8 -*-
"""
Created on Fri Aug  4 17:23:09 2017

@author: Choon
"""

#P_evo test

#import numpy as np

#print('Hello World!')

K_b=1/200
q_nc=1
q_c=7.5

comp_list=[1,1,1,0,0,1,1,0,1,1]
#comp_list=[1,1,1,1,0,1,0,1,1,1,0,0,1,1,0,1,0,1,1,1]

Q={0:q_nc,1:q_c}

P_tx=[[1,0,0,0,0,0,0,0,0,0]]
P=P_tx[0]


P_temp=P

P_temp[0]=P_temp[0]+K_b*(Q[comp_list[0]]*(-P[0])+P[1]-P[0])

"""
for i in range(1,len(P)-1):    
    P_temp[i]+=K_b*(Q[comp_list[i]]*(P[i-1]-P[i])+P[i+1]-P[i]) #f(P)
P_temp[len(P)-1]+=K_b*(Q[comp_list[len(P)-1]]*(P[len(P)-2]-P[len(P)-1])-P[len(P)-1])
    
for j in range(0,len(P)):
    P[j]=P_temp[j]/sum(P_temp)

P_tx.append(P)
""" 

#P_tx=np.array(P_tx)
#WP_xt=np.transpose(P_tx)