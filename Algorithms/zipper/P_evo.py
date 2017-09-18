import numpy as np
import matplotlib.pyplot as plt

K_b=0.005
q_nc=1
q_c=7.5

comp_list=[1,1,1,0,0,1,1,0,1,1]
ts_len=len(comp_list)

Q={0:q_nc,1:q_c}

t_max=10
P_tx=np.array([[0]*ts_len]*t_max)
P_tx[0][0]=1

t=0
P_tx[1][0]=P_tx[0][0]+K_b*(7.5*(-1*P_tx[0][0])+P_tx[0][1]-P_tx[0][0])
"""
for i in range(1,ts_len-1):    
    P_tx[t+1][i]=P_tx[t][i]+K_b*(Q[comp_list[i]]*(P_tx[t][i-1]-P_tx[t][i])+P_tx[t][i+1]-P_tx[t][i])
P_tx[t+1][ts_len-1]=P_tx[t][ts_len-1]+K_b*(Q[comp_list[ts_len-1]]*(P_tx[t][ts_len-2]-P_tx[t][ts_len-1])-P_tx[t][ts_len-1])
""" 
print(P_tx)
print(P_tx[t][0]+K_b*(Q[comp_list[0]]*(-P_tx[t][0])+P_tx[t][1]-P_tx[t][0]))

#python P_evo.py

"""
for t in range(0,t_max-1):
    P_tx[t+1][0]=P_tx[t][0]+K_b*(Q[comp_list[0]]*(-P_tx[t][0])+P_tx[t][1]-P_tx[t][0])
    for i in range(1,ts_len-1):    
        P_tx[t+1][i]=P_tx[t][i]+K_b*(Q[comp_list[i]]*(P_tx[t][i-1]-P_tx[t][i])+P_tx[t][i+1]-P_tx[t][i])
    P_tx[t+1][ts_len-1]=P_tx[t][ts_len-1]+K_b*(Q[comp_list[ts_len-1]]*(P_tx[t][ts_len-2]-P_tx[t][ts_len-1])-P_tx[t][ts_len-1])
    
    
P_tx=np.array(P_tx)
"""

"""
P_xt=np.transpose(P_tx)

t_plot=np.arange(0,t_max)
for k in range(0,ts_len):
    plt.plot(t_plot,P_xt[k])
    
plt.show()
"""