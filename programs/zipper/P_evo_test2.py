import numpy as np
import matplotlib.pyplot as plt

K_b=0.005
q_nc=1
q_c=7.5

comp_list=[1,1,1,0,0,1,1,0,1,1]
#comp_list=[1,1,1,1,0,1,0,1,1,1,0,0,1,1,0,1,0,1,1,1]
ts_len=len(comp_list)

Q={0:q_nc,1:q_c}

#P=[1,0,0,0,0,0,0,0,0,0]
#P_tx=[P]

t_max=10
P_tx=[[0]*ts_len]*t_max
P_tx[0][0]=1

#P_temp=P
    
#P_temp[0]+=K_b*(Q[comp_list[0]]*(-P[0])+P[1]-P[0])

t=1

P_tx[t][0]=P_tx[t-1][0]+K_b*(Q[comp_list[0]]*(-P_tx[t-1][0])+P_tx[t-1][1]-P_tx[t-1][0])
for i in range(1,ts_len-1):    
    P_tx[t][i]=P_tx[t-1][i]+K_b*(Q[comp_list[i]]*(P_tx[t-1][i-1]-P_tx[t-1][i])+P_tx[t-1][i+1]-P_tx[t-1][i])
    
P_tx[t][ts_len-1]=P_tx[t-1][ts_len-1]+K_b*(Q[comp_list[ts_len-1]]*(P_tx[t-1][ts_len-2]-P_tx[t-1][ts_len-1])-P_tx[t-1][ts_len-1])


"""
for t in range(1,t_max):
    #P=P_tx[t-1]
    P_temp[0]=P_tx[t-1][0]+K_b*(Q[comp_list[0]]*(-P_tx[t-1][0])+P_tx[t-1][1]-P_tx[t-1][0])
    for i in range(1,ts_len-1):    
        P_temp[i]=P_tx[t-1][i]+K_b*(Q[comp_list[i]]*(P_tx[t-1][i-1]-P_tx[t-1][i])+P_tx[t-1][i+1]-P_tx[t-1][i]) #f(P)
    P_temp[ts_len-1]=P_tx[t-1][ts_len-1]+K_b*(Q[comp_list[ts_len-1]]*(P_tx[t-1][ts_len-2]-P_tx[t-1][ts_len-1])-P_tx[t-1][ts_len-1])
    
    P_tx[t]=P_temp
    #P_tx.append(P_temp)
    
P_tx=np.array(P_tx)

P_xt=np.transpose(P_tx)

t_plot=np.arange(0,t_max)
for k in range(0,ts_len):
    plt.plot(t_plot,P_xt[k])
    
plt.show()
"""