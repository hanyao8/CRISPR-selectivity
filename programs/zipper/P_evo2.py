import numpy as np
import matplotlib.pyplot as plt

t_max=10000

K_b=10/float(t_max)
#K_b=0.0001
q_nc=1
q_c=7.5

#comp_list=[0,0,0,0,0,0,0,0,0,0]
comp_list=[1,1,1,1,1,1,1,1,1,1]
ts_len=len(comp_list)

Q={0:q_nc,1:q_c}


P_tx=[]
for i in range(0,t_max):
    P_tx.append([])
P_tx[0]=[1,0,0,0,0,0,0,0,0,0]

for t in range(0,t_max-1):
    P_tx[t+1].append(P_tx[t][0]+K_b*(Q[comp_list[0]]*(-P_tx[t][0])+P_tx[t][1]-P_tx[t][0]))
    for i in range(1,ts_len-1):    
        P_tx[t+1].append(P_tx[t][i]+K_b*(Q[comp_list[i]]*(P_tx[t][i-1]-P_tx[t][i])+P_tx[t][i+1]-P_tx[t][i]))
    P_tx[t+1].append(P_tx[t][ts_len-1]+K_b*(Q[comp_list[ts_len-1]]*(P_tx[t][ts_len-2]-P_tx[t][ts_len-1])-P_tx[t][ts_len-1]))

    local_sum=sum(P_tx[t+1])
    for j in range(0,ts_len):
        P_tx[t+1][j]=P_tx[t+1][j]/local_sum

P_tx=np.array(P_tx)





P_xt=np.transpose(P_tx)

t_plot=np.linspace(0,1,t_max)
#for k in [9]:
#    plt.plot(t_plot,P_xt[k],label="K_b=%f"%(K_b))
    
for k in range(0,ts_len):
    plt.plot(t_plot,P_xt[k],label="%d"%(k+1))
    
plt.grid(True)
plt.legend(loc=1)
plt.ylabel("Normalized Population")
plt.xlabel("Time")
plt.title("Numerical Simulation Smoothness")
plt.show()

print(P_tx[-1])