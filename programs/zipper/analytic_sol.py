#linalg
#explore analytic solution for the case of l=5

from numpy import linalg as la
import numpy as np
import matplotlib.pyplot as plt


l=10
K_b=0.001
K_f=7.5*K_b

P_0=np.zeros(l)
P_0[0]=1


t_list=np.linspace(0,100000,100)

#initialisation
#A=np.array([[1,2],[3,4]])
A=np.array([np.zeros(l) for i in range(0,l)])
A[0][0]=-(K_f+K_b)
A[0][1]=K_b
for j in range(1,l-1):
    A[j][j-1]=K_f
    A[j][j]=-(K_f+K_b)
    A[j][j+1]=K_b
A[l-1][l-2]=K_f
A[l-1][l-1]=-(K_f+K_b)

#A2=np.array([[-17,2],[15,-17]])

eigen=la.eig(A)
e_vals=eigen[0]
S=eigen[1]
S_inv=la.inv(S)

P_tx=[]
for t in t_list:
    e_tL=np.array([np.zeros(l) for i in range(0,l)])
    for i in range(0,l):
        e_tL[i][i]=np.exp(t*e_vals[i])
        
    e_tA=np.matmul(S,np.matmul(e_tL,S_inv))
    
    P_tx.append(np.matmul(e_tA,P_0))


P_xt=np.transpose(P_tx)

for k in range(0,l):
    plt.plot(t_list,P_xt[k],label="%d"%(k+1))
    
plt.grid(True)
plt.legend(loc=1)
plt.ylabel("Unnormalized Population")
plt.xlabel("Time")
plt.title("Analytic Solution")
plt.ylim(0,1)
plt.show()