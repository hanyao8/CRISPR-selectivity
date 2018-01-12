import numpy as np

A=np.array([[1,0,0,0],[0.25,0.25,0.25,0.25],[0.25,0.25,0.25,0.25],[0.25,0.25,0.25,0.25]])
A_tp=np.transpose(A)
b=np.array([[0],[0],[0],[0]])

#x=np.linalg.solve(A,b)

def cm_2_p(cm):
    cm_tp=np.transpose(cm)
    local_soln=np.linalg.eig(cm_tp)
    eigenvals=local_soln[0]
    n_eigens=len(eigenvals)
    for i in range(0,n_eigens):
        if eigenvals[i]==1:
            return(local_soln[1][:,i])
    

