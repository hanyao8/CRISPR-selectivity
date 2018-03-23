import numpy as np
import os
import matplotlib.pyplot as plt

#########################
#SantaLucia Parameters
#########################

#Initialisation of energies involved in the 'Santa Lucia rules'

#The free energy of a duplex appears in the following notation:
#G_dpx[n_W][n_X][n_Y][n_Z] is equivalent to the free energy of the duplex
#5' WX 3'
#3' YZ 5'


G_dpx=np.empty((4,4,4,4))

G_dpx[0][0][1][1]=G_dpx[1][1][0][0]=-1.00 #AA/TT      
G_dpx[0][1][1][0]=-0.88 #AT/TA
G_dpx[1][0][0][1]=-0.58 #TA/AT
G_dpx[2][0][3][1]=G_dpx[1][3][0][2]=-1.45 #CA/GT
G_dpx[3][1][2][0]=G_dpx[0][2][1][3]=-1.44 #GT/CA
G_dpx[2][1][3][0]=G_dpx[0][3][1][2]=-1.28 #CT/GA
G_dpx[3][0][2][1]=G_dpx[1][2][0][3]=-1.30 #GA/CT
G_dpx[2][3][3][2]=-2.17 #CG/GC
G_dpx[3][2][2][3]=-2.24 #GC/CG
G_dpx[3][3][2][2]=G_dpx[2][2][3][3]=-1.84 #GG/CC

ND_uniq=np.array([-1.00,-0.88,-0.58,-1.45,-1.44,-1.28,-1.30,-2.17,-2.24,-1.84])

#single defects
G_dpx[0][0][1][0]=G_dpx[0][1][0][0]=0.61
G_dpx[0][0][1][2]=G_dpx[2][1][0][0]=0.88
G_dpx[0][0][1][3]=G_dpx[3][1][0][0]=0.14
G_dpx[0][1][1][1]=G_dpx[1][1][1][0]=0.69
G_dpx[0][1][1][2]=G_dpx[2][1][1][0]=0.73
G_dpx[0][1][1][3]=G_dpx[3][1][1][0]=0.07
G_dpx[0][2][1][0]=G_dpx[0][1][2][0]=0.77
G_dpx[0][2][1][1]=G_dpx[1][1][2][0]=0.64
G_dpx[0][2][1][2]=G_dpx[2][1][2][0]=1.33
G_dpx[0][3][1][0]=G_dpx[0][1][3][0]=0.02
G_dpx[0][3][1][1]=G_dpx[1][1][3][0]=0.71
G_dpx[0][3][1][3]=G_dpx[3][1][3][0]=-0.13

G_dpx[1][0][0][0]=G_dpx[0][0][0][1]=0.69
G_dpx[1][0][0][2]=G_dpx[2][0][0][1]=0.92
G_dpx[1][0][0][3]=G_dpx[3][0][0][1]=0.42
G_dpx[1][1][0][1]=G_dpx[1][0][1][1]=0.68
G_dpx[1][1][0][2]=G_dpx[2][0][1][1]=0.75
G_dpx[1][1][0][3]=G_dpx[3][0][1][1]=0.34
G_dpx[1][2][0][0]=G_dpx[0][0][2][1]=1.33
G_dpx[1][2][0][1]=G_dpx[1][0][2][1]=0.97
G_dpx[1][2][0][2]=G_dpx[2][0][2][1]=1.05
G_dpx[1][3][0][0]=G_dpx[0][0][3][1]=0.74
G_dpx[1][3][0][1]=G_dpx[1][0][3][1]=0.43
G_dpx[1][3][0][3]=G_dpx[3][0][3][1]=0.44

G_dpx[2][0][3][0]=G_dpx[0][3][0][2]=0.43
G_dpx[2][0][3][2]=G_dpx[2][3][0][2]=0.75
G_dpx[2][0][3][3]=G_dpx[3][3][0][2]=0.03
G_dpx[2][1][3][1]=G_dpx[1][3][1][2]=-0.12
G_dpx[2][1][3][2]=G_dpx[2][3][1][2]=0.40
G_dpx[2][1][3][3]=G_dpx[3][3][1][2]=-0.32
G_dpx[2][2][3][0]=G_dpx[0][3][2][2]=0.79
G_dpx[2][2][3][1]=G_dpx[1][3][2][2]=0.62
G_dpx[2][2][3][2]=G_dpx[2][3][2][2]=0.70
G_dpx[2][3][3][0]=G_dpx[0][3][3][2]=0.11
G_dpx[2][3][3][1]=G_dpx[1][3][3][2]=-0.47
G_dpx[2][3][3][3]=G_dpx[3][3][3][2]=-0.11

G_dpx[3][0][2][0]=G_dpx[0][2][0][3]=0.17
G_dpx[3][0][2][2]=G_dpx[2][2][0][3]=0.81
G_dpx[3][0][2][3]=G_dpx[3][2][0][3]=-0.25
G_dpx[3][1][2][1]=G_dpx[1][2][1][3]=0.45
G_dpx[3][1][2][2]=G_dpx[2][2][1][3]=0.98
G_dpx[3][1][2][3]=G_dpx[3][2][1][3]=-0.59
G_dpx[3][2][2][0]=G_dpx[0][2][2][3]=0.47
G_dpx[3][2][2][1]=G_dpx[1][2][2][3]=0.62
G_dpx[3][2][2][2]=G_dpx[2][2][2][3]=0.79
G_dpx[3][3][2][0]=G_dpx[0][2][3][3]=-0.52
G_dpx[3][3][2][1]=G_dpx[1][2][3][3]=0.08
G_dpx[3][3][2][3]=G_dpx[3][2][3][3]=-1.11

#Double defects
#Using calculated average 'isolated defect energies'
G_avg=np.zeros(8)
G_avg[0]=G_avg_AA=1.20
G_avg[1]=G_avg_AC=1.56
G_avg[2]=G_avg_AG=0.81
G_avg[3]=G_avg_TT=1.15
G_avg[4]=G_avg_TC=1.44
G_avg[5]=G_avg_TG=0.75
G_avg[6]=G_avg_CC=1.69
G_avg[7]=G_avg_GG=0.50

G_dpx[0][0][0][0]=(G_avg_AA+G_avg_AA)
G_dpx[0][0][0][2]=G_dpx[0][0][2][0]=G_dpx[0][2][0][0]=G_dpx[2][0][0][0]=(G_avg_AA+G_avg_AC)
G_dpx[0][0][0][3]=G_dpx[0][0][3][0]=G_dpx[0][3][0][0]=G_dpx[3][0][0][0]=(G_avg_AA+G_avg_AG)
G_dpx[0][1][0][1]=G_dpx[1][0][1][0]=(G_avg_AA+G_avg_TT)
G_dpx[0][1][0][2]=G_dpx[2][0][1][0]=G_dpx[1][0][2][0]=G_dpx[0][2][0][1]=(G_avg_AA+G_avg_TC)
G_dpx[0][1][0][3]=G_dpx[3][0][1][0]=G_dpx[1][0][3][0]=G_dpx[0][3][0][1]=(G_avg_AA+G_avg_TG)
G_dpx[0][2][0][2]=G_dpx[2][0][2][0]=(G_avg_AA+G_avg_CC)
G_dpx[0][3][0][3]=G_dpx[3][0][3][0]=(G_avg_AA+G_avg_GG)

G_dpx[0][0][2][2]=G_dpx[0][2][2][0]=G_dpx[2][0][0][2]=G_dpx[2][2][0][0]=(G_avg_AC+G_avg_AC)
G_dpx[0][0][2][3]=G_dpx[0][0][3][2]=G_dpx[2][3][0][0]=G_dpx[3][2][0][0]=G_dpx[2][0][0][3]=G_dpx[3][0][0][2]=G_dpx[0][3][2][0]=G_dpx[0][2][3][0]=(G_avg_AC+G_avg_AG)
G_dpx[0][1][2][1]=G_dpx[1][2][1][0]=G_dpx[2][1][0][1]=G_dpx[1][0][1][2]=(G_avg_AC+G_avg_TT)
G_dpx[0][1][2][2]=G_dpx[2][2][1][0]=G_dpx[2][1][0][2]=G_dpx[2][0][1][2]=G_dpx[0][2][2][1]=G_dpx[1][2][2][0]=G_dpx[1][0][2][2]=G_dpx[2][2][0][1]=(G_avg_AC+G_avg_TC)
G_dpx[0][1][2][3]=G_dpx[3][2][1][0]=G_dpx[0][3][2][1]=G_dpx[1][2][3][0]=G_dpx[1][0][3][2]=G_dpx[2][3][0][1]=G_dpx[2][1][0][3]=G_dpx[3][0][1][2]=(G_avg_AC+G_avg_TG)
G_dpx[0][2][2][2]=G_dpx[2][2][2][0]=G_dpx[2][0][2][2]=G_dpx[2][2][0][2]=(G_avg_AC+G_avg_CC)
G_dpx[0][3][2][3]=G_dpx[3][2][3][0]=G_dpx[3][0][3][2]=G_dpx[2][3][0][3]=(G_avg_AC+G_avg_GG)

G_dpx[0][0][3][3]=G_dpx[0][3][3][0]=G_dpx[3][0][0][3]=G_dpx[3][3][0][0]=(G_avg_AG+G_avg_AG)
G_dpx[0][1][3][1]=G_dpx[1][3][1][0]=G_dpx[1][0][1][3]=G_dpx[3][1][0][1]=(G_avg_AG+G_avg_TT)
G_dpx[0][1][3][2]=G_dpx[2][3][1][0]=G_dpx[0][2][3][1]=G_dpx[1][3][2][0]=G_dpx[1][0][2][3]=G_dpx[3][2][0][1]=G_dpx[3][1][0][2]=G_dpx[2][0][1][3]=(G_avg_AG+G_avg_TC)
G_dpx[0][1][3][3]=G_dpx[3][3][1][0]=G_dpx[3][1][0][3]=G_dpx[3][0][1][3]=G_dpx[0][3][3][1]=G_dpx[1][3][3][0]=G_dpx[1][0][3][3]=G_dpx[3][3][0][1]=(G_avg_AG+G_avg_TG)
G_dpx[0][2][3][2]=G_dpx[2][3][2][0]=G_dpx[2][0][2][3]=G_dpx[3][2][0][2]=(G_avg_AG+G_avg_CC)
G_dpx[0][3][3][3]=G_dpx[3][3][3][0]=G_dpx[3][0][3][3]=G_dpx[3][3][0][3]=(G_avg_AG+G_avg_GG)

G_dpx[1][1][1][1]=(G_avg_TT+G_avg_TT)
G_dpx[1][1][1][2]=G_dpx[1][1][2][1]=G_dpx[1][2][1][1]=G_dpx[2][1][1][1]=(G_avg_TT+G_avg_TC)
G_dpx[1][1][1][3]=G_dpx[1][1][3][1]=G_dpx[1][3][1][1]=G_dpx[3][1][1][1]=(G_avg_TT+G_avg_TG)
G_dpx[2][1][2][1]=G_dpx[1][2][1][2]=(G_avg_TT+G_avg_CC)
G_dpx[3][1][3][1]=G_dpx[1][3][1][3]=(G_avg_TT+G_avg_GG)

G_dpx[1][1][2][2]=G_dpx[1][2][2][1]=G_dpx[2][1][1][2]=G_dpx[2][2][1][1]=(G_avg_TC+G_avg_TC)
G_dpx[1][1][2][3]=G_dpx[1][1][3][2]=G_dpx[2][3][1][1]=G_dpx[3][2][1][1]=G_dpx[2][1][1][3]=G_dpx[3][1][1][2]=G_dpx[1][3][2][1]=G_dpx[1][2][3][1]=(G_avg_TC+G_avg_TG)
G_dpx[1][2][2][2]=G_dpx[2][2][2][1]=G_dpx[2][1][2][2]=G_dpx[2][2][1][2]=(G_avg_TC+G_avg_CC)
G_dpx[1][3][2][3]=G_dpx[3][2][3][1]=G_dpx[3][1][3][2]=G_dpx[2][3][1][3]=(G_avg_TC+G_avg_GG)

G_dpx[1][1][3][3]=G_dpx[1][3][3][1]=G_dpx[3][1][1][3]=G_dpx[3][3][1][1]=(G_avg_TG+G_avg_TG)
G_dpx[1][2][3][2]=G_dpx[2][3][2][1]=G_dpx[2][1][2][3]=G_dpx[3][2][1][2]=(G_avg_TG+G_avg_CC)
G_dpx[1][3][3][3]=G_dpx[3][3][3][1]=G_dpx[3][1][3][3]=G_dpx[3][3][1][3]=(G_avg_TG+G_avg_GG)

G_dpx[2][2][2][2]=(G_avg_CC+G_avg_CC)
G_dpx[2][3][2][3]=G_dpx[3][2][3][2]=(G_avg_CC+G_avg_GG)

G_dpx[3][3][3][3]=(G_avg_GG+G_avg_GG)

#16x16 matrix
#The matrix created will be such that ROW nt identity 'WX' corresponds to the 'upper' strand 5'-WX-3'
#To illustrate this:

#         "AG"             "CA"             "CT"        ...     ...     ...
#"AG" deltaG(AG/AG)    deltaG(AG/CA)    deltaG(AG/CT)  
#"CA" deltaG(CA/AG)    deltaG(CA/CA)    deltaG(CA/CT) 
#"CT" deltaG(CT/AG)    deltaG(CT/CA)    deltaG(CT/CT)    
#...
#...

#Where deltaG(WX/YZ) is the deltaG of:
#  (5')-(W)-(X)-(3')
#        |   |
#  (3')-(Y)-(Z)-(5')
#I.e. the conventional Santa Lucia definition

rows=["AA","AT","AC","AG","CA","CT","CC","CG","GC","GG","GA","GT","TC","TG","TA","TT"]
nt_2_num={"A":0,"T":1,"C":2,"G":3}
G_dpx_16x16=np.empty((16,16))
for i in range(0,16):
    for j in range(0,16):
        G_dpx_16x16[i][j]=G_dpx[nt_2_num[ rows[i][0] ]] [nt_2_num[ rows[i][1] ]][nt_2_num[ rows[j][0] ]][nt_2_num[ rows[j][1] ]]
        
print(G_dpx)
print(G_dpx_16x16)

energies=G_dpx_16x16.flatten()

pw_comp_dict={"A":"T","T":"A","C":"G","G":"C"}
comp_dict={"TA":"AT","AT":"TA","AA":"TT","TT":"AA","CT":"GA","AG":"TC","GA":"CT","TC":"AG","GT":"CA","AC":"TG","CA":"GT","TG":"AC","GG":"CC","CC":"GG","CG":"GC","GC":"CG"}
def score(a,b):
    num=0
    for i in range(0,2):
        if not(b[i] == pw_comp_dict[a[i]]):
            num+=1
    return(num)
    
ND_array=np.array([])
SD_array=np.array([])
DD_array=np.array([])

ND_uniq=np.array([])
#SD_uniq=np.array([])
DD_uniq=np.array([])
for i in range(0,16):
    for j in range(0,16):
        if score(rows[i],rows[j]) == 2:
            DD_array=np.append(DD_array,G_dpx_16x16[i][j])
            DD_uniq=np.append(DD_uniq,G_dpx_16x16[i][j])
            if (rows[i]+rows[j] == rows[j][1]+rows[j][0]+rows[i][1]+rows[i][0]): #if symmetric, then do it again
                DD_uniq=np.append(DD_uniq,G_dpx_16x16[i][j])
        elif score(rows[i],rows[j]) == 1:
            SD_array=np.append(SD_array,G_dpx_16x16[i][j])   
            #There are no symmetric single defect dimers
        else:
            ND_array=np.append(ND_array,G_dpx_16x16[i][j])
            ND_uniq=np.append(ND_uniq,G_dpx_16x16[i][j])
            if (rows[i]+rows[j] == rows[j][1]+rows[j][0]+rows[i][1]+rows[i][0]):
                ND_uniq=np.append(ND_uniq,G_dpx_16x16[i][j])    

SD_uniq=SD_array

ND_uniq.sort()
SD_uniq.sort()
DD_uniq.sort()

ND_uniq=ND_uniq[0::2]
SD_uniq=SD_uniq[0::2]
DD_uniq=DD_uniq[0::2]

print(np.mean(ND_array),"+/-",np.std(ND_array))
print(np.mean(SD_array),"+/-",np.std(SD_array))
print(np.mean(DD_array),"+/-",np.std(DD_array))

print(np.mean(ND_uniq),"+/-",np.std(ND_uniq))
print(np.mean(SD_uniq),"+/-",np.std(SD_uniq))
print(np.mean(DD_uniq),"+/-",np.std(DD_uniq))

"""
plt.hist(DD_array,bins=100)
plt.hist(SD_array,bins=100)
plt.hist(ND_array,bins=100)
"""

f1=plt.figure()
f2=plt.figure()
f3=plt.figure()
f4=plt.figure()

ax1=f1.add_subplot(111)
ax2=f2.add_subplot(111)
ax3=f3.add_subplot(111)
ax4=f4.add_subplot(111)

ax1.hist(DD_uniq,bins=500)
ax1.hist(SD_uniq,bins=500)
ax1.hist(ND_uniq,bins=500)

ax2.hist(DD_uniq,bins=20)
ax3.hist(SD_uniq,bins=20)
ax4.hist(ND_uniq,bins=20)

plt.show()