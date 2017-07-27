#Translation binding test of model 3
#Algorithm will be similar to the parameter extraction algorithm



import numpy as np
import math
import matplotlib.pyplot as plt


def F_cs(s):
    s=list(s)
    comp={'A':'T','T':'A','C':'G','G':'C'}
    for i in range(0,len(s)):
        s[i]=comp[s[i]]
    s=''.join(s)
    return s
    

nt_2_n={'A':0,'T':1,'C':2,'G':3}


k_B=1.38e-23
N_avo=6.022e23
T=310
beta=1/k_B/T
BE=-0.038*(1.6e-19) #Average binding energy corresponding to the below sequence used, derived from Santa Lucia 1998


#Initialisation of energies involved in the 'Santa Lucia rules'

#The free energy of a duplex appears in the following notation:
#G_dpx[n_W][n_X][n_Y][n_Z] is equivalent to the free energy of the duplex
#3' WX 5'
#5' YZ 3'


G_dpx=[]
for i in range(0,4):
    G_dpx.append([])
    for j in range(0,4):
        G_dpx[i].append([])
        for k in range(0,4):
            G_dpx[i][j].append([0,0,0,0])

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

#Double defects
#Using calculated average 'isolated defect energies'
G_avg_AA=1.20*4184/N_avo
G_avg_AC=1.56*4184/N_avo
G_avg_AG=0.81*4184/N_avo
G_avg_TT=1.15*4184/N_avo
G_avg_TC=1.44*4184/N_avo
G_avg_TG=0.75*4184/N_avo
G_avg_CC=1.69*4184/N_avo
G_avg_GG=0.50*4184/N_avo

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

#Total number of energies does add up to 256

#initiation energies
G_init=[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
G_init[2][3]=G_init[3][2]= +0.98*4184/N_avo
G_init[0][1]=G_init[1][0]= +1.03*4184/N_avo

#(+)Defect initiation energies are missing!!


#ts_min=20
#ts_max=20
ts_min=ts_max=20

with open(r'C:\Users\Choon\Desktop\CRISPR\Code\nt_data\chromosome_20_complete.txt', 'r') as myfile:
    #data=myfile.read()
    data=myfile.read().replace('\n', '')

PAM='CT' #NCT
"""
sites=[]
for i in range(ts_max+(1+len(PAM))-1,len(data)):
    if data[i-len(PAM)+1:i+1]==PAM:
        sites.append(i)
"""

i_list=[]
j_list=[]
S_list=[]
Z_nc_list=[]
#site1=min(sites)
site1=60045
#bindings=64444000
final=64444000
for i in range(ts_min,ts_max+1):
    print i
    cs=data[site1-(1+len(PAM))+1-i:site1-(1+len(PAM))+1]
    ts=F_cs(cs)

    Z=0
    q_comp=0
    S=0
    
    #for j in range(site1,site1+bindings):
    for j in range(site1,final):
        tgds=data[j-(1+len(PAM))+1-i:j-(1+len(PAM))+1]    
        if tgds.count('N')==0:
            G_binding=0
            for k in range(0,i-1):
                G_binding+=G_dpx[ nt_2_n[tgds[k]] ][ nt_2_n[tgds[k+1]] ][ nt_2_n[ts[k]] ][ nt_2_n[ts[k+1]] ]
                
            G_binding+= (G_init[ nt_2_n[tgds[0]] ][ nt_2_n[ts[0]] ] + G_init[ nt_2_n[tgds[i-1]] ][ nt_2_n[ts[i-1]] ])
            Z+=np.exp(-beta*G_binding)
            
            if tgds==cs:
                q_comp+=np.exp(-beta*G_binding)
                
            if (j%1000000)==0:
                print j
                print "running..."
                
            j_list.append(j)
            Z_nc_list.append(Z-q_comp)
                
    S=1/((Z/q_comp)-1)
    i_list.append(i)
    S_list.append(S)
    #Z_nc_list.append(Z-q_comp)
    
x=np.array(j_list)
y=np.array(Z_nc_list)
plt.plot(x,y)
plt.title("nc Component of Z")
plt.ylabel("Z_nc")
plt.xlabel("Binding Events")
plt.show()

    
    


