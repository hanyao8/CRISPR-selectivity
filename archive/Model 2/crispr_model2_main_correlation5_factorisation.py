#crispr model 2: incorporating nucleotide correlation - attempt utilising partition function factorisation
#WORKING!

import numpy as np
import math

k_B=1.38e-23
T=310
beta=1/k_B/T
BE=-0.038*(1.6e-19) #Average binding energy corresponding to the below sequence used, derived from Santa Lucia 1998

p_A=0.25
p_T=0.25
p_C=0.25
p_G=0.25

"""
p_A=0.1
p_T=0.2
p_C=0.3
p_G=0.4
"""

p=[0,0,0,0]
p[0]=p_A
p[1]=p_T
p[2]=p_C
p[3]=p_G


n_A=0
n_T=0
n_C=0
n_G=0

ts=["T","T","T","A","T","A","T","A","C","T","T","T","T","T","G","T","T","T","T","G"]
ts_len=len(ts) #Length of the targeting sequence

cs=[]
p_cs=1

ts_n=[] #A numerical targeting sequence array can be constructed, with A=0, T=1, C=2, G=3
cs=[]
cs_n=[]

A_pos=[]
T_pos=[]
C_pos=[]
G_pos=[]

for a in range(0,ts_len):
    if ts[a]=='A':
        n_A+=1
        ts_n.append(0)
        cs.append('T')
        cs_n.append(1)
        A_pos.append(a)

    if ts[a]=='T':
        n_T+=1
        ts_n.append(1)
        cs.append('A')
        cs_n.append(0)
        T_pos.append(a)

    if ts[a]=='C':
        n_C+=1
        ts_n.append(2)
        cs.append('G')
        cs_n.append(3)
        C_pos.append(a)

    if ts[a]=='G':
        n_G+=1
        ts_n.append(3)
        cs.append('C')
        cs_n.append(2)
        G_pos.append(a)

    if ts[a]!='A' and  ts[a]!='T' and ts[a]!='C' and ts[a]!='G':
        raise Exception("Invalid Nucleotide")
        
        
#Note that the correlated probabilities are only defined for a specific direction
#Either 3' --> 5' or 5' --> 3'

#Taking the 3' --> 5' case:
#and using the convention XY = P(Y/X) where X is closer to the 3' and Y is closer to the 5' end
#i.e. the probability that the adjacent nucleotide in the 5' end direction is Y, given
#that the original one is X
#Simple consider that the 'left' end of the strand is the 3' end, and the right end is the 5' end.

cm=[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
#correlation matrix
#AA AT AC AG
#TA TT TC TG
#CA CT CC CG
#GA GT GC GG

#ATCG, A=0, T=1, C=2, G=3

cm[0][0]=0.25 #AA
cm[0][1]=0.25 #AT
cm[0][2]=0.25 #AC
cm[0][3]=0.25 #AG

cm[1][0]=0.25 #TA
cm[1][1]=0.25 #TT
cm[1][2]=0.25 #TC
cm[1][3]=0.25 #TG

cm[2][0]=0.25 #CA
cm[2][1]=0.25 #CT
cm[2][2]=0.25 #CC
cm[2][3]=0.25 #CG

cm[3][0]=0.25 #GA
cm[3][1]=0.25 #GT
cm[3][2]=0.25 #GC
cm[3][3]=0.25 #GG

"""
cm[0][0]=0.23 #AA
cm[0][1]=0.24 #AT
cm[0][2]=0.26 #AC
cm[0][3]=0.27 #AG

cm[1][0]=0.22 #TA
cm[1][1]=0.25 #TT
cm[1][2]=0.25 #TC
cm[1][3]=0.28 #TG

cm[2][0]=0.1 #CA
cm[2][1]=0.2 #CT
cm[2][2]=0.3 #CC
cm[2][3]=0.4 #CG

cm[3][0]=0.2 #GA
cm[3][1]=0.21 #GT
cm[3][2]=0.29 #GC
cm[3][3]=0.3 #GG
"""

if cs[0]=='A':
    p_cs=p_A
if cs[0]=='T':
    p_cs=p_T
if cs[0]=='C':
    p_cs=p_C
if cs[0]=='G':
    p_cs=p_G

for a in range(0,ts_len-1):
    if cs[a]=='A':
        if cs[a+1]=='A':
            p_cs*=cm[0][0]
        if cs[a+1]=='T':
            p_cs*=cm[0][1]
        if cs[a+1]=='C':
            p_cs*=cm[0][2]
        if cs[a+1]=='G':
            p_cs*=cm[0][3]

    if cs[a]=='T':
        if cs[a+1]=='A':
            p_cs*=cm[1][0]
        if cs[a+1]=='T':
            p_cs*=cm[1][1]
        if cs[a+1]=='C':
            p_cs*=cm[1][2]
        if cs[a+1]=='G':
            p_cs*=cm[1][3]          
            
    if cs[a]=='C':
        if cs[a+1]=='A':
            p_cs*=cm[2][0]
        if cs[a+1]=='T':
            p_cs*=cm[2][1]
        if cs[a+1]=='C':
            p_cs*=cm[2][2]
        if cs[a+1]=='G':
            p_cs*=cm[2][3]
            
    if cs[a]=='G':
        if cs[a+1]=='A':
            p_cs*=cm[3][0]
        if cs[a+1]=='T':
            p_cs*=cm[3][1]
        if cs[a+1]=='C':
            p_cs*=cm[3][2]
        if cs[a+1]=='G':
            p_cs*=cm[3][3]
            
#Incorporating non-degenerate defect energies
#Expressing this information in a matrix
#ATCG 0123

DE=[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]

DE[0][0]=0 #0.1*BE #AA

DE[1][0]=DE[0][1]=BE #TA
DE[1][1]=0 #0.1*BE #TT

DE[2][0]=DE[0][2]=0 #CA
DE[2][1]=DE[1][2]=0 #CT
DE[2][2]=0 #0.1*BE #CC

DE[3][0]=DE[0][3]=0 #GA
DE[3][1]=DE[1][3]=0 #GT
DE[3][2]=DE[2][3]=BE #GC
DE[3][3]=0 #0.1*BE #GG
        

Z=1
#Initiation of the partition function factorisation- i.e. the first nucleotide position of the sequence
#Defining a vector to aid the calculation of Z
#ATCG 0123
v_Z=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
re_ord_list=[0,4,8,12,1,5,9,13,2,6,10,14,3,7,11,15]
re_ord_list_A=re_ord_list[0:4]
re_ord_list_T=re_ord_list[4:8]
re_ord_list_C=re_ord_list[8:12]
re_ord_list_G=re_ord_list[12:16]
for i in range(0,16):
    v_Z[i]=p[int(math.floor(i/4))]*np.exp(-beta*DE[int(math.floor(i/4))][ ts_n[0] ])

for pos in range(1,ts_len):
    v_Z_temp=[]
    #v_Z_temp=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    for i in range(0,len(v_Z)):
        #v_Z_temp[i]=v_Z[i] * cm[int(math.floor(i/4))][i%4] * np.exp(-beta*DE[][])
        v_Z_temp.append(v_Z[i] * cm[int(math.floor(i/4))][int(i%4)] * np.exp(-beta*DE[int(i%4)][ ts_n[pos] ]))
    v_Z=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    v_Z_temp2=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    for i in range(0,len(v_Z_temp)):
        v_Z_temp2[i]=v_Z_temp[re_ord_list[i]]
    v_Z[0]=v_Z[1]=v_Z[2]=v_Z[3]=sum(v_Z_temp2[0:4])
    v_Z[4]=v_Z[5]=v_Z[6]=v_Z[7]=sum(v_Z_temp2[4:8])
    v_Z[8]=v_Z[9]=v_Z[10]=v_Z[11]=sum(v_Z_temp2[8:12])
    v_Z[12]=v_Z[13]=v_Z[14]=v_Z[15]=sum(v_Z_temp2[12:16])
    
#Division by four because the last summation should not take place
Z=sum(v_Z)/4
    
q_compl=p_cs * np.exp(-beta*ts_len*BE)
S=1/((Z/q_compl)-1) 