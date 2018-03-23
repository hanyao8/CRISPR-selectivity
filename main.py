#CRISPR CAS9 Investigating Specificity: Model 1 

#MODEL USING PARTITION FUNCTION FACTORISATION

"""The aim is to create a model to express the specificity S of a user-defined
targeting sequence (usually 20 nucleotides in length) to a targeted section of a
DNA strand in terms of the length of the targeting sequence, epsilon and the 
relative abundance of nucleotides A, T, C and G

Settings available:
11: Model 1 MFT with Z factorisation
12: Model 1 MFT without Z factorisation
21: Model 2 MFT with Z factorisation
31: Model 3 MFT with Z factorisation
32: Model 3 MFT ts composition and entropy vs complementary exponent, Z, S- data generation
33: Model 3 MFT ts_len vs complementary exponent, Z, S- data generation
34: Model 3 MFT GENOME (not ts) entropy vs Z,S- data generation (using a set of 4 standard, unbiased targeting seqs)
35: Model 3 MFT randomised binding energy matrix (fixed ts) (using genome with 'uniform statistics') (unbiased ts and unbiased genome)
36: Model 3 MFT randomised binding energy matrix (fixed matrix) (using genome with 'uniform statistics')
41: Chromo-walk with fixed ts length with ts taken from different sites
42: Chromo-walk at fixed site with varying length
43: Chromo-walk on random genome for model 1
44: Chromo-walk on random genome for model 2
45: Chromo-walk on random genome for model 3
"""

import numpy as np
import os
import io
import time

run_model=35
#run_model=11
#run_model=1,2,3

k_B=1.38e-23
N_avo=6.022e23

T=310 #(High temp. approximation?)
#BE=-0.0545*(1.6e-19) #epsilon- given by hydrogen bond interactions between complementary NTs
BE=-0.038*(1.6e-19)
N_G=(64444167-499910) #Length of chromosome 20 minus 'N' sites
#63944000

#ts='ACATTAGGGT'
ts='CGGTCCGT'
#ts='TTTATATACT'
#ts='TTTATATACTTTTTGTTTTG'
#ts='GCGCGCGCGCGCGCGCGCGC'

########################
#run_model=32 settings
ts_len_32=20

#run_model=36 settings
ts_len_36=17
########################



#########################
#Chromosome-walk settings
#########################
n_iters=1000000
#n_iters=30000
#n_iters=5
#n_iters=int(1e6)
PAM='CT' #NCT
ts_len41=8
ts_len_start42=15




#######################
#Statistical Parameters
#######################



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
"""
p_A=0.2794
p_T=0.2826
p_C=0.2176
p_G=0.2204
"""
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

"""
cm[0][0]=0.1 #AA
cm[0][1]=0.2 #AT
cm[0][2]=0.3 #AC
cm[0][3]=0.4 #AG
        
cm[1][0]=0.1 #TA
cm[1][1]=0.2 #TT
cm[1][2]=0.3 #TC
cm[1][3]=0.4 #TG

cm[2][0]=0.1 #CA
cm[2][1]=0.2 #CT
cm[2][2]=0.3 #CC
cm[2][3]=0.4 #CG

cm[3][0]=0.1 #GA
cm[3][1]=0.2 #GT
cm[3][2]=0.3 #GC
cm[3][3]=0.4 #GG
"""


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
cm[0][0]=0.3148 #AA
cm[0][1]=0.2445 #AT
cm[0][2]=0.1805 #AC
cm[0][3]=0.2602 #AG
        
cm[1][0]=0.1953 #TA
cm[1][1]=0.3208 #TT
cm[1][2]=0.2151 #TC
cm[1][3]=0.2688 #TG

cm[2][0]=0.3425 #CA
cm[2][1]=0.3332 #CT
cm[2][2]=0.2687 #CC
cm[2][3]=0.0556 #CG

cm[3][0]=0.2800 #GA
cm[3][1]=0.2317 #GT
cm[3][2]=0.2175 #GC
cm[3][3]=0.2708 #GG
"""



#########################
#Thermodynamic Parameters
#########################

#Incorporating non-degenerate defect energies
#Expressing this information in a matrix
#ATCG 0123

DE=[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]



#DE[1][0]=DE[0][1]=BE #TA
#DE[3][2]=DE[2][3]=BE #GC



DE[0][0]=2*1.20*4184/N_avo #0.1*BE #AA

DE[1][0]=DE[0][1]=-2*0.41*4184/N_avo #TA
DE[1][1]=2*1.15*4184/N_avo #0.1*BE #TT

DE[2][0]=DE[0][2]=2*1.56*4184/N_avo #CA
DE[2][1]=DE[1][2]=2*1.44*4184/N_avo #CT
DE[2][2]=2*1.69*4184/N_avo #0.1*BE #CC

DE[3][0]=DE[0][3]=2*0.81*4184/N_avo #GA
DE[3][1]=DE[1][3]=2*0.75*4184/N_avo #GT
DE[3][2]=DE[2][3]=-2*1.04*4184/N_avo #GC
DE[3][3]=2*0.50*4184/N_avo #0.1*BE #GG



#########################
#SantaLucia Parameters
#########################

#Initialisation of energies involved in the 'Santa Lucia rules'

#The free energy of a duplex appears in the following notation:
#G_dpx[n_W][n_X][n_Y][n_Z] is equivalent to the free energy of the duplex
#5' WX 3'
#3' YZ 5'


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
G_avg=np.zeros(8)
G_avg[0]=G_avg_AA=1.20*4184/N_avo
G_avg[1]=G_avg_AC=1.56*4184/N_avo
G_avg[2]=G_avg_AG=0.81*4184/N_avo
G_avg[3]=G_avg_TT=1.15*4184/N_avo
G_avg[4]=G_avg_TC=1.44*4184/N_avo
G_avg[5]=G_avg_TG=0.75*4184/N_avo
G_avg[6]=G_avg_CC=1.69*4184/N_avo
G_avg[7]=G_avg_GG=0.50*4184/N_avo

"""
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
"""

#Total number of energies does add up to 256






#Initiation Energies
"""
G_init=[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
G_init[2][3]=G_init[3][2]=0# +0.98*4184/N_avo
G_init[0][1]=G_init[1][0]=0# +1.03*4184/N_avo#
"""    

#initiation energies
#WC Pairs
G_init=[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
G_init[2][3]=G_init[3][2]= +0.98*4184/N_avo
G_init[0][1]=G_init[1][0]= +1.03*4184/N_avo

#non-WC Pairs, energy values due for review
G_init[0][0]= +1.005*4184/N_avo
G_init[0][2]=G_init[2][0]= +1.005*4184/N_avo
G_init[0][3]=G_init[3][0]= +1.005*4184/N_avo
G_init[1][1]= +1.005*4184/N_avo
G_init[1][2]=G_init[2][1]= +1.005*4184/N_avo
G_init[1][3]=G_init[3][1]= +1.005*4184/N_avo
G_init[2][2]= +1.005*4184/N_avo
G_init[3][3]= +1.005*4184/N_avo

   
         



















#########################

#End of user interface

#########################



def F_cs(s): #Find complementary sequence
    s=list(s)
    comp={'A':'T','T':'A','C':'G','G':'C'}
    for i in range(0,len(s)):
        s[i]=comp[s[i]]
    s=''.join(s)
    return s
    

nt_2_n={'A':0,'T':1,'C':2,'G':3}


def cm_2_p(cm):
    cm_tp=np.transpose(cm)
    eigen_soln=np.linalg.eig(cm_tp)
    for i in range(0,len(eigen_soln[0])):
        if abs(np.real(eigen_soln[0][i])-1)<1e-5:
            cm_2_p_soln=np.real(eigen_soln[1][:,i])
            return(cm_2_p_soln/sum(cm_2_p_soln))
    





constants=[k_B,N_avo]
params1=[T,BE,N_G]
p=[p_A,p_T,p_C,p_G]
model1_params=np.array([constants,params1,p])
model2_params=np.array([constants,params1,p,cm,DE])
model3_params=np.array([constants,params1,p,cm,G_dpx,G_avg,G_init])

#ts='TTATCTGTTCTGGTGTTCGT'
ts_len=len(ts)

cs=[]
        
for a in range(0,ts_len):
    if ts[a]=='A':
        cs.append('T')
    if ts[a]=='T':
        cs.append('A')
    if ts[a]=='C':
        cs.append('G')
    if ts[a]=='G':
        cs.append('C')
    if ts[a]!='A' and  ts[a]!='T' and ts[a]!='C' and ts[a]!='G':
        raise Exception("Invalid Nucleotide")

if run_model==11:
    import model1
    sim11=model1.PF(model1_params,ts)

    print("S_real=",sim11._PF__S_real,"q_exp=",sim11._PF__q_exp,"Z=",sim11._PF__Z)

    print(sim11._PF__S_real,sim11._PF__q_compl*(4**20),sim11._PF__Z*sim11._PF__N_G)

if run_model==12:
    import model1
    sim12=model1.NPF(model1_params,ts)
    print(sim12._NPF__S)
    
if run_model==21:
    import model2
    sim21=model2.PF(model2_params,ts)
    print("S_real=",sim21._PF__S_real,"q_comp_exp=",sim21._PF__q_exp,"Z=",sim21._PF__Z)

if run_model==31:
    import model3
    sim31=model3.PF(model3_params,ts)
    print("S_real=",sim31._PF__S_real,"q_comp_exp=",sim31._PF__exp_comp,"Z=",sim31._PF__Z)
    
if run_model==32:
    import model3
    data_sheet=io.open("C:\\Users\\Choon\\Desktop\\CRISPR\\Code\\results\\result41.txt",'a') #w indicates that the file is writable
    #data_sheet=os.open("C:\\Users\\Choon\\Desktop\\CRISPR\\Code\\results\\result1.txt",'a') #w indicates that the file is writable
    #data_sheet=open("C:\Users\Choon\Desktop\CRISPR\Code\results\result1.txt","w+")
    
    for i in range(0,n_iters):
        
        if i%5000==0:
            print("running... i=",i)
    
        ts=np.random.choice(['A','T','C','G'])
        for i in range(1,ts_len_32):
            ts+=np.random.choice(['A','T','C','G'])
        
        cs=F_cs(ts)
        ts_n=[nt_2_n[ts[i]] for i in range(0,len(ts))]
        cs_n=[nt_2_n[cs[i]] for i in range(0,len(cs))]
        
        sim32=model3.PF(model3_params,ts)
        
        A_frac=ts.count('A')/ts_len_32
        T_frac=ts.count('T')/ts_len_32
        C_frac=ts.count('C')/ts_len_32
        G_frac=ts.count('G')/ts_len_32
        frac=np.array([A_frac,T_frac,C_frac,G_frac])
        
        """
        info=0.0
        for i in range(0,4):
            for j in range(0,4):
                info+= -(cm[i][j]*p[i])*np.log(cm[i][j]*p[i])
        """

        info1=0.0
        for i in range(0,4):
            if frac[i] > 1/(ts_len_32+1):
                info1+= -frac[i]*np.log(frac[i])
        
        #info2=0.0
        
        
        GC_frac=(ts.count('G')+ts.count('C'))/ts_len_32
        q_comp_exp=sim32._PF__exp_comp
        Z=sim32._PF__Z
        N_G_32=sim32._PF__N_G
        S_real=sim32._PF__S_real
        Boltz_prob=sim32._PF__Boltz_prob

        row= ts +","+ str(ts_len_32) +","+ str(info1)  +","+ str(GC_frac) +","+ str(G_frac) +","+ str(q_comp_exp) +","+ str(Z) +","+ str(N_G_32) +","+ str(S_real) +","+ str(Boltz_prob) + "\n"        
        #row= ts +","+ str(ts_len_32) +","+ str(info1) +","+ str(info2) +","+ str(GC_frac) +","+ str(G_frac) +","+ str(q_comp_exp) +","+ str(Z) +","+ str(N_G_32) +","+ str(S_real) +","+ str(Boltz_prob) + "\n"
        #row= ts +","+ ts_len_32 +","+ info +","+ GC_compo +","+ G_compo +","+ q_comp_exp +","+ Z +","+ N_G_32 +","+ S_real +","+ Boltz_prob + "\n"
        data_sheet.write(row)

if run_model==33: #Looping over sequences with different ts_len
    import model3
    data_sheet=io.open("C:\\Users\\Choon\\Desktop\\CRISPR\\Code\\results\\result5.txt",'a') #w indicates that the file is writable
    
    ts_len_test_min=5
    ts_len_test_max=35
    
    for i in range(0,n_iters): 
        if i%500==0:
            print("running... i=",i) 
            
        for j in range(ts_len_test_min,ts_len_test_max+1):
            ts=np.random.choice(['A','T','C','G'])
            for i in range(1,j):
                ts+=np.random.choice(['A','T','C','G'])
            
            cs=F_cs(ts)
            ts_n=[nt_2_n[ts[i]] for i in range(0,len(ts))]
            cs_n=[nt_2_n[cs[i]] for i in range(0,len(cs))]
            
            sim33=model3.PF(model3_params,ts)
            
            A_frac=ts.count('A')/j
            T_frac=ts.count('T')/j
            C_frac=ts.count('C')/j
            G_frac=ts.count('G')/j
            frac=np.array([A_frac,T_frac,C_frac,G_frac])
            
            """
            info=0.0
            for i in range(0,4):
                for j in range(0,4):
                    info+= -(cm[i][j]*p[i])*np.log(cm[i][j]*p[i])
            """
    
            info1=0.0
            for i in range(0,4):
                if frac[i] > 1/(j+1):
                    info1+= -frac[i]*np.log(frac[i])
            
            #info2=0.0
            
            
            GC_frac=(ts.count('G')+ts.count('C'))/j
            q_comp_exp=sim33._PF__exp_comp
            Z=sim33._PF__Z
            N_G_32=sim33._PF__N_G
            S_real=sim33._PF__S_real
            Boltz_prob=sim33._PF__Boltz_prob
    
            row= ts +","+ str(j) +","+ str(info1)  +","+ str(GC_frac) +","+ str(G_frac) +","+ str(q_comp_exp) +","+ str(Z) +","+ str(N_G_32) +","+ str(S_real) +","+ str(Boltz_prob) + "\n"        
            #row= ts +","+ str(ts_len_32) +","+ str(info1) +","+ str(info2) +","+ str(GC_frac) +","+ str(G_frac) +","+ str(q_comp_exp) +","+ str(Z) +","+ str(N_G_32) +","+ str(S_real) +","+ str(Boltz_prob) + "\n"
            #row= ts +","+ ts_len_32 +","+ info +","+ GC_compo +","+ G_compo +","+ q_comp_exp +","+ Z +","+ N_G_32 +","+ S_real +","+ Boltz_prob + "\n"
            data_sheet.write(row)            
            
if run_model==34:
    import model3
    data_sheet=io.open("C:\\Users\\Choon\\Desktop\\CRISPR\\Code\\results\\genome_entropy3.txt",'a')    
    
    model34_ts_array=['AATACAGTTCCTGCGGA','TTATCTGAACCAGCGGT','CCACTCGAATTAGTGGC','GGAGTGCAATTACTCCCG']
    

    
    #model34_params=model3_params
    #DE34=np.ones((4,4))*2.275*4184/N_avo
    #DE34[1][0]=DE34[0][1]=DE34[3][2]=DE34[2][3]=-1.45*4184/N_avo
    
    #0.4125*4184/N_avo
    
    G_dpx34=np.zeros((4,4,4,4))
    
    #no defects
    G_dpx34[0][0][1][1]=G_dpx34[1][1][0][0]=-1.45*4184/N_avo #AA/TT      
    G_dpx34[0][1][1][0]=-1.45*4184/N_avo #AT/TA
    G_dpx34[1][0][0][1]=-1.45*4184/N_avo #TA/AT
    G_dpx34[2][0][3][1]=G_dpx34[1][3][0][2]=-1.45*4184/N_avo #CA/GT
    G_dpx34[3][1][2][0]=G_dpx34[0][2][1][3]=-1.45*4184/N_avo #GT/CA
    G_dpx34[2][1][3][0]=G_dpx34[0][3][1][2]=-1.45*4184/N_avo #CT/GA
    G_dpx34[3][0][2][1]=G_dpx34[1][2][0][3]=-1.45*4184/N_avo #GA/CT
    G_dpx34[2][3][3][2]=-1.45*4184/N_avo #CG/GC
    G_dpx34[3][2][2][3]=-1.45*4184/N_avo #GC/CG
    G_dpx34[3][3][2][2]=G_dpx34[2][2][3][3]=-1.45*4184/N_avo #GG/CC
    
    #single defects
    G_dpx34[0][0][1][0]=G_dpx34[0][1][0][0]=0.4125*4184/N_avo
    G_dpx34[0][0][1][2]=G_dpx34[2][1][0][0]=0.4125*4184/N_avo
    G_dpx34[0][0][1][3]=G_dpx34[3][1][0][0]=0.4125*4184/N_avo
    G_dpx34[0][1][1][1]=G_dpx34[1][1][1][0]=0.4125*4184/N_avo
    G_dpx34[0][1][1][2]=G_dpx34[2][1][1][0]=0.4125*4184/N_avo
    G_dpx34[0][1][1][3]=G_dpx34[3][1][1][0]=0.4125*4184/N_avo
    G_dpx34[0][2][1][0]=G_dpx34[0][1][2][0]=0.4125*4184/N_avo
    G_dpx34[0][2][1][1]=G_dpx34[1][1][2][0]=0.4125*4184/N_avo
    G_dpx34[0][2][1][2]=G_dpx34[2][1][2][0]=0.4125*4184/N_avo
    G_dpx34[0][3][1][0]=G_dpx34[0][1][3][0]=0.4125*4184/N_avo
    G_dpx34[0][3][1][1]=G_dpx34[1][1][3][0]=0.4125*4184/N_avo
    G_dpx34[0][3][1][3]=G_dpx34[3][1][3][0]=0.4125*4184/N_avo
    
    G_dpx34[1][0][0][0]=G_dpx34[0][0][0][1]=0.4125*4184/N_avo
    G_dpx34[1][0][0][2]=G_dpx34[2][0][0][1]=0.4125*4184/N_avo
    G_dpx34[1][0][0][3]=G_dpx34[3][0][0][1]=0.4125*4184/N_avo
    G_dpx34[1][1][0][1]=G_dpx34[1][0][1][1]=0.4125*4184/N_avo
    G_dpx34[1][1][0][2]=G_dpx34[2][0][1][1]=0.4125*4184/N_avo
    G_dpx34[1][1][0][3]=G_dpx34[3][0][1][1]=0.4125*4184/N_avo
    G_dpx34[1][2][0][0]=G_dpx34[0][0][2][1]=0.4125*4184/N_avo
    G_dpx34[1][2][0][1]=G_dpx34[1][0][2][1]=0.4125*4184/N_avo
    G_dpx34[1][2][0][2]=G_dpx34[2][0][2][1]=0.4125*4184/N_avo
    G_dpx34[1][3][0][0]=G_dpx34[0][0][3][1]=0.4125*4184/N_avo
    G_dpx34[1][3][0][1]=G_dpx34[1][0][3][1]=0.4125*4184/N_avo
    G_dpx34[1][3][0][3]=G_dpx34[3][0][3][1]=0.4125*4184/N_avo
    
    G_dpx34[2][0][3][0]=G_dpx34[0][3][0][2]=0.4125*4184/N_avo
    G_dpx34[2][0][3][2]=G_dpx34[2][3][0][2]=0.4125*4184/N_avo
    G_dpx34[2][0][3][3]=G_dpx34[3][3][0][2]=0.4125*4184/N_avo
    G_dpx34[2][1][3][1]=G_dpx34[1][3][1][2]=0.4125*4184/N_avo
    G_dpx34[2][1][3][2]=G_dpx34[2][3][1][2]=0.4125*4184/N_avo
    G_dpx34[2][1][3][3]=G_dpx34[3][3][1][2]=0.4125*4184/N_avo
    G_dpx34[2][2][3][0]=G_dpx34[0][3][2][2]=0.4125*4184/N_avo
    G_dpx34[2][2][3][1]=G_dpx34[1][3][2][2]=0.4125*4184/N_avo
    G_dpx34[2][2][3][2]=G_dpx34[2][3][2][2]=0.4125*4184/N_avo
    G_dpx34[2][3][3][0]=G_dpx34[0][3][3][2]=0.4125*4184/N_avo
    G_dpx34[2][3][3][1]=G_dpx34[1][3][3][2]=0.4125*4184/N_avo
    G_dpx34[2][3][3][3]=G_dpx34[3][3][3][2]=0.4125*4184/N_avo
    
    G_dpx34[3][0][2][0]=G_dpx34[0][2][0][3]=0.4125*4184/N_avo
    G_dpx34[3][0][2][2]=G_dpx34[2][2][0][3]=0.4125*4184/N_avo
    G_dpx34[3][0][2][3]=G_dpx34[3][2][0][3]=0.4125*4184/N_avo
    G_dpx34[3][1][2][1]=G_dpx34[1][2][1][3]=0.4125*4184/N_avo
    G_dpx34[3][1][2][2]=G_dpx34[2][2][1][3]=0.4125*4184/N_avo
    G_dpx34[3][1][2][3]=G_dpx34[3][2][1][3]=0.4125*4184/N_avo
    G_dpx34[3][2][2][0]=G_dpx34[0][2][2][3]=0.4125*4184/N_avo
    G_dpx34[3][2][2][1]=G_dpx34[1][2][2][3]=0.4125*4184/N_avo
    G_dpx34[3][2][2][2]=G_dpx34[2][2][2][3]=0.4125*4184/N_avo
    G_dpx34[3][3][2][0]=G_dpx34[0][2][3][3]=0.4125*4184/N_avo
    G_dpx34[3][3][2][1]=G_dpx34[1][2][3][3]=0.4125*4184/N_avo
    G_dpx34[3][3][2][3]=G_dpx34[3][2][3][3]=0.4125*4184/N_avo
    
    G_avg34=np.ones(8)*2.275/2*4184/N_avo
    G_init34=np.ones((4,4))*1.005*4184/N_avo
    
    model34_params=np.array([constants,params1,p,cm,G_dpx34,G_avg34,G_init34])
    
    for i in range(0,n_iters): 
        cm_gnm=np.array([np.random.uniform(0,1,4) for i in range(0,4)])
        cm_gnm=np.array([cm_gnm[i]/sum(cm_gnm[i]) for i in range(0,4)])
        cm_gnm_tp=np.transpose(cm_gnm)
        
        p_gnm=cm_2_p(cm_gnm)
        check=np.matmul(cm_gnm_tp,p_gnm)
        if abs(p_gnm[0]-check[0])>1e-3:
            raise Exception()
        #print(cm_gnm)
        #print(p_gnm)
        #print("check:",check)
        
        info_gnm=0.0
        for j in range(0,4):
            for k in range(0,4):
                info_gnm+= -(cm_gnm[j][k]*p_gnm[j])*np.log(cm_gnm[j][k]*p_gnm[j])
            
        model34_params[2]=p_gnm
        model34_params[3]=cm_gnm
            
        if i%500==0:
            print("running... i=",i) 
            
        result=np.zeros(5)
        for j in range(0,4):
            sim34=model3.PF(model34_params,model34_ts_array[j])
            result[0]+=sim34._PF__Z
            result[1]+=sim34._PF__S_real
            result[2]+=sim34._PF__Boltz_prob
            result[3]+=sim34._PF__exp_comp
            result[4]+=sim34._PF__N_G
            #print(result)
        result/=4
        
            
        Z=result[0]
        S_real=result[1]
        Boltz_prob=result[2]
        q_comp_exp=result[3]
        N_G34=result[4]

        row= str(i) +","+ str(info_gnm) +","+ str(Z)  +","+ str(S_real) +","+ str(Boltz_prob)  +","+ str(q_comp_exp) +","+ str(N_G34) + "\n"        
        data_sheet.write(row)            
    
if run_model==35:
    import model3
    import programs.stats_2_BEtensor as s2B
    
    pw_comp_dict={"A":"T","T":"A","C":"G","G":"C"}
    num_comp_dict={0:1,1:0,2:3,3:2}
    comp_dict={"TA":"AT","AT":"TA","AA":"TT","TT":"AA","CT":"GA","AG":"TC","GA":"CT","TC":"AG","GT":"CA","AC":"TG","CA":"GT","TG":"AC","GG":"CC","CC":"GG","CG":"GC","GC":"CG"}
    nt_2_num={"A":0,"T":1,"C":2,"G":3}
    num_2_nt={0:"A",1:"T",2:"C",3:"G"}
    
    """
    def ntscore(a,b):
        num=0
        for i in range(0,2):
            if not(b[i] == pw_comp_dict[a[i]]):
                num+=1
        return(num)
        
    def numscore(a,b,c,d):
        num=0
        if not(a == num_comp_dict[c]):
            num+=1
        if not(b == num_comp_dict[d]):
            num+=1            
        return(num)
    """
    
    cur_time=time.time()
    data_sheet=io.open(os.getcwd()+"\\main35result_%d.txt"%cur_time,'a')    
    data_sheet.write("i" +","+ "info_gnm" +","+ "Z"  +","+ "S_real" +","+ "Boltz_prob"  +","+ "q_comp_exp" +","+ "N_G35" + "\n"     )    
    
    model35_ts_array=['AATACAGTTCCTGCGGA','TTATCTGAACCAGCGGT','CCACTCGAATTAGTGGC','GGAGTGCAATTACTCCCG']
    
    #pairs=["AA","AT","AC","AG","CA","CT","CC","CG","GC","GG","GA","GT","TC","TG","TA","TT"]
    
    BE_stats=np.array([[-1.4180,0.4152,2.2767],[0.512,0.507,0.564]]) #unique
    #BE_stats2=np.array([[-1.4056,0.4152,2.2767],[0.433,0.507,0.543]])

    G_avg35=np.ones(8)*2.275/2*4184/N_avo
    G_init35=np.ones((4,4))*1.005*4184/N_avo    

    #model34_params=model3_params
    #DE34=np.ones((4,4))*2.275*4184/N_avo
    #DE34[1][0]=DE34[0][1]=DE34[3][2]=DE34[2][3]=-1.45*4184/N_avo
    
    #0.4125*4184/N_avo
    
    info_gnm=0.0
    for j in range(0,4):
        for k in range(0,4):
            info_gnm+= -(cm[j][k]*p[j])*np.log(cm[j][k]*p[j])
    
    for i in range(0,n_iters): 
        G_dpx35=s2B.convert(BE_stats)*4184/N_avo


        
        model35_params=np.array([constants,params1,p,cm,G_dpx35,G_avg35,G_init35])
        

            
            
        if i%500==0:
            print("running... i=",i) 
            
        result=np.zeros(5)
        for j in range(0,4):
            sim35=model3.PF(model35_params,model35_ts_array[j],DD_useavg=False)
            result[0]+=sim35._PF__Z
            result[1]+=sim35._PF__S_real
            result[2]+=sim35._PF__Boltz_prob
            result[3]+=sim35._PF__exp_comp
            result[4]+=sim35._PF__N_G
            #print(result)
        result/=4
        
            
        Z=result[0]
        S_real=result[1]
        Boltz_prob=result[2]
        q_comp_exp=result[3]
        N_G35=result[4]

        row= str(i) +","+ str(info_gnm) +","+ str(Z)  +","+ str(S_real) +","+ str(Boltz_prob)  +","+ str(q_comp_exp) +","+ str(N_G35) + "\n"        
        data_sheet.write(row)                
    
if run_model==36:

    import model3
    import programs.stats_2_BEtensor as s2B
    
    pw_comp_dict={"A":"T","T":"A","C":"G","G":"C"}
    num_comp_dict={0:1,1:0,2:3,3:2}
    comp_dict={"TA":"AT","AT":"TA","AA":"TT","TT":"AA","CT":"GA","AG":"TC","GA":"CT","TC":"AG","GT":"CA","AC":"TG","CA":"GT","TG":"AC","GG":"CC","CC":"GG","CG":"GC","GC":"CG"}
    nt_2_num={"A":0,"T":1,"C":2,"G":3}
    num_2_nt={0:"A",1:"T",2:"C",3:"G"}


    cur_time=time.time()
    data_sheet=io.open(os.getcwd()+"\\main36result_%d.txt"%cur_time,'a')    
    #ts +","+ str(ts_len_32) +","+ str(info1)  +","+ str(GC_frac) +","+ str(G_frac) +","+ str(q_comp_exp) +","+ str(Z) +","+ str(N_G_32) +","+ str(S_real) +","+ str(Boltz_prob) +","+ str(info_gnm) + "\n"   
    data_sheet.write("ts" +","+ "ts_len_36" +","+ "ts_info"  +","+ "GC_frac" +","+ "G_frac" +","+ "q_comp_exp" +","+ "Z" +","+ "N_G_36" +","+ "S_real" +","+ "Boltz_prob" +","+ "info_gnm" + "\n"     )    

    #pairs=["AA","AT","AC","AG","CA","CT","CC","CG","GC","GG","GA","GT","TC","TG","TA","TT"]
    
    BE_stats=np.array([[-1.4180,0.4152,2.2767],[0.512,0.507,0.564]]) #unique
    #BE_stats2=np.array([[-1.4056,0.4152,2.2767],[0.433,0.507,0.543]])

    G_dpx36=s2B.convert(BE_stats)*4184/N_avo
    G_avg36=np.ones(8)*2.275/2*4184/N_avo
    G_init36=np.ones((4,4))*1.005*4184/N_avo    
    
    model36_params=np.array([constants,params1,p,cm,G_dpx36,G_avg36,G_init36])

    info_gnm=0.0
    for j in range(0,4):
        for k in range(0,4):
            info_gnm+= -(cm[j][k]*p[j])*np.log(cm[j][k]*p[j])
            
    for i in range(0,n_iters):
        
        if i%5000==0:
            print("running... i=",i)
    
        ts=np.random.choice(['A','T','C','G'])
        for i in range(1,ts_len_36):
            ts+=np.random.choice(['A','T','C','G'])
        
        #cs=F_cs(ts)
        #ts_n=[nt_2_n[ts[i]] for i in range(0,len(ts))]
        #cs_n=[nt_2_n[cs[i]] for i in range(0,len(cs))]
        
        sim36=model3.PF(model36_params,ts,DD_useavg=False)
        
        A_frac=ts.count('A')/ts_len_36
        T_frac=ts.count('T')/ts_len_36
        C_frac=ts.count('C')/ts_len_36
        G_frac=ts.count('G')/ts_len_36
        frac=np.array([A_frac,T_frac,C_frac,G_frac])
        

        info1=0.0
        for i in range(0,4):
            if frac[i] > 1/(ts_len_36+1):
                info1+= -frac[i]*np.log(frac[i])
        
        #info2=0.0
        
        
        GC_frac=(ts.count('G')+ts.count('C'))/ts_len_36
        q_comp_exp=sim36._PF__exp_comp
        Z=sim36._PF__Z
        N_G_36=sim36._PF__N_G
        S_real=sim36._PF__S_real
        Boltz_prob=sim36._PF__Boltz_prob

        row= ts +","+ str(ts_len_36) +","+ str(info1)  +","+ str(GC_frac) +","+ str(G_frac) +","+ str(q_comp_exp) +","+ str(Z) +","+ str(N_G_36) +","+ str(S_real) +","+ str(Boltz_prob) +","+ str(info_gnm) + "\n"   
        #row= ts +","+ str(ts_len_32) +","+ str(info1) +","+ str(info2) +","+ str(GC_frac) +","+ str(G_frac) +","+ str(q_comp_exp) +","+ str(Z) +","+ str(N_G_32) +","+ str(S_real) +","+ str(Boltz_prob) + "\n"
        #row= ts +","+ ts_len_32 +","+ info +","+ GC_compo +","+ G_compo +","+ q_comp_exp +","+ Z +","+ N_G_32 +","+ S_real +","+ Boltz_prob + "\n"
        data_sheet.write(row)

    
if run_model==41:
    import chwalk.cmain
    sim41=chwalk.cmain.sim(model3_params,ts,ts_len41,run_model,n_iters,PAM)

if run_model==42:
    import chwalk.cmain
    sim42=chwalk.cmain.sim(model3_params,ts,ts_len_start42,run_model,n_iters,PAM)

if run_model==43:
    import chwalk.cmain
    sim43=chwalk.cmain.sim(model1_params,ts,ts_len,run_model,n_iters,PAM)

if run_model==44:
    import chwalk.cmain
    sim44=chwalk.cmain.sim(model2_params,ts,ts_len,run_model,n_iters,PAM)
    
if run_model==45:
    import chwalk.cmain
    sim45=chwalk.cmain.sim(model3_params,ts,ts_len,run_model,n_iters,PAM)
	
#C:\\Users\\Choon\\Documents\\GitHub\\CRISPR-selectivity-OO\\chwalk