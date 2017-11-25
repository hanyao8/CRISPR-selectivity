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
32: Model 3 MFT Composition and entropy vs complementary exponent, Z, S data generation
33: Model 3 MFT ts_len vs complementary exponent, Z, S data generation
41: Chromo-walk with fixed ts length with ts taken from different sites
42: Chromo-walk at fixed site with varying length
43: Chromo-walk on random genome for model 1
44: Chromo-walk on random genome for model 2
45: Chromo-walk on random genome for model 3
"""

import numpy as np
#import os
import io

<<<<<<< HEAD
run_model=32
=======
run_model=11
>>>>>>> f7cdfc3b85b03d755910bb538f9631b652bf4e75
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


"""
p_A=0.25
p_T=0.25
p_C=0.25
p_G=0.25


p_A=0.1
p_T=0.2
p_C=0.3
p_G=0.4
"""
p_A=0.2794
p_T=0.2826
p_C=0.2176
p_G=0.2204

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




#########################
#Thermodynamic Parameters
#########################

#Incorporating non-degenerate defect energies
#Expressing this information in a matrix
#ATCG 0123

DE=[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]



DE[1][0]=DE[0][1]=BE #TA
DE[3][2]=DE[2][3]=BE #GC



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