#CRISPR CAS9 Investigating Specificity: Model 1 

#MODEL USING PARTITION FUNCTION FACTORISATION

"""The aim is to create a model to express the specificity S of a user-defined
targeting sequence (usually 20 nucleotides in length) to a targeted section of a
DNA strand in terms of the length of the targeting sequence, epsilon and the 
relative abundance of nucleotides A, T, C and G

Key assumptions used in the basic 'Model 1' will be that each complementary pair
is associated with a binding energy epsilon and non-complementary pairs have 
zero binding energy
"""

"""
Alternative approach using partition function factorisation
"""


import numpy as np
import math

def nCr(n,r):
    return math.factorial(n)/math.factorial(n-r)/math.factorial(r)

k_B=1.38e-23
T=310 #(High temp. approximation?)
beta=1/k_B/T
BE=-0.0545*(1.6e-19) #epsilon- given by hydrogen bond interactions between complementary NTs
N_G=(64444167-499910)

p_A=0.2794
p_T=0.2826
p_C=0.2176
p_G=0.2204


"""
p_A=0.25
p_T=0.25
p_C=0.25
p_G=0.25

p_A=0.1
p_T=0.2
p_C=0.3
p_G=0.4 #Note none of the probabilities can be set to zero
"""
#if p_A+p_T+p_C+p_G !=1 :
#    raise Exception("Invalid probability sum")

n_A=0
n_T=0
n_C=0
n_G=0

ts=["T","T","A","T","C","T","G","T","T","C","T","G","G","T","G","T","T","C","G","T"]
#ts=["T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T"]
ts_len=len(ts) #Length of the targeting sequence,

cs=[]
p_cs=1
Z=1

for a in range(0,ts_len):
    if ts[a]=='A':
        n_A+=1
        cs.append('T')
        p_cs*=p_T
        Z*=( 1-p_T + p_T*np.exp(-beta*BE))
    if ts[a]=='T':
        n_T+=1
        cs.append('A')
        p_cs*=p_A
        Z*=( 1-p_A + p_A*np.exp(-beta*BE))
    if ts[a]=='C':
        n_C+=1
        cs.append('G')
        p_cs*=p_G
        Z*=( 1-p_G + p_G*np.exp(-beta*BE))
    if ts[a]=='G':
        n_G+=1
        cs.append('C')
        p_cs*=p_C
        Z*=( 1-p_C + p_C*np.exp(-beta*BE))
    if ts[a]!='A' and  ts[a]!='T' and ts[a]!='C' and ts[a]!='G':
        raise Exception("Invalid Nucleotide")


q_compl=p_cs * np.exp(-beta*ts_len*BE)
S=1/((Z/q_compl)-1)
S_real=S/(N_G-ts_len)/p_cs #Assuming correct binding to to one unique site in the genome
print ts_len,S,S_real

#ln_S=np.log(q_compl)-np.log(Z-q_compl)