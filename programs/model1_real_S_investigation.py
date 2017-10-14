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

switch=1
#the switch variable defines the analysis to be carried out
#0,1,2: S(l),S(p),S(BE)

import numpy as np
import math
import matplotlib.pyplot as plt

def nCr(n,r):
    return math.factorial(n)/math.factorial(n-r)/math.factorial(r)

k_B=1.38e-23
T=310
beta=1/k_B/T
BE=-0.054*(1.6e-19) #epsilon- given by hydrogen bond interactions between complementary NTs

p_A=0.2794
p_T=0.2826
p_C=0.2176
p_G=0.2204
"""
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

overall_ts=["T","T","A","T","C","T","G","T","T","C","T","G","G","T","G","T","T","C","G","T"]#,"T","T","T"]#,"T","T","T","T","T","T","T","T","T","T","T","T"]
#overall_ts=["T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T"]
overall_ts_len=len(overall_ts) #Length of the targeting sequence,

#Selectivity as a funciton of ts (or equivalently l)
if switch == 0:
    i_list=[]
    S_list=[]
    S_real_list=[]
    S_3=[]
    
    for i in range(1,overall_ts_len+1):
        ts=overall_ts[:i]
    
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
        S_real=S/p_cs/(3e9)+S
        
        i_list.append(i)
        S_list.append(S)
        S_real_list.append(S_real)
        S_3.append(S*S_real)
        
        print ts_len,i,p_cs,Z
        
        #ln_S=np.log(q_compl)-np.log(Z-q_compl) 
        
    x=np.array(i_list)
    y1=np.array(S_list)
    y2=np.array(S_real_list)
    y3=np.array(S_3)

    #plt.plot(x,y1,label="S_ideal")
    plt.plot(x,y2,label="S_real")
    #plt.plot(x,y3)
    plt.title("S(l)")
    plt.ylabel("Selectivity")
    plt.xlabel("Targeting sequence length")
    plt.legend()
    plt.show()

    
#Selectivity as a function of probability
if switch == 1:
    N_G=3e9
    ts_len=overall_ts_len
    p_list=np.linspace(0,0.999999999999,200000+1)
    
    
    
    def S_ideal(dummy_p):
        return ((dummy_p*np.exp(-beta*BE))**ts_len)/(( 1+dummy_p* (np.exp(-beta*BE)-1) )**ts_len - (dummy_p*np.exp(-beta*BE))**ts_len )
        
    def S_real(dummy_p):
        return (np.exp(-beta*BE*ts_len))/N_G/( ( 1+ dummy_p*(np.exp(-beta*BE)-1) )**ts_len - (dummy_p*np.exp(-beta*BE))**ts_len )    
    
    plt.plot(p_list,S_ideal(p_list),label="S_ideal")
    plt.plot(p_list,S_real(p_list),label="S_real")
    plt.ylim(0,10)
    
    plt.title("S(p_A)")
    plt.ylabel("Selectivity")
    plt.xlabel("p_A")
    plt.legend()
    plt.show()

#Selectivity as a function of Binding energy

if switch == 2:
    N_G=3e9
    ts_len=overall_ts_len
    arg_list=np.linspace(0,10,10000+1)
    
    
    def S_ideal_2(dummy_arg):
        return ((p_A*np.exp(dummy_arg))**ts_len)/(( 1+p_A* (np.exp(dummy_arg)-1) )**ts_len - (p_A*np.exp(dummy_arg))**ts_len )
        
    def S_real_2(dummy_arg):
        return (np.exp(dummy_arg*ts_len))/N_G/( ( 1+ p_A*(np.exp(dummy_arg)-1) )**ts_len - (p_A*np.exp(dummy_arg))**ts_len )    
    
    plt.plot(arg_list,S_ideal_2(arg_list),label="S_ideal")
    plt.plot(arg_list,S_real_2(arg_list),label="S_real")
    plt.ylim(0,10)
    
    plt.title("S(-beta*epsilon)")
    plt.ylabel("Selectivity")
    plt.xlabel("-beta*epsilon")
    plt.legend()
    plt.show()
