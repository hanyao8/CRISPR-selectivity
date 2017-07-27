#Working model with different defect energies, results consistent with partition function factorisation

import numpy as np
import math


k_B=1.38e-23
T=300
beta=1/k_B/T
BE=-0.1*(1.6e-19) #0.1eV

p_A=0.1
p_T=0.2
p_C=0.3
p_G=0.4 


n_A=0
n_T=0
n_C=0
n_G=0

ts=["T","T","T","A","T","A","T","A","C","T","T","T","T","T","G","T","T","T","T","G"]
ts_len=len(ts)

cs=[]
p_cs=1

for a in range(0,ts_len):
    if ts[a]=='A':
        n_A+=1
        cs.append('T')
        p_cs*=p_T

    if ts[a]=='T':
        n_T+=1
        cs.append('A')
        p_cs*=p_A

    if ts[a]=='C':
        n_C+=1
        cs.append('G')
        p_cs*=p_G

    if ts[a]=='G':
        n_G+=1
        cs.append('C')
        p_cs*=p_C

    if ts[a]!='A' and  ts[a]!='T' and ts[a]!='C' and ts[a]!='G':
        raise Exception("Invalid Nucleotide")
        

#Incorporating non-degenerate defect energies
#Expressing this information in a matrix
#ATCG 0123

DE=[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]

DE[0][0]=0 #AA

DE[1][0]=DE[0][1]=BE #TA
DE[1][1]=0 #TT

DE[2][0]=DE[0][2]=0 #CA
DE[2][1]=DE[1][2]=0 #CT
DE[2][2]=0 #CC

DE[3][0]=DE[0][3]=0 #GA
DE[3][1]=DE[1][3]=0 #GT
DE[3][2]=DE[2][3]=BE #GC
DE[3][3]=0 #GG
        
        
"""Big Loop
"""


Z=0
Z_comp_list=[]
for i in range(0,ts_len+1):
    Z_comp_list.append(0)

for i in range(0,ts_len+1):
    print Z
    #print Z_comp_list
    for j in range(0,ts_len-i+1):
        for k in range(0,ts_len-i-j+1):
            l=ts_len-i-j-k
            

            for m_A_T in range(max(0,i-n_A-n_C-n_G),min(i,n_T)+1):
                for m_A_C in range(max(0,i-m_A_T-n_A-n_G),min(i-m_A_T,n_C)+1):
                    for m_A_G in range(max(0,i-m_A_C-m_A_T-n_A),min(i-m_A_T-m_A_C,n_G)+1):
                        m_A_A=i-m_A_T-m_A_C-m_A_G
                                        
                        for m_T_A in range(max(0,j-(n_T-m_A_T)-(n_C-m_A_C)-(n_G-m_A_G)),min(j,n_A-m_A_A)+1): #Here- you need to consider the n_ncp spaces taken up by A nucleotides in the targeted sequence
                            for m_T_C in range(max(0,j-m_T_A-(n_T-m_A_T)-(n_G-m_A_G)),min(j-m_T_A,n_C-m_A_C)+1):
                                for m_T_G in range(max(0,j-m_T_A-m_T_C-(n_T-m_A_T)),min(j-m_T_A-m_T_C,n_G-m_A_G)+1):
                                    m_T_T=j-m_T_A-m_T_C-m_T_G
                                                    
                                    for m_C_G in range(max(0,k-(n_T-m_A_T-m_T_T)-(n_A-m_A_A-m_T_A)-(n_C-m_A_C-m_T_C)),min(k,n_G-m_A_G-m_T_G)+1):  
                                        for m_C_A in range(max(0,k-m_C_G-(n_T-m_A_T-m_T_T)-(n_C-m_A_C-m_T_C)),min(k-m_C_G,n_A-m_T_A-m_A_A)+1):
                                            for m_C_T in range(max(0,k-m_C_G-m_C_A-(n_C-m_A_C-m_T_C)),min(k-m_C_A-m_C_G,n_T-m_A_T-m_T_T)+1):
                                                m_C_C=k-m_C_G-m_C_T-m_C_A      
                                                                
                                                m_G_C=n_C-m_A_C-m_T_C-m_C_C
                                                m_G_A=n_A-m_A_A-m_T_A-m_C_A
                                                m_G_T=n_T-m_A_T-m_T_T-m_C_T
                                                m_G_G=n_G-m_A_G-m_T_G-m_C_G
                                                
                                                n_cp=m_A_T+m_T_A+m_C_G+m_G_C
                                                n_ncp=ts_len-n_cp
                                                                

                                                D = math.factorial(n_T)/math.factorial(m_A_T)/math.factorial(m_T_T)/math.factorial(m_C_T)/math.factorial(m_G_T)*math.factorial(n_A)/math.factorial(m_T_A)/math.factorial(m_A_A)/math.factorial(m_C_A)/math.factorial(m_G_A)*math.factorial(n_G)/math.factorial(m_C_G)/math.factorial(m_A_G)/math.factorial(m_T_G)/math.factorial(m_G_G)*math.factorial(n_C)/math.factorial(m_G_C)/math.factorial(m_A_C)/math.factorial(m_T_C)/math.factorial(m_C_C)
                                                                
                                                
                                                E_tot=m_A_A*DE[0][0] + (m_T_A+m_A_T)*DE[0][1] + m_T_T*DE[1][1] + (m_C_A+m_A_C)*DE[2][0] + (m_C_T+m_T_C)*DE[2][1] + m_C_C*DE[2][2] + (m_G_A+m_A_G)*DE[3][0] + (m_G_T+m_T_G)*DE[3][1] + (m_G_C+m_C_G)*DE[3][2] + m_G_G*DE[3][3]
                                                                    
                                                Z += D*(p_A**i)*(p_T**j)*(p_C**k)*(p_G**l)*np.exp(-beta*E_tot)
                                                Z_comp_list[n_cp] += D*(p_A**i)*(p_T**j)*(p_C**k)*(p_G**l)*np.exp(-beta*E_tot)
                            
                            
    
E_compl=(n_A+n_T)*DE[1][0]+(n_C+n_G)*DE[3][2]

q_compl=p_cs * np.exp(-beta*E_compl)
S=1/(Z/q_compl-1)

print Z


"""

            for n_cp_T in range(max(0,i-(ts_len-n_T)),min(i,n_T)+1): #The number of complementary A that are attached to T nucleotides in the targeting sequence 
                for n_cp_A in range(max(0,j-(ts_len-n_A)),min(j,n_A)+1): #The number of complementary T that are attached to A in the ts
                    for n_cp_G in range(max(0,k-(ts_len-n_G)),min(k,n_G)+1): #The number of complementary C that are attached to G in the ts
                        for n_cp_C in range(max(0,l-(ts_len-n_C)),min(l,n_C)+1): #The number of complementary G that are attached to C in the ts
                            #print i
                            print i,j#,k,l,n_cp_T,n_cp_A,n_cp_G,n_cp_C
                            m_A=i-n_cp_T
                            m_T=j-n_cp_A
                            m_C=k-n_cp_G
                            m_G=l-n_cp_C
                            n_cp=n_cp_T+n_cp_A+n_cp_G+n_cp_C
                            n_ncp=ts_len-n_cp
                            n_ncp_A=n_A-n_cp_A
                            n_ncp_T=n_T-n_cp_T
                            n_ncp_C=n_C-n_cp_C
                            n_ncp_G=n_G-n_cp_G
                            
                            D=nCr(n_T,n_cp_T)*nCr(n_A,n_cp_A)*nCr(n_G,n_cp_G)*nCr(n_C,n_cp_C) * math.factorial(n_ncp)/math.factorial(m_A)/math.factorial(m_T)/math.factorial(m_C)/math.factorial(m_G)
                            
                            Z+=D*(p_A**i)*(p_T**j)*(p_C**k)*(p_G**l)*np.exp(-beta*n_cp*BE)

                                 """   