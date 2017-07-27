#The main feature of main5 is that the loop relating to the targeted sequence begins from the complementary limit (i.e. n_ncp=0)
#Note that n_, n_cp_, n_ncp_ refer to the targeting sequence; while m_, s_ refer to the targeted sequence

import numpy as np
import math

def nCr(n,r):
    return math.factorial(n)/math.factorial(n-r)/math.factorial(r)

k_B=1.38e-23
T=300
beta=1/k_B/T
BE=-0.1*(1.6e-19) #epsilon- given by hydrogen bond interactions between complementary NTs

p_A=0.1
p_T=0.2
p_C=0.3
p_G=0.4 #Note none of the probabilities can be set to zero

#if p_A+p_T+p_C+p_G !=1 :
#    raise Exception("Invalid probability sum")

n_A=0
n_T=0
n_C=0
n_G=0

ts=["T","T","T","A","T","A","T","A","C","T","T","T","T","T","G","T","T","T","T","G"]
ts_len=len(ts) #Length of the targeting sequence

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


"""Big Loop
"""

 #Contains the product of the degeneracy and probability of each possible combination within a specific n_cp

Z=0
#Define function instead?
#n_cp=0 #number of complementary pairs

Z_comp_list=[] #Contains 
Dp_list=[]
for i in range(0,ts_len+1):
    Dp_list.append(0)
    Z_comp_list.append(0)


for i in range(0,ts_len+1):
    print Z
    print Z_comp_list
    for j in range(0,ts_len-i+1):
        for k in range(0,ts_len-i-j+1):
            l=ts_len-i-j-k
            
            #for n_cp_T in range(max(0,i-(ts_len-n_T)),min(i,n_T)+1): #The number of complementary A that are attached to T nucleotides in the targeting sequence 
            for n_ncp_T in range(max(0,n_T-i),min(n_T,ts_len-i)+1):
                #for n_cp_A in range(max(0,j-(ts_len-n_A)),min(j,n_A)+1): #The number of complementary T that are attached to A in the ts
                for n_ncp_A in range(max(0,n_A-j),min(n_A,ts_len-j)+1):
                    #for n_cp_G in range(max(0,k-(ts_len-n_G)),min(k,n_G)+1): #The number of complementary C that are attached to G in the ts
                    for n_ncp_G in range(max(0,n_G-k),min(n_G,ts_len-k)+1):
                        #for n_cp_C in range(max(0,l-(ts_len-n_C)),min(l,n_C)+1): #The number of complementary G that are attached to C in the ts
                        for n_ncp_C in range(max(0,n_C-l),min(n_C,ts_len-l)+1):
                            #print i
                            #print i,j#,k,l,n_cp_T,n_cp_A,n_cp_G,n_cp_C

                            n_cp_A=n_A-n_ncp_A
                            n_cp_T=n_T-n_ncp_T
                            n_cp_C=n_C-n_ncp_C
                            n_cp_G=n_G-n_ncp_G

                            m_A=i-n_cp_T
                            m_T=j-n_cp_A
                            m_C=k-n_cp_G
                            m_G=l-n_cp_C
                            n_cp=n_cp_T+n_cp_A+n_cp_G+n_cp_C
                            n_ncp=n_ncp_T+n_ncp_A+n_ncp_G+n_ncp_C

                            
                            D=nCr(n_T,n_cp_T)*nCr(n_A,n_cp_A)*nCr(n_G,n_cp_G)*nCr(n_C,n_cp_C) * math.factorial(n_ncp)/math.factorial(m_A)/math.factorial(m_T)/math.factorial(m_C)/math.factorial(m_G)
                            
                            Z+=D*(p_A**i)*(p_T**j)*(p_C**k)*(p_G**l)*np.exp(-beta*n_cp*BE)
                            Z_comp_list[n_cp]+=D*(p_A**i)*(p_T**j)*(p_C**k)*(p_G**l)*np.exp(-beta*n_cp*BE)
                            Dp_list[n_cp]+=D*(p_A**i)*(p_T**j)*(p_C**k)*(p_G**l)

                                    
                            for m_A_T in range(max(0,m_A-n_ncp_A-n_ncp_G-n_ncp_C),min(m_A,n_ncp_T)+1): #s_ stands for surplus- above the 
                                m_A_nT=i-n_cp_T-m_A_T
                                for m_A_C in range(max(0,(m_A-m_A_T)-n_ncp_A-n_ncp_G),min(m_A-m_A_T,n_ncp_C)+1): #<--*
                                    for m_A_G in range(max(0,(m_A-m_A_T-m_A_C)-n_ncp_A),min(m_A-m_A_T-m_A_C,n_ncp_G)+1):
                                        m_A_A=m_A-m_A_T-m_A_C-m_A_G
                                        
                                        for m_T_A in range(max(0,m_T-(n_ncp_T-m_A_T)-(n_ncp_C-m_A_C)-(n_ncp_G-m_A_G)),min(m_T,n_ncp_A-m_A_A)+1): #Here- you need to consider the n_ncp spaces taken up by A nucleotides in the targeted sequence

                                            s_T_nA=j-n_cp_A-m_T_A
                                            for m_T_C in range(max(0,m_T-m_T_A-(n_ncp_T-m_A_T)-(n_ncp_G-m_A_G)),min(m_T-m_T_A,n_ncp_C-m_A_C)+1):
                                                for m_T_G in range(max(0,m_T-m_T_A-m_T_C-(n_ncp_T-m_A_T)),min(m_T-m_T_A-m_T_C,n_ncp_G-m_A_G)+1):
                                                    s_T_T=m_T-m_T_A-m_T_C-m_T_G
                                                    
                                                    for m_C_G in range(max(0,m_C-(n_ncp_T-m_A_T-s_T_T)-(n_ncp_A-m_A_A-m_T_A)-(n_ncp_C-m_A_C-m_T_C)),min(m_C,n_ncp_G-m_A_G-m_T_G)+1):
                                                    
                                                        s_C_nG=k-n_cp_G-m_C_G
                                                        for m_C_A in range(max(0,m_C-m_C_G-(n_ncp_T-m_A_T-s_T_T)-(n_ncp_C-m_A_C-m_T_C)),min(n_ncp_A-m_T_A-m_A_A,m_C-m_C_G)+1):
                                                            for m_C_T in range(max(0,m_C-m_C_G-m_C_A-(n_ncp_C-m_A_C-m_T_C)),min(n_ncp_T-m_A_T-s_T_T,m_C-m_C_A-m_C_G)+1):
                                                                m_C_C=m_C-m_C_G-m_C_T-m_C_A      
                                                                
                                                                m_G_C=n_ncp_C-m_A_C-m_T_C-m_C_C
                                                                m_G_A=n_ncp_A-m_A_A-m_T_A-m_C_A
                                                                m_G_T=n_ncp_T-m_A_T-s_T_T-m_C_T
                                                                m_G_G=n_ncp_G-m_A_G-m_T_G-m_C_G
                                                                

                                                                    
                                                                if (m_A_T != 0) or (m_T_A != 0) or (m_G_C != 0) or (m_C_G != 0): 
                                                                    d = math.factorial(n_T)/math.factorial(n_cp_T+m_A_T)/math.factorial(s_T_T)/math.factorial(m_C_T)/math.factorial(m_G_T)*math.factorial(n_A)/math.factorial(n_cp_A+m_T_A)/math.factorial(m_A_A)/math.factorial(m_C_A)/math.factorial(m_G_A)*math.factorial(n_G)/math.factorial(n_cp_G+m_C_G)/math.factorial(m_A_G)/math.factorial(m_T_G)/math.factorial(m_G_G)*math.factorial(n_C)/math.factorial(n_cp_C+m_G_C)/math.factorial(m_A_C)/math.factorial(m_T_C)/math.factorial(m_C_C)
                                                                else:
                                                                    d = 0
                                                                    
                                                                Z -= d*(p_A**i)*(p_T**j)*(p_C**k)*(p_G**l)*np.exp(-beta*n_cp*BE)
                                                                Z_comp_list[n_cp] -= d*(p_A**i)*(p_T**j)*(p_C**k)*(p_G**l)*np.exp(-beta*n_cp*BE)
                                                                Dp_list[n_cp]-=d*(p_A**i)*(p_T**j)*(p_C**k)*(p_G**l)
                                                                    
    
print Z_comp_list

q_compl=p_cs * np.exp(-beta*ts_len*BE)
S=1/(Z/q_compl-1)

print S