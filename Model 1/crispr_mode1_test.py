#CRISPR CAS9 Investigating Specificity: Model 1 (Basic)
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
T=300 #(High temp. approximation?)
beta=1/k_B/T
BE=-0.1*(1.6e-19) #epsilon- given by hydrogen bond interactions between complementary NTs

p_A=0.1
p_T=0.5
p_C=0.2
p_G=0.2 #Note none of the probabilities can be set to zero

#if p_A+p_T+p_C+p_G !=1 :
#    raise Exception("Invalid probability sum")

n_A=0
n_T=0
n_C=0
n_G=0

ts=["C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","T"]
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
        
#An alternative approach could be the utilisation of a dictionary
#The calculation of the complementary sequence probability should also be done via the big loop


"""Big Loop
"""

Dp_list=[] #Contains the product of the degeneracy and probability of each possible combination within a specific n_cp
Z_comp_list=[] #Contains 
Z=0
#Define function instead?
#n_cp=0 #number of complementary pairs
for n_cp in range(0,ts_len+1):
    n_ncp=ts_len-n_cp #Number of non-complementary pairs defined
    Dp_list.append(0) 
    
    if n_T-n_ncp >= 0:
        i_start=n_T-n_ncp
    if n_T-n_ncp < 0: #notice that this naively corresponds to an unphysical situation of negative i_start
        i_start=0
    if n_T+n_ncp <= ts_len:
        i_end=n_T+n_ncp
    if n_T+n_ncp > ts_len: #notice that this naively corresponds to an unphysical situation of i_end larger than ts_len
        i_end=ts_len-n_T+n_cp
        
    if i_start<0: #These should not really happen
        raise Exception()
    if i_end>ts_len:
        raise Exception()
    for i in range(i_start,i_end+1): #Number of A nucleotides in targeted sequence
        
        if i<n_T:
            j_start=n_A-(n_ncp-(n_T-i))
            j_end=n_A+n_ncp
        if i==n_T:
            j_start=n_A-n_ncp
            j_end=n_A+n_ncp

        if i>n_T:
            j_start=n_A-n_ncp
            j_end=n_A+(n_ncp-(i-n_T))
            
        if j_start<0:
            j_start=0
        if i+j_end>ts_len or j_end>ts_len: #Need more thought here
            j_end=min(ts_len-i,ts_len-n_ncp)
         
        #if n_A-(ts_len-n_cp) >= 0:
         #   j_start=n_T-(ts_len-n_cp)
        #if n_A-(ts_len-n_cp) < 0: #notice that this naively corresponds to an unphysical situation of negative i_start
         #   j_start=0
        #if n_A+(ts_len-n_cp) <= ts_len:
         #   j_end=n_T+(ts_len-n_cp)
        #if n_A+(ts_len-n_cp) > ts_len: #notice that this naively corresponds to an unphysical situation of i_end larger than ts_len
         #   j_end=ts_len-n_T+n_cp
         
        #j_start=0
        #j_end=ts_len-n_T      
          
        for j in range(j_start,j_end+1): #Number of T nucleotides in targeted sequence
            if i<n_T:
                if j<n_A:
                    k_start=n_G-(n_ncp-(n_T-i+n_A-j))
                    k_end=n_G+n_ncp
                if j==n_A:
                    k_start=n_G-(n_ncp-(n_T-i)) #Note that this is just a particular case of the condition above
                    k_end=n_G+n_ncp
                if j>n_A:
                    k_start=n_G-(n_ncp-(n_T-i))
                    k_end=n_G+(n_ncp-(j-n_A))
                    
            if i==n_T:
                if j<n_A:
                    k_start=n_G-(n_ncp-(n_A-j))
                    k_end=n_G+n_ncp
                if j==n_A:
                    k_start=n_G-n_ncp
                    k_end=n_G+n_ncp
                if j>n_A:
                    k_start=n_G-n_ncp
                    k_end=n_G+(n_ncp-(j-n_A))
                
            if i>n_T:
                if j<n_A:
                    k_start=n_G-(n_ncp-(n_A-j)) #The explanation for this is that the deficit (i.e. j<n_A) is what matters rather than the surplus i>n_T
                    k_end=n_G+(n_ncp-(i-n_T))
                if j==n_A:
                    k_start=n_G-n_ncp
                    k_end=n_G+(n_ncp-(i-n_T))
                if j>n_A:
                    k_start=n_G-n_ncp
                    k_end=n_G+(n_ncp-(i-n_T+j-n_A))
                    
            if k_start<0 or n_ncp-i-j-k_start>0: #Need more thought here
                k_start=max(0,n_ncp-i-j)
            if i+j+k_end>ts_len or k_end>ts_len: #Need more thought here
                k_end=min(ts_len-i-j,ts_len-n_ncp)           
            for k in range(k_start,k_end+1): #Number of C nucleotides in targeted sequence
                l=ts_len-i-j-k #Number of G nucleotides in targeted sequence
                print n_cp,i,j,k,l
                
                omega=0
                
#                i-n_T
#                j-n_A
#                k-n_G
#                l-n_C
                #omega=omega(ts_len,n_cp,i-n_T,j-n_A,k-n_G,l-n_C)
                #defining m_A, m_T, m_C, m_G- as 'floating nucleotides' which are extra, in addition to the fixed complementary pairs
                
                #if :
                
                n_ncp_T_start=n_ncp-n_T-n_C-n_G
                if n_ncp_T_start<0:
                    n_ncp_T_start=0
                    
                if n_ncp_T_start>n_T:
                    n_ncp_T_start=n_T
                
                if n_ncp%2==0:
                    n_ncp_T_end=n_ncp/2
                if n_ncp%2==1:
                    n_ncp_T_end=(n_ncp-1)/2
                    
                
                    
                if n_ncp_T_end>n_T:
                    n_ncp_T_end=n_T
                        
                for n_ncp_T in range(n_ncp_T_start,n_ncp_T_end+1): #all possible combinations of 
                    m_A=i-(n_T-n_ncp_T)
                    if m_A>i:
                        m_A=i
                    if m_A<0:
                        m_A=0
                    
                    n_ncp_A_start=n_ncp-n_A-n_C-n_G
                    if n_ncp_A_start<0:
                        n_ncp_A_start=0
                    
                    if n_ncp_A_start>n_A:
                        n_ncp_A_start=n_A
                        
                    if n_ncp%2==0:
                        n_ncp_A_end=n_ncp/2 #This should be independent of the value of n_ncp_T
                    if n_ncp%2==1:
                        n_ncp_A_end=(n_ncp-1)/2 #This should be independent of the value of n_ncp_T
                        
                    if n_ncp>n_A:
                            n_ncp_A_end=n_A
                            
                    for n_ncp_A in range(n_ncp_A_start,n_ncp_A_end+1):
                        m_T=j-(n_A-n_ncp_A)
                        if m_T>j:
                            m_T=j
                        if m_T<0:
                            m_T=0
                        #if :
                        n_ncp_G_start=n_ncp-n_A-n_C-n_T
                        if n_ncp_G_start<0:
                            n_ncp_G_start=0
                    
                        if n_ncp_G_start>n_G:
                            n_ncp_G_start=n_G
                            
                        n_ncp_G_end=n_ncp-n_ncp_T-n_ncp_A
                            
                        if n_ncp>n_G:
                            n_ncp_G_end=n_G
                            
                        if n_ncp_A+n_ncp_T+n_ncp_G_end > n_ncp:
                            n_ncp_G_end = n_ncp - n_ncp_A - n_ncp_T
                                
                        for n_ncp_G in range(n_ncp_G_start,n_ncp_G_end+1):
                            n_ncp_C=n_ncp-n_ncp_T-n_ncp_A-n_ncp_G
                            if n_ncp_A+n_ncp_T+n_ncp_C+n_ncp_G==n_ncp: #It is actually possible from the above that the sum does not add up to the total number of non-complementary pairs
                                
                                m_C=k-(n_G-n_ncp_G)
                                if m_C>k:
                                    m_C=k
                                if m_C<0:
                                    m_C=0
                            
                                
                                m_G=l-(n_C-n_ncp_C)
                                if m_G>l:
                                    m_G=l
                                if m_G<0:
                                    m_G=0
                                
                                omega += math.factorial(n_ncp) / math.factorial(n_ncp_A) / math.factorial(n_ncp_T) / math.factorial(n_ncp_C) / math.factorial(n_ncp_G)
                                ###Negative values appear here for factorial arguments
                                
                                for a in range(1,i-(n_T-n_ncp_T)+1):
                                    omega -= nCr(n_ncp_T,a)*nCr(n_ncp-n_ncp_T,m_A-a) * math.factorial(n_ncp-m_A)/math.factorial(m_T)/math.factorial(m_C)/math.factorial(m_G) #Check this again, make sure it is right
                                for b in range(1,j-(n_A-n_ncp_A)+1):
                                    omega -= nCr(n_ncp_A,b)*nCr(n_ncp-n_ncp_A,m_T-b) * math.factorial(n_ncp-m_T)/math.factorial(m_C)/math.factorial(m_G)/math.factorial(m_A) #Check this again, make sure it is right
                                for c in range(1,k-(n_G-n_ncp_G)+1):
                                    omega -= nCr(n_ncp_G,c)*nCr(n_ncp-n_ncp_G,m_C-c) * math.factorial(n_ncp-m_C)/math.factorial(m_G)/math.factorial(m_A)/math.factorial(m_T) #Check this again, make sure it is right
                                for d in range(1,l-(n_C-n_ncp_C)+1):
                                    omega -= nCr(n_ncp_C,d)*nCr(n_ncp-n_ncp_C,m_G-d) * math.factorial(n_ncp-m_G)/math.factorial(m_A)/math.factorial(m_T)/math.factorial(m_C) #Check this again, make sure it is right
                                
                            else:
                                omega += 0
                                                            
                            #Subtracting off the sequences which violate the fixed n_cp number within the loop
                            #The subtracted off sequences correspond only to ones with a greater number of complementary pairs
                            
                            for a in range(1,i-(n_T-n_ncp_T)+1):
                                omega -= nCr(n_ncp_T,a)*nCr(n_ncp-n_ncp_T,m_A-a) * math.factorial(n_ncp-m_A)/math.factorial(m_T)/math.factorial(m_C)/math.factorial(m_G) #Check this again, make sure it is right
                            for b in range(1,j-(n_A-n_ncp_A)+1):
                                omega -= nCr(n_ncp_A,b)*nCr(n_ncp-n_ncp_A,m_T-b) * math.factorial(n_ncp-m_T)/math.factorial(m_C)/math.factorial(m_G)/math.factorial(m_A) #Check this again, make sure it is right
                            for c in range(1,k-(n_G-n_ncp_G)+1):
                                omega -= nCr(n_ncp_G,c)*nCr(n_ncp-n_ncp_G,m_C-c) * math.factorial(n_ncp-m_C)/math.factorial(m_G)/math.factorial(m_A)/math.factorial(m_T) #Check this again, make sure it is right
                            for d in range(1,l-(n_C-n_ncp_C)+1):
                                omega -= nCr(n_ncp_C,d)*nCr(n_ncp-n_ncp_C,m_G-d) * math.factorial(n_ncp-m_G)/math.factorial(m_A)/math.factorial(m_T)/math.factorial(m_C) #Check this again, make sure it is right
                                    
                    
                    #if : #if the total number of sub-n_ncp's sum up to n_ncp
                     #   omega += n_ncp_A
                    
                #omega= math.factorial(ts_len) / (math.factorial(i))/(math.factorial(j))/(math.factorial(j))/(math.factorial(l))
                Dp_list[n_cp] += ( omega * (p_A**i)*(p_T**j)*(p_C**k)*(p_G**l) ) #The number of possible combinations is the same as the distinguishable particles statistical weight formula
                
    Z += (Dp_list[n_cp]) * np.exp(-beta*n_cp*BE) 
"""


q_compl=p_cs * np.exp(-beta*ts_len*BE)
S=1/((Z/q_compl)-1)

print S
 
                #sub-degeneracy times by sequence probability here
            

"""
D_list=[]
q_list=[]
Dq_list=[]

for j in range(0,ts_len+1):
    D_list.append( math.factorial(ts_len)/math.factorial(j)/math.factorial(ts_len-j) * 
    (1**j) * (3**(ts_len-j)) 
    )

for i in range(0,ts_len+1):
    Dq_list.append(D_list[i]*np.exp(-beta*i*BE))

#for j in range(0,ts_len+1):
#    q_list.append()

#Dq_list=np.array([])

Z=sum(Dq_list)
q_compl=np.exp(-beta*ts_len*BE)

#Still: what is the thing that has a statistical weight of 1?
S=1/((Z/q_compl)-1)
"""