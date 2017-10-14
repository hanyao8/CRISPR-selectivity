"""Model to investigate selectivity in a CRISPR-CAS9 system, considering the correlated probability of
nucleotides occurring and non-degenerate defect energies
"""

#Correlation Model
#The correlation model is such that it defines the correlated probability between nucleotides
#given a specific direction (i.e. the fact that DNA sequences end in 3-prime and 5-prime ends)

import numpy as np
import math



def nCr(n,r):
    return math.factorial(n)/math.factorial(n-r)/math.factorial(r)

def all_perms(elements):
    if len(elements) <=1:
        yield elements
    else:
        for perm in all_perms(elements[1:]):
            for i in range(len(elements)):
                # nb elements[0:1] works in both string and list contexts
                yield perm[:i] + elements[0:1] + perm[i:]
                

k_B=1.38e-23
T=300
beta=1/k_B/T
BE=-0.1*(1.6e-19) #epsilon- given by hydrogen bond interactions between complementary NTs

p_A=0.25
p_T=0.25
p_C=0.25
p_G=0.25 #Note none of the probabilities can be set to zero

p=[0,0,0,0]
p[0]=p_A
p[1]=p_T
p[2]=p_C
p[3]=p_G

rho=[0,0,0,0]
rho[0]=p_T
rho[1]=p_A
rho[2]=p_G
rho[3]=p_C

#if p_A+p_T+p_C+p_G !=1:
#    raise Exception("Invalid probability sum")

n_A=0
n_T=0
n_C=0
n_G=0

ts=["T","T","T","A","T","A","T","A","C","T","T","T","T","T","G","T","T","T","T","G"]
ts_len=len(ts) #Length of the targeting sequence

ts_n=[] #A numerical targeting sequence array can be constructed
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
            
"""
#ATCG 0123
for n_ncp in range(0,ts_len+1):
    n_cp=ts_len-n_ncp
    
    for i in range(0,ts_len):
        r=n_ncp
        local_p_prod=1
        if r==0:
            local_p_prod*=p[ cs_n[i] ]
            Dp+=local_p_prod
        if r>0:
            #blah
            r-=1    
            
            Dp+=local_p_prod


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
                                                
                                                m_matrix=[[m_A_A,m_A_T,m_A_C,m_A_G],[m_T_A,m_T_T,m_T_C,m_T_G],[m_C_A,m_C_T,m_C_C,m_C_G],[m_G_A,m_G_T,m_G_C,m_G_G]]
                                                
                                                n_cp=m_A_T+m_T_A+m_C_G+m_G_C
                                                n_ncp=ts_len-n_cp
                                                
                                                for a in all_perms(A_pos):
                                                    for b in all_perms(T_pos):
                                                        for c in all_perms(C_pos):
                                                            for d in all_perms(G_pos):
                                                                abcd=a+b+c+d
                                                                ntn=0
                                                                
                                                                DG={}
                                                                
                                                                #ATCG 0123
                                                                #Targeted sequence   m_A_A, m_T_A, m_C_A, m_G_A     m_A_T, m_T_T,...    m_A_C, ...    m_A_G, ...
                                                                #Targeting sequence       A A A A A ...                T T T T ...       C C ...       G G ...
                                                                
                                                                for x in range(ntn,ntn+m_A_A):
                                                                    DG[abcd[x]]=[0,0]
                                                                    #This can be called as DG[n][0] and DG[n][1]
                                                                    ntn=x
                                                                for x in range(ntn,ntn+m_T_A):
                                                                    DG[abcd[x]]=[1,0]
                                                                    ntn=x
                                                                for x in range(ntn,ntn+m_C_A):
                                                                    DG[abcd[x]]=[2,0]
                                                                    ntn=x
                                                                for x in range(ntn,ntn+m_G_A):
                                                                    DG[abcd[x]]=[3,0]
                                                                    ntn=x
                                                                for x in range(ntn,ntn+m_A_T):
                                                                    DG[abcd[x]]=[0,1]
                                                                    ntn=x
                                                                for x in range(ntn,ntn+m_T_T):
                                                                    DG[abcd[x]]=[1,1]
                                                                    ntn=x
                                                                for x in range(ntn,ntn+m_C_T):
                                                                    DG[abcd[x]]=[2,1]
                                                                    ntn=x
                                                                for x in range(ntn,ntn+m_G_T):
                                                                    DG[abcd[x]]=[3,1]
                                                                    ntn=x
                                                                for x in range(ntn,ntn+m_A_C):
                                                                    DG[abcd[x]]=[0,2]
                                                                    ntn=x
                                                                for x in range(ntn,ntn+m_T_C):
                                                                    DG[abcd[x]]=[1,2]
                                                                    ntn=x
                                                                for x in range(ntn,ntn+m_C_C):
                                                                    DG[abcd[x]]=[2,2]
                                                                    ntn=x
                                                                for x in range(ntn,ntn+m_G_C):
                                                                    DG[abcd[x]]=[3,2]
                                                                    ntn=x
                                                                for x in range(ntn,ntn+m_A_G):
                                                                    DG[abcd[x]]=[0,3]
                                                                    ntn=x
                                                                for x in range(ntn,ntn+m_T_G):
                                                                    DG[abcd[x]]=[1,3]
                                                                    ntn=x
                                                                for x in range(ntn,ntn+m_C_G):
                                                                    DG[abcd[x]]=[2,3]
                                                                    ntn=x
                                                                for x in range(ntn,ntn+m_G_G):
                                                                    DG[abcd[x]]=[3,3]
                                                                    ntn=x
                                                                    
                                                                local_list=[]
                                                                for y in range(0,ts_len):
                                                                    local_list.append(DG[y][0])
                                                                    
                                                                P=p[ local_list[0] ]
                                                                for z in range(0,ts_len-1):
                                                                    P*=cm[local_list[z]][local_list[z+1]]
                                                    
                                                                #For non-degenerate defect energies: use ( E_tot=___ )
                                                    
                                                                Z+=(1/math.factorial(m_A_A)/math.factorial(m_A_T)/math.factorial(m_A_C)/math.factorial(m_A_G) /math.factorial(m_T_A)/math.factorial(m_T_T)/math.factorial(m_T_G)/math.factorial(m_T_C) /math.factorial(m_C_A)/math.factorial(m_C_T)/math.factorial(m_C_C)/math.factorial(m_C_G) /math.factorial(m_G_A)/math.factorial(m_G_T)/math.factorial(m_G_C)/math.factorial(m_G_G))*P*np.exp(-beta*n_cp*BE)
                                                                    
                                                            
                                                                
          
 

q_compl=p_cs * np.exp(-beta*ts_len*BE)

S=1/(Z/q_compl-1)
