import numpy as np
import math





                

class PF:
    


    
    def __init__(self,param,targeting_seq):

        k_B=param[0][0]
        self.__T=param[1][0] #(High temp. approximation?)
        self.__beta=1/k_B/self.__T
        self.__BE=param[1][1] #epsilon- given by hydrogen bond interactions between complementary NTs
        self.__N_G=param[1][2]
        N_avo=param[0][1]
        
        self.__p_A=param[2][0]
        self.__p_T=param[2][1]
        self.__p_C=param[2][2]
        self.__p_G=param[2][3]

        p=[0,0,0,0]
        p[0]=self.__p_A
        p[1]=self.__p_T
        p[2]=self.__p_C
        p[3]=self.__p_G
        
        
        n_A=0
        n_T=0
        n_C=0
        n_G=0
        
        ts=list(targeting_seq)
        #ts=["T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T"]
        ts_len=len(ts) #Length of the targeting sequence,

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
        
        cm[0][0]=param[3][0][0] #AA
        cm[0][1]=param[3][0][1] #AT
        cm[0][2]=param[3][0][2] #AC
        cm[0][3]=param[3][0][3] #AG
        
        cm[1][0]=param[3][1][0] #TA
        cm[1][1]=param[3][1][1] #TT
        cm[1][2]=param[3][1][2] #TC
        cm[1][3]=param[3][1][3] #TG
        
        cm[2][0]=param[3][2][0] #CA
        cm[2][1]=param[3][2][1] #CT
        cm[2][2]=param[3][2][2] #CC
        cm[2][3]=param[3][2][3] #CG
        
        cm[3][0]=param[3][3][0] #GA
        cm[3][1]=param[3][3][1] #GT
        cm[3][2]=param[3][3][2] #GC
        cm[3][3]=param[3][3][3] #GG
        
        
        if cs[0]=='A':
            p_cs=self.__p_A
        if cs[0]=='T':
            p_cs=self.__p_T
        if cs[0]=='C':
            p_cs=self.__p_C
        if cs[0]=='G':
            p_cs=self.__p_G
        
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
        
        DE[0][0]=param[4][0][0] #0.1*BE #AA
        
        DE[1][0]=DE[0][1]=param[4][1][0] #TA
        DE[1][1]=param[4][1][1] #0.1*BE #TT
        
        DE[2][0]=DE[0][2]=param[4][2][0] #CA
        DE[2][1]=DE[1][2]=param[4][2][1] #CT
        DE[2][2]=param[4][2][2] #0.1*BE #CC
        
        DE[3][0]=DE[0][3]=param[4][3][0] #GA
        DE[3][1]=DE[1][3]=param[4][3][1] #GT
        DE[3][2]=DE[2][3]=param[4][3][2] #GC
        DE[3][3]=param[4][3][3] #0.1*BE #GG
                
        
        #Initiation Energy
        #At each end:
        E_init=1.005*4184/N_avo
        
        #Z=1
        
        #Initiation of the partition function factorisation- i.e. the first nucleotide position of the sequence
        #Defining a vector to aid the calculation of Z
        #ATCG 0123
        v_Z=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
        re_ord_list=[0,4,8,12,1,5,9,13,2,6,10,14,3,7,11,15]
        #re_ord_list_A=re_ord_list[0:4]
        #re_ord_list_T=re_ord_list[4:8]
        #re_ord_list_C=re_ord_list[8:12]
        #re_ord_list_G=re_ord_list[12:16]
        for i in range(0,16):
            v_Z[i]=p[int(math.floor(i/4))]*np.exp(-self.__beta*(DE[int(math.floor(i/4))][ ts_n[0] ] +2*E_init) )
        
        for pos in range(1,ts_len):
            v_Z_temp=[]
            #v_Z_temp=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
            for i in range(0,len(v_Z)):
                #v_Z_temp[i]=v_Z[i] * cm[int(math.floor(i/4))][i%4] * np.exp(-beta*DE[][])
                v_Z_temp.append(v_Z[i] * cm[int(math.floor(i/4))][int(i%4)] * np.exp(-self.__beta*DE[int(i%4)][ ts_n[pos] ]))
            v_Z=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
            v_Z_temp2=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
            for i in range(0,len(v_Z_temp)):
                v_Z_temp2[i]=v_Z_temp[re_ord_list[i]]
            v_Z[0]=v_Z[1]=v_Z[2]=v_Z[3]=sum(v_Z_temp2[0:4])
            v_Z[4]=v_Z[5]=v_Z[6]=v_Z[7]=sum(v_Z_temp2[4:8])
            v_Z[8]=v_Z[9]=v_Z[10]=v_Z[11]=sum(v_Z_temp2[8:12])
            v_Z[12]=v_Z[13]=v_Z[14]=v_Z[15]=sum(v_Z_temp2[12:16])
            
        #Division by four because the last summation should not take place
        self.__Z=sum(v_Z)/4
            
        G_compl=2*E_init
        for i in range(0,ts_len):
            G_compl += DE[cs_n[i]][ts_n[i]]
    
            
        self.__q_compl=p_cs * np.exp(-self.__beta*G_compl)
        self.__q_exp=np.exp(-self.__beta*G_compl)
        self.__S=1/((self.__Z/self.__q_compl)-1) 
        self.__S_real=np.exp(-self.__beta*G_compl)/self.__N_G/self.__Z
        
        #print(self.__S)
        
        
        
        
"""     
class NPF:
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
                    
    
    #'Algorithm L' from website: http://blog.bjrn.se/2008/04/lexicographic-permutations-using.html
    #This algorithm is used to generate distinct permutations from a list with repeated elements
    #For a list with unique elements, it is quite easy to use functions in the itertools library, but when the list contains repeated elements it seems more difficult to implement an algorithm
    def next_permutation(seq, pred=cmp):

    
        def reverse(seq, start, end):
            # seq = seq[:start] + reversed(seq[start:end]) + \
            #       seq[end:]
            end -= 1
            if end <= start:
                return
            while True:
                seq[start], seq[end] = seq[end], seq[start]
                if start == end or start+1 == end:
                    return
                start += 1
                end -= 1
        
        if not seq:
            raise StopIteration
    
        try:
            seq[0]
        except TypeError:
            raise TypeError("seq must allow random access.")
    
        first = 0
        last = len(seq)
        seq = seq[:]
    
        # Yield input sequence as the STL version is often
        # used inside do {} while.
        yield seq
        
        if last == 1:
            raise StopIteration
    
        while True:
            next = last - 1
    
            while True:
                # Step 1.
                next1 = next
                next -= 1
                
                if pred(seq[next], seq[next1]) < 0:
                    # Step 2.
                    mid = last - 1
                    while not (pred(seq[next], seq[mid]) < 0):
                        mid -= 1
                    seq[next], seq[mid] = seq[mid], seq[next]
                    
                    # Step 3.
                    reverse(seq, next1, last)
    
                    # Change to yield references to get rid of
                    # (at worst) |seq|! copy operations.
                    yield seq[:]
                    break
                if next == first:
                    raise StopIteration
        raise StopIteration
        
    def __init__(self,param,targeting_seq):
        

        k_B=param[0]
        self.__T=param[1][0] #(High temp. approximation?)
        self.__beta=1/k_B/self.__T
        self.__BE=param[1][1] #epsilon- given by hydrogen bond interactions between complementary NTs
        self.__N_G=param[1][2]
        
        self.__p_A=param[2][0]
        self.__p_T=param[2][1]
        self.__p_C=param[2][2]
        self.__p_G=param[2][3]

        p=[0,0,0,0]
        p[0]=self.__p_A
        p[1]=self.__p_T
        p[2]=self.__p_C
        p[3]=self.__p_G
        
        
        n_A=0
        n_T=0
        n_C=0
        n_G=0
        
        ts=list(targeting_seq)
        #ts=["T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T"]
        ts_len=len(ts) #Length of the targeting sequence,

        
                      

        
        
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
                
              
        #Incorporating non-degenerate defect energies
        #Expressing this information in a matrix
        #ATCG 0123
        
        DE=[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
        
        DE[0][0]=param[4][0][0] #0.1*BE #AA
        
        DE[1][0]=DE[0][1]=param[4][1][0] #TA
        DE[1][1]=param[4][1][1] #0.1*BE #TT
        
        DE[2][0]=DE[0][2]=param[4][2][0] #CA
        DE[2][1]=DE[1][2]=param[4][2][1] #CT
        DE[2][2]=param[4][2][2] #0.1*BE #CC
        
        DE[3][0]=DE[0][3]=param[4][3][0] #GA
        DE[3][1]=DE[1][3]=param[4][3][1] #GT
        DE[3][2]=DE[2][3]=param[4][3][2] #GC
        DE[3][3]=param[4][3][3] #0.1*BE #GG
                
        
        
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
        
        cm[0][0]=param[3][0][0] #AA
        cm[0][1]=param[3][0][1] #AT
        cm[0][2]=param[3][0][2] #AC
        cm[0][3]=param[3][0][3] #AG
        
        cm[1][0]=param[3][1][0] #TA
        cm[1][1]=param[3][1][1] #TT
        cm[1][2]=param[3][1][2] #TC
        cm[1][3]=param[3][1][3] #TG
        
        cm[2][0]=param[3][2][0] #CA
        cm[2][1]=param[3][2][1] #CT
        cm[2][2]=param[3][2][2] #CC
        cm[2][3]=param[3][2][3] #CG
        
        cm[3][0]=param[3][3][0] #GA
        cm[3][1]=param[3][3][1] #GT
        cm[3][2]=param[3][3][2] #GC
        cm[3][3]=param[3][3][3] #GG
        
        
        
        if cs[0]=='A':
            p_cs=self.__p_A
        if cs[0]=='T':
            p_cs=self.__p_T
        if cs[0]=='C':
            p_cs=self.__p_C
        if cs[0]=='G':
            p_cs=self.__p_G
        
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

        
        Z=0
        Z_comp_list=[]
        for i in range(0,ts_len+1):
            Z_comp_list.append(0)
        
        
        for n_ncp in range(0,ts_len+1):
            print(Z)
            n_cp=ts_len-n_ncp
            for m_A_A in range(0,min(n_ncp,n_A)+1):
                #print m_A_A
                for m_T_A in range(max(0,n_cp-n_T-n_C-n_G),min(n_A-m_A_A,n_cp)+1):
                    for m_C_A in range(0,min(n_ncp-m_A_A,n_A-m_A_A-m_T_A)+1):
                        m_G_A = n_A-m_A_A-m_T_A-m_C_A
                        n_ncp_A=m_A_A+m_C_A+m_G_A
        
                        
                        for m_A_T in range(max(0,n_cp-m_T_A-n_C-n_G),min(n_cp-m_T_A,n_T)+1):
                            for m_T_T in range(0,min(n_ncp-n_ncp_A,n_T-m_A_T)+1):
                                for m_C_T in range(0,min(n_ncp-n_ncp_A-m_T_T,n_T-m_A_T-m_T_T)+1):
                                    m_G_T = n_T-m_A_T-m_T_T-m_C_T
                                    n_ncp_T=m_T_T+m_C_T+m_G_T
                                    
                                    for m_A_C in range(0,min(n_ncp-n_ncp_A-n_ncp_T,n_C)+1):
                                        #print m_A_C
                                        for m_T_C in range(0,min(n_ncp-n_ncp_A-n_ncp_T-m_A_C,n_C-m_A_C)+1):
                                            for m_C_C in range(0,min(n_ncp-n_ncp_A-n_ncp_T-m_A_C-m_T_C,n_C-m_A_C-m_T_C)+1):
                                                m_G_C = n_C-m_A_C-m_T_C-m_C_C
                                                n_ncp_C=m_A_C+m_T_C+m_C_C
                                                
                                                for m_A_G in range(0,min(n_ncp-n_ncp_A-n_ncp_T-n_ncp_C,n_G)+1):
                                                    for m_T_G in range(0,min(n_ncp-n_ncp_A-n_ncp_T-n_ncp_C-m_A_G,n_G-m_A_G)+1):
                                                        m_C_G = n_cp-m_A_T-m_T_A-m_G_C
                                                        m_G_G = n_G-m_A_G-m_T_G-m_C_G
                                                        n_ncp_G=m_A_G+m_T_G+m_G_G
                                                        
                                                        if n_ncp_A+n_ncp_T+n_ncp_C+n_ncp_G != n_ncp:
                                                            raise Exception()
                                                        
                                                        m_A=m_A_A+m_A_T+m_A_C+m_A_G #equivalent to i
                                                        m_T=m_T_A+m_T_T+m_T_C+m_T_G #equivalent to j
                                                        m_C=m_C_A+m_C_T+m_C_C+m_C_G #equivalent to k
                                                        m_G=m_G_A+m_G_T+m_G_C+m_G_G #equivalent to l
                                                    
                                                        seq_A=[0]*m_A_A+[1]*m_T_A+[2]*m_C_A+[3]*m_G_A
                                                        seq_T=[0]*m_A_T+[1]*m_T_T+[2]*m_C_T+[3]*m_G_T
                                                        seq_C=[0]*m_A_C+[1]*m_T_C+[2]*m_C_C+[3]*m_G_C
                                                        seq_G=[0]*m_A_G+[1]*m_T_G+[2]*m_C_G+[3]*m_G_G
                                                        
                                                        
                                                        
                                                        for perm_A in next_permutation(seq_A):
                                                            for perm_T in next_permutation(seq_T):
                                                                for perm_C in next_permutation(seq_C):
                                                                    for perm_G in next_permutation(seq_G):
                                                                        
                                                                        pos_nt={}
                                                                        
                                                                        for a in range(0,n_A):
                                                                            pos_nt[A_pos[a]]=perm_A[a]
                                                                            
                                                                        for b in range(0,n_T):
                                                                            pos_nt[T_pos[b]]=perm_T[b]
                                                                        
                                                                        for c in range(0,n_C):
                                                                            pos_nt[C_pos[c]]=perm_C[c]
                                                                        
                                                                        for d in range(0,n_G):
                                                                            pos_nt[G_pos[d]]=perm_G[d]
                                                                        
                                                                        P=p[ pos_nt[0] ]
                                                                        
                                                                        for x in range(0,ts_len-1):
                                                                            P*=cm[pos_nt[x]][pos_nt[x+1]]
                                                                            
                                                                        E_tot=m_A_A*DE[0][0] + (m_T_A+m_A_T)*DE[0][1] + m_T_T*DE[1][1] + (m_C_A+m_A_C)*DE[2][0] + (m_C_T+m_T_C)*DE[2][1] + m_C_C*DE[2][2] + (m_G_A+m_A_G)*DE[3][0] + (m_G_T+m_T_G)*DE[3][1] + (m_G_C+m_C_G)*DE[3][2] + m_G_G*DE[3][3]
                                                                            
                                                                        Z+=P*np.exp(-self.__beta*E_tot)
                                                                        
                                                            
                                                        
                                                        
        q_compl=p_cs * np.exp(-self.__beta*ts_len*self.__BE)
        
        S=1/(Z/q_compl-1)                       
                                                        
                                                        
                                                        

            
        
        for i in range(0,ts_len+1):
            
            #print Z_comp_list
            for j in range(0,ts_len-i+1):
                for k in range(0,ts_len-i-j+1):
                    
                    l=ts_len-i-j-k
                    print i,l,Z
                    
                    ijkl_seq=[0]*i+[1]*j+[2]*k+[3]*l
                    
                    for perm in next_permutation(ijkl_seq):
                        
                        P=p[ perm[0] ]
                        n_cp=0
                        
                        for a in range(0,ts_len-1):
                            P*=cm[perm[a]][perm[a+1]]
                            
                            if perm[a] == cs[a]:
                                n_cp+=1
                                
                        n_ncp=ts_len-n_cp
                        
                        Z+=P*np.exp(-beta*n_cp*BE)
        
                        #if (change of Z is small):
                            #raise Exception
                            
                        
        q_compl=p_cs * np.exp(-beta*ts_len*BE)
        
        S=1/(Z/q_compl-1)
        """