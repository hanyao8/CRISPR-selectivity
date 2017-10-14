#crispr model 3: obtains binding energy given by duplex free energies, obtained from Santa Lucia rules
#Incorporates nucleotide correlation
#Utilises partition function factorisation

import numpy as np
import math

class PF:
    def __init__(self,param,targeting_seq):

        k_B=param[0][0]
        #N_avo=param[0][1]
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
        ts_len=len(ts) #Length of the targeting sequence,



        
        #Complementary sequence and its probability of occuring in the genome to be computed
        cs=[]
        p_cs=1
        
        ts_n=[] #A numerical targeting sequence array can be constructed, with A=0, T=1, C=2, G=3
        cs=[]
        cs_n=[]
        
        #An array of the positions of the ATCG nucleotides in the targeting sequence
        #The position should be defined with respect to 3' and 5' ends of the strand
        A_pos=[]
        T_pos=[]
        C_pos=[]
        G_pos=[]
        
        #Loop to calculate the variables and sequences initialised above
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
                
        
        #Section to initialise nucleotide correlation matrix
        
              
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
        

        
        #Calculation of the probability of the complementary sequence given the information from the correlation matrix
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
        G_dpx[0][0][1][1]=G_dpx[1][1][0][0]=param[4][0][0][1][1] #AA/TT      
        G_dpx[0][1][1][0]=param[4][0][1][1][0] #AT/TA
        G_dpx[1][0][0][1]=param[4][1][0][0][1] #TA/AT
        G_dpx[2][0][3][1]=G_dpx[1][3][0][2]=param[4][2][0][3][1] #CA/GT
        G_dpx[3][1][2][0]=G_dpx[0][2][1][3]=param[4][3][1][2][0] #GT/CA
        G_dpx[2][1][3][0]=G_dpx[0][3][1][2]=param[4][2][1][3][0] #CT/GA
        G_dpx[3][0][2][1]=G_dpx[1][2][0][3]=param[4][3][0][2][1] #GA/CT
        G_dpx[2][3][3][2]=param[4][2][3][3][2] #CG/GC
        G_dpx[3][2][2][3]=param[4][3][2][2][3] #GC/CG
        G_dpx[3][3][2][2]=G_dpx[2][2][3][3]=param[4][3][3][2][2] #GG/CC
        
        #single defects
        G_dpx[0][0][1][0]=G_dpx[0][1][0][0]=param[4][0][0][1][0]
        G_dpx[0][0][1][2]=G_dpx[2][1][0][0]=param[4][0][0][1][2]
        G_dpx[0][0][1][3]=G_dpx[3][1][0][0]=param[4][0][0][1][3]
        G_dpx[0][1][1][1]=G_dpx[1][1][1][0]=param[4][0][1][1][1]
        G_dpx[0][1][1][2]=G_dpx[2][1][1][0]=param[4][0][1][1][2]
        G_dpx[0][1][1][3]=G_dpx[3][1][1][0]=param[4][0][1][1][3]
        G_dpx[0][2][1][0]=G_dpx[0][1][2][0]=param[4][0][2][1][0]
        G_dpx[0][2][1][1]=G_dpx[1][1][2][0]=param[4][0][2][1][1]
        G_dpx[0][2][1][2]=G_dpx[2][1][2][0]=param[4][0][2][1][2]
        G_dpx[0][3][1][0]=G_dpx[0][1][3][0]=param[4][0][3][1][0]
        G_dpx[0][3][1][1]=G_dpx[1][1][3][0]=param[4][0][3][1][1]
        G_dpx[0][3][1][3]=G_dpx[3][1][3][0]=param[4][0][3][1][3]
        
        G_dpx[1][0][0][0]=G_dpx[0][0][0][1]=param[4][1][0][0][0]
        G_dpx[1][0][0][2]=G_dpx[2][0][0][1]=param[4][1][0][0][2]
        G_dpx[1][0][0][3]=G_dpx[3][0][0][1]=param[4][1][0][0][3]
        G_dpx[1][1][0][1]=G_dpx[1][0][1][1]=param[4][1][1][0][1]
        G_dpx[1][1][0][2]=G_dpx[2][0][1][1]=param[4][1][1][0][2]
        G_dpx[1][1][0][3]=G_dpx[3][0][1][1]=param[4][1][1][0][3]
        G_dpx[1][2][0][0]=G_dpx[0][0][2][1]=param[4][1][2][0][0]
        G_dpx[1][2][0][1]=G_dpx[1][0][2][1]=param[4][1][2][0][1]
        G_dpx[1][2][0][2]=G_dpx[2][0][2][1]=param[4][1][2][0][2]
        G_dpx[1][3][0][0]=G_dpx[0][0][3][1]=param[4][1][3][0][0]
        G_dpx[1][3][0][1]=G_dpx[1][0][3][1]=param[4][1][3][0][1]
        G_dpx[1][3][0][3]=G_dpx[3][0][3][1]=param[4][1][3][0][3]
        
        G_dpx[2][0][3][0]=G_dpx[0][3][0][2]=param[4][2][0][3][0]
        G_dpx[2][0][3][2]=G_dpx[2][3][0][2]=param[4][2][0][3][2]
        G_dpx[2][0][3][3]=G_dpx[3][3][0][2]=param[4][2][0][3][3]
        G_dpx[2][1][3][1]=G_dpx[1][3][1][2]=param[4][2][1][3][1]
        G_dpx[2][1][3][2]=G_dpx[2][3][1][2]=param[4][2][1][3][2]
        G_dpx[2][1][3][3]=G_dpx[3][3][1][2]=param[4][2][1][3][3]
        G_dpx[2][2][3][0]=G_dpx[0][3][2][2]=param[4][2][2][3][0]
        G_dpx[2][2][3][1]=G_dpx[1][3][2][2]=param[4][2][2][3][1]
        G_dpx[2][2][3][2]=G_dpx[2][3][2][2]=param[4][2][2][3][2]
        G_dpx[2][3][3][0]=G_dpx[0][3][3][2]=param[4][2][3][3][0]
        G_dpx[2][3][3][1]=G_dpx[1][3][3][2]=param[4][2][3][3][1]
        G_dpx[2][3][3][3]=G_dpx[3][3][3][2]=param[4][2][3][3][3]
        
        G_dpx[3][0][2][0]=G_dpx[0][2][0][3]=param[4][3][0][2][0]
        G_dpx[3][0][2][2]=G_dpx[2][2][0][3]=param[4][3][0][2][2]
        G_dpx[3][0][2][3]=G_dpx[3][2][0][3]=param[4][3][0][2][3]
        G_dpx[3][1][2][1]=G_dpx[1][2][1][3]=param[4][3][1][2][1]
        G_dpx[3][1][2][2]=G_dpx[2][2][1][3]=param[4][3][1][2][2]
        G_dpx[3][1][2][3]=G_dpx[3][2][1][3]=param[4][3][1][2][3]
        G_dpx[3][2][2][0]=G_dpx[0][2][2][3]=param[4][3][2][2][0]
        G_dpx[3][2][2][1]=G_dpx[1][2][2][3]=param[4][3][2][2][1]
        G_dpx[3][2][2][2]=G_dpx[2][2][2][3]=param[4][3][2][2][2]
        G_dpx[3][3][2][0]=G_dpx[0][2][3][3]=param[4][3][3][2][0]
        G_dpx[3][3][2][1]=G_dpx[1][2][3][3]=param[4][3][3][2][1]
        G_dpx[3][3][2][3]=G_dpx[3][2][3][3]=param[4][3][3][2][3]

        #Double defects
        #Using calculated average 'isolated defect energies'
        
        G_avg_AA=param[5][0]
        G_avg_AC=param[5][1]
        G_avg_AG=param[5][2]
        G_avg_TT=param[5][3]
        G_avg_TC=param[5][4]
        G_avg_TG=param[5][5]
        G_avg_CC=param[5][6]
        G_avg_GG=param[5][7]
        
        
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
        G_init[2][3]=G_init[3][2]=param[6][2][3]
        G_init[0][1]=G_init[1][0]=param[6][0][1]
        
        #non-WC Pairs, energy values due for review
        G_init[0][0]=param[6][0][0]
        G_init[0][2]=G_init[2][0]=param[6][0][2]
        G_init[0][3]=G_init[3][0]=param[6][0][3]
        G_init[1][1]=param[6][1][1]
        G_init[1][2]=G_init[2][1]=param[6][1][2]
        G_init[1][3]=G_init[3][1]=param[6][1][3]
        G_init[2][2]=param[6][2][2]
        G_init[3][3]=param[6][3][3]
        
           
         
         
         
           
        
        #Initiation of the partition function factorisation- i.e. the first nucleotide position of the sequence
        #Defining a vector to aid the calculation of Z
        #ATCG 0123
        v_Z=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
        re_ord_list=[0,4,8,12,1,5,9,13,2,6,10,14,3,7,11,15]
        
        for i in range(0,16):
            v_Z[i]=p[int(math.floor(i/4))]*np.exp(-self.__beta*G_init[int(math.floor(i/4))][ts_n[0]])
        
        for pos in range(1,ts_len):
            v_Z_temp=[]
            for i in range(0,len(v_Z)):
                v_Z_temp.append(v_Z[i] * cm[int(math.floor(i/4))][int(i%4)] * np.exp(-self.__beta*G_dpx[int(math.floor(i/4))][int(i%4)][ts_n[pos-1]][ts_n[pos]]))
            v_Z=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
            v_Z_temp2=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
            for i in range(0,len(v_Z_temp)):
                v_Z_temp2[i]=v_Z_temp[re_ord_list[i]]
            v_Z[0]=v_Z[1]=v_Z[2]=v_Z[3]=sum(v_Z_temp2[0:4])
            v_Z[4]=v_Z[5]=v_Z[6]=v_Z[7]=sum(v_Z_temp2[4:8])
            v_Z[8]=v_Z[9]=v_Z[10]=v_Z[11]=sum(v_Z_temp2[8:12])
            v_Z[12]=v_Z[13]=v_Z[14]=v_Z[15]=sum(v_Z_temp2[12:16])
            
        #Collapsing the v_Z array
        v_Z=[v_Z[0]*np.exp(-self.__beta*G_init[0][ts_n[ts_len-1]]),v_Z[4]*np.exp(-self.__beta*G_init[1][ts_n[ts_len-1]]),v_Z[8]*np.exp(-self.__beta*G_init[2][ts_n[ts_len-1]]),v_Z[12]*np.exp(-self.__beta*G_init[3][ts_n[ts_len-1]])]
            
        Z=sum(v_Z)
        
        
        G_compl=G_init[cs_n[0]][ts_n[0]]+G_init[cs_n[ts_len-1]][ts_n[ts_len-1]]
        for i in range(1,ts_len):
            G_compl += G_dpx[cs_n[i-1]][cs_n[i]][ts_n[i-1]][ts_n[i]]
        
        a=np.exp(-self.__beta*G_compl)
        self.__q_compl=p_cs * np.exp(-self.__beta*G_compl)
        self.__S=1/((Z/self.__q_compl)-1) 
        #S_real=S/(N_G-ts_len)/p_cs
        self.__S_real=np.exp(-self.__beta*G_compl)/self.__N_G/Z

        
        #The model requires defect energies to be known. At this moment it gives a somewhat erroneous value for the partition function