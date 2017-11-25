#Translation binding test of model 3
#Algorithm will be similar to the parameter extraction algorithm



import numpy as np
import matplotlib.pyplot as plt


def F_cs(s): #Find complementary sequence
    s=list(s)
    comp={'A':'T','T':'A','C':'G','G':'C'}
    for i in range(0,len(s)):
        s[i]=comp[s[i]]
    s=''.join(s)
    return s
    

nt_2_n={'A':0,'T':1,'C':2,'G':3}



class iter:
    def __repr__(self):
        return "A simulation which corresponds to a single walk of a targeting\
        sequence along a chromosome or any strand of DNA"
        #ts:%s, \
        #first position: %g" #% (ts,fp)
        
    def __init__(self,param,data,site0,PAM,ts_len):
        k_B=param[0][0]
        self.__T=param[1][0] #(High temp. approximation?)
        self.__beta=1/k_B/self.__T
        self.__BE=param[1][1] #epsilon- given by hydrogen bond interactions between complementary NTs
        self.__N_G=param[1][2]
        

             
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
        
        
        self.__data=data
        self.__site0=site0
        self.__PAM=PAM        
        self.__ts_len=ts_len
        
        self.__data_len=len(self.__data)

  
        self.__j_list=[]
        self.__Z_nc_list=[]
        self.__Z_j_list=[]
        self.__S_j_list=[]


        cs=self.__data[self.__site0-(1+len(self.__PAM))+1-self.__ts_len:self.__site0-(1+len(self.__PAM))+1]
        ts=F_cs(cs)
        
        print(ts)
        print(cs)
        
        self.__Z=0
        self.__q_comp=0
        self.__comp_count=0
        self.__S=0
    
        #for j in range(site1,site1+bindings):
        for j in range(0,self.__data_len-self.__ts_len+1):
            tgds=data[j:j+self.__ts_len]    
            if tgds.count('N')==0:
                G_binding=0
                for k in range(0,self.__ts_len-1):
                    G_binding+=G_dpx[ nt_2_n[tgds[k]] ][ nt_2_n[tgds[k+1]] ][ nt_2_n[ts[k]] ][ nt_2_n[ts[k+1]] ]
                    
                G_binding+= (G_init[ nt_2_n[tgds[0]] ][ nt_2_n[ts[0]] ] + G_init[ nt_2_n[tgds[self.__ts_len-1]] ][ nt_2_n[ts[self.__ts_len-1]] ])
                self.__Z+=np.exp(-self.__beta*G_binding)
                self.__Z_nc=self.__Z-self.__q_comp
                S_j=self.__q_comp/self.__Z_nc
                
                if tgds==cs:
                    self.__q_comp=np.exp(-self.__beta*G_binding)
                    self.__comp_count+=1
                    #print("comp site! count=%d" %(self.__comp_count))
            
                if (j%5000000)==0:
                    print("running... %d" % (j))
                
                if (j%50000)==(50000-1):
                    self.__j_list.append(j)
                    self.__Z_nc_list.append(self.__Z_nc)
                    self.__Z_j_list.append(self.__Z_nc/j)
                    self.__S_j_list.append(S_j)
            
                if (j%10000000)==(10000000-1):
                    #print(tgds)
                    print(ts,S_j,self.__q_comp,self.__Z_nc,self.__N_G,j)
            
        self.__S=1/((self.__Z/self.__q_comp)-1)

class random_genome1:
    
    def __repr__(self):
        return "A simulation which corresponds to a single walk of a targeting\
        sequence along a random of DNA defined by statistical parameters input"
        #ts:%s, \
        #first position: %g" #% (ts,fp)
        
    def __init__(self,param,ts):
        k_B=param[0][0]
        T=param[1][0] #(High temp. approximation?)
        self.__beta=1/k_B/T
        self.__BE=param[1][1] #epsilon- given by hydrogen bond interactions between complementary NTs
        self.__N_G=param[1][2]
        #N_avo=param[0][1]
        
             
        
      
        #Incorporating non-degenerate defect energies
        #Expressing this information in a matrix
        #ATCG 0123
        
        DE=[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
    
        
        DE[1][0]=DE[0][1]=self.__BE #TA

        DE[3][2]=DE[2][3]=self.__BE #GC

    
        

        
        

        ts_len=len(ts)
        
        ts_n=[nt_2_n[ts[i]] for i in range(0,ts_len)]
        


  
        self.__j_list=[]
        self.__S_j_list=[]
        self.__Z_j_list=[]


        cs=F_cs(ts)
        
        print(ts)
        print(cs)
        
        #self.__Z=0
        #self.__q_comp=0
        self.__comp_count=0
        #self.__S=0
        self.__Z_nc=0.0
        
        G_binding=0
        for k in range(0,ts_len):
            G_binding+=DE[ nt_2_n[cs[k]] ][ nt_2_n[ts[k]] ]
                    
        self.__q_comp=self.__Z=np.exp(-self.__beta*G_binding)
        #print(self.__q_comp)
        
        #for j in range(0,self.__N_G-ts_len+1):
        for j in range(0,1000000):
            
            if (j%100000)==0:
                mini_genome=  np.random.choice([0,1,2,3],size=100000+ts_len-1,p=param[2]) 

            
            x=j%100000
            
            tgds=mini_genome[x:x+ts_len]
            
            G_binding=0
            for k in range(0,ts_len):
                G_binding+=DE[ tgds[k] ][ ts_n[k] ]
                    

            #print(np.exp(-self.__beta*G_binding))
            self.__Z_nc+=np.exp(-self.__beta*G_binding)
            #print(self.__Z_nc)
            
            if (j%50000)==(50000-1):
                self.__j_list.append(j)
                self.__Z_j_list.append(self.__Z_nc/j)
                S_j=self.__q_comp/(self.__Z_nc*self.__N_G/j)
                self.__S_j_list.append(S_j)
                
            if (j%100000)==(100000-1):
                #print(tgds)
                print(S_j,self.__q_comp,self.__Z_nc,self.__N_G,j)
            

            
        self.__S=self.__q_comp/self.__Z_nc


class random_genome2:
    
    def __repr__(self):
        return "A simulation which corresponds to a single walk of a targeting\
        sequence along a random of DNA defined by statistical parameters input"
        #ts:%s, \
        #first position: %g" #% (ts,fp)
        
    def __init__(self,param,ts):
        k_B=param[0][0]
        T=param[1][0] #(High temp. approximation?)
        self.__beta=1/k_B/T
        self.__BE=param[1][1] #epsilon- given by hydrogen bond interactions between complementary NTs
        self.__N_G=param[1][2]
        N_avo=param[0][1]
        
             
        
      
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

        
        

        ts_len=len(ts)
        
        ts_n=[nt_2_n[ts[i]] for i in range(0,ts_len)]
        


  
        self.__j_list=[]
        self.__S_j_list=[]
        self.__Z_j_list=[]

        cs=F_cs(ts)
        
        print(ts)
        print(cs)
        
        #self.__Z=0
        #self.__q_comp=0
        self.__comp_count=0
        #self.__S=0
        self.__Z_nc=0.0
        
        G_binding=0
        for k in range(0,ts_len):
            G_binding+=DE[ nt_2_n[cs[k]] ][ nt_2_n[ts[k]] ]
                    
        G_binding+= 2*E_init
        self.__q_comp=self.__Z=np.exp(-self.__beta*G_binding)
        #print(self.__q_comp)
        
        for j in range(0,self.__N_G-ts_len+1):
        #for j in range(0,1000000):
        #for j in range(0,10):
            
            if (j%100000)==0:
                mini_genome= [ np.random.choice([0,1,2,3],p=param[2]) ] 
                for k in range(0,100000+ts_len-1):
                #for k in range(0,ts_len-1):
                    #tgds=[np.random.choice(["A","T","C","G"],p=param[2])]
                    #tgds.append( np.random.choice( ["A","T","C","G"],p=np.array(param[3][0]) ) )
                    mini_genome.append( np.random.choice( [0,1,2,3],p=param[3][mini_genome[k]] ) ) 
            
            x=j%100000
            
            tgds=mini_genome[x:x+ts_len]
            
            G_binding=0
            for k in range(0,ts_len):
                G_binding+=DE[ tgds[k] ][ ts_n[k] ]
                    
            G_binding+= 2*E_init
            #print(np.exp(-self.__beta*G_binding))
            self.__Z_nc+=np.exp(-self.__beta*G_binding)
            #print(self.__Z_nc)
            
            if (j%50000)==(50000-1):
                self.__j_list.append(j)
                self.__Z_j_list.append(self.__Z_nc/j)
                S_j=self.__q_comp/(self.__Z_nc*self.__N_G/j)
                self.__S_j_list.append(S_j)
                
            if (j%100000)==(100000-1):
                #print(tgds)
                print(S_j,self.__q_comp,self.__Z_nc,self.__N_G,j)

            #print(tgds,self.__q_comp,np.exp(-self.__beta*G_binding),self.__Z_nc,j)
            
            """
            if (j%5000)==(5000-1):
                print("running... %d" % (j))
                print(self.__Z)
                S_j=1/((self.__N_G/j*self.__Z/self.__q_comp)-1)
                print(S_j)
            
            """
            """
            if (j%1000)==(1000-1):
                #print("running... %d" % (j))
                self.__j_list.append(j)
                S_j=self.__q_comp/(self.__Z_nc*self.__N_G/j)
                #S_j=1/((self.__N_G/j*self.__Z/self.__q_comp)-1)
                self.__S_j_list.append(S_j)
                print(self.__Z)
                print(S_j)
                #plt.scatter(j,self.__Z,s=0.2)
                #plt.scatter(j,S_j,s=0.2)
                #plt.show()
            """
            
        self.__S=self.__q_comp/self.__Z_nc

class random_genome3:
    
    def __repr__(self):
        return "A simulation which corresponds to a single walk of a targeting\
        sequence along a random of DNA defined by statistical parameters input"
        #ts:%s, \
        #first position: %g" #% (ts,fp)
        
    def __init__(self,param,ts):
        k_B=param[0][0]
        self.__T=param[1][0] #(High temp. approximation?)
        self.__beta=1/k_B/self.__T
        self.__BE=param[1][1] #epsilon- given by hydrogen bond interactions between complementary NTs
        self.__N_G=param[1][2]
        
             
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
        
        

        ts_len=len(ts)
        
        ts_n=[nt_2_n[ts[i]] for i in range(0,ts_len)]
  
        self.__j_list=[]
        self.__S_j_list=[]
        self.__Z_j_list=[]

        cs=F_cs(ts)
        
        print(ts)
        print(cs)
        
        #self.__Z=0
        #self.__q_comp=0
        self.__comp_count=0
        #self.__S=0
        self.__Z_nc=0.0
        
        G_binding=0
        for k in range(0,ts_len-1):
            G_binding+=G_dpx[ nt_2_n[cs[k]] ][ nt_2_n[cs[k+1]] ][ nt_2_n[ts[k]] ][ nt_2_n[ts[k+1]] ]
                    
        G_binding+= (G_init[ nt_2_n[cs[0]] ][ nt_2_n[ts[0]] ] + G_init[ nt_2_n[cs[ts_len-1]] ][ nt_2_n[ts[ts_len-1]] ])
        self.__q_comp=self.__Z=np.exp(-self.__beta*G_binding)
        #print(self.__q_comp)
        
        for j in range(0,self.__N_G-ts_len+1):
        #for j in range(0,1000000):
            
            if (j%100000)==0:
                mini_genome= [ np.random.choice([0,1,2,3],p=param[2]) ] 
                for k in range(0,100000+ts_len-1):
                #for k in range(0,ts_len-1):
                    #tgds=[np.random.choice(["A","T","C","G"],p=param[2])]
                    #tgds.append( np.random.choice( ["A","T","C","G"],p=np.array(param[3][0]) ) )
                    mini_genome.append( np.random.choice( [0,1,2,3],p=param[3][mini_genome[k]] ) ) 
            
            x=j%100000
            
            tgds=mini_genome[x:x+ts_len]
            
            G_binding=0
            for k in range(0,ts_len-1):
                G_binding+=G_dpx[ tgds[k] ][ tgds[k+1] ][ ts_n[k] ][ ts_n[k+1] ]
                    
            G_binding+= (G_init[ tgds[0] ][ ts_n[0] ] + G_init[ tgds[ts_len-1] ][ ts_n[ts_len-1] ])
            #print(np.exp(-self.__beta*G_binding))
            self.__Z_nc+=np.exp(-self.__beta*G_binding)
            #print(self.__Z_nc)
            
            if (j%50000)==(50000-1):
                self.__j_list.append(j)
                self.__Z_j_list.append(self.__Z_nc/j)
                S_j=self.__q_comp/(self.__Z_nc*self.__N_G/j)
                self.__S_j_list.append(S_j)
                
            if (j%200000)==(200000-1):
                #print(tgds)
                print(S_j,self.__q_comp,self.__Z_nc,self.__N_G,j)
            
            """
            if (j%5000)==(5000-1):
                print("running... %d" % (j))
                print(self.__Z)
                S_j=1/((self.__N_G/j*self.__Z/self.__q_comp)-1)
                print(S_j)
            
            """
            """
            if (j%1000)==(1000-1):
                #print("running... %d" % (j))
                self.__j_list.append(j)
                S_j=self.__q_comp/(self.__Z_nc*self.__N_G/j)
                #S_j=1/((self.__N_G/j*self.__Z/self.__q_comp)-1)
                self.__S_j_list.append(S_j)
                print(self.__Z)
                print(S_j)
                #plt.scatter(j,self.__Z,s=0.2)
                #plt.scatter(j,S_j,s=0.2)
                #plt.show()
            """
            
        self.__S=self.__q_comp/self.__Z_nc

"""                 
    x=np.array(j_list)
    y=np.array(Z_nc_list)
    plt.plot(x,y)
    plt.title("nc Component of Z")
    plt.ylabel("Z_nc")
    plt.xlabel("Binding Events")
    plt.show()

    
        

#site1=min(sites)
site1=60045
#bindings=64444000
final=64444000
"""
    