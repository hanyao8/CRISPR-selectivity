#Translation binding test of model 3
#Algorithm will be similar to the parameter extraction algorithm



import numpy as np
#import matplotlib.pyplot as plt


def F_cs(s):
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
                
                if tgds==cs:
                    self.__q_comp=np.exp(-self.__beta*G_binding)
                    self.__comp_count+=1
                    print("comp site! count=%d" %(self.__comp_count))
            
                if (j%5000000)==0:
                    print("running... %d" % (j))
                
                if (j%10000)==0:
                    self.__j_list.append(j)
                    self.__Z_nc_list.append(self.__Z-self.__q_comp)
            
        self.__S=1/((self.__Z/self.__q_comp)-1)


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
    