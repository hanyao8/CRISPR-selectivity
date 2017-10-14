#iteration

import numpy as np

k_B=1.38e-23
T=310 #(High temp. approximation?)
beta=1/k_B/T
BE=-0.038*(1.6e-19) #0.1eV 
"""
The value of the binding requires revision. Recall that the BE value of 0.038eV only correponds to a particular sequence.
A better value for the average for a random sequence can be found in Santa Lucia et al 1998
"""

ts_string='TGATTTAGAACCTGAAAGCA'#'ACTAAATCTTGGACTTTCGT'
ts_array=["T","G","A","T","T","T","A","G","A","A","C","C","T","G","A","A","A","G","C","A"]#["A","C","T","A","A","A","T","C","T","T","G","G","A","C","T","T","T","C","G","T"]

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

pair_2_G={'AA':DE[0][0],'AT':DE[0][1],'AC':DE[0][2],'AG':DE[0][3],'TA':DE[1][0],'TT':DE[1][1],'TC':DE[1][2],'TG':DE[1][3],'CA':DE[2][0],'CT':DE[2][1],'CC':DE[2][2],'CG':DE[2][3],'GA':DE[3][0],'GT':DE[3][1],'GC':DE[3][2],'GG':DE[3][3]} 

cs_string=''
cs_array=[]
p_cs=1

ts_len=len(ts_array)

for a in range(0,ts_len):
    if ts_array[a]=='A':
        cs_array.append('T')
        cs_string+='T'
    if ts_array[a]=='T':
        cs_array.append('A')
        cs_string+='A'
    if ts_array[a]=='C':
        cs_array.append('G')
        cs_string+='G'
    if ts_array[a]=='G':
        cs_array.append('C')
        cs_string+='C'
"""
snt=12888830
data_len=30
N_binds=(data_len-ts_len)*100
"""
#np.random.choice(elements, 10, p=probabilities)


class iter_object:
    
    def __init__ (self,starting_nt,data_len,Po_factor):
        #self.__N_iters=N_iters
        self.__snt=starting_nt
        self.__data_len=data_len 
        
        self.__ts=ts_string
        self.__cs=cs_string
        self.__ts_len=len(self.__ts)
        self.__N_binds=(self.__data_len-self.__ts_len)*Po_factor
        self.__S_list=[0.0]
        self.__binding_events=0 #Z=0
        self.__compl_events=0
        
        
        with open(r'C:\Users\Choon\Desktop\CRISPR\Code\nt_data\chromosome_20_complete.txt', 'r') as myfile:
            data=myfile.read().replace('\n', '')
        self.__data=data[self.__snt:(self.__snt+self.__data_len)]        
        
        
        i=1
        while i < self.__N_binds+1:  

            r=int(np.random.uniform(0,self.__data_len-self.__ts_len))
            tgds=data[r:r+self.__ts_len]
        
            G_20=20*BE
            G_19=19*BE
            A=np.exp(-0.5*(G_20+G_19)*beta)
                
            G_binding=0
            for j in range(0,self.__ts_len):
                G_binding+=pair_2_G[tgds[j]+ts_array[j]]
                    
                
            p=1/(1+A*np.exp(G_binding*beta))
            self.__binding_events+=int((np.random.choice([0,1],1,[1-p,p]))[0]) #np.random.choice(elements, 10, p=probabilities)
            if tgds==self.__cs:
                self.__compl_events+=1
                    

                    
            if (self.__compl_events==0) or (self.__binding_events==self.__compl_events):
                self.__S_list.append(0.0)
                    
            else:
                self.__S_list.append(1/((float(self.__binding_events)/float(self.__compl_events))-1))
               
            print (tgds)
            print (i)
                                       
            i+=1
                
