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




#np.random.choice(elements, 10, p=probabilities)


class iter_object:
    
    def __init__ (self,starting_nt,data_len,N_binds,ts,cs):
        #self.__N_iters=N_iters
        self.__snt=starting_nt
        self.__data_len=data_len 
        
        self.__ts=ts
        self.__cs=cs
        self.__ts_len=len(self.__ts)
        self.__N_binds=N_binds
        self.__S_list=[0]
        self.__binding_events=0 #Z=0
        self.__compl_events=0
        
        
        with open(r'C:\Users\Choon\Desktop\CRISPR\Code\nt_data\chromosome_20_complete.txt', 'r') as myfile:
            data=myfile.read().replace('\n', '')
        self.__data=data[self.__snt:self.__snt+self.__data_len]        
        
        
        i=1
        while i < self.__N_binds+1:  

            r=int(np.random.uniform(0,self.__data_len-self.__ts_len))
            tgds=data[r:r+self.__ts_len]
            
            if tgds.count('N')==0:
                G_20=20*BE
                G_19=19*BE
                A=np.exp(-0.5*(G_20+G_19)*beta)
                
                G_binding=0
                for j in range(0,self.__ts_len):
                    G_binding+=pair_2_G[tgds[j]+ts[j]]
                    
                """
                for j in range(0,ts_len-1):
                    G_binding+=duplex_2_G[]
                """
                
                p=1/(1+A*np.exp(G_binding*beta))
                self.__binding_events+=int((np.random.choice([0,1],1,[1-p,p]))[0]) #np.random.choice(elements, 10, p=probabilities)
                if tgds==self.__cs:
                    self.__compl_events+=1
                    
                """
                if self.__compl_events>0:
                    self.__S_list.append(1/((self.__binding_events/self.__compl_events)-1))

                else:
                    self.__S_list.append(0) 
                """
                    
                if (self.__compl_events==0) or (self.__binding_events==self.__compl_events):
                    self.__S_list.append(0.0)
                    
                else:
                    self.__S_list.append(1/((float(self.__binding_events)/float(self.__compl_events))-1))
               
                print (i)
                                       
                i+=1
                
