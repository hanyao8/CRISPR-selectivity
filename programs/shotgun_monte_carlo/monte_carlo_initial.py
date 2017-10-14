#Investigating the convergence of a Monte Carlo simulation to mean field models
#Initial attempt at writing the Monte Carlo program for moddel 1:
#Using OOP

from iteration import iter_object
import matplotlib.pyplot as plt
import numpy as np


ts_string='TGATTTAGAACCTGAAAGCA' #'ACTAAATCTTGGACTTTCGT'
ts_array=["T","G","A","T","T","T","A","G","A","A","C","C","T","G","A","A","A","G","C","A"] #["A","C","T","A","A","A","T","C","T","T","G","G","A","C","T","T","T","C","G","T"]
#complementary sequence: 'TGATTTAGAACCTGAAAGCA'
ts_len=len(ts_string)

N_iters=1 #Number of iterations
data_len=30
starting_nt=12888830
N_binds=(data_len-ts_len)*3

cs_string=''
cs_array=[]
p_cs=1

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
    if ts[a]!='A' and  ts[a]!='T' and ts[a]!='C' and ts[a]!='G':
        raise Exception("Invalid Nucleotide")
"""
#Extraction of statistical parameters from genome nucleotide sequence data
with open(r'C:\Users\Choon\Desktop\CRISPR\Code\nt_data\chromosome_20_complete.txt', 'r') as myfile:
    #data=myfile.read()
    data=myfile.read().replace('\n', '')
    data=data[starting_nt:starting_nt+data_len]

len_data=len(data)

N_A_G=data[0:len_data-1].count('A')
N_T_G=data[0:len_data-1].count('T')
N_C_G=data[0:len_data-1].count('C')
N_G_G=data[0:len_data-1].count('G')

p_A=N_A_G/(N_A_G+N_T_G+N_C_G+N_G_G)
p_T=N_T_G/(N_A_G+N_T_G+N_C_G+N_G_G)
p_C=N_C_G/(N_A_G+N_T_G+N_C_G+N_G_G)
p_G=N_G_G/(N_A_G+N_T_G+N_C_G+N_G_G)

p_N=(len_data-(N_A_G+N_T_G+N_C_G+N_G_G))/len_data

p=[p_A,p_T,p_C,p_G]

for i in range(0,N_iters):
    it_instance=iter_object(starting_nt,data_len,N_binds)
    x=np.array(np.linspace(0,N_binds,N_binds+1))
    y=it_instance._iter_object__S_list
    #plt.plot(np.linspace(0,5,6),[1,2,3,4,5,6])
    plt.plot(x,y)
    
    
plt.show()

