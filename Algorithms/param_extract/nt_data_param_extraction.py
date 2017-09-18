#Algorithm to compute parameters of nucleotide frequency and correlation in a genome


#txt_data = np.genfromtxt(r'C:\Users\Choon\Desktop\CRISPR\Code\nt_data\chromosome_1_1_1000000.txt')
#print txt_data

with open(r'C:\Users\Choon\Desktop\CRISPR\Code\nt_data\chromosome_20_complete.txt', 'r') as myfile:
    #data=myfile.read()
    data=myfile.read().replace('\n', '')
    #data=data[12888830:12888860]

len_data=len(data)

N_A_G=float(data[0:len_data-1].count('A'))
N_T_G=float(data[0:len_data-1].count('T'))
N_C_G=float(data[0:len_data-1].count('C'))
N_G_G=float(data[0:len_data-1].count('G'))

p_A=N_A_G/(N_A_G+N_T_G+N_C_G+N_G_G)
p_T=N_T_G/(N_A_G+N_T_G+N_C_G+N_G_G)
p_C=N_C_G/(N_A_G+N_T_G+N_C_G+N_G_G)
p_G=N_G_G/(N_A_G+N_T_G+N_C_G+N_G_G)

p_N=(len_data-(N_A_G+N_T_G+N_C_G+N_G_G))/len_data

p=[p_A,p_T,p_C,p_G]


        
pair_2_count={'AA':0,'AT':0,'AC':0,'AG':0,'TA':0,'TT':0,'TC':0,'TG':0,'CA':0,'CT':0,'CC':0,'CG':0,'GA':0,'GT':0,'GC':0,'GG':0,'AN':0,'TN':0,'CN':0,'GN':0}

#correlation matrix
#AA AT AC AG
#TA TT TC TG
#CA CT CC CG
#GA GT GC GG

#ATCG, A=0, T=1, C=2, G=3

for i in range(0,len(data)-1):
    if data[i] !='N': 
        pair_2_count[data[i:i+2]]+=1
        
cm=[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]

cm[0][0]=pair_2_count['AA']/(N_A_G-pair_2_count['AN']) #AA
cm[0][1]=pair_2_count['AT']/(N_A_G-pair_2_count['AN']) #AT
cm[0][2]=pair_2_count['AC']/(N_A_G-pair_2_count['AN']) #AC
cm[0][3]=pair_2_count['AG']/(N_A_G-pair_2_count['AN']) #AG

cm[1][0]=pair_2_count['TA']/(N_T_G-pair_2_count['TN']) #TA
cm[1][1]=pair_2_count['TT']/(N_T_G-pair_2_count['TN']) #TT
cm[1][2]=pair_2_count['TC']/(N_T_G-pair_2_count['TN']) #TC
cm[1][3]=pair_2_count['TG']/(N_T_G-pair_2_count['TN']) #TG

cm[2][0]=pair_2_count['CA']/(N_C_G-pair_2_count['CN']) #CA
cm[2][1]=pair_2_count['CT']/(N_C_G-pair_2_count['CN']) #CT
cm[2][2]=pair_2_count['CC']/(N_C_G-pair_2_count['CN']) #CC
cm[2][3]=pair_2_count['CG']/(N_C_G-pair_2_count['CN']) #CG

cm[3][0]=pair_2_count['GA']/(N_G_G-pair_2_count['GN']) #GA
cm[3][1]=pair_2_count['GT']/(N_G_G-pair_2_count['GN']) #GT
cm[3][2]=pair_2_count['GC']/(N_G_G-pair_2_count['GN']) #GC
cm[3][3]=pair_2_count['GG']/(N_G_G-pair_2_count['GN']) #GG

print (p)
print (cm)