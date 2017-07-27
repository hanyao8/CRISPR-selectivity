#Algorithm to compute parameters of nucleotide frequency and correlation in a genome


#txt_data = np.genfromtxt(r'C:\Users\Choon\Desktop\CRISPR\Code\nt_data\chromosome_1_1_1000000.txt')
#print txt_data

with open(r'C:\Users\Choon\Desktop\CRISPR\Code\nt_data\chromosome_1_1_1000000.txt', 'r') as myfile:
    #data=myfile.read()
    data=myfile.read().replace('\n', '')

len_data=len(data)
N_samples=10
len_sample=len_data/N_samples

"""
N_A_G=data.count('A')
N_T_G=data.count('T')
N_C_G=data.count('C')
N_G_G=data.count('G')


"""
for i in range(0,N_samples):
    
    N_A_G=0
    N_T_G=0
    N_C_G=0
    N_G_G=0


    AA_count=0
    AT_count=0
    AC_count=0
    AG_count=0
    TA_count=0
    TT_count=0
    TC_count=0
    TG_count=0
    CA_count=0
    CT_count=0
    CC_count=0
    CG_count=0
    GA_count=0
    GT_count=0
    GC_count=0
    GG_count=0
    
    cm=[[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
    
    for j in range(int(i*len_sample),int((i+1)*len_sample)):
        if j < len_data-1:

            if data[j]=='A':
                N_A_G+=1
                if data[j+1]=='A':
                    AA_count+=1
                if data[j+1]=='T':
                    AT_count+=1
                if data[j+1]=='C':
                    AC_count+=1
                if data[j+1]=='G':
                    AG_count+=1
            if data[j]=='T':
                N_T_G+=1
                if data[j+1]=='A':
                    TA_count+=1
                if data[j+1]=='T':
                    TT_count+=1
                if data[j+1]=='C':
                    TC_count+=1
                if data[j+1]=='G':
                    TG_count+=1
            if data[j]=='C':
                N_C_G+=1
                if data[j+1]=='A':
                    CA_count+=1
                if data[j+1]=='T':
                    CT_count+=1
                if data[j+1]=='C':
                    CC_count+=1
                if data[j+1]=='G':
                    CG_count+=1
            if data[j]=='G':
                N_G_G+=1
                if data[j+1]=='A':
                    GA_count+=1
                if data[j+1]=='T':
                    GT_count+=1
                if data[j+1]=='C':
                    GC_count+=1
                if data[j+1]=='G':
                    GG_count+=1
            
    p_A=N_A_G/(N_A_G+N_T_G+N_C_G+N_G_G)
    p_T=N_T_G/(N_A_G+N_T_G+N_C_G+N_G_G)
    p_C=N_C_G/(N_A_G+N_T_G+N_C_G+N_G_G)
    p_G=N_G_G/(N_A_G+N_T_G+N_C_G+N_G_G)
        
    p=[p_A,p_T,p_C,p_G]
    
    #correlation matrix
    #AA AT AC AG
    #TA TT TC TG
    #CA CT CC CG
    #GA GT GC GG

    #ATCG, A=0, T=1, C=2, G=3
    
    cm[0][0]=AA_count/N_A_G #AA
    cm[0][1]=AT_count/N_A_G #AT
    cm[0][2]=AC_count/N_A_G #AC
    cm[0][3]=AG_count/N_A_G #AG

    cm[1][0]=TA_count/N_T_G #TA
    cm[1][1]=TT_count/N_T_G #TT
    cm[1][2]=TC_count/N_T_G #TC
    cm[1][3]=TG_count/N_T_G #TG

    cm[2][0]=CA_count/N_C_G #CA
    cm[2][1]=CT_count/N_C_G #CT
    cm[2][2]=CC_count/N_C_G #CC
    cm[2][3]=CG_count/N_C_G #CG

    cm[3][0]=GA_count/N_G_G #GA
    cm[3][1]=GT_count/N_G_G #GT
    cm[3][2]=GC_count/N_G_G #GC
    cm[3][3]=GG_count/N_G_G #GG
        
    print ('p',p)
    print ('cm',cm)
        




