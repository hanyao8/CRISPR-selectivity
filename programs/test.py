# -*- coding: utf-8 -*-
"""
Created on Sat Feb 10 20:55:36 2018

@author: Choon
"""

import numpy as np

def stats_2_BEtensor(stats):
    pw_comp_dict={"A":"T","T":"A","C":"G","G":"C"}
    num_comp_dict={0:1,1:0,2:3,3:2}
    #comp_dict={"TA":"AT","AT":"TA","AA":"TT","TT":"AA","CT":"GA","AG":"TC","GA":"CT","TC":"AG","GT":"CA","AC":"TG","CA":"GT","TG":"AC","GG":"CC","CC":"GG","CG":"GC","GC":"CG"}
    #nt_2_num={"A":0,"T":1,"C":2,"G":3}
    #num_2_nt={0:"A",1:"T",2:"C",3:"G"}
    def ntscore(a,b):
        num=0
        for i in range(0,2):
            if not(b[i] == pw_comp_dict[a[i]]):
                num+=1
        return(num)
        
    def numscore(a,b,c,d):
        num=0
        if not(a == num_comp_dict[c]):
            num+=1
        if not(b == num_comp_dict[d]):
            num+=1            
        return(num)
        
    #cur_time=time.time()
    #data_sheet=io.open(os.getcwd()+"\\main35result_%d.txt"%cur_time,'a')    

    #BE_stats=np.array([[-1.4180,0.4152,2.2767],[0.512,0.507,0.564]]) #unique
    #BE_stats2=np.array([[-1.4056,0.4152,2.2767],[0.433,0.507,0.543]])
    
    BE_stats=
    E_uniq=[np.random.normal(BE_stats[0][0],BE_stats[1][0],10),np.random.normal(BE_stats[0][1],BE_stats[1][1],48),np.random.normal(BE_stats[0][2],BE_stats[1][2],78)]  
    
    count=[0,0,0]
    
    G_dpx35=np.zeros((4,4,4,4))
    """
    for i in range(0,4): #populating the symmetric terms
        for j in range(0,4):
            if numscore(i,j,j,i) == 2:
                G_dpx35[i][j][j][i]=E_uniq[2][count[2]]
                count[2]+=1
            else:
                G_dpx35[i][j][j][i]=E_uniq[0][count[0]]
                count[0]+=1
    """
    for i in range(0,4):
        for l in range(i,4):
            for j in range(0,4):
                if i!=l:
                    k_start=0
                else:
                    k_start=j
                for k in range(k_start,4):
                    if not(i==l and j==k):
                        if numscore(i,j,k,l) == 2: 
                            G_dpx35[i][j][k][l]=E_uniq[2][count[2]]
                            G_dpx35[l][k][j][i]=E_uniq[2][count[2]]
                            count[2]+=1
                        elif numscore(i,j,k,l) == 1: 
                            G_dpx35[i][j][k][l]=E_uniq[1][count[1]]
                            G_dpx35[l][k][j][i]=E_uniq[1][count[1]]
                            count[1]+=1
                        else:
                            G_dpx35[i][j][k][l]=E_uniq[0][count[0]]
                            G_dpx35[l][k][j][i]=E_uniq[0][count[0]]
                            count[0]+=1                            
                    else:
                        if numscore(i,j,k,l) == 2:
                            G_dpx35[i][j][k][l]=E_uniq[2][count[2]]
                            count[2]+=1
                        else:
                            G_dpx35[i][j][k][l]=E_uniq[0][count[0]]
                            count[0]+=1

    return(G_dpx35)