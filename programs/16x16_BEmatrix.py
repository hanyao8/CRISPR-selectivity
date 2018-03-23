# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 12:10:00 2018

@author: Choon
"""
import numpy as np
import os

pw_comp_dict={"A":"T","T":"A","C":"G","G":"C"}
comp_dict={"TA":"AT","AT":"TA","AA":"TT","TT":"AA","CT":"GA","AG":"TC","GA":"CT","TC":"AG","GT":"CA","AC":"TG","CA":"GT","TG":"AC","GG":"CC","CC":"GG","CG":"GC","GC":"CG"}

#cols=["TA","AT","AA","TT","CT","AG","GA","TC","GT","AC","CA","TG","GG","CC","CG","GC"]
cols=["AA","AT","TA","AG",  "CT","TT","TG","TC",  "GT","GC","CG","CC",  "GG","GA","AC","CA"]
rows=[comp_dict[cols[i]] for i in range(0,16)]



def score(a,b):
    num=0
    for i in range(0,2):
        if not(b[i] == pw_comp_dict[a[i]]):
            num+=1
    return(num)

M_BE=[[] for i in range(0,16)]
#M_BE=np.empty((16,16))
for i in range(0,16):
    for j in range(0,16):
        #M_BE[i][j]=score(rows[i],cols[j])
        M_BE[i].append(score(rows[i],cols[j]))
        
for i in range(0,16):
    print(M_BE[i])
#M_BE.astype(int)
#np.savetxt(os.getcwd()+"\\M_BE.txt",M_BE)
    
    