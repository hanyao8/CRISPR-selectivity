#chromo_walk main

import os
os.chdir("C:\\Users\\Choon\\Documents\\GitHub\\CRISPR-selectivity-OO\\chwalk")

import numpy as np
import numpy.random as nr
import matplotlib.pyplot as plt
import single_walk as sw


#setting=0: Computation of selectivity of (n_iters) random sequences of the same length
#setting=1: Computation of selectivity of binding between ts and cs of different l at the same site in the genome

class sim:

    def __init__(self,param,ts_number,switch,n_iters,PAM):

        with open(r'C:\Users\Choon\Desktop\CRISPR\Code\nt_data\chromosome_20_complete.txt', 'r') as myfile1:
            #data=myfile.read()
            data=myfile1.read().replace('\n', '')
            
        sites = np.genfromtxt(r'C:\Users\Choon\Desktop\CRISPR\Code\Algorithms\chromo_walk\sites.txt', delimiter='\n')
        sites = np.array(sites)
            
        print (len(sites))
        
        S_list=[]  
        q_comp_list=[]
        comp_count_list=[]
        Z_list=[]
        
        if switch==41:
        
            #ts_min=20
            #ts_max=20
            #ts_len_min=ts_len_max=20
            ts_len=ts_number
        
            random_sites=nr.choice(sites,n_iters)
            #random_sites=sites[0:1]
            random_sites = [int(x) for x in random_sites]
            
            print(random_sites)
            
        
            for i in random_sites:
                print("i=%d" %(i))
                a=sw.iter(param,data,i,PAM,ts_len)
                
                S_list.append(a._iter__S)
                q_comp_list.append(a._iter__q_comp)
                comp_count_list.append(a._iter__comp_count)
                Z_list.append(a._iter__Z)
                
                x=np.array(a._iter__j_list)
                y=np.array(a._iter__Z_nc_list)
                y=np.log(y)/np.log(10)
                plt.plot(x,y)
                
                print("len(x)=%d, len(y)=%d" %(len(x),len(y)))
                
        if switch==42:
            random_site=int(nr.choice(sites))
            print("random site=%d" %random_site)
            
            ts_start=ts_number
            ts_end=ts_start+n_iters-1
            
            for ts_len in range(ts_start,ts_end+1):
                print("l=%d" %(ts_len))
                a=sw.iter(param,data,random_site,PAM,ts_len)
                
                S_list.append(a._iter__S)
                q_comp_list.append(a._iter__q_comp)
                comp_count_list.append(a._iter__comp_count)
                Z_list.append(a._iter__Z)
                
                x=np.array(a._iter__j_list)
                y=np.array(a._iter__Z_nc_list)
                y=np.log(y)/np.log(10)
                plt.plot(x,y,label="l=%d"%(ts_len))
                plt.legend(loc=4)
                
                print("len(x)=%d, len(y)=%d" %(len(x),len(y)))        
            
                
        plt.title("nc Component of Z")
        plt.ylabel("log(Z_nc)")
        plt.xlabel("Binding Events")
        plt.show()