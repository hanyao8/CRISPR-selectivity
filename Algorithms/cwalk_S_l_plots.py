import numpy as np
import matplotlib.pyplot as plt

switch=1

txt_data1=np.genfromtxt(r'C:\Users\Choon\Desktop\CRISPR\Code\console_sim_data\cwalk_S_l_main.csv',delimiter=',')
l1=txt_data1[1:,1]
S1=txt_data1[1:,2]
log_S1=np.log(S1)/np.log(10)
S1_pred=txt_data1[1:,3]
log_S1_pred=np.log(S1_pred)/np.log(10)
f_GC_1=txt_data1[1:,4]

txt_data2=np.genfromtxt(r'C:\Users\Choon\Desktop\CRISPR\Code\console_sim_data\cwalk_S_l_extra.csv',delimiter=',')
l2=txt_data2[1:,1]
S2=txt_data2[1:,2]
log_S2=np.log(S2)/np.log(10)
S2_pred=txt_data2[1:,3]
log_S2_pred=np.log(S2_pred)/np.log(10)
f_GC_2=txt_data2[1:,4]

l=np.append(l1,l2)
S=np.append(S1,S2)
log_S=np.log(S)/np.log(10)
log_S_pred=np.append(log_S1_pred,log_S2_pred)
f_GC=np.append(f_GC_1,f_GC_2)

#cd_pod is a list containing the columns at which a sequence of consecutive data starts
#'consecutive data' refers to data derived from targeting sequences taken next to a fixed point in the genome
cd_pos=[0]
for i in range(1,len(l1)):
    if (l1[i]-l1[i-1]) != 1:
        cd_pos.append(i)
cd_pos.append(len(l1))

if switch==0:
    plt.scatter(l,log_S)
    for i, txt in enumerate(f_GC_2):
        plt.annotate(round(txt,2),(l2[i],log_S2[i]))
        
    for j in range(0,len(cd_pos)-1):
        plt.plot(l1[cd_pos[j]:cd_pos[j+1]],log_S1[cd_pos[j]:cd_pos[j+1]])
    plt.title("log(S_ch20) against l")
    plt.xlabel("targeting sequence length")
    plt.ylabel("log_10 (S_ch20)")

#Plotting the deviation of computed selectivities from predicted
if switch==1:
    dev=(log_S-log_S_pred)/log_S_pred
    
    plt.scatter(l,dev)

    plt.title("Deviation against l")
    plt.xlabel("targeting sequence length")
    plt.ylabel("Deviation (Delta)")


plt.grid(True)
plt.show()