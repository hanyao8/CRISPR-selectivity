
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
from mpl_toolkits.mplot3d import Axes3D

plotswitch=2
plot_array=[11,13,15,17,19]

txt_data=np.genfromtxt(r'C:\Users\Choon\Desktop\CRISPR\Code\results\result4.txt',delimiter=',')
#txt_data=np.delete(txt_data,(385208),axis=0)

#txt_data=txt_data[0:100000]

data_len=len(txt_data)

ts=txt_data[:,0]
ts_len=txt_data[1][1]
info=txt_data[:,2]
GC_frac=txt_data[:,3]
G_frac=txt_data[:,4]
q_comp_exp=txt_data[:,5]
Z=txt_data[:,6]
N_G=txt_data[1][7]
S_real=txt_data[:,8]
Boltz_prob=txt_data[:,9]

n_S_hists=20 #number of selectivity histograms
n_bins=100

info_max=max(info)*1.001
info_min=min(info)
info_interval=info_max-info_min
info_d=info_interval/n_S_hists


info_index=((info-info_min)/info_interval*n_S_hists ).astype(int)

#Generate histograms


S_populations=[[] for i in range(0,n_S_hists)]
logS_populations=[[] for i in range(0,n_S_hists)]
info_pos=np.linspace(info_min+0.5*info_d,info_max-0.5*info_d,n_S_hists)


for i in range(0,data_len):
    S_populations[info_index[i]].append(S_real[i])
    logS_populations[info_index[i]].append( np.log(S_real[i])/np.log(10) )


if plotswitch==1:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    
    """
    ax1=f1.add_subplot(111,projection='3d')
    for i in range(0,n_S_hists):
        hist,bins=np.histogram(logS_populations[i],bins=n_bins)
        xs=(bins[:-1]+bins[1:])/2
        ax1.bar(xs,hist,zs=info_pos,zdir='y',alpha=0.5)
    """
    
    for i in range(0,n_S_hists):
    #for z in info_pos:
        ys = logS_populations[i]
    
        hist, bins = np.histogram(ys, bins=n_bins)
        xs = (bins[:-1] + bins[1:])/2
    
        ax.bar(xs, hist, zs=info_pos[i], zdir='y', alpha=0.3)
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    
if plotswitch==2:
    f1=plt.figure()
    ax1=f1.add_subplot(111)
    ax1.hist(logS_populations[plot_array[0]],n_bins)
    ax1.set_title("log(S_real) Histogram, ts entropy range: %f ~ %f"%(\
                  info_pos[plot_array[0]]-info_d/2,info_pos[plot_array[0]]+info_d/2))
    ax1.set_ylabel("Frequency")
    ax1.set_xlabel("log_10 (S_real)")
    
    f2=plt.figure()
    ax2=f2.add_subplot(111)
    ax2.hist(logS_populations[plot_array[1]],n_bins)
    ax2.set_title("log(S_real) Histogram, ts entropy range: %f ~ %f"%(\
                  info_pos[plot_array[1]]-info_d/2,info_pos[plot_array[1]]+info_d/2))
    ax2.set_ylabel("Frequency")
    ax2.set_xlabel("log_10 (S_real)")

    f3=plt.figure()
    ax3=f3.add_subplot(111)
    ax3.hist(logS_populations[plot_array[2]],n_bins)
    ax3.set_title("log(S_real) Histogram, ts entropy range: %f ~ %f"%(\
                  info_pos[plot_array[2]]-info_d/2,info_pos[plot_array[2]]+info_d/2))
    ax3.set_ylabel("Frequency")
    ax3.set_xlabel("log_10 (S_real)")
    
    f4=plt.figure()
    ax4=f4.add_subplot(111)
    ax4.hist(logS_populations[plot_array[3]],n_bins)
    ax4.set_title("log(S_real) Histogram, ts entropy range: %f ~ %f"%(\
                  info_pos[plot_array[3]]-info_d/2,info_pos[plot_array[3]]+info_d/2))
    ax4.set_ylabel("Frequency")
    ax4.set_xlabel("log_10 (S_real)")
    
    f5=plt.figure()
    ax5=f5.add_subplot(111)
    ax5.hist(logS_populations[plot_array[4]],n_bins)    
    ax5.set_title("log(S_real) Histogram, ts entropy range: %f ~ %f"%(\
                  info_pos[plot_array[4]]-info_d/2,info_pos[plot_array[4]]+info_d/2))
    ax5.set_ylabel("Frequency")
    ax5.set_xlabel("log_10 (S_real)")    

plt.show()