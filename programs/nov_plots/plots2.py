import numpy as np
import matplotlib.pyplot as plt
import scipy.stats

txt_data=np.genfromtxt(r'C:\Users\Choon\Desktop\CRISPR\Code\results\result5.txt',delimiter=',')
#txt_data=np.delete(txt_data,(385208),axis=0)

ts=txt_data[:,0]
ts_len=txt_data[:,1]
info=txt_data[:,2]
GC_frac=txt_data[:,3]
G_frac=txt_data[:,4]
q_comp_exp=txt_data[:,5]
Z=txt_data[:,6]
N_G=txt_data[1][7]
S_real=txt_data[:,8]
Boltz_prob=txt_data[:,9]

PR_ts_len_S_real=scipy.stats.pearsonr(ts_len,S_real) #Note that the linear fit is measured here.
#PR_GC_S_real=scipy.stats.pearsonr(GC_frac,S_real)
#PR_G_S_real=scipy.stats.pearsonr(G_frac,S_real)

SR_ts_len_S_real=scipy.stats.spearmanr(ts_len,S_real) #Spearman's rank assesses monotonic, not necessarily linear relationships
#SR_GC_S_real=scipy.stats.spearmanr(GC_frac,S_real)
#SR_G_S_real=scipy.stats.spearmanr(G_frac,S_real)

print("info-S_real Pearson correlation=",PR_ts_len_S_real)
#print("GC-S_real Pearson correlation=",PR_GC_S_real)
#print("G-S_real Pearson correlation=",PR_G_S_real)

print("info-S_real Spearman correlation=",SR_ts_len_S_real)
#print("GC-S_real Spearman correlation=",SR_GC_S_real)
#print("G-S_real Spearman correlation=",SR_G_S_real)

print("info-q_comp_exp Spearman correlation=",scipy.stats.spearmanr(ts_len,q_comp_exp))
#print("GC-q_comp_exp Spearman correlation=",scipy.stats.spearmanr(GC_frac,q_comp_exp))
#print("G-q_comp_exp Spearman correlation=",scipy.stats.spearmanr(G_frac,q_comp_exp))

print("info-Z Spearman correlation=",scipy.stats.spearmanr(ts_len,Z))
#print("GC-Z Spearman correlation=",scipy.stats.spearmanr(GC_frac,Z))
#print("G-Z Spearman correlation=",scipy.stats.spearmanr(G_frac,Z))

print("info-Off_bind_prob Spearman correlation=",scipy.stats.spearmanr(ts_len,1-Boltz_prob))
#print("GC-Off_bind_prob Spearman correlation=",scipy.stats.spearmanr(GC_frac,1-Boltz_prob))
#print("G-Off_bind_prob Spearman correlation=",scipy.stats.spearmanr(G_frac,1-Boltz_prob))

f1=plt.figure()
f2=plt.figure()
f3=plt.figure()
f4=plt.figure()
f5=plt.figure()
f6=plt.figure()


ax1=f1.add_subplot(111)
ax2=f2.add_subplot(111)
ax3=f3.add_subplot(111)
ax4=f4.add_subplot(111)
ax5=f5.add_subplot(111)
ax6=f6.add_subplot(111)
#ax7=f7.add_subplot(111)



ax1.scatter(ts_len,S_real,s=0.1)
ax1.set_title("ts len vs S_real")
ax1.set_xlabel("ts length")
ax1.set_ylabel("S_real")



ax2.scatter(ts_len,np.log(S_real)/np.log(10),s=0.1)
ax2.set_title("ts len vs log(S_real)")
ax2.set_xlabel("ts length")
ax2.set_ylabel("log_10 (S_real)")



ax3.scatter(ts_len,np.log(q_comp_exp)/np.log(10),s=0.1)
ax3.set_title("ts len vs log(q_comp_exp)")
ax3.set_xlabel("ts length")
ax3.set_ylabel("log_10 (q_comp_exp)")



ax4.scatter(ts_len,np.log(Z)/np.log(10),s=0.1)
ax4.set_title("ts len vs log(Z)")
ax4.set_xlabel("ts length")
ax4.set_ylabel("log_10 (Z)")



ax5.scatter(ts_len,np.log(1-Boltz_prob)/np.log(10),s=0.1)
ax5.set_title("ts len vs log(Off-bind prob)")
ax5.set_xlabel("ts length")
ax5.set_ylabel("log_10 (Off-bind prob)")


ax6.scatter(ts_len,1-Boltz_prob,s=0.1)
ax6.set_title("ts len vs Off-bind prob")
ax6.set_xlabel("ts length")
ax6.set_ylabel("Off-bind prob")

plt.show()