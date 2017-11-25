# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 00:18:36 2017

@author: Choon
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats

txt_data=np.genfromtxt(r'C:\Users\Choon\Desktop\CRISPR\Code\results\result4.txt',delimiter=',')
#txt_data=np.delete(txt_data,(385208),axis=0)

#txt_data=txt_data[0:100000]

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

info_min=min(info)
info_max=max(info)

PR_info_S_real=scipy.stats.pearsonr(info,S_real) #Note that the linear fit is measured here.
PR_GC_S_real=scipy.stats.pearsonr(GC_frac,S_real)
PR_G_S_real=scipy.stats.pearsonr(G_frac,S_real)

SR_info_S_real=scipy.stats.spearmanr(info,S_real) #Spearman's rank assesses monotonic, not necessarily linear relationships
SR_GC_S_real=scipy.stats.spearmanr(GC_frac,S_real)
SR_G_S_real=scipy.stats.spearmanr(G_frac,S_real)

print("info-S_real Pearson correlation=",PR_info_S_real)
print("GC-S_real Pearson correlation=",PR_GC_S_real)
print("G-S_real Pearson correlation=",PR_G_S_real)

print("info-S_real Spearman correlation=",SR_info_S_real)
print("GC-S_real Spearman correlation=",SR_GC_S_real)
print("G-S_real Spearman correlation=",SR_G_S_real)

print("info-q_comp_exp Spearman correlation=",scipy.stats.spearmanr(info,q_comp_exp))
print("GC-q_comp_exp Spearman correlation=",scipy.stats.spearmanr(GC_frac,q_comp_exp))
print("G-q_comp_exp Spearman correlation=",scipy.stats.spearmanr(G_frac,q_comp_exp))

print("info-Z Spearman correlation=",scipy.stats.spearmanr(info,Z))
print("GC-Z Spearman correlation=",scipy.stats.spearmanr(GC_frac,Z))
print("G-Z Spearman correlation=",scipy.stats.spearmanr(G_frac,Z))

print("info-Off_bind_prob Spearman correlation=",scipy.stats.spearmanr(info,1-Boltz_prob))
print("GC-Off_bind_prob Spearman correlation=",scipy.stats.spearmanr(GC_frac,1-Boltz_prob))
print("G-Off_bind_prob Spearman correlation=",scipy.stats.spearmanr(G_frac,1-Boltz_prob))

f1=plt.figure()
f2=plt.figure()
f3=plt.figure()
"""
f4=plt.figure()
f5=plt.figure()
f6=plt.figure()
f7=plt.figure()
f8=plt.figure()
f9=plt.figure()
f10=plt.figure()
f11=plt.figure()
f12=plt.figure()
"""
"""
ax1=f1.add_subplot(211)
ax2=f2.add_subplot(211)
ax3=f3.add_subplot(211)


#ax1=f1.add_subplot(211)
ax1hist=f1.add_subplot(212)
#ax1hist=ax1.twinx()
#ax2=f2.add_subplot(211)
ax2hist=f2.add_subplot(212)
#ax2hist=ax2.twinx()
#ax3=f3.add_subplot(211)
ax3hist=f3.add_subplot(212)
#ax3hist=ax3.twinx()
"""

ax1hist=f1.add_subplot(111)
ax2hist=f2.add_subplot(111)
ax3hist=f3.add_subplot(111)

"""

ax4=f4.add_subplot(211)
ax4hist=f4.add_subplot(212)
#ax4hist=ax4.twinx()
ax5=f5.add_subplot(211)
ax5hist=f5.add_subplot(212)
#ax5hist=ax5.twinx()
ax6=f6.add_subplot(211)
ax6hist=f6.add_subplot(212)
#ax6hist=ax6.twinx()
ax7=f7.add_subplot(211)
ax7hist=f7.add_subplot(212)
#ax7hist=ax7.twinx()
ax8=f8.add_subplot(211)
ax8hist=f8.add_subplot(212)
#ax8hist=ax8.twinx()
ax9=f9.add_subplot(211)
ax9hist=f9.add_subplot(212)
#ax9hist=ax9.twinx()
ax10=f10.add_subplot(211)
ax10hist=f10.add_subplot(212)
#ax10hist=ax10.twinx()
ax11=f11.add_subplot(211)
ax11hist=f11.add_subplot(212)
#ax11hist=ax11.twinx()
ax12=f12.add_subplot(211)
ax12hist=f12.add_subplot(212)
#ax12hist=ax12.twinx()

"""

"""
ax1.scatter(info,S_real,s=0.1)
ax1.set_title("ts entropy vs S_real")
ax1.set_xlabel("ts Shannon Entropy")
ax1.set_ylabel("S_real")

ax2.scatter(GC_frac,S_real,s=0.1)
ax2.set_title("ts GC composition vs S_real")
ax2.set_xlabel("ts GC composition")
ax2.set_ylabel("S_real")

ax3.scatter(G_frac,S_real,s=0.1)
ax3.set_title("ts G composition vs S_real")
ax3.set_xlabel("ts G composition")
ax3.set_ylabel("S_real")

"""


"""
ax1sc=ax1.scatter(info,np.log(S_real)/np.log(10),s=0.1,c=GC_frac,cmap=plt.cm.coolwarm)
#ax1.scatter(info,np.log(S_real)/np.log(10),s=0.1,alpha=0.7,c=1-GC_frac,cmap=plt.cm.coolwarm)
ax1.set_title("ts entropy vs log(S_real)")
ax1.set_xlabel("ts Shannon Entropy")
ax1.set_ylabel("log_10 (S_real)")
cbar1=f1.colorbar(ax1sc)
cbar1.ax.set_yticklabels(["0.0 GC","0.2 GC","0.4 GC","0.6 GC","0.8 GC","1.0 GC"])
"""

ax1hist.hist(info,bins=100,edgecolor='black',linewidth=0.6,fill=False)
#ax1hist.hist(info,bins=100,edgecolor='black',linewidth=0.35,fill=False)
#ax1hist.hist(info,log=True,bins=100,alpha=0.5)
ax1hist.set_ylabel("Frequency")
#ax1hist.yscale('log')
ax1hist.set_xlabel("ts Shannon Entropy")

"""
ax2sc=ax2.scatter(GC_frac,np.log(S_real)/np.log(10),s=0.1,c=info,cmap=plt.cm.coolwarm)
ax2.set_title("ts GC composition vs log(S_real)")
ax2.set_xlabel("ts GC composition")
ax2.set_ylabel("log_10 (S_real)")
cbar2=f2.colorbar(ax2sc)
cbar2.ax.set_yticklabels(np.round(np.linspace(info_min,info_max,8),2))
#cbar2.ax.set_yticklabels(["entropy=%f"%info_min,"entropy=%f"%info_max])
"""

ax2hist.hist(GC_frac,bins=100,edgecolor='black',linewidth=0.6,fill=False)
#ax2hist.hist(GC_frac,bins=100,edgecolor='black',linewidth=0.35,fill=False)
ax2hist.set_ylabel("Frequency")
ax2hist.set_xlabel("ts GC composition")

"""
ax3sc=ax3.scatter(G_frac,np.log(S_real)/np.log(10),s=0.1,c=info,cmap=plt.cm.coolwarm)
ax3.set_title("ts G composition vs log(S_real)")
ax3.set_xlabel("ts G composition")
ax3.set_ylabel("log_10 (S_real)")
cbar3=f3.colorbar(ax3sc)
cbar3.ax.set_yticklabels(np.round(np.linspace(info_min,info_max,8),2))
#cbar3.ax.set_yticklabels(["entropy=%f"%info_min,"entropy=%f"%info_max])

"""
ax3hist.hist(G_frac,bins=100,edgecolor='black',linewidth=0.6,fill=False)
#ax3hist.hist(G_frac,bins=100,edgecolor='black',linewidth=0.35,fill=False)
ax3hist.set_xlabel("ts G composition")
ax3hist.set_ylabel("Frequency")

"""
ax4.scatter(info,np.log(q_comp_exp)/np.log(10),s=0.1,c=1-GC_frac,cmap=plt.cm.coolwarm)
ax4.set_title("ts entropy vs log(q_comp_exp)")
ax4.set_xlabel("ts Shannon Entropy")
ax4.set_ylabel("log_10 (q_comp_exp)")

ax4hist.hist(info,bins=100,edgecolor='black',linewidth=0.6,fill=False)
#ax4hist.hist(info,bins=100,edgecolor='black',linewidth=0.35,fill=False)
ax4hist.set_ylabel("Frequency")
ax4hist.set_xlabel("ts Shannon Entropy")

ax5.scatter(GC_frac,np.log(q_comp_exp)/np.log(10),s=0.1,c=1.4-info,cmap=plt.cm.coolwarm)
ax5.set_title("ts GC composition vs log(q_comp_exp)")
ax5.set_xlabel("ts GC composition")
ax5.set_ylabel("log_10 (q_comp_exp)")

ax5hist.hist(GC_frac,bins=100,edgecolor='black',linewidth=0.6,fill=False)
#ax5hist.hist(GC_frac,bins=100,edgecolor='black',linewidth=0.35,fill=False)
ax5hist.set_ylabel("Frequency")
ax5hist.set_xlabel("ts GC composition")

ax6.scatter(G_frac,np.log(q_comp_exp)/np.log(10),s=0.1,c=1.4-info,cmap=plt.cm.coolwarm)
ax6.set_title("ts G composition vs log(q_comp_exp)")
ax6.set_xlabel("ts G composition")
ax6.set_ylabel("log_10 (q_comp_exp)")

ax6hist.hist(G_frac,bins=100,edgecolor='black',linewidth=0.6,fill=False)
#ax6hist.hist(G_frac,bins=100,edgecolor='black',linewidth=0.35,fill=False)
ax6hist.set_ylabel("Frequency")
ax6hist.set_xlabel("ts G composition")

ax7.scatter(info,np.log(Z)/np.log(10),s=0.1,c=1-GC_frac,cmap=plt.cm.coolwarm)
ax7.set_title("ts entropy vs log(Z)")
ax7.set_xlabel("ts Shannon Entropy")
ax7.set_ylabel("log_10 (Z)")

ax7hist.hist(info,bins=100,edgecolor='black',linewidth=0.6,fill=False)
#ax7hist.hist(info,bins=100,edgecolor='black',linewidth=0.35,fill=False)
ax7hist.set_ylabel("Frequency")
ax7hist.set_xlabel("ts Shannon Entropy")

ax8.scatter(GC_frac,np.log(Z)/np.log(10),s=0.1,c=1.4-info,cmap=plt.cm.coolwarm)
ax8.set_title("ts GC composition vs log(Z)")
ax8.set_xlabel("ts GC composition")
ax8.set_ylabel("log_10 (Z)")

ax8hist.hist(GC_frac,bins=100,edgecolor='black',linewidth=0.6,fill=False)
#ax8hist.hist(GC_frac,bins=100,edgecolor='black',linewidth=0.35,fill=False)
ax8hist.set_ylabel("Frequency")
ax8hist.set_xlabel("ts GC composition")

ax9.scatter(G_frac,np.log(Z)/np.log(10),s=0.1,c=1.4-info,cmap=plt.cm.coolwarm)
ax9.set_title("ts G composition vs log(Z)")
ax9.set_xlabel("ts G composition")
ax9.set_ylabel("log_10 (Z)")

ax9hist.hist(G_frac,bins=100,edgecolor='black',linewidth=0.6,fill=False)
#ax9hist.hist(G_frac,bins=100,edgecolor='black',linewidth=0.35,fill=False)
ax9hist.set_ylabel("Frequency")
ax9hist.set_xlabel("ts G composition")

ax10.scatter(info,np.log(1-Boltz_prob)/np.log(10),s=0.1,c=1-GC_frac,cmap=plt.cm.coolwarm)
ax10.set_title("ts entropy vs log(Off-bind prob)")
ax10.set_xlabel("ts Shannon entropy")
ax10.set_ylabel("log_10 (Off-bind prob)")

ax10hist.hist(info,bins=100,edgecolor='black',linewidth=0.6,fill=False)
#ax10hist.hist(info,bins=100,edgecolor='black',linewidth=0.35,fill=False)
ax10hist.set_ylabel("Frequency")
ax10hist.set_xlabel("ts Shannon Entropy")

ax11.scatter(GC_frac,np.log(1-Boltz_prob)/np.log(10),s=0.1,c=1.4-info,cmap=plt.cm.coolwarm)
ax11.set_title("ts GC composition vs log(Off-bind prob)")
ax11.set_xlabel("ts GC composition")
ax11.set_ylabel("log_10 (Off-bind prob)")

ax11hist.hist(GC_frac,bins=100,edgecolor='black',linewidth=0.6,fill=False)
#ax11hist.hist(GC_frac,bins=100,edgecolor='black',linewidth=0.35,fill=False)
ax11hist.set_ylabel("Frequency")
ax11hist.set_xlabel("ts GC composition")

ax12.scatter(G_frac,np.log(1-Boltz_prob)/np.log(10),s=0.1,c=1.4-info,cmap=plt.cm.coolwarm)
ax12.set_title("ts G composition vs log(Off-bind prob)")
ax12.set_xlabel("ts G composition")
ax12.set_ylabel("log_10 (Off-bind prob)")

ax12hist.hist(G_frac,bins=100,edgecolor='black',linewidth=0.6,fill=False)
#ax12hist.hist(G_frac,bins=100,edgecolor='black',linewidth=0.35,fill=False)
ax12hist.set_ylabel("Frequency")
ax12hist.set_xlabel("ts G composition")
"""
plt.show()