#GC_composition counting script

import numpy as np

txt_data=np.genfromtxt(r'C:\Users\Choon\Desktop\CRISPR\Code\console_sim_data\cwalk_S_l_main.csv',delimiter=',',dtype=None)
seqs=txt_data[1:,0]

ts='TGAAAGAATGGGGACAACAG'
for ts in seqs:
    ts=list(ts)
    a=float(ts.count('G')+ts.count('C'))
    GC_frac=a/float(len(ts))
    print(GC_frac)