import numpy as np
import matplotlib.pyplot as plt




"""
Z_cml=[3.36864913908e+21,7.60677889983e+21,1.01394939395e+22,1.10960139377e+22,1.13522111908e+22,1.14039736024e+22,1.14121638707e+22,1.14132037831e+22]

x=[]
for i in range(0,len(Z_cml)):
    x.append(i)

#fig = plt.figure(1, figsize=(9,8))
plt.scatter(x,Z_cml)
#plt.plot(phi,fit_func(phi,po[0]), linewidth=2.5, label='Best Fit')
plt.title("Cumulative Value of Partition Function")
#plt.ylabel("")
plt.xlabel("Number of Defects")
plt.show()
"""

def f(n):
    for i in range(1,n+1):
        yield np.array([1,2,3,4,5])*i
    
F=f(5)
#f=[[1,2,3,4,5],[4,3,5,3,1],[6,4,2,6,1]]
x=[1,2,3,4,5]

for j in F:
    plt.plot(x,j, linewidth=2.5)
    
plt.show()



"""
ts_len=20
x=[]
for i in range(0,ts_len+1):
    x.append(i)


Z_comp_list=[0.034435978146603409,
 6.2888607189803221,
 528.34949409492219,
 27221.955809008527,
 966947.3764696985,
 25225426.163360886,
 502483638.912067,
 7840569151.0436554,
 97493494531.505829,
 977082644408.30054,
 7946707949608.8369,
 52608221544996.234,
 283314605189383.25,
 1235714389855250.7,
 4326640162339753.5,
 11984160478879270.0,
 25664252757704416.0,
 40982089489876960.0,
 45936969432537944.0,
 32246320827622364.0,
 10667142427152504.0]

#fig = plt.figure(1, figsize=(9,8))
plt.scatter(x,Z_comp_list)
#plt.plot(phi,fit_func(phi,po[0]), linewidth=2.5, label='Best Fit')
plt.title("Parition Function Components")
#plt.ylabel("")
plt.xlabel("Number of Complementary Pairs")
plt.show()
"""