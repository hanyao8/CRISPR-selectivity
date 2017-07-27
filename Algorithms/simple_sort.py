#Simple sort

a=[5,7,3,7,4,54,2,465,54,13,57,8765,236,4]
temp=0
X=False
X_list=[]
for j in range(0,len(a)-1):
    X_list.append(False)

while not(X):
    X=True
    for i in range(0,len(a)-1):
        if a[i] > a[i+1]:
            temp=a[i]
            a[i]=a[i+1]
            a[i+1]=temp  
    
    for j in range(0,len(a)-1):
        X_list[j]=(a[j]<=a[j+1]) 
        X=X and X_list[j]
        
        
        