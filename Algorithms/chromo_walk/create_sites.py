with open(r'C:\Users\Choon\Desktop\CRISPR\Code\nt_data\chromosome_20_complete.txt', 'r') as myfile:
    #data=myfile.read()
    data=myfile.read().replace('\n', '')

f=open("C:\Users\Choon\Desktop\CRISPR\Code\Algorithms\chromo_walk\sites.txt","w+")

ts_len=20
PAM='CT' #NCT

sites=[]
for i in range(ts_len+(1+len(PAM))-1,len(data)):
    #for i in range(ts_len_max+(1+len(PAM))-1,len(data)):
    if data[i-len(PAM)+1:i+1]==PAM:
        f.write("%d\n" %i)