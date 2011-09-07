def tester(maxi, Lwin, Lts, Lsoi):
    for i in range(maxi):
        print 'soi['+str(max(0,i-Lts+Lwin))+':'+str(min(Lsoi-1,i+Lwin-1))+'] == ts[' + str(max(0,Lts-Lwin-i)) + ':' +str(Lts-1-max(0,Lwin+i-Lsoi)) + ']'