import numpy
import pickle
import random
import time
import os
import matplotlib.pylab as plt
import sys

def BuildSeqsID(soiLen,lenStd,subLen,subStd,nSub,minLenSeq,minLenSub,iRate,dRate):
    t0 = time.time()
    soi = []
    ts = []
    tsInfo = []
    AA = ['G','P','A','V','L','I','M','C','F','Y','W','H','K','R','Q','N','E','D','S','T']
    # make soi
    for i in range(soiLen):
        soi.append(random.choice(AA))
    # make test sequences
    for i in range(10): # for identity 1.0 to 0.1 by 0.1 increments
        for j in range(nSub): # make nSub of each
            seq = []
            # the seqlen is a N(subLen,subStd) > minLen
            seqLen = int(max(minLenSeq,random.gauss(soiLen,lenStd)))
            # the sub seq len is a soiLen > N(subLen,subStd) > minLenSub
            seqSubLen = int(min(min(soiLen,seqLen),max(minLenSub,random.gauss(subLen,subStd))))
            # the relative placements are uniform random
            seqSubPlace = int(random.uniform(0,seqLen-seqSubLen-1))
            relSoiPlace = int(random.uniform(0,soiLen-seqSubLen-1))
            #print seqLen,seqSubLen,seqSubPlace,relSoiPlace
            # build a random sequence
            for k in range(seqLen):
                seq.append(random.choice(AA))
            # add in identity within the subregion
            m = 0
            for k in range(seqSubLen):
                if random.uniform(0,1) > i*0.1:
                    #print k,seqSubPlace,seqLen,m,relSoiPlace,soiLen
                    seq[k + seqSubPlace] = soi[m + relSoiPlace]
                else:
                    while(seq[k + seqSubPlace] == soi[m + relSoiPlace]):
                        seq[k + seqSubPlace] = random.choice(AA)
                m += 1
            # insertion and deletion time!
            tmpLen = seqSubLen
            k = 0
            while(k < seqSubLen):
                if random.uniform(0,1) < iRate:
                    seq.insert(k + seqSubPlace,random.choice(AA))
                    seqSubLen += 1
                    seqLen += 1
                if random.uniform(0,1) < dRate:
                    seq.pop(k + seqSubPlace)
                    seqSubLen -= 1
                    seqLen -= 1
                    k -= 1
                k += 1
            seqSubLen = tmpLen
            # add info to the lists
            ts.append(seq)
            tsInfo.append([i*nSub+j,seqLen,seqSubLen,seqSubPlace,relSoiPlace])
    print 'Generated seqs in ' + str(time.time() - t0) + ' seconds.'
    t0 = time.time()
    outdir = 'benchmark' + str(soiLen) + 'ID'
    if outdir not in os.listdir('.'):
        os.mkdir(outdir)
    # output the soi
    fsoi = open(outdir + '/soiSeqID.faa','w')
    fts = open(outdir + '/testSeqsID.faa','w')
    fvr = open(outdir + '/benchVRID.faa','w')
    soiStr = ''
    for i in range(soiLen):
        soiStr = soiStr + soi[i]
    fsoi.write('>SeqOfInt_'+str(soiLen)+','+str(lenStd)+','+str(subLen)+','+str(subStd)+','+str(nSub)+','+str(minLenSeq)+','+str(minLenSub)+'\n')
    fsoi.write(soiStr + '\n')
    fvr.write('>SeqOfInt_'+str(soiLen)+','+str(lenStd)+','+str(subLen)+','+str(subStd)+','+str(nSub)+','+str(minLenSeq)+','+str(minLenSub)+'\n')
    fvr.write(soiStr + '\n')
    for i in range(10*nSub):
        tsStr = ''
        for j in range(tsInfo[i][1]):
            tsStr = tsStr + ts[i][j]
        fts.write('>ts'+str(i)+'nsub'+str(nSub)+'\n')
        fts.write(tsStr + '\n')
        fvr.write('>ts'+str(i)+'nsub'+str(nSub)+'\n')
        fvr.write(tsStr + '\n')
    fsoi.close()
    fts.close()
    fvr.close()
    # pickle out
    output = open(outdir + '/tsInfo.pkl', 'wb')
    pickle.dump([[soiLen,lenStd,subLen,subStd,nSub,minLenSeq,minLenSub],tsInfo],output)
    print 'Outputed seqs in ' + str(time.time() - t0) + ' seconds.'

def BuildSeqs(soiLen,lenStd,subLen,subStd,nSub,minLenSeq,minLenSub):
    t0 = time.time()
    soi = []
    ts = []
    tsInfo = []
    AA = ['G','P','A','V','L','I','M','C','F','Y','W','H','K','R','Q','N','E','D','S','T']
    # make soi
    for i in range(soiLen):
        soi.append(random.choice(AA))
    # make test sequences
    for i in range(10): # for identity 1.0 to 0.1 by 0.1 increments
        for j in range(nSub): # make nSub of each
            seq = []
            # the seqlen is a N(subLen,subStd) > minLen
            seqLen = int(max(minLenSeq,random.gauss(soiLen,lenStd)))
            # the sub seq len is a soiLen > N(subLen,subStd) > minLenSub
            seqSubLen = int(min(min(soiLen,seqLen),max(minLenSub,random.gauss(subLen,subStd))))
            # the relative placements are uniform random
            seqSubPlace = int(random.uniform(0,seqLen-seqSubLen-1))
            relSoiPlace = int(random.uniform(0,soiLen-seqSubLen-1))
            #print seqLen,seqSubLen,seqSubPlace,relSoiPlace
            # build a random sequence
            for k in range(seqLen):
                seq.append(random.choice(AA))
            # add in identity within the subregion
            m = 0
            
            for k in range(seqSubLen):
                if random.uniform(0,1) > i*0.1:
                    #print k,seqSubPlace,seqLen,m,relSoiPlace,soiLen
                    seq[k + seqSubPlace] = soi[m + relSoiPlace]
                else:
                    while(seq[k + seqSubPlace] == soi[m + relSoiPlace]):
                        seq[k + seqSubPlace] = random.choice(AA)
                m += 1
            
            # add info to the lists
            ts.append(seq)
            tsInfo.append([i*nSub+j,seqLen,seqSubLen,seqSubPlace,relSoiPlace])
    print 'Generated seqs in ' + str(time.time() - t0) + ' seconds.'
    t0 = time.time()
    outdir = 'benchmarks/benchmark' + str(soiLen)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    # output the soi
    fsoi = open(outdir + '/soiSeq.faa','w')
    fts = open(outdir + '/testSeqs.faa','w')
    fvr = open(outdir + '/benchVR.faa','w')
    soiStr = ''
    for i in range(soiLen):
        soiStr = soiStr + soi[i]
    fsoi.write('>SeqOfInt_'+str(soiLen)+','+str(lenStd)+','+str(subLen)+','+str(subStd)+','+str(nSub)+','+str(minLenSeq)+','+str(minLenSub)+'\n')
    fsoi.write(soiStr + '\n')
    fvr.write('>SeqOfInt_'+str(soiLen)+','+str(lenStd)+','+str(subLen)+','+str(subStd)+','+str(nSub)+','+str(minLenSeq)+','+str(minLenSub)+'\n')
    fvr.write(soiStr + '\n')
    for i in range(10*nSub):
        tsStr = ''
        for j in range(tsInfo[i][1]):
            tsStr = tsStr + ts[i][j]
        fts.write('>ts'+str(i)+'nsub'+str(nSub)+'\n')
        fts.write(tsStr + '\n')
        fvr.write('>ts'+str(i)+'nsub'+str(nSub)+'\n')
        fvr.write(tsStr + '\n')
    fsoi.close()
    fts.close()
    fvr.close()
    # pickle out
    output = open(outdir + '/tsInfo.pkl', 'wb')
    pickle.dump([[soiLen,lenStd,subLen,subStd,nSub,minLenSeq,minLenSub],tsInfo],output)
    print 'Outputed seqs in ' + str(time.time() - t0) + ' seconds.'
    
def findProposedOverlaps(proposed,length):
    prop = numpy.zeros(length)
    for i in range(len(proposed)):
        prop[proposed[i][0]:proposed[i][1]] += numpy.ones(proposed[i][1]-proposed[i][0])
    return len(numpy.where(prop > 2)[0])
    
def getHitOut(actual,proposed,length):
    act = numpy.zeros(length)
    prop = numpy.zeros(length)
    for i in range(len(actual)):
        act[actual[i][0]:actual[i][1]] = numpy.ones(actual[i][1]-actual[i][0])
    for i in range(len(proposed)):
        prop[proposed[i][0]:proposed[i][1]] = numpy.ones(proposed[i][1]-proposed[i][0])
        
    propOverlap = findProposedOverlaps(proposed,length)
    
    overlap = len(numpy.where((act == prop)*(act==1))[0])
    outer = length - len(numpy.where(act)[0]) + propOverlap
    totalprop = len(numpy.where((prop==1))[0])
    
    ratioIn = overlap/float(length - outer)
    ratioOut = (totalprop-overlap)/float(outer)
    #print actual, proposed, ratioIn, ratioOut
    return ratioIn,ratioOut
    
def plotBenchMulti(aligns,seqDat,blaster,hmmer):
    # read in data from VR
    pklAlg = open(aligns,'rb')
    [setOfAligns] = pickle.load(pklAlg)
    pklAlg.close()
    pklDat = open(seqDat,'rb')
    [[soiLen,lenStd,subLen,subStd,nSub,minLenSeq,minLenSub],tsInfo] = pickle.load(pklDat)
    pklDat.close()
    
    outdir = 'benchmarks/benchmark' + str(soiLen)
    
    # read in data from HMMer
    hmmin = open(hmmer,'rb')
    setter = False
    start = -1
    stop = -1
    prop = []
    hin = [0]*10
    hout = [0]*10
    numTrack = numpy.zeros(nSub*10)
    for line in hmmin:
        if line.strip() == 'Internal pipeline statistics summary:':
            if setter:
                prop.append([start,stop])
                actual = []
                if len(numpy.shape(tsInfo[number][4])) > 0:
                    for i in range(len(tsInfo[number][4])):
                        actual.append([tsInfo[number][4][i],tsInfo[number][4][i] + tsInfo[number][2][i]])
                else:
                    actual = [[tsInfo[number][4],tsInfo[number][4] + tsInfo[number][2]]]
                ratioIn,ratioOut = getHitOut(actual,prop,soiLen)
                #if number < 10:
                    #print number, start, stop, (min(tsInfo[number][4] + tsInfo[number][2],stop) - max(tsInfo[number][4],start))/float(tsInfo[number][2]), ((stop-start)-(min(tsInfo[number][4] + tsInfo[number][2],stop) - max(tsInfo[number][4],start)))/float(soiLen-tsInfo[number][2])
                #print 'posted phmmer'
                hin[number/nSub] += ratioIn
                hout[number/nSub] += ratioOut
            break
        if line[0:2] == '>>':
            if setter:
                prop.append([start,stop])
                actual = []
                if len(numpy.shape(tsInfo[number][4])) > 0:
                    for i in range(len(tsInfo[number][4])):
                        actual.append([tsInfo[number][4][i],tsInfo[number][4][i] + tsInfo[number][2][i]])
                else:
                    actual = [[tsInfo[number][4],tsInfo[number][4] + tsInfo[number][2]]]
                ratioIn,ratioOut = getHitOut(actual,prop,soiLen)
                #if number < 10:
                    #print number, start, stop, (min(tsInfo[number][4] + tsInfo[number][2],stop) - max(tsInfo[number][4],start))/float(tsInfo[number][2]), ((stop-start)-(min(tsInfo[number][4] + tsInfo[number][2],stop) - max(tsInfo[number][4],start)))/float(soiLen-tsInfo[number][2])
                #print 'posted phmmer'
                hin[number/nSub] += ratioIn
                hout[number/nSub] += ratioOut
            if line[0:5] == '>>END':
                break
            number = line[3:].strip()
            number = number.split('ts')[1]
            number = int(number.split('nsub')[0])
            numTrack[number] += 1
            if numTrack[number] > 1:
                print 'WARNING'
            start = 0
            stop = 0
            prop = []
            setter = False
        if start != 0:
            if line[2:11] == '== domain':
                prop.append([start,stop])
                start = 0
                stop = 0
        if line.strip().split('_')[0] == 'SeqOfInt':
            #print line.split(' ')
            #seqStrLen = len(line.split(' ')[2])
            #print seqStrLen
            #print line[2+seqStrLen:2+seqStrLen+4].strip()
            #print int(line[2+seqStrLen:2+seqStrLen+4].strip())
            if start == 0:
                liner = line.split(' ')
                #print liner
                for i in range(3,len(liner)):
                    if liner[i] != '':
                        start = int(liner[i])
                        break
            #print line
            #print line.strip().split(' ')[-1]
            stop = int(line.strip().split(' ')[-1])
            setter = True
                
    for i in range(10):
        hin[i] = hin[i]/float(nSub)
        hout[i] = 1 - hout[i]/float(nSub)
    
    # read in data from BLAST
    blastin = open(blaster,'rb')
    bin = [0]*10
    bout = [0]*10
    start = -1
    prop = []
    setter = False
    
    for line in blastin:
        if line.strip() == 'Lambda     K      H':
            break
        if line[0] == '>':
            if setter:
                prop.append([start,stop])
                actual = []
                if len(numpy.shape(tsInfo[number][4])) > 0:
                    for i in range(len(tsInfo[number][4])):
                        actual.append([tsInfo[number][4][i],tsInfo[number][4][i] + tsInfo[number][2][i]])
                else:
                    actual = [[tsInfo[number][4],tsInfo[number][4] + tsInfo[number][2]]]
                ratioIn,ratioOut = getHitOut(actual,prop,soiLen)
                #if number < 100:
                #print number
                    #print start,stop,tsInfo[number][4],tsInfo[number][4] + tsInfo[number][2]
                    #print (min(tsInfo[number][4] + tsInfo[number][2],stop) - max(tsInfo[number][4],start))/float(tsInfo[number][2])
                    #print ((stop-start)-(min(tsInfo[number][4] + tsInfo[number][2],stop) - max(tsInfo[number][4],start)))/float(soiLen-tsInfo[number][2])
                #print 'posted BLASTp'
                bin[number/nSub] += ratioIn
                bout[number/nSub] += ratioOut
            if line[0:4] == '>END':
                break
            start = 0
            stop = 0
            prop = []
            setter = False
            number = line[2:-1].split('ts')[1]
            number = int(number.split('nsub')[0])
            #numTrack[number] += 1
            #print line
        if start != 0:
            if line[0:8] == ' Score =':
                prop.append([start,stop])
                start = 0
                stop = 0
        if line[0:5] == 'Query' and start == 0:
            #print line
            start = int(line[6:11].strip())
        if line[0:5] == 'Query' and start != -1:
            #print line
            stop = int(line[-5:].strip())
            setter = True
    for i in range(10):
        bin[i] = bin[i]/float(nSub)
        bout[i] = 1 - bout[i]/float(nSub)
        
    #plt.figure()
    #plt.plot(numTrack)
    
    # make VR data
    score = []
    maxer = []
    hit = []
    out = []
    for i in range(10):
        score.append(0)
        maxer.append(0)
        hit.append(0)
        out.append(0)
        for j in range(nSub):
            penalty = 0
            maxpenalty = tsInfo[nSub*i + j][2]
            if len(setOfAligns.aligns[nSub*i + j].localAligns) > 0:
                #if len(setOfAligns.aligns[nSub*i + j].localAligns) > 1:
                    #print i,j,len(setOfAligns.aligns[nSub*i + j].localAligns)
                percentCovered = 0
                prop = []
                for k in range(len(setOfAligns.aligns[nSub*i + j].localAligns)):
                    number = i*nSub + j
                    actual = []
                    if len(numpy.shape(tsInfo[number][4])) > 0:
                        for i in range(len(tsInfo[number][4])):
                            actual.append([tsInfo[number][4][i],tsInfo[number][4][i] + tsInfo[number][2][i]])
                    else:
                        actual = [[tsInfo[number][4],tsInfo[number][4] + tsInfo[number][2]]]
                    for l in range(len(setOfAligns.aligns[nSub*i + j].localAligns[k].SOIrange)):
                        start = setOfAligns.aligns[nSub*i + j].localAligns[k].SOIrange[l][0]
                        stop = setOfAligns.aligns[nSub*i + j].localAligns[k].SOIrange[l][1]
                        prop.append([start,stop])
                ratioIn,ratioOut = getHitOut(actual,prop,soiLen)
                #print start,stop,actual,ratioIn,ratioOut
                #print 'posted vrope'
                out[-1] += ratioOut
                hit[-1] += ratioIn
            if penalty == 0:
                penalty = tsInfo[nSub*i + j][2]
            score[-1] += penalty
            maxer[-1] += maxpenalty
        score[-1] = score[-1]/float(maxer[-1])
        hit[-1] = hit[-1]/float(nSub)
        out[-1] = 1 - out[-1]/float(nSub)
    
        
    # plot figures
    fig = plt.figure(figsize=(17,7))
    ax = fig.add_subplot(121)
    ax.set_title('Specificity')
    ax.plot(out,label="Vrope")
    ax.plot(bout,label="BLAST")
    ax.plot(hout,label="HMMer")
    ax.set_xticklabels([1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1])
    ax.set_xlabel('Expected percent identity')
    ax.set_ylim(0,1)
    ax.set_yticks([0.0,0.25,0.5,0.75,1.0])
    ax.set_yticklabels([0.0,0.25,0.5,0.75,1.0])
    ax.set_ylabel('Average Percent of Non-Alignment Correctly Identified')
    ax.legend(loc = 4)
    ax = fig.add_subplot(122)
    ax.set_title('Sensitivity')
    ax.plot(hit,label="Vrope")
    ax.plot(bin,label="BLAST")
    ax.plot(hin,label="HMMer")
    ax.set_xticklabels([1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1])
    ax.set_xlabel('Expected percent identity')
    ax.set_ylim(0,1)
    ax.set_yticks([0.0,0.25,0.5,0.75,1.0])
    ax.set_yticklabels([0.0,0.25,0.5,0.75,1.0])
    ax.set_ylabel('Average Percent of Alignment Correctly Identified')
    ax.legend(loc = 3)
    fig.savefig(outdir + '/STATS_OUT.pdf')
    fig.savefig(outdir + '/STATS_OUT.png', bbox_inches='tight', pad_inches=0.03)
    plt.show()
    
# add error checking
if __name__ == '__main__':
    soiLen = int(sys.argv[1])
    lenStd = numpy.sqrt(soiLen)
    subLen = int(sys.argv[2])
    subStd = numpy.sqrt(subLen)
    nSub = int(sys.argv[3])
    minLenSeq = soiLen/2
    minLenSub = subLen/2
    BuildSeqs(soiLen,lenStd,subLen,subStd,nSub,minLenSeq,minLenSub)