#import matplotlib
#matplotlib.rcParams['font.size'] = 6
import matplotlib.pylab as plt
import numpy, logging
import cookb_signalsmooth as smooth
from PARAMS import *
from VelvetropeClasses import *

# Given a vector find out how many of each amino acid are in it
def GetVecComp(vec):
    if len(vec) > 1:
        type1 = numpy.zeros(max(max(vec)+1,21),dtype=numpy.int16)
    else:
        type1 = numpy.zeros(21)
    for i in range(len(vec)):
        type1[vec[i]] += 1
    return type1
    
# Given two compositions find the expected number of matches
def GetNumMatches(type1,type2):
    total = 0
    for i in range(min(len(type1),len(type2))):
        total += type1[i]*type2[i]
    return total

def generatePlotData(setOfAligns,parameters):
    offset = -1
    for align in setOfAligns.aligns:
        lensoi = align.seqOfInt.length
        lents = align.testSeq.length
        origMatchVec = numpy.zeros(lensoi)
        
        matchVec = numpy.zeros(lensoi)
        
        for locA in align.localAligns:
            
            for i in range(len(locA.SOIrange)):
                soiStart = locA.fSOIrange[0]
                soiEnd = locA.fSOIrange[1]
                tsStart = locA.fTSrange[0]
                tsEnd = locA.fTSrange[1]
                soiLstop = locA.SOIrange[i][0]
                soiRstop = locA.SOIrange[i][1]
                tsLstop = locA.TSrange[i][0]
                tsRstop = locA.TSrange[i][1]
            
                if align.testSeq.readingFrames == 1:
                    soi = numpy.array(align.seqOfInt.seq,dtype = numpy.int32)
                    ts = numpy.array(align.testSeq.seq,dtype = numpy.int32)
                else:
                    soi = numpy.array(align.seqOfInt.seq[0],dtype = numpy.int32)
                    ts = numpy.array(align.testSeq.seq[locA.offset],dtype = numpy.int32)
            
                fullScoreVec = numpy.zeros(lensoi)
                matchVecIn = numpy.zeros(lensoi)
            
                if align.seqOfInt.dataDepths[0] != 0:
                    fullScoreVec[soiStart:soiEnd] = numpy.asarray(soi[soiStart:soiEnd]%2**align.seqOfInt.dataDepths[0] == ts[tsStart:tsEnd]%2**align.seqOfInt.dataDepths[0],dtype = numpy.int32)
                
                    matchVecIn[soiLstop:soiRstop] = numpy.asarray(soi[soiLstop:soiRstop]%2**align.seqOfInt.dataDepths[0] == ts[tsLstop:tsRstop]%2**align.seqOfInt.dataDepths[0],dtype = numpy.int32)
                
                    soiVecComp = GetVecComp(soi[soiStart:soiEnd]%2**align.seqOfInt.dataDepths[0])
                    tsVecComp = GetVecComp(ts[tsStart:tsEnd]%2**align.seqOfInt.dataDepths[0])
                else:
                    fullScoreVec[soiStart:soiEnd] = numpy.asarray(soi[soiStart:soiEnd] == ts[tsStart:tsEnd],dtype = numpy.int32)
                
                    matchVecIn[soiLstop:soiRstop] = numpy.asarray(soi[soiLstop:soiRstop] == ts[tsLstop:tsRstop],dtype = numpy.int32)
                
                    soiVecComp = GetVecComp(soi[soiStart:soiEnd])
                    tsVecComp = GetVecComp(ts[tsStart:tsEnd])
            
                if align.seqOfInt.dataDepths[0] != 0:
                    colScore = numpy.zeros((len(align.seqOfInt.dataDepths),lensoi))
                else:
                    colScore = numpy.zeros(lensoi)
            
                # increments club
                if align.seqOfInt.dataDepths[0] == 0:
                    locA.score += matchVecIn
                else:
                    colScore[0] = matchVecIn
                    for j in range(len(align.seqOfInt.dataDepths)-1):
                        colScore[1+j][soiLstop:soiRstop] = matchVecIn[soiLstop:soiRstop]*numpy.asarray(soi[soiLstop:soiRstop]%2**align.seqOfInt.dataDepths[j+1] == ts[tsLstop:tsRstop]%2**align.seqOfInt.dataDepths[j+1],dtype = numpy.int32)
                    locA.score += colScore
                
            
                
            
                # set cdfL/cdfR, only technically needs to be done once
                if i == 0:
                    locA.fullScoreVec += fullScoreVec
                    
                    cdfL = numpy.zeros(lensoi)
                    cdfR = numpy.zeros(lensoi)
            
                    numMatches = GetNumMatches(soiVecComp,tsVecComp)
            
                    singleexp = numMatches/(float(soiEnd-soiStart)*(soiEnd-soiStart))
            
                    correctionL = numpy.cumsum(numpy.ones(soiEnd-soiStart)*singleexp)
                    cdfL[soiStart:soiEnd] = numpy.cumsum(fullScoreVec[soiStart:soiEnd]) - correctionL
            
                    cdfR[soiStart:soiEnd] = numpy.asarray(numpy.cumsum(fullScoreVec[soiStart:soiEnd][::-1]) - correctionL)[::-1]
                
                    locA.cdfL = cdfL
                    locA.cdfR = cdfR
            
                #print numpy.shape(align.score[:,locA.SOIrange[0]:locA.SOIrange[1]])
                #print numpy.shape(locA.score)
                
                align.inClub[locA.SOIrange[i][0]:locA.SOIrange[i][1]] += numpy.ones((locA.SOIrange[i][1]-locA.SOIrange[i][0]))
                
            if align.seqOfInt.dataDepths[0] != 0:
                align.score += locA.score
            else:
                align.score += locA.score
                
            
            
        setOfAligns.inClub += align.inClub
        setOfAligns.scoreSum += align.score
            
    return setOfAligns

def RegExpOut(alignments, parameters):
    numGenes = len(alignments.testSeqs)/alignments.seqOfInt.readingFrames
    aligns = []
    if len(numpy.shape(alignments.scoreSum)) != 1:
        if alignments.seqOfInt.readingFrames > 1:
            length = len(alignments.seqOfInt.seq[0])
        else:
            length = len(alignments.seqOfInt.seq)
        for i in range(length):
            if alignments.scoreSum[0][i]/float(numGenes) > parameters['REG_EXP_RATIO']:
                aligns.append(1)
            else:
                aligns.append(0)
    else:
        if alignments.seqOfInt.readingFrames > 1:
            length = len(alignments.seqOfInt.seq[0])
        else:
            length = len(alignments.seqOfInt.seq)
        for i in range(length):
            if alignments.scoreSum[i]/float(numGenes) > parameters['REG_EXP_RATIO']:
                aligns.append(1)
            else:
                aligns.append(0)
    
    start = False
    newReg = ''
    regs = []
    for i in range(length):
        if start == False:
            newReg = str(i) + ' .*'
            gapVal = 0
            if aligns[i] == 1:
                if alignments.seqOfInt.readingFrames == 1:
                    if alignments.seqOfInt.dataDepths[0] != 0:
                        newReg = newReg + FIVEBIT_TO_AA[alignments.seqOfInt.seq[i]%2**5]
                    else:
                        newReg = newReg + FIVEBIT_TO_AA[alignments.seqOfInt.seq[i]]
                else:
                    newReg = newReg + FIVEBIT_TO_AA[alignments.seqOfInt.seq[0][i]%2**5]
                start = True
            else:
                pass
        else:
            if aligns[i] == 1:
                if alignments.seqOfInt.readingFrames == 1:
                    if alignments.seqOfInt.dataDepths[0] != 0:
                        newReg = newReg + FIVEBIT_TO_AA[alignments.seqOfInt.seq[i]%2**5]
                    else:
                        newReg = newReg + FIVEBIT_TO_AA[alignments.seqOfInt.seq[i]]
                else:
                    newReg = newReg + FIVEBIT_TO_AA[alignments.seqOfInt.seq[0][i]%2**5]
                gapVal = 0
            else:
                gapVal += 1
                if gapVal < parameters['REG_EXP_GAP']:
                    newReg = newReg + '.'
                else:
                    if len(newReg.split('*')[1]) > parameters['LOCAL_WINDOW']/2:
                        newReg = newReg[:-4]
                        newReg = newReg + '.* ' + str(i-parameters['REG_EXP_GAP'])
                        regs.append(newReg)
                    start = False
    if start == True:
        if len(newReg.split('*')[1]) > parameters['LOCAL_WINDOW']/2:
            newReg = newReg + '.*'
            regs.append(newReg)
                    
    logging.debug('Regular expressions of alignments...')
    
    for i in range(len(regs)):
        logging.debug(regs[i])
        
    return regs
                    

def GenOutput(alignments, parameters):
    # Read in the data from file
    
    #pkl_file = open(inputFile,'rb')
    #[colVecs,colVecsIn,offset,offVal,dist,bps,sigDif,locfilter,filterParam,data,decoder,lstopers,rstopers,runLengthDist] = pickle.load(pkl_file)
    
    colVecs = alignments.inClubSum
    colVecsIn = alignments.scoreSum

    sigDif = parameters['GLOBAL_SIG_LEVEL']
    locfilter = 'CDF'
    
    soiInfo = alignments.seqOfInt.name
    if alignments.seqOfInt.readingFrames == 1:
        soi = alignments.seqOfInt.seq
    else:
        soi = alignments.seqOfInt.seq[0]
    testSeq = []
    testSeqInfo = []
    for ts in alignments.testSeqs:
        testSeq.append(ts.seq)
        testSeqInfo.append(ts.name)
    soi = numpy.array(soi,dtype = numpy.int32)
    
    lstopers = numpy.zeros(alignments.seqOfInt.length)
    rstopers = numpy.zeros(alignments.seqOfInt.length)
    #for alignment in alignments.aligns:
        #if len(alignment.localAligns)>0:
            #for locA in alignment.localAligns:
                #lstopers[locA.SOIrange[0]] += 1
                #rstopers[locA.SOIrange[1]] += 1
    
    ratio = parameters['REG_EXP_RATIO']
    numGenes = len(alignments.testSeqs)/alignments.seqOfInt.readingFrames
    #print numGenes
        
    conservedFrac = smooth.smooth(colVecsIn[0]/numpy.asarray(colVecs+0.01,dtype=numpy.float),window_len=11,window='flat')
    if len(numpy.shape(colVecsIn)) != 1:
        pass
        #silentMutations = smooth.smooth(colVecsIn[1]/numpy.asarray(colVecsIn[0]+0.01,dtype=numpy.float),window_len=11,window='flat')
        
        ##silentMutations2 = smooth.smooth(numpy.arcsin(numpy.sqrt(colVecsIn[1]/numpy.asarray(colVecsIn[0]+0.01,dtype=numpy.float))),window_len=11,window='flat') # arcsin transform
        #print numpy.mean(silentMutations)
        #print numpy.mean(silentMutations[numpy.where(silentMutations > 0.01)])
        #print numpy.mean(conservedFrac)
        #print numpy.mean(conservedFrac[numpy.where(conservedFrac > 0.01)])
    
    #print max(colVecs)
    #print max(colVecsIn)
    
    # decode the soi
    genomelist = []
    for i in range(len(soi)):
        if alignments.seqOfInt.sequenceType == 'AA':
            genomelist.append(FIVEBIT_TO_AA[soi[i]])
        elif alignments.seqOfInt.sequenceType == 'DNA':
            genomelist.append(FIVEBIT_TO_AA[soi[i]%2**5])
    #print 'Making run length distribution...'
    
    #fig = plt.figure()
    #y,x = GetRunLengthDistribution(numpy.array(offVal))
    #plt.loglog(x,y,align='center')
    #plt.title('Run Length Distribution after global (sig > ' + str(sigDif) + ') and local ' + locfilter + ' filtering')
    #plt.xlabel('Run Lengths')
    #plt.ylabel('Instances')
    
    #maxRun = max(runLengthDist)
    #plt.figure()
    #runLenLog = numpy.zeros(maxRun+1)
    #for i in range(len(runLengthDist)):
        #runLenLog[runLengthDist[i]] += 1
    #plt.loglog(range(0,maxRun+1),runLenLog)
    #t1 = plt.loglog(range(1,maxRun+1),runLenLog[1]*numpy.arange(1.0,maxRun+1)**-1.5)
    #t2 = plt.loglog(range(1,maxRun+1),runLenLog[1]*numpy.arange(1.0,maxRun+1)**-2.0)
    #plt.title('Run Length Distribution ' + str(numGenes) + ' Bact RpoBs')
    #plt.legend((t1,t2),('l^{-1.5}','l^{-2.0}'))
    
    logging.info('Making histograms...')

    
    
    residuesPerGraph = 400
    residuesPerLine = 100
    numGraphs = len(soi)/residuesPerGraph + 1
    for i in range(numGraphs):
        fig = plt.figure(figsize=(22,17))
        #fig.text(0.5,0.95,'InClub and InSetInClub matches across all alignments.\nGlobal filter of sigma > ' + str(sigDif) + ' then local filter of 5 in 8.',horizontalalignment='center')
        fig.text(0.5,0.92,'InClub and Score matches across all alignments.\nGlobal filter of sigma > ' + str(sigDif) + ' then local filter of ' + str(locfilter) + '.\nTotal genes compared: '+ str(numGenes)+'. Reference Gene: ' + str(soiInfo),horizontalalignment='center')
        fig.text(0.5,0.05,'Sequence compared against, red text indicates at least ratio=' + str(ratio) + ' InSetInClub\nSequence from ' + str(i*residuesPerGraph+1) + ' to ' + str(min((i+1)*residuesPerGraph,len(soi))) + ' of ' + str(len(soi)) + '.',horizontalalignment='center')
        fig.text(0.05,0.5,'Frequency of matches across all other genes, Score = blue and InClub = green and further in shades of red',rotation=90,verticalalignment='center')
        ax = []
        for j in numpy.arange(0,8,2):
            if len(soi) > residuesPerGraph*i + (residuesPerLine/2)*(j+2):
                ax.append(fig.add_subplot(int('81'+str(j+1))))
                fig.text(0.95,0.82-j*0.1,'Residues:\n'+str(i*residuesPerGraph+(residuesPerLine/2)*j+1)+'-'+str(i*residuesPerGraph+(residuesPerLine/2)*j+residuesPerLine),horizontalalignment='center')
                fig.text(0.95,0.82-(j+1)*0.1,'Cons. Ident.\n Start(m) and End(c)',horizontalalignment='center')
                ax[j].bar(numpy.arange(residuesPerLine),colVecs[(residuesPerLine/2)*j+residuesPerGraph*i:(residuesPerLine/2)*(j+2)+residuesPerGraph*i]/float(numGenes),color='green',alpha=0.25,align='center')
                #print colVecs[(residuesPerLine/2)*j+residuesPerGraph*i:(residuesPerLine/2)*(j+2)+residuesPerGraph*i]/float(numGenes)
                if len(numpy.shape(colVecsIn)) == 1:
                    ax[j].bar(numpy.arange(residuesPerLine),colVecsIn[(residuesPerLine/2)*j+residuesPerGraph*i:(residuesPerLine/2)*(j+2)+residuesPerGraph*i]/float(numGenes),align='center',color='blue')
                else:
                    #print numpy.shape(colVecsIn)
                    ax[j].bar(numpy.arange(residuesPerLine),colVecsIn[0][(residuesPerLine/2)*j+residuesPerGraph*i:(residuesPerLine/2)*(j+2)+residuesPerGraph*i]/float(numGenes),align='center',color='blue')
                    for k in range(1,numpy.shape(colVecsIn)[0]):
                        ax[j].bar(numpy.arange(residuesPerLine),colVecsIn[k][(residuesPerLine/2)*j+residuesPerGraph*i:(residuesPerLine/2)*(j+2)+residuesPerGraph*i]/float(numGenes),align='center',color='red')
                #ax[j].bar(numpy.arange(residuesPerLine),colVecs[(residuesPerLine/2)*j+residuesPerGraph*i:(residuesPerLine/2)*(j+2)+residuesPerGraph*i]/float(numGenes),color='green',alpha=0.25,align='center')
                # NOT LOG
                #ax[j].bar(numpy.arange(residuesPerLine),colVecs[(residuesPerLine/2)*j+residuesPerGraph*i:(residuesPerLine/2)*(j+2)+residuesPerGraph*i]/float(numGenes),color='green',alpha=0.25,align='center')
                #if len(numpy.shape(colVecsIn)) == 1:
                    #ax[j].bar(numpy.arange(residuesPerLine),colVecsIn[(residuesPerLine/2)*j+residuesPerGraph*i:(residuesPerLine/2)*(j+2)+residuesPerGraph*i]/float(numGenes),align='center',color='blue')
                #else:
                    ##print numpy.shape(colVecsIn)
                    #ax[j].bar(numpy.arange(residuesPerLine),colVecsIn[0][(residuesPerLine/2)*j+residuesPerGraph*i:(residuesPerLine/2)*(j+2)+residuesPerGraph*i]/float(numGenes),align='center',color='blue')
                    #for k in range(1,numpy.shape(colVecsIn)[0]):
                        #ax[j].bar(numpy.arange(residuesPerLine),colVecsIn[k][(residuesPerLine/2)*j+residuesPerGraph*i:(residuesPerLine/2)*(j+2)+residuesPerGraph*i]/float(numGenes),align='center',color='red')
                
                    #ax[j].plot(numpy.arange(100),silentMutations2[residuesPerLine*j+800*i:residuesPerLine*(j+1)+800*i],color='magenta') # arcsin transform
                #ax[j].bar(numpy.arange(100),lstopers[residuesPerLine*j+800*i:residuesPerLine*(j+1)+800*i]/float(numGenes),color='magenta',alpha=0.5,align='center')
                #ax[j].bar(numpy.arange(100),rstopers[residuesPerLine*j+800*i:residuesPerLine*(j+1)+800*i]/float(numGenes),color='black',alpha=0.5,align='center')
                ax[j].set_ylim(0,1)
                ax[j].set_xlim(0,residuesPerLine-1)
                ax[j].set_xticks(numpy.arange(residuesPerLine))
                ax[j].set_xticklabels(genomelist[(residuesPerLine/2)*j+residuesPerGraph*i:(residuesPerLine/2)*(j+2)+residuesPerGraph*i])
                ax[j].set_yticks([0.0,0.5,1.0])
                ax[j].set_yticklabels([0.0,0.5,1.0])
                k = 0
                for lab in ax[j].xaxis.get_ticklabels():
                    if len(numpy.shape(colVecsIn)) != 1:
                        if colVecsIn[0][(residuesPerLine/2)*j+residuesPerGraph*i+k]/float(numGenes) > ratio:
                            lab.set_color('red')
                        k += 1
                    else:
                        if colVecsIn[(residuesPerLine/2)*j+residuesPerGraph*i+k]/float(numGenes) > ratio:
                            lab.set_color('red')
                        k += 1
                ax.append(fig.add_subplot(int('81'+str(j+2))))
                ax[j+1].plot(numpy.arange(residuesPerLine),lstopers[(residuesPerLine/2)*j+residuesPerGraph*i:(residuesPerLine/2)*(j+2)+residuesPerGraph*i]/float(numGenes),color='magenta')
                ax[j+1].plot(numpy.arange(residuesPerLine),rstopers[(residuesPerLine/2)*j+residuesPerGraph*i:(residuesPerLine/2)*(j+2)+residuesPerGraph*i]/float(numGenes),color='cyan')
                ax[j+1].plot(numpy.arange(residuesPerLine),conservedFrac[(residuesPerLine/2)*j+residuesPerGraph*i:(residuesPerLine/2)*(j+2)+residuesPerGraph*i],color='black')
                if len(numpy.shape(colVecsIn)) != 1:
                    pass
                    #ax[j+1].plot(numpy.arange(residuesPerLine),silentMutations[(residuesPerLine/2)*j+residuesPerGraph*i:(residuesPerLine/2)*(j+2)+residuesPerGraph*i],color='red')
                ax[j+1].set_ylim(0,1)
                ax[j+1].set_xlim(0,residuesPerLine-1)
                ax[j+1].set_xticks(numpy.arange(residuesPerLine))
                ax[j+1].set_xticklabels(genomelist[(residuesPerLine/2)*j+residuesPerGraph*i:(residuesPerLine/2)*(j+2)+residuesPerGraph*i])
                ax[j+1].set_yticks([0.0,0.5,1.0])
                ax[j+1].set_yticklabels([0.0,0.5,1.0])
                k = 0
                for lab in ax[j+1].xaxis.get_ticklabels():
                    if lstopers[(residuesPerLine/2)*j+residuesPerGraph*i+k] + rstopers[(residuesPerLine/2)*j+residuesPerGraph*i+k] > 0:
                        lab.set_color('red')
                    k += 1
            else:
                ax.append(fig.add_subplot(int('81'+str(j+1))))
                tmplist1 = []
                tmplist2 = colVecsIn.transpose()[(residuesPerLine/2)*j+residuesPerGraph*i:len(soi)].transpose()
                tmplist3 = colVecs[(residuesPerLine/2)*j+residuesPerGraph*i:len(soi)]
                tmplist1.extend(genomelist[(residuesPerLine/2)*j+residuesPerGraph*i:len(soi)])
                #print tmplist5
                #tmplist2.extend(list(colVecs[residuesPerLine*j+800*i:len(soi)]))
                #tmplist3.extend(list(colVecsIn[residuesPerLine*j+800*i:len(soi)]))
                fig.text(0.95,0.82-j*0.1,'Residue Numbers:\n'+str(i*residuesPerGraph+(residuesPerLine/2)*j+1)+'-'+str(len(soi)),horizontalalignment='center')
                for k in range(residuesPerLine-len(tmplist1)):
                    tmplist1.append(' ')
                    #print numpy.shape(colVecsIn)
                    #print numpy.shape(tmplist2)
                    if len(numpy.shape(colVecsIn)) != 1:
                        #tmplist2 = numpy.append(tmplist2,numpy.zeros((numpy.shape(colVecsIn)[0],1)),1)
                        tmplist2 = numpy.append(tmplist2,numpy.zeros((numpy.shape(colVecsIn)[0],1)),1)
                    else:
                        tmplist2 = numpy.append(tmplist2,[0],1)
                    tmplist3 = numpy.append(tmplist3,[0],1)
                #print tmplist2[0][0:residuesPerLine]
                #print tmplist5
                ax[j].bar(numpy.arange(residuesPerLine),(numpy.array(tmplist3[0:residuesPerLine])/float(numGenes)),color='green',alpha=0.25,align='center')
                if len(numpy.shape(colVecsIn)) == 1:
                    ax[j].bar(numpy.arange(residuesPerLine),(numpy.array(tmplist2[0:residuesPerLine])/float(numGenes)),align='center', color='blue')
                else:
                    ax[j].bar(numpy.arange(residuesPerLine),(numpy.array(tmplist2[0][0:residuesPerLine])/float(numGenes)),align='center', color='blue')
                    for k in range(1,numpy.shape(colVecsIn)[0]):
                        ax[j].bar(numpy.arange(residuesPerLine),(numpy.array(tmplist2[k][0:residuesPerLine])/float(numGenes)),align='center', color='red')
                ## NOT LOG
                #ax[j].bar(numpy.arange(residuesPerLine),numpy.array(tmplist3[0:residuesPerLine])/float(numGenes),color='green',alpha=0.25,align='center')
                #if len(numpy.shape(colVecsIn)) == 1:
                    #ax[j].bar(numpy.arange(residuesPerLine),numpy.array(tmplist2[0:residuesPerLine])/float(numGenes),align='center', color='blue')
                #else:
                    #ax[j].bar(numpy.arange(residuesPerLine),numpy.array(tmplist2[0][0:residuesPerLine])/float(numGenes),align='center', color='blue')
                    #for k in range(1,numpy.shape(colVecsIn)[0]):
                        #ax[j].bar(numpy.arange(residuesPerLine),numpy.array(tmplist2[k][0:residuesPerLine])/float(numGenes),align='center', color='red')
                ax[j].set_ylim(0,1)
                ax[j].set_xlim(0,residuesPerLine-1)
                ax[j].set_xticks(numpy.arange(residuesPerLine))
                ax[j].set_xticklabels(tmplist1[0:residuesPerLine])
                ax[j].set_yticks([0.0,0.5,1.0])
                ax[j].set_yticklabels([0.0,0.5,1.0])
                k = 0
                for lab in ax[j].xaxis.get_ticklabels():
                    if len(numpy.shape(colVecsIn)) != 1:
                        if tmplist2[0][k]/float(numGenes) > ratio:
                            lab.set_color('red')
                        k += 1
                    else:
                        if tmplist2[k]/float(numGenes) > ratio:
                            lab.set_color('red')
                        k += 1
                ax.append(fig.add_subplot(int('81'+str(j+2))))
                ax[j+1].plot(numpy.arange(len(lstopers[(residuesPerLine/2)*j+residuesPerGraph*i:])),lstopers[(residuesPerLine/2)*j+residuesPerGraph*i:]/float(numGenes),color='magenta')
                ax[j+1].plot(numpy.arange(len(rstopers[(residuesPerLine/2)*j+residuesPerGraph*i:])),rstopers[(residuesPerLine/2)*j+residuesPerGraph*i:]/float(numGenes),color='cyan')
                ax[j+1].plot(numpy.arange(len(conservedFrac[(residuesPerLine/2)*j+residuesPerGraph*i:])),conservedFrac[(residuesPerLine/2)*j+residuesPerGraph*i:],color='black')
                if len(numpy.shape(colVecsIn)) != 1:
                    pass
                    #ax[j+1].plot(numpy.arange(len(silentMutations[(residuesPerLine/2)*j+residuesPerGraph*i:])),silentMutations[(residuesPerLine/2)*j+residuesPerGraph*i:],color='red')
                ax[j+1].set_ylim(0,1)
                ax[j+1].set_xlim(0,residuesPerLine-1)
                ax[j+1].set_xticks(numpy.arange(residuesPerLine))
                ax[j+1].set_xticklabels(genomelist[(residuesPerLine/2)*j+residuesPerGraph*i:(residuesPerLine/2)*(j+2)+residuesPerGraph*i])
                ax[j+1].set_yticks([0.0,0.5,1.0])
                ax[j+1].set_yticklabels([0.0,0.5,1.0])
                k = 0
                for lab in ax[j+1].xaxis.get_ticklabels():
                    if k < len(lstopers[(residuesPerLine/2)*j+residuesPerGraph*i:]):
                        if lstopers[(residuesPerLine/2)*j+residuesPerGraph*i+k] + rstopers[(residuesPerLine/2)*j+residuesPerGraph*i+k] > 0:
                            lab.set_color('red')
                    k += 1
                j = 8
                break
        fig.savefig(parameters['OUTFILE'] + '/images/COMBO' + str(i+1) + '.pdf')
        fig.savefig(parameters['OUTFILE'] + '/images/COMBO' + str(i+1) + '.png', bbox_inches='tight', pad_inches=0.03)
    #plt.show()
    #plt.figure()
    #n, bins, patches = plt.hist(rowList,100,normed=1)
    #rowMean = numpy.mean(rowList)
    #rowStd = numpy.std(rowList)
    #plt.axvspan(rowMean-rowStd,rowMean+rowStd,facecolor='g',alpha=0.1)
    #plt.axvspan(rowMean-expStd,rowMean+expStd,facecolor='r',alpha=0.1)
    #plt.axvline(x=expMatch,color='r')
    #plt.axvline(x=rowMean,color='b')
    #plt.title(r'$\mathrm{Histogram\ of\ Row\ Matches:}\ \mu=\mathrm{Red\ line.},\ \sigma_{true}=\mathrm{Green\ bar.},\ \sigma_{approx}=\mathrm{Red\ bar.}$')
    ## add a 'best fit' line
    #y = mlab.normpdf( bins, rowMean, rowStd)
    #l = plt.plot(bins, y, 'r--', linewidth=1)
    #plt.xlabel('Matches per row')
    #plt.ylabel('Frequency')
    #plt.grid(True)

def GetRunLengthDistribution(rows):
    # add a zero column to the end
    runLens = []
    #print rows
    runMat = numpy.append(rows,numpy.reshape(numpy.zeros(numpy.shape(rows)[0]),(numpy.shape(rows)[0],1)),axis=1)
    #print runMat
    #print sum(runMat)
    # Start looking for runs of length 3 only
    runMat = runMat*numpy.roll(runMat,1)
    while(len(numpy.where(runMat)[0]) > 0):
        runMat = runMat*numpy.roll(runMat,1)
        runs = len(numpy.where(runMat)[0])
        runLens.append(runs)
    #print runLens
    for i in range(len(runLens)-1):
        for j in range(len(runLens)-1-i):
            runLens[-j-2-i] -= runLens[-i-1]*(j+2)
            #print -j-2-i
            #print runLens[-j-2-i]
            #print -i-1, (j+1)
            #print runLens[-i-1]*(j+1)
            #print runLens
    runLensX = numpy.arange(3,len(runLens)+2)
    return runLens[0:-1], runLensX


def makeIt():
    # makes the pretty latex tabular and prints matrix, emphasises edges, make sure to put a '-' at end
    #A = ['S','A','T','A','N','S','S','P','A','W','N','-','S','A','T','A','N','S','P','R','A','W','N','-']
    A = ['S','A','T','A','N','-','S','P','A','W','N','-','P','A','T','A','N','S','-','P','R','A','W','N','-']
    theMat, temp = MakeMatrixFull(A)
    string = "\\begin{tabular}{|"
    for i in range(len(A)-1):
        string = string + 'c'
    string = string + '|}'
    print string
    print '\hline'
    string = ""
    for i in range(len(A)-2):
        string = string + A[i] + " & "
    string = string + A[-2] + " \\\\"
    print string
    print '\hline'
    for i in range(len(A)-1):
        string = ""
        for j in range(len(A)):
            if theMat[i+1][j] == 0:
                if A[(j+i+1)%len(A)] == '-':
                    if j < len(A) - 1:
                        string = string + " \\textcolor{blue}{$\star$} "
                    theMat[i+1][j] = -1
                else:
                    if A[j] == '-':
                        if j < len(A) - 1:
                            string = string + " \\textcolor{blue}{$\star$} "
                        theMat[i+1][j] = -1
                    else:
                        string = string + " " + A[(j+i+1)%len(A)] + " "
            else:
                string = string + " \\textcolor{red}{" + A[(j+i+1)%len(A)] + "} "
            if j < len(A) - 2:
                string = string + "&"
            else:
                string = string + "\\"
        print string
    print '\hline'
    print '\end{tabular}'
    plt.matshow(theMat)
    
def makeItt2():
    # makes the pretty latex tabular and prints matrix, emphasises edges, make sure to put a '-' at end
    #A = ['S','A','T','A','N','S','S','P','A','W','N','-','S','A','T','A','N','S','P','R','A','W','N','-']
    #A = ['S','A','T','A','N','-','S','P','A','W','N','-','P','A','T','A','N','S','-','P','R','A','W','N','-']
    A = ['V','E','L','V','E','T','R','Q','P','E','A','L','G','Q','R','I','T','H','M','R','Q','C','K','S']
    B = ['P','E','A','T','C','Q','A','T','V','E','L','V','E','T','R','Q','C','K','S','A','L','L','Q','R','I','T','H','A','M','Q','C','K','S','L','V','A','R','P','Q','E','V','T','Y','M','E','T']
    theMat, temp = MakeMatrixFull(A)
    string = "\\begin{tabular}{|"
    for i in range(len(A)-1):
        string = string + 'c'
    string = string + '|}'
    print string
    print '\hline'
    string = ""
    for i in range(len(A)-2):
        string = string + A[i] + " & "
    string = string + A[-2] + " \\\\"
    print string
    print '\hline'
    for i in range(len(A)-1):
        string = ""
        for j in range(len(A)):
            if theMat[i+1][j] == 0:
                if A[(j+i+1)%len(A)] == '-':
                    if j < len(A) - 1:
                        string = string + " \\textcolor{blue}{$\star$} "
                    theMat[i+1][j] = -1
                else:
                    if A[j] == '-':
                        if j < len(A) - 1:
                            string = string + " \\textcolor{blue}{$\star$} "
                        theMat[i+1][j] = -1
                    else:
                        string = string + " " + A[(j+i+1)%len(A)] + " "
            else:
                string = string + " \\textcolor{red}{" + A[(j+i+1)%len(A)] + "} "
            if j < len(A) - 2:
                string = string + "&"
            else:
                string = string + "\\"
        print string
    print '\hline'
    print '\end{tabular}'
    plt.matshow(theMat)
        
def BitCount(num):
    count = 0
    while(num.any()):
        count += num&1
        num = num>>1
    return count
    
def makeIt3():
    """
    makeIt3(): Makes a latex picture showing what regions of the matrix are needed using len(geneLens) genes of relative sizes given
    """
    # size of the entire picture
    frameWidth = 500
    # border for colors
    eps = 5
    # Text info for G_1xG_2
    textWidth = 20
    textHeight = 10
    # Must add up to 1!
    geneLens = [.25,.25,.25,.25]
    # Start picture
    print '\\begin{picture}(' + str(frameWidth) + ',' + str(frameWidth) + ')'
    # draw outline
    print '% Draw outline'
    print '\put(0,0){\line(1,0){' + str(frameWidth) + '}}'
    print '\put(0,0){\line(0,1){' + str(frameWidth) + '}}'
    print '\put(' + str(frameWidth) + ',' + str(frameWidth) + '){\line(0,-1){' + str(frameWidth) + '}}'
    print '\put(' + str(frameWidth) + ',' + str(frameWidth) + '){\line(-1,0){' + str(frameWidth) + '}}'
    # draw vertical cuts
    print '% Draw vertical cuts'
    for i in range(len(geneLens)-1):
        print '\put(' + str(int(sum(geneLens[0:i+1])*frameWidth)) + ',0){\line(0,1){' + str(frameWidth) + '}}'
    # draw diagonal cuts
    print '% Draw diagonal cuts'
    for i in range(len(geneLens)):
        print '\put(0,'  + str(frameWidth - int(sum(geneLens[0:i+1])*frameWidth)) + '){\line(1,1){' + str(int(sum(geneLens[0:i+1])*frameWidth)) + '}}'
    for i in range(len(geneLens)-1):
        print '\put(' + str(int(sum(geneLens[0:i+1])*frameWidth)) + ',0){\line(1,1){' + str(frameWidth - int(sum(geneLens[0:i+1])*frameWidth)) + '}}'
    # draw interiors
    print '% draw interiors'
    # draw int upper triangles
    print '% draw int upper (blue) triangles'
    for i in range(len(geneLens)):
        print '\put(' + str(int(sum(geneLens[0:i])*frameWidth) + eps) + ',' + str(frameWidth - eps) + '){\\textcolor{blue}{\line(1,0){' + str(int(geneLens[i]*frameWidth) - 3*eps) + '}}}'
        print '\put(' + str(int(sum(geneLens[0:i])*frameWidth) + eps) + ',' + str(frameWidth - eps) + '){\\textcolor{blue}{\line(0,-1){' + str(int(geneLens[i]*frameWidth) - 3*eps) + '}}}'
        print '\put(' + str(int(sum(geneLens[0:i])*frameWidth) + eps) + ',' + str(frameWidth + 2*eps - int(geneLens[i]*frameWidth)) + '){\\textcolor{blue}{\line(1,1){' + str(int(geneLens[i]*frameWidth) - 3*eps) + '}}}'
        # put info in the box
        print '\put(' + str(int(.333*geneLens[i]*frameWidth+sum(geneLens[0:i])*frameWidth)-textWidth) + ',' + str(int((frameWidth - .333*int(geneLens[i]*frameWidth)))) + '){\\textcolor{blue}{$G_{' + str(i+1) + '} \\times G_{' + str(i+1) + '}$}}'
    # draw lower triangles
    print '% draw int lower (red) triangles'
    for i in range(len(geneLens)):
        # bottom
        print '\put(' + str(int(sum(geneLens[0:i])*frameWidth) + 2*eps) + ',' + str(eps) + '){\\textcolor{red}{\line(1,0){' + str(int(geneLens[i]*frameWidth) - 3*eps) + '}}}'
        # hypot
        print '\put(' + str(int(sum(geneLens[0:i])*frameWidth) + 2*eps) + ',' + str(eps) + '){\\textcolor{red}{\line(1,1){' + str(int(geneLens[i]*frameWidth) - 3*eps) + '}}}'
        # right
        print '\put(' + str(int(sum(geneLens[0:i+1])*frameWidth) - eps) + ',' + str(eps) + '){\\textcolor{red}{\line(0,1){' + str(int(geneLens[i]*frameWidth) - 3*eps) + '}}}'
        # put info in the box
        print '\put(' + str(int(.666*geneLens[i]*frameWidth+sum(geneLens[0:i])*frameWidth)-textWidth) + ',' + str(eps + .333*geneLens[i]*frameWidth - textHeight) + '){\\textcolor{red}{$G_{' + str(i+1) + '} \\times G_{' + str(i+1) + '}$}}'
    # draw green parrallels
    print '% draw upper (green) parrallels'
    for j in range(len(geneLens)-1):
        for i in range(j,len(geneLens)-1):
            #right
            print '\put(' + str(int(sum(geneLens[0:j])*frameWidth) + eps) + ',' + str(frameWidth - eps - sum(geneLens[j:i+1])*frameWidth) + '){\\textcolor{green}{\line(0,-1){' + str(int(geneLens[i+1]*frameWidth) - 4*eps) + '}}}'
            # left
            print '\put(' + str(int(sum(geneLens[0:j+1])*frameWidth) - eps) + ',' + str(frameWidth - 3*eps - sum(geneLens[j+1:i+1])*frameWidth) + '){\\textcolor{green}{\line(0,-1){' + str(int(geneLens[i+1]*frameWidth) - 4*eps) + '}}}'
            #up
            print '\put(' + str(int(sum(geneLens[0:j])*frameWidth) + eps) + ',' + str(frameWidth - eps - sum(geneLens[j:i+1])*frameWidth) + '){\\textcolor{green}{\line(1,1){' + str(int(geneLens[j]*frameWidth) - 2*eps) + '}}}'
            #down
            print '\put(' + str(int(sum(geneLens[0:j])*frameWidth) + eps) + ',' + str(frameWidth + 3*eps - sum(geneLens[j:i+1])*frameWidth - int(geneLens[i+1]*frameWidth)) + '){\\textcolor{green}{\line(1,1){' + str(int(geneLens[j]*frameWidth) - 2*eps) + '}}}'
            #put info in the box
            print '\put(' + str(int(sum(geneLens[0:j])*frameWidth) + 0.5*geneLens[j]*frameWidth - textWidth) + ',' + str(((frameWidth + 3*eps - sum(geneLens[j:i+1])*frameWidth - int(geneLens[i+1]*frameWidth))+(frameWidth - 3*eps - sum(geneLens[j+1:i+1])*frameWidth))/2 - textHeight) + '){\\textcolor{green}{$G_{' + str(j+1) + '} \\times G_{' + str(i+2) + '}$}}'
    # draw red parrellels
    print '% draw lower (red) parrallels'
    for j in range(len(geneLens)-1):
        for i in range(j,len(geneLens)-1):
            #left
            print '\put(' + str(int(sum(geneLens[0:i+1])*frameWidth) + eps) + ',' + str(sum(geneLens[j:i+1])*frameWidth - eps) + '){\\textcolor{red}{\line(0,-1){' + str(int(geneLens[j]*frameWidth) - 4*eps) + '}}}'
            # top
            print '\put(' + str(int(sum(geneLens[0:i+1])*frameWidth) + eps) + ',' + str(sum(geneLens[j:i+1])*frameWidth - eps) + '){\\textcolor{red}{\line(1,1){' + str(int(geneLens[i+1]*frameWidth) - 2*eps) + '}}}'
            # bottom
            print '\put(' + str(int(sum(geneLens[0:i+1])*frameWidth) + eps) + ',' + str(sum(geneLens[j:i+1])*frameWidth - eps - (int(geneLens[j]*frameWidth) - 4*eps)) + '){\\textcolor{red}{\line(1,1){' + str(int(geneLens[i+1]*frameWidth) - 2*eps) + '}}}'
            #right
            print '\put(' + str(int(sum(geneLens[0:i+2])*frameWidth) - eps) + ',' + str(sum(geneLens[j:i+2])*frameWidth - 3*eps) + '){\\textcolor{red}{\line(0,-1){' + str(int(geneLens[j]*frameWidth) - 4*eps) + '}}}'
            #put info in the box
            print '\put(' + str(int(sum(geneLens[0:i+1])*frameWidth + 0.5*geneLens[i+1]*frameWidth) - textWidth) + ',' + str(((sum(geneLens[j:i+1])*frameWidth - geneLens[j]*frameWidth) + (sum(geneLens[j:i+2])*frameWidth))/2 - textHeight) + '){\\textcolor{red}{$G_{' + str(i+2) + '} \\times G_{' + str(j+1) + '}$}}'
    print '\end{picture}'
    
def makeIt4(numRows):
    """
    makeIt3(): Makes a latex picture showing the wraparound scheme for fixed row length for numRows rows (try 5)
    """
    # size of the entire picture
    frameWidth = 300
    # border for rows
    eps = 10
    nudge = 4
    # number of rows to show
    if numRows == 0:
        numRows = 5
    #Start picture
    print '\\begin{picture}(' + str(frameWidth) + ',' + str(frameWidth) + ')'
    # draw outline
    print '% Draw outline'
    print '\put(0,0){\line(1,1){' + str(frameWidth/2) + '}}'
    print '\put(0,0){\line(0,1){' + str(frameWidth/2) + '}}'
    print '\put(0,' + str(frameWidth/2) + '){\line(1,1){' + str(frameWidth/2) + '}}'
    print '\put(' + str(frameWidth/2) + ',' + str(frameWidth/2) + '){\line(0,1){' + str(frameWidth/2) + '}}'
    # draw upper lines
    print '% draw upper lines'
    for i in range(numRows):
        print '\put(' + str(frameWidth/2) + ',' + str(frameWidth-eps*(i+1)) + '){\line(-1,0){' + str(eps*(i+1)) + '}}'
        print '\put(' + str(frameWidth/2-eps*(i+2)-2) + ',' + str(frameWidth-eps*(i+1)-nudge) + '){$' + str(i+1) + '$}'
    # draw lower lines
    print '% draw lower lines'
    for i in range(numRows):
        print '\put(0,' + str(frameWidth/2-eps*(i+1)) + '){\line(1,0){' + str(frameWidth/2 - eps*(i+1)) + '}}'
        print '\put(' + str(frameWidth/2 - eps*(i)-2) + ',' + str(frameWidth/2-eps*(i+1)-nudge) + '){$' + str(i+1) + '$}'
    print '\end{picture}'