import numpy #for arrays, vector operations
import time #timing...
import os
import logging
import intronFinder
from VelvetropeClasses import *
from PARAMS import *
import time
import math

###################
# Clubs Generator #
###################
"""
    This script makes all of the alignments using either python, C or CUDA
    it requires the seqOfInt and tesSeqs from the Clubs Reader
    it also requires the parameters dictionary
    
    Speed in order of fastest to slowest is CUDA, C, python
"""

def GenClubs(seqOfInt, testSeqs, parameters):
    if parameters['COMPUTE_TYPE'] == 'python':
        logging.info('Computing clubs using python (slow, try C or CUDA)...')
        logging.info('Currently outdated, uses old algorithm...')
        return GenClubsPython(seqOfInt, testSeqs, parameters)
    elif parameters['COMPUTE_TYPE'] == 'CUDA':
        logging.info('Computing clubs using CUDA...')
        return GenClubsCUDA(seqOfInt, testSeqs, parameters)
    elif parameters['COMPUTE_TYPE'] == 'C':
        logging.info('Computing clubs using C (for large datasets try CUDA)...')
        # check to make sure install was done correctly
        directoryContents = os.listdir('src')
        foundCfile = False
        for files in directoryContents:
            if files == 'clubgen_c.py':
                foundCfile = True
        if foundCfile == False:
            logging.critical('Could not find clubgen_c.py!!!')
            logging.critical('You need to run make in the src/ directory to compile the C club finder!')
            logging.critical('Please read the installation directions in the ReadMe file!')
            return 0
        else:
            return GenClubsC(seqOfInt, testSeqs, parameters)
    logging.critical('Could not recognize COMPUTE_TYPE!')
    return 0

# Find clubs using CUDA
def GenClubsCUDA(seqOfInt, testSeqs, parameters):
    import clubgen_pycuda
    t0 = time.time()
    alignmentInfo, soiLen, tsLen, counter = clubgen_pycuda.GenClubs(seqOfInt, testSeqs, parameters)
    logging.info('Time in cuda function: ' + str(time.time() - t0))
    
    #import matplotlib.pylab as plt
    #plt.plot(alignmentInfo)
    #plt.show()
    #return 0

    MAX_MATCH = 20
    
    alignments = SetOfAlignments(seqOfInt)
            
    i = 0
    soiOn = 0
    tsOn = 0
    seqOn = 0
    rframe = 0
    
    locWin = parameters['LOCAL_WINDOW']
    theAlignment = Alignment(seqOfInt, testSeqs[seqOn])
    
    
    
    while(i < MAX_MATCH*3):
        #print alignmentInfo[soiOn*MAX_MATCH*3+tsOn*3*MAX_MATCH*soiLen+i:soiOn*MAX_MATCH*3+tsOn*3*MAX_MATCH*soiLen+i+3]
        start = alignmentInfo[soiOn*MAX_MATCH*3+tsOn*3*MAX_MATCH*soiLen+i]
        stop = alignmentInfo[soiOn*MAX_MATCH*3+tsOn*3*MAX_MATCH*soiLen+i+1]
        offset = alignmentInfo[soiOn*MAX_MATCH*3+tsOn*3*MAX_MATCH*soiLen+i+2]
        if soiOn*MAX_MATCH*3+tsOn*3*MAX_MATCH*soiLen+i+5 < MAX_MATCH*3*tsLen*soiLen:
            nextOffset = alignmentInfo[soiOn*MAX_MATCH*3+tsOn*3*MAX_MATCH*soiLen+i+5]
        else:
            nextOffset = -1
        #print start,stop,offset,nextOffset
        if start != stop:
            # this part of the alignmentInfo contains useful information!
            # logging.debug('alignmentInfo: ' + str(start) + ' ' + str(stop) + ' ' + str(offset))
            lensoi = min(1000,seqOfInt.length-1000*soiOn)
            lents = min(1000,testSeqs[seqOn%len(testSeqs)].length-1000*counter[tsOn])
            soiShift = 1000*soiOn
            tsShift = 1000*counter[tsOn]
            
            SOIrange = []
            TSrange = []
            
            while(offset == nextOffset):
                tSrange, tTrange = cFindRanges(start,stop,offset,lensoi,lents,locWin,soiShift,tsShift)
                SOIrange.append(tSrange)
                TSrange.append(tTrange)
                i+=3
                #print alignmentInfo[soiOn*MAX_MATCH*3+tsOn*3*MAX_MATCH*soiLen+i:soiOn*MAX_MATCH*3+tsOn*3*MAX_MATCH*soiLen+i+3]
                start = alignmentInfo[soiOn*MAX_MATCH*3+tsOn*3*MAX_MATCH*soiLen+i]
                stop = alignmentInfo[soiOn*MAX_MATCH*3+tsOn*3*MAX_MATCH*soiLen+i+1]
                offset = alignmentInfo[soiOn*MAX_MATCH*3+tsOn*3*MAX_MATCH*soiLen+i+2]
                if soiOn*MAX_MATCH*3+tsOn*3*MAX_MATCH*soiLen+i+5 < MAX_MATCH*3*tsLen*soiLen:
                    nextOffset = alignmentInfo[soiOn*MAX_MATCH*3+tsOn*3*MAX_MATCH*soiLen+i+5]
                else:
                    nextOffset = -1
            tSrange, tTrange = cFindRanges(start,stop,offset,lensoi,lents,locWin,soiShift,tsShift)
            SOIrange.append(tSrange)
            TSrange.append(tTrange)
            
            soiStart = max(0,offset-lents+locWin) + soiShift
            soiEnd = min(lensoi-1,offset+locWin-1) + soiShift
            tsStart = max(0,lents-locWin-offset) + tsShift
            tsEnd = lents-1 - max(0,locWin+offset-lensoi) + tsShift
            
            logging.debug('Added align! ' + str(SOIrange) + ' to ' + str(TSrange) + ' of ' + testSeqs[seqOn%len(testSeqs)].name + ' rf ' + str(rframe))
            locAlign = LocalAlignment(SOIrange,TSrange,[soiStart,soiEnd],[tsStart,tsEnd],rframe)
            theAlignment.AddLocalAlign(locAlign) # adds the local alignment to the total alignment
            i += 3
        else:
            # useful information has passed... look at next block.
            if tsOn < tsLen-1 or soiOn < soiLen-1:
                # not done yet
                if soiOn < soiLen-1:
                    #look at next block
                    soiOn += 1
                    i = 0
                else:
                    # look at next section
                    tsOn += 1
                    soiOn = 0
                    i = 0
                    if counter[tsOn] == 0:
                        # reach the end of a test sequence
                        seqOn += 1
                        alignments.AddAlign(theAlignment)
                        if seqOn%len(testSeqs) == 0:
                            rframe += 1
                        theAlignment = Alignment(seqOfInt, testSeqs[seqOn%len(testSeqs)])
            else:
                #done
                alignments.AddAlign(theAlignment)
                i = MAX_MATCH*3
    return alignments

# Find clubs using C
def GenClubsC(seqOfInt, testSeqs, parameters):
    import clubgen_c
    # build up sequence of interest
    soiPartTotal = seqOfInt.length/1000+1
    soiParts = clubgen_c.int_array_1000(soiPartTotal)
    if seqOfInt.readingFrames > 1:
        soiseq = seqOfInt.seq[0]
    else:
        soiseq = seqOfInt.seq
    for i in range(soiPartTotal):
        if seqOfInt.sequenceType == 'AA':
            for j in range(min(1000,seqOfInt.length-1000*i)):
                clubgen_c.a_set(i,j,numpy.int(soiseq[j+1000*i]),soiParts)
            for j in range(min(1000,seqOfInt.length-1000*i),1000):
                clubgen_c.a_set(i,j,0,soiParts)
        else:
            for j in range(min(1000,seqOfInt.length-1000*i)):
                clubgen_c.a_set(i,j,numpy.int(soiseq[j+1000*i]%2**seqOfInt.dataDepths[0]),soiParts)
            for j in range(min(1000,seqOfInt.length-1000*i),1000):
                clubgen_c.a_set(i,j,0,soiParts)
            
    tsPartTotal = 0
    tsPartInfo = []
    tsPartNum = 0
    if seqOfInt.readingFrames == 1:
        for j in range(len(testSeqs)):
            tmpPart = testSeqs[j].length/1000+1
            tsPartTotal += tmpPart
            for i in range(tmpPart):
                tsPartInfo.append(tsPartNum)
            tsPartNum += 1
    else:
        for k in range(seqOfInt.readingFrames):
            for j in range(len(testSeqs)):
                tmpPart = len(testSeqs[j].seq[k])/1000+1
                tsPartTotal += tmpPart
                for i in range(tmpPart):
                    tsPartInfo.append(tsPartNum)
                tsPartNum += 1
                
    #print tsPartInfo
       
    thePart = -1
    counter = 0
    tsPartNum = 0
    tsParts = clubgen_c.int_array_1000(tsPartTotal)
    
    if seqOfInt.readingFrames == 1:
        if seqOfInt.sequenceType == 'AA':
            for i in range(tsPartTotal):
                if thePart == tsPartInfo[i]:
                    counter += 1
                else:
                    thePart = tsPartInfo[i]
                    counter = 0
                for j in range(min(1000,testSeqs[tsPartInfo[i]].length-1000*counter)):
                    clubgen_c.a_set(i,j,numpy.int(testSeqs[tsPartInfo[i]].seq[j+1000*counter]),tsParts)
                for j in range(min(1000,testSeqs[tsPartInfo[i]].length-1000*counter),1000):
                    clubgen_c.a_set(i,j,0,tsParts)
        else:
            for i in range(tsPartTotal):
                if thePart == tsPartInfo[i]:
                    counter += 1
                else:
                    thePart = tsPartInfo[i]
                    counter = 0
                for j in range(min(1000,testSeqs[tsPartInfo[i]].length-1000*counter)):
                    clubgen_c.a_set(i,j,numpy.int(testSeqs[tsPartInfo[i]].seq[j+1000*counter]%2**seqOfInt.dataDepths[0]),tsParts)
                for j in range(min(1000,testSeqs[tsPartInfo[i]].length-1000*counter),1000):
                    clubgen_c.a_set(i,j,0,tsParts)
    else:
        #ts = numpy.zeros(tsPartTotal*1000,dtype = numpy.uint8)
        for k in range(seqOfInt.readingFrames):
            for i in range(tsPartTotal/seqOfInt.readingFrames):
                infoLen = len(testSeqs[tsPartInfo[i]%len(testSeqs)].seq[k])
                #print thePart,tsPartInfo[i+k*tsPartTotal/seqOfInt.readingFrames]
                if thePart == tsPartInfo[i+k*tsPartTotal/seqOfInt.readingFrames]:
                    counter += 1
                else:
                    thePart = tsPartInfo[i+k*tsPartTotal/seqOfInt.readingFrames]
                    counter = 0
                    
                #ts[1000*i+k*tsPartTotal/seqOfInt.readingFrames*1000:1000*i+k*tsPartTotal/seqOfInt.readingFrames*1000+min(1000,infoLen-1000*counter)] = testSeqs[tsPartInfo[i]%len(testSeqs)].seq[k][1000*counter:1000*counter+min(1000,infoLen-1000*counter)]%2**seqOfInt.dataDepths[0]
                
                #ts[1000*i+k*tsPartTotal/seqOfInt.readingFrames*1000+min(1000,infoLen-1000*counter):1000*i+k*tsPartTotal/seqOfInt.readingFrames*1000+1000] = numpy.zeros(1000-min(1000,infoLen-1000*counter),dtype = numpy.uint8)
                #print counter
                #print infoLen-1000*counter
                for j in range(min(1000,infoLen-1000*counter)):
                    clubgen_c.a_set(i+k*tsPartTotal/seqOfInt.readingFrames,j,numpy.int(testSeqs[tsPartInfo[i]%len(testSeqs)].seq[k][j+1000*counter]%2**seqOfInt.dataDepths[0]),tsParts)
                    #print 'set', i+k*tsPartTotal/seqOfInt.readingFrames,j,numpy.int(testSeqs[tsPartInfo[i]%len(testSeqs)].seq[k][j+1000*counter]%2**seqOfInt.dataDepths[0])
                    #ts[j+1000*i+k*tsPartTotal/seqOfInt.readingFrames*1000] = testSeqs[tsPartInfo[i]%len(testSeqs)].seq[k][j+1000*counter]%2**seqOfInt.dataDepths[0]
                for j in range(min(1000,infoLen-1000*counter),1000):
                    clubgen_c.a_set(i+k*tsPartTotal/seqOfInt.readingFrames,j,0,tsParts)
                    #ts[j+1000*i+k*tsPartTotal/seqOfInt.readingFrames*1000] = 0
    
    cTsPartInfo = clubgen_c.intArray(tsPartTotal)
    for i in range(tsPartTotal):
        cTsPartInfo[i] = tsPartInfo[i]
        
    # make vector that contains all of the info
    MAX_MATCH = 20
    MAX_INLINE = 5
    alignmentInfo = clubgen_c.intArray(MAX_MATCH*5*len(testSeqs)*seqOfInt.readingFrames)
    for i in range(MAX_MATCH*5*len(testSeqs)*seqOfInt.readingFrames):
        alignmentInfo[i] = 0;
        
    if seqOfInt.readingFrames == 1:
        tsLens = clubgen_c.intArray(len(testSeqs))
        for i in range(len(testSeqs)):
            tsLens[i] = testSeqs[i].length
    else:
        tsLens = clubgen_c.intArray(len(testSeqs)*seqOfInt.readingFrames)
        for k in range(seqOfInt.readingFrames):
            for i in range(len(testSeqs)):
                tsLens[i+k*len(testSeqs)] = len(testSeqs[i].seq[k])
    
    #tsLengths = clubgen_c.intArray(len(testSeqs))
    #for i in range(len(testSeqs)):
        #tsLengths[i] = testSeqs[i].length
    
    
    
    #print 'about to go into c'
    clubgen_c.TesterAlign_c(soiPartTotal, soiParts, tsPartTotal, tsParts, cTsPartInfo, alignmentInfo, seqOfInt.length, tsLens, parameters['GLOBAL_SIG_LEVEL'], parameters['LOCAL_SIG_LEVEL'], parameters['LOCAL_WINDOW'], MAX_MATCH, MAX_INLINE, parameters['LOCAL_BRIDGE_WIDTH'], parameters['COMB_SPEED'])

    alignments = SetOfAlignments(seqOfInt)
            
    i = 0
    seqOn = 0
    rframe = 0
    
    locWin = parameters['LOCAL_WINDOW']
    theAlignment = Alignment(seqOfInt, testSeqs[seqOn])
    while(i < MAX_MATCH*5):
        
        if alignmentInfo[seqOn*MAX_MATCH*5+i] != alignmentInfo[seqOn*MAX_MATCH*5+i+1]:
            #logging.debug('alignmentInfo: ' + str(alignmentInfo[seqOn*MAX_MATCH*5+i]) + ',' + str(alignmentInfo[seqOn*MAX_MATCH*5+i+1]) + ',' + str(alignmentInfo[seqOn*MAX_MATCH*5+i+2]) + ',' + str(alignmentInfo[seqOn*MAX_MATCH*5+i+3]) + ',' + str(alignmentInfo[seqOn*MAX_MATCH*5+i+4]))
            lensoi = min(1000,seqOfInt.length-1000*alignmentInfo[seqOn*MAX_MATCH*5+i+3])
            lents = min(1000,testSeqs[seqOn%len(testSeqs)].length-1000*alignmentInfo[seqOn*MAX_MATCH*5+i+4])
            soiShift = 1000*alignmentInfo[seqOn*MAX_MATCH*5+i+3]
            tsShift = 1000*alignmentInfo[seqOn*MAX_MATCH*5+i+4]
                
            SOIrange = []
            TSrange = []
                
            while(alignmentInfo[seqOn*MAX_MATCH*5+i+2] == alignmentInfo[seqOn*MAX_MATCH*5+i+7] and alignmentInfo[seqOn*MAX_MATCH*5+i] != alignmentInfo[seqOn*MAX_MATCH*5+i+1]):
                t3,t4 = cFindRanges(alignmentInfo[seqOn*MAX_MATCH*5+i], alignmentInfo[seqOn*MAX_MATCH*5+i+1], alignmentInfo[seqOn*MAX_MATCH*5+i+2], lensoi, lents, locWin, soiShift, tsShift)
                SOIrange.append(t3)
                TSrange.append(t4)
                i += 5
                logging.debug('alignmentInfo: ' + str(alignmentInfo[seqOn*MAX_MATCH*5+i]) + ',' + str(alignmentInfo[seqOn*MAX_MATCH*5+i+1]) + ',' + str(alignmentInfo[seqOn*MAX_MATCH*5+i+2]) + ',' + str(alignmentInfo[seqOn*MAX_MATCH*5+i+3]) + ',' + str(alignmentInfo[seqOn*MAX_MATCH*5+i+4]))
            t3,t4 = cFindRanges(alignmentInfo[seqOn*MAX_MATCH*5+i], alignmentInfo[seqOn*MAX_MATCH*5+i+1], alignmentInfo[seqOn*MAX_MATCH*5+i+2], lensoi, lents, locWin, soiShift, tsShift)
            SOIrange.append(t3)
            TSrange.append(t4)
            
            soiStart = max(0,alignmentInfo[seqOn*MAX_MATCH*5+i+2]-lents+locWin) + soiShift
            soiEnd = min(lensoi-1,alignmentInfo[seqOn*MAX_MATCH*5+i+2]+locWin-1) + soiShift
            tsStart = max(0,lents-locWin-alignmentInfo[seqOn*MAX_MATCH*5+i+2]) + tsShift
            tsEnd = lents-1 - max(0,locWin+alignmentInfo[seqOn*MAX_MATCH*5+i+2]-lensoi) + tsShift
            
            #locAlign = LocalAlignment(SOIrange,TSrange,colScore,cdfL,cdfR,fullScoreVec)
            logging.debug('Added align! ' + str(SOIrange) + ' to ' + str(TSrange) + ' of ' + testSeqs[seqOn%len(testSeqs)].name + ' rf ' + str(rframe))
            locAlign = LocalAlignment(SOIrange,TSrange,[soiStart,soiEnd],[tsStart,tsEnd],rframe)
            theAlignment.AddLocalAlign(locAlign) # adds the local alignment to the total alignment
            i += 5
        else:
            #print 'added the alignment'
            alignments.AddAlign(theAlignment)
            if seqOn < len(testSeqs)*seqOfInt.readingFrames - 1:
                i = 0
                seqOn += 1
                if seqOn%len(testSeqs) == 0:
                    rframe += 1
                theAlignment = Alignment(seqOfInt, testSeqs[seqOn%len(testSeqs)])
            else:
                i = MAX_MATCH*5

    return alignments

def GenClubsPython(seqOfInt, testSeqs, parameters):
    """
        GenClubs(seqOfInt, testSeqs):
        Generates all alignment objects from the initial seqOfInt[0] and all sequences in the testSeqs objects
        
        uses GenAlignment()
    """
    alignments = SetOfAlignments(seqOfInt)
    if seqOfInt.readingFrames == 1:
        i = 1
        for testSeq in testSeqs:
            logging.debug('(' + str(i) + '/' + str(len(testSeqs)) + ') Comparing ' + seqOfInt.name + ' to ' + testSeq.name)
            alignments.AddAlign(GenAlignmentPy(seqOfInt,testSeq,seqOfInt.seq,testSeq.seq,parameters))
            i += 1
    elif seqOfInt.readingFrames == 6:
        i = 1
        for testSeq in testSeqs:
            j = 1
            for ts in testSeq.seq:
                logging.debug('(' + str((i-1)/6+1) + '+' + str(j) + '/' + str(len(testSeqs)) + '+' + str(seqOfInt.readingFrames) + ') Comparing ' + seqOfInt.name + ' to ' + testSeq.name)
                if len(ts) > 5:
                    alignments.AddAlign(GenAlignmentPy(seqOfInt,testSeq,seqOfInt.seq[0],ts,parameters))
                i += 1
                j += 1
    elif seqOfInt.readingFrames == 3:
        i = 1
        for testSeq in testSeqs:
            j = 1
            for ts in testSeq.seq:
                logging.debug('(' + str((i-1)/3+1) + '+' + str(j) + '/' + str(len(testSeqs)) + '+' + str(seqOfInt.readingFrames) + ') Comparing ' + seqOfInt.name + ' to ' + testSeq.name)
                if len(ts) > 2:
                    alignments.AddAlign(GenAlignmentPy(seqOfInt,testSeq,seqOfInt.seq[0],ts,parameters))
                i += 1
                j += 1
    #plt.show()
    return alignments
    
def GenAlignmentPy(soi_seq, ts_seq, soi, ts,parameters):
    pass
    
def GenAlignment(soi_seq, ts_seq, soi, ts,parameters):
    lensoi = int(soi_seq.length)
    lents = int(len(ts))
    
    breakpoints = soi_seq.dataDepths
    
    c_soi = clubgen_c.intArray(lensoi)
    if breakpoints[0] == 0:
        for i in range(lensoi):
            c_soi[i] = soi[i]
    else:
        for i in range(lensoi):
            c_soi[i] = soi[i]%2**breakpoints[0]
        
    c_ts = clubgen_c.intArray(lents)
    if breakpoints[0] == 0:
        for i in range(lents):
            c_ts[i] = ts[i]
    else:
        for i in range(lents):
            c_ts[i] = ts[i]%2**breakpoints[0]
    
    MAX_MATCH = 100
    MAX_INLINE = 3
    c_localStart = clubgen_c.intArray(MAX_MATCH)
    c_localEnd = clubgen_c.intArray(MAX_MATCH)
    c_localShift = clubgen_c.intArray(MAX_MATCH)
    for i in range(MAX_MATCH):
        c_localStart[i] = 0
        c_localEnd[i] = 0
        c_localShift[i] = 0
    
    #print soiVecComp
    #print tsVecComp
    #print numMatches
    #print 'lol?'
    #print 'local window', parameters['LOCAL_WINDOW']

    t0 = time.time()
    a = clubgen_c.GenAlignment_c(c_soi,c_ts,lensoi,lents, parameters['GLOBAL_SIG_LEVEL'], parameters['LOCAL_SIG_LEVEL'], parameters['LOCAL_WINDOW'], MAX_MATCH, MAX_INLINE, parameters['LOCAL_BRIDGE_WIDTH'], parameters['COMB_SPEED'], c_localStart,c_localEnd,c_localShift)
    TIMER = time.time() - t0

    theAlignment = Alignment(soi_seq, ts_seq)
    
    ##############################
    ##### C does local ###########
    ##############################
    
    if c_localStart[0] == c_localEnd[0]:
        return theAlignment
    
    cdfL = numpy.zeros(lensoi)
    cdfR = numpy.zeros(lensoi)
    
    if breakpoints[0] == 0:
        colScore = numpy.zeros(lensoi)
    else:
        colScore = numpy.zeros((len(breakpoints),lensoi))
        
    locWin = parameters['LOCAL_WINDOW']
    
    i = 0
    while(i < MAX_MATCH):
        if c_localStart[i] != c_localEnd[i]:

            fullScoreVec = numpy.zeros(lensoi)
            if soi_seq.dataDepths[0] != 0:
                colScore = numpy.zeros((len(soi_seq.dataDepths),lensoi))
            else:
                colScore = numpy.zeros(lensoi)
                
            SOIrange = []
            TSrange = []
                
            while(c_localShift[i] == c_localShift[i+1] and c_localStart[i] != c_localEnd[i]):
                #t1,t2,t3,t4,t5,t6 = cToPyDataFixer(c_localStart[i], c_localEnd[i], c_localShift[i], lensoi, lents, locWin, theAlignment)
                #colScore += t1
                #fullScoreVec += t2
                t3,t4 = cFindRanges(c_localStart[i], c_localEnd[i], c_localShift[i], lensoi, lents, locWin)
                SOIrange.append(t3)
                TSrange.append(t4)
                i += 1
            #t1,t2,t3,t4,cdfL,cdfR = cToPyDataFixer(c_localStart[i], c_localEnd[i], c_localShift[i], lensoi, lents, locWin, theAlignment)
            #colScore += t1
            #fullScoreVec += t2
            t3,t4 = cFindRanges(c_localStart[i], c_localEnd[i], c_localShift[i], lensoi, lents, locWin)
            SOIrange.append(t3)
            TSrange.append(t4)
            
            soiStart = max(0,c_localShift[i]-lents+locWin)
            soiEnd = min(lensoi-1,c_localShift[i]+locWin-1)
            tsStart = max(0,lents-locWin-c_localShift[i])
            tsEnd = lents-1 - max(0,locWin+c_localShift[i]-lensoi)
            
            #locAlign = LocalAlignment(SOIrange,TSrange,colScore,cdfL,cdfR,fullScoreVec)
            locAlign = LocalAlignment(SOIrange,TSrange,[soiStart,soiEnd],[tsStart,tsEnd])
            theAlignment.AddLocalAlign(locAlign) # adds the local alignment to the total alignment
            logging.debug('Added align! ' + str(SOIrange))
            i += 1
        else:
            break
            i = MAX_MATCH
    return theAlignment
    
def cFindRanges(localStart, localEnd, localShift, lensoi, lents, locWin, soiShift, tsShift):
    soiStart = max(0,localShift-lents+locWin)
    soiEnd = min(lensoi-1,localShift+locWin-1)
    tsStart = max(0,lents-locWin-localShift)
    tsEnd = lents-1 - max(0,locWin+localShift-lensoi)
    soiLstop = soiStart + localStart + soiShift
    soiRstop = min(soiStart + localEnd,soiEnd) + soiShift
    tsLstop = tsStart + localStart + tsShift
    tsRstop = tsStart + localEnd + tsShift
    return [soiLstop,soiRstop],[tsLstop,tsRstop]
    
def cToPyDataFixer(localStart, localEnd, localShift, lensoi, lents, locWin, align):
    soiStart = max(0,localShift-lents+locWin)
    soiEnd = min(lensoi-1,localShift+locWin-1)
    tsStart = max(0,lents-locWin-localShift)
    tsEnd = lents-1 - max(0,locWin+localShift-lensoi)
    soiLstop = soiStart + localStart
    soiRstop = min(soiStart + localEnd,soiEnd)
    tsLstop = tsStart + localStart
    tsRstop = tsStart + localEnd
    
    if align.testSeq.readingFrames == 1:
        soi = numpy.array(align.seqOfInt.seq,dtype = numpy.int32)
        ts = numpy.array(align.testSeq.seq,dtype = numpy.int32)
    else:
        soi = numpy.array(align.seqOfInt.seq[0],dtype = numpy.int32)
        ########### CHEAT FOR MULTIPLE READING FRAMES ################
        ts = numpy.array(align.testSeq.seq[offset%align.testSeq.readingFrames],dtype = numpy.int32)
            
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
        locAscore = matchVecIn
    else:
        colScore[0] = matchVecIn
        for j in range(len(align.seqOfInt.dataDepths)-1):
            colScore[1+j][soiLstop:soiRstop] = matchVecIn[soiLstop:soiRstop]*numpy.asarray(soi[soiLstop:soiRstop]%2**align.seqOfInt.dataDepths[j+1] == ts[tsLstop:tsRstop]%2**align.seqOfInt.dataDepths[j+1],dtype = numpy.int32)
        locAscore = colScore # [:,soiLstop:soiRstop]
            
    # set cdfL/cdfR
            
    cdfL = numpy.zeros(lensoi)
    cdfR = numpy.zeros(lensoi)
            
    numMatches = GetNumMatches(soiVecComp,tsVecComp)
            
    singleexp = numMatches/(float(soiEnd-soiStart)*(soiEnd-soiStart))
            
    correctionL = numpy.cumsum(numpy.ones(soiEnd-soiStart)*singleexp)
    cdfL[soiStart:soiEnd] = numpy.cumsum(fullScoreVec[soiStart:soiEnd]) - correctionL
            
    cdfR[soiStart:soiEnd] = numpy.asarray(numpy.cumsum(fullScoreVec[soiStart:soiEnd][::-1]) - correctionL)[::-1]
            
            
    #locA.fullScoreVec = fullScoreVec
    #locA.SOIrange = [soiLstop,soiRstop]
    #locA.TSrange = [tsLstop,tsRstop]
    
    return locAscore,fullScoreVec,[soiLstop,soiRstop],[tsLstop,tsRstop],cdfL,cdfR
    
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
    
def GenAlignmentPy(soi_seq, ts_seq, soi, ts,parameters):
    """
        GenAlignment(soi_seq, ts_seq, soi, ts):
        Generates an alignment object between the sequence objects soi and ts
        These are contained in the objects soi_seq and ts_seq respectively
        
        First performs a global filter determining if there are a significant number of matches beyond the expectation
        in a particular offset, quantified in multiples of standard deviations above the expectation
        
        The offsets that pass this filter are then filtered locally looking for areas of high local density beyond
        the expected density, again quantified in multiples of standard deviations above the expectation
        
        If an offset passes both filters it is considered a local alignment and is added to the alignment object
        All offsets are considered and then the alignment is passed back
    """
    # initialize variables
    theAlignment = Alignment(soi_seq, ts_seq)
    
    sigDif = parameters['GLOBAL_SIG_LEVEL']
    
    lensoi = soi_seq.length
    #lents = ts_seq.length
    breakpoints = soi_seq.dataDepths
    # extends test sequence
    tmpts = list(ts)
    ts = list(ts)
    lents = len(tmpts)
    for i in range(lensoi/lents+1):
        ts.extend(tmpts)
    ts = numpy.array(ts,dtype=numpy.int32) # switches it to an array
    soi = numpy.array(soi,dtype=numpy.int32) # switches it to an array

    expMatch = 0
    expStd = 0
    
    totaltime = 0
    
    # if we are looking at more than one data depth generate the column score array
    if breakpoints[0] != 0:
        colScore = numpy.zeros((len(breakpoints),lensoi))
    else:
        colScore = numpy.zeros(lensoi)
    colClub = numpy.zeros(lensoi)
    
    # gets the composition of the soi and ts
    if breakpoints[0] == 0:
        soiVecComp = GetVecComp(soi)
        tsVecComp = GetVecComp(ts[0:lensoi])
    else:
        soiVecComp = GetVecComp(soi%2**breakpoints[0])
        tsVecComp = GetVecComp(ts[0:lensoi]%2**breakpoints[0])
    
    # gets the expected number of matches given the composition of the soi
    numMatches = GetNumMatches(soiVecComp,tsVecComp)
    
    # searches through all offsets
    for i in range(0,lents):

        #finds the actual number of matches
        if breakpoints[0] == 0:
            matchVec = numpy.asarray(soi == ts[i:lensoi+i],dtype = numpy.int32)
        else:
            matchVec = numpy.asarray(soi%2**breakpoints[0] == ts[i:lensoi+i]%2**breakpoints[0],dtype = numpy.int32)

        # Calculates the actual number of matches
        if numpy.shape(numpy.where(matchVec))[0] == 0:
            actMatches = 0
        else:
            actMatches = numpy.shape(numpy.where(matchVec)[0])[0]
        # calculates what we would expect the number of matches and the standard deviation to be
        #expMatch += numMatches
        #expStd += math.sqrt(numMatches/float(lensoi))
        if numMatches < 0:
            numMatches = 0

        # Global filter, must be at least sigDif standard deviations away to go to the next round of filtering
        if actMatches > numMatches/float(lensoi) + sigDif*math.sqrt(numMatches/float(lensoi)):
            
            #local filter
            #if PLOT_CDF_ON == True: # captures a bit more data for visualization, default on
                #origMatchVec = matchVec[:]
                #matchVecIn, matchVec, lstop, rstop, cdfL, cdfR = ClubCullCDFVerbose(matchVec, numMatches/float(lensoi*lensoi))
            #else:
                #matchVecIn, matchVec, lstop, rstop = ClubCullCDF(matchVec, numMatches/float(lensoi*lensoi))
            origMatchVec = matchVec[:]
            
            
            
            ###################################
            ##### experimental section!!! #####
            ###################################
            
            t0 = time.time()
            
            #matchVecIn, matchVec, lstop, rstop, cdfL, cdfR = ClubCullCDF(matchVec, numMatches/float(lensoi*lensoi))
            #matchVecIn, matchVec, lstop, rstop, cdfL, cdfR = ClubCullCDFlocal(matchVec, soi, ts, breakpoints)
            #matchVecIn, matchVec, lstop, rstop, cdfL, cdfR = ClubCullRatio(matchVec)
            
            matchVecIn, matchVec, lstop, rstop, cdfL, cdfR = ClubCullCDFmulti(matchVec, numMatches/float(lensoi*lensoi),parameters)
            
                    
            
            totaltime += time.time() - t0
            
            ###################################
            ###################################
            ###################################
            
            # finds the remaining matches after both rounds of filtering
            actMatches2 = len(numpy.where(matchVecIn)[0])
            
            # increments club
            if breakpoints[0] == 0:
                colScore = matchVecIn
            else:
                colScore[0] = matchVecIn
                for j in range(len(breakpoints)-1):
                    colScore[1+j] = matchVecIn*numpy.asarray(soi%2**breakpoints[j+1] == ts[i:lensoi+i]%2**breakpoints[j+1],dtype = numpy.int32)
            
            for k in range(len(lstop)):
                if k < len(rstop):
                    if lstop[k] != rstop[k]:
                    # graphical output option (default on), adds more data about the stop and start positions
                        if parameters['PLOT_CDF'] == True:
                            if breakpoints[0] != 0:
                                locAlign = LocalAlignment([[lstop[k],rstop[k]]],[[(lstop[k]+i)%lents,(lstop[k]+i)%lents+rstop[k]-lstop[k]]],[0,lensoi-1],[0,lents-1])
                            else:
                                locAlign = LocalAlignment([[lstop[k],rstop[k]]],[[(lstop[k]+i)%lents,(lstop[k]+i)%lents+rstop[k]-lstop[k]]],[0,lensoi-1],[0,lents-1])
                        else:
                            if breakpoints[0] != 0:
                                locAlign = LocalAlignment([[lstop[k],rstop[k]]],[[(lstop[k]+i)%lents,(lstop[k]+i)%lents+rstop[k]-lstop[k]]],[0,lensoi-1],[0,lents-1])
                            else:
                                locAlign = LocalAlignment([[lstop[k],rstop[k]]],[[(lstop[k]+i)%lents,(lstop[k]+i)%lents+rstop[k]-lstop[k]]],[0,lensoi-1],[0,lents-1])
                        theAlignment.AddLocalAlign(locAlign) # adds the local alignment to the total alignment
                        logging.debug('Added align! ' + str(lstop[k]) + ' - ' + str(rstop[k]) + ' -- ' + str((actMatches,numMatches/float(lensoi) + sigDif*math.sqrt(numMatches/float(lensoi)))) + ' - off: ' + str(i))
    
        # offsets the number of matches by shifting the window of interest
        if breakpoints[0] == 0:
            numMatches += soiVecComp[ts[i+lensoi-1]]
            numMatches -= soiVecComp[ts[i-1]]
        else:
            numMatches += soiVecComp[ts[i+lensoi-1]%2**breakpoints[0]]
            numMatches -= soiVecComp[ts[i-1]%2**breakpoints[0]]
    return theAlignment
    
def ClubCullCDFmulti(match,singleexp,parameters):
    """
        ClubCullCDF(match,singleexp):
    """
    # initialize variables
    width = parameters['LOCAL_WINDOW']
    sigAbove = parameters['LOCAL_SIG_LEVEL']

    needed = math.sqrt(singleexp*width)*sigAbove
    sig = singleexp*5

    correctionL = numpy.cumsum(numpy.ones(len(match))*singleexp)
    cdfL = numpy.cumsum(match) - correctionL
    #print cdfL
    #Lstop = numpy.max(numpy.where(cdfL < sig))
    diffLt = cdfL[width:] - cdfL[:-width]
    diffL = numpy.where(numpy.asarray(diffLt  > needed,dtype=numpy.int32)[:-1]*numpy.asarray(diffLt[1:]<=diffLt[:-1],dtype=numpy.int32))[0]
    #plt.figure()
    #plt.plot(diffLt)
    #plt.show()
    #print diffL
    #print len(diffL)
    if len(diffL) == 0:
        Lstop = [0]
    else:
        Lstop = []
        stop = False
        while stop == False:
            Lstop2 = numpy.min(diffL)
            stop2 = False
            while(stop2 == False):
                if match[Lstop2] == 0:
                    Lstop2 += 1
                else:
                    stop2 = True
            Lstop.append(Lstop2)
            startAgain = numpy.where(diffL[1:] - diffL[:-1] > width*parameters['LOCAL_BRIDGE_WIDTH'])[0]
            if len(startAgain) == 0:
                stop = True
            else:
                diffL = diffL[numpy.min(startAgain)+1:]
    #print Lstop
    #correctionR = correctionL[::-1]
    cdfR = numpy.asarray(numpy.cumsum(match[::-1]) - correctionL)[::-1]
    diffRt = cdfR[:-width] - cdfR[width:]
    diffR = numpy.where(numpy.asarray(diffRt> needed,dtype=numpy.int32)[1:]*numpy.asarray(diffRt[1:]<=diffRt[:-1],dtype=numpy.int32))[0]
    #print diffR
    if len(diffR) == 0:
        if Lstop[0] != 0:
            Rstop = [len(match)]
        else:
            Rstop = [0]
    else:
        Rstop = []
        stop = False
        while stop == False:
            Rstop2 = numpy.max(diffR) + width
            stop2 = False
            while(stop2 == False):
                if match[Rstop2] == 0:
                    Rstop2 -= 1
                else:
                    stop2 = True
            Rstop.append(Rstop2)
            startAgain = numpy.where(diffR[1:] - diffR[:-1] > width*parameters['LOCAL_BRIDGE_WIDTH'])[0]
            if len(startAgain) == 0:
                stop = True
            else:
                diffR = diffR[:numpy.max(startAgain)-1]
            if len(diffR) == 0:
                stop = True
    Rstop.reverse()
    
    ############ NEEDS WORK ##########
    ## hack culling of short aligns
    removal = []
    for i in range(min(len(Lstop),len(Rstop))):
        if Rstop[i] - Lstop[i] < width:
            removal.append(i)
    Rstop2 = []
    Lstop2 = []
    for i in range(min(len(Lstop),len(Rstop))):
        if i not in removal:
            Rstop2.append(Rstop[i])
            Lstop2.append(Lstop[i])
    Rstop = Rstop2[:]
    Lstop = Lstop2[:]
    
    
    #print cdfR
    #Rstop = numpy.min(numpy.where(cdfR < sig),0)
    #print Rstop
    
    # secondary refinement

    cdfClub = numpy.zeros(len(match))
    for i in range(len(Lstop)):
        if i < len(Rstop):
            cdfClub[Lstop[i]:Rstop[i]] = 1
    return match*cdfClub, cdfClub, Lstop, Rstop, cdfL, cdfR