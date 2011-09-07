import numpy, logging
from PARAMS import *
from VelvetropeClasses import *
import matplotlib.pylab as plt
import clubs_generator

def IntronFinder(alignments,parameters):
    logging.info('Finding introns...')
    ratioSOI = parameters['INTRON_OVERLAP']
    ratioTS = parameters['INTRON_OVERLAP']
    ratioTrue = parameters['INTRON_TRUTH_RATIO']
    for indAligns in alignments.aligns:
        intronArea = numpy.zeros(indAligns.seqOfInt.length)
        removal = []
        i = -1
        for locAligns in indAligns.localAligns:
            i += 1
            for locAlignsComp in indAligns.localAligns:
                if locAligns != locAlignsComp:
                    if ratioOverlap(locAligns.SOIrange,locAligns.SOIrange) > ratioSOI:
                        if ratioOverlap(locAligns.TSrange,locAlignsComp.TSrange) > ratioTS:
                            for k in range(locAligns.SOIrange[0],locAligns.SOIrange[1]):
                                intronArea[k] += 1
                            break
                            
        #plt.figure()
        #plt.plot(intronArea)
        #plt.plot([0,indAligns.seqOfInt.length],[ratioTrue*len(indAligns.localAligns),ratioTrue*len(indAligns.localAligns)])
        #plt.show()
        
        # find intron areas
        intronInfos = numpy.asarray(numpy.asarray(intronArea,dtype=numpy.int32) > ratioTrue*len(indAligns.localAligns),dtype=numpy.int32)
        
            
            
        intronRanges = rangeFinder(intronInfos)
        
        if len(intronRanges) < 1:
            logging.info('Found no introns... ' + str(intronRanges))
            break
            
        logging.info('Found introns... ' + str(intronRanges))
        
        logging.debug('Old seq length: ' + str(len(indAligns.seqOfInt.seq)))
        # cut out intron areas
        for i in range(len(intronRanges)):
            if i == 0:
                newSOIseq = indAligns.seqOfInt.seq[0:intronRanges[0][0]]
                logging.debug('Saved area from 0 to ' + str(intronRanges[0][0]))
            else:
                newSOIseq.extend(indAligns.seqOfInt.seq[intronRanges[i-1][1]:intronRanges[i][0]])
                logging.debug('Saved area from ' + str(intronRanges[i-1][1]) + ' to ' + str(intronRanges[i][0]))
            if len(intronRanges) == i+1:
                newSOIseq.extend(indAligns.seqOfInt.seq[intronRanges[i][1]:])
                logging.debug('Saved area from ' + str(intronRanges[i][1]) + ' to ' + str(len(indAligns.seqOfInt.seq)))
        logging.debug('New seq length: ' + str(len(newSOIseq)))
        
        newSOI = Sequence(indAligns.seqOfInt.sequenceType, indAligns.seqOfInt.dataDepths, indAligns.seqOfInt.readingFrames, indAligns.seqOfInt.name, newSOIseq)
        
        # realign
        subalign = clubs_generator.GenClubs(newSOI, [indAligns.testSeq],parameters)
        
        # correct alignments
        indAligns.localAligns = []
        for newAlign in subalign.aligns:
            toAdd = Alignment(indAligns.seqOfInt,indAligns.testSeq)
            for newLoc in newAlign.localAligns:
                newLoc = recoverAlignSOI(newLoc,intronRanges)
                if parameters['PLOT_CDF'] == True:
                    newLoc = recoverAlignPlot(newLoc,intronRanges)
                    #tmpvec = list(newLoc.fullScoreVec[0:intronInfos[0]])
                    #tmpvec.extend(-1*numpy.ones(intronInfos[-1] - intronInfos[0]))
                    #tmpvec.extend(newLoc.fullScoreVec[intronInfos[0]:])
                    #newLoc.fullScoreVec = numpy.asarray(tmpvec)
                    #tmpCDFL = list(newLoc.cdfL[0:intronInfos[0]])
                    #tmpCDFL.extend(numpy.zeros(intronInfos[-1] - intronInfos[0]))
                    #tmpCDFL.extend(newLoc.cdfL[intronInfos[0]:])
                    #newLoc.cdfL = numpy.asarray(tmpCDFL)
                    #tmpCDFR = list(newLoc.cdfR[0:intronInfos[0]])
                    #tmpCDFR.extend(numpy.zeros(intronInfos[-1] - intronInfos[0]))
                    #tmpCDFR.extend(newLoc.cdfR[intronInfos[0]:])
                    #newLoc.cdfR = numpy.asarray(tmpCDFR)
                indAligns.AddLocalAlign(newLoc)
                logging.debug('Added corrected alignment inbetween intron')
                #else:
                    #if PLOT_CDF_ON == False:
                        #toAddLoc = LocalAlignment([newLoc.SOIrange[0] + intronInfos[-1] - intronInfos[0], newLoc.SOIrange[1] + intronInfos[-1] - intronInfos[0]],newLoc.TSrange,newLoc.score)
                        #indAligns.AddLocalAlign(toAddLoc)
                    #else:
                        #toAddLoc = LocalAlignment([newLoc.SOIrange[0] + intronInfos[-1] - intronInfos[0], newLoc.SOIrange[1] + intronInfos[-1] - intronInfos[0]],newLoc.TSrange,newLoc.score,newLoc.cdfL,newLoc.cdfR,newLoc.fullScoreVec)
                    #indAligns.AddLocalAlign(toAddLoc)
                    #logging.info('Added corrected alignment after intron')
                    
    return alignments
    
def recoverAlignPlot(newLoc,intronRanges):
    tmpvec = list(newLoc.fullScoreVec[0:intronRanges[0][0]])
    #logging.debug('Added region of length ' + str(len(newLoc.fullScoreVec[0:intronRanges[0][0]])))
    #logging.debug('New length = ' + str(len(tmpvec)))
    tmpvec.extend(-1*numpy.ones(intronRanges[0][1] - intronRanges[0][0]))
    #logging.debug('Added intron of length ' + str(intronRanges[0][1] - intronRanges[0][0]))
    #logging.debug('New length = ' + str(len(tmpvec)))
    for i in range(1,len(intronRanges)):
        tmpvec.extend(newLoc.fullScoreVec[intronRanges[i-1][0]:intronRanges[i][0] - PrevInts(intronRanges,i)])
        #logging.debug('Added region of length ' + str(len(newLoc.fullScoreVec[intronRanges[i-1][0]:intronRanges[i][0] - PrevInts(intronRanges,i)])))
        #logging.debug('Added region from ' + str(intronRanges[i-1][0]) + ' to ' + str(intronRanges[i][0] - PrevInts(intronRanges,i)))
        #logging.debug('New length = ' + str(len(tmpvec)))
        tmpvec.extend(-1*numpy.ones(intronRanges[i][1] - intronRanges[i][0]))
        #logging.debug('Added intron of length ' + str(intronRanges[i][1] - intronRanges[i][0]))
        #logging.debug('New length = ' + str(len(tmpvec)))
    tmpvec.extend(newLoc.fullScoreVec[intronRanges[len(intronRanges)-1][0]-PrevInts(intronRanges,len(intronRanges)-1):])
    #logging.debug('Added region of length ' + str(len(newLoc.fullScoreVec[intronRanges[len(intronRanges)-1][0]-PrevInts(intronRanges,len(intronRanges)-1):])))
    #logging.debug('Added region from ' + str(intronRanges[len(intronRanges)-1][0]-PrevInts(intronRanges,len(intronRanges)-1)))
    #logging.debug('Assembled length = ' + str(len(tmpvec)))
    
    tmpCDFL = list(newLoc.cdfL[0:intronRanges[0][0]])
    tmpCDFL.extend(-1*numpy.ones(intronRanges[0][1] - intronRanges[0][0]))
    for i in range(1,len(intronRanges)):
        tmpCDFL.extend(newLoc.cdfL[intronRanges[i-1][0]:intronRanges[i][0] - PrevInts(intronRanges,i)])
        tmpCDFL.extend(-1*numpy.ones(intronRanges[i][1] - intronRanges[i][0]))
    tmpCDFL.extend(newLoc.cdfL[intronRanges[len(intronRanges)-1][0]-PrevInts(intronRanges,len(intronRanges)-1):])
    
    
    tmpCDFR = list(newLoc.cdfR[0:intronRanges[0][0]])
    tmpCDFR.extend(-1*numpy.ones(intronRanges[0][1] - intronRanges[0][0]))
    for i in range(1,len(intronRanges)):
        tmpCDFR.extend(newLoc.cdfR[intronRanges[i-1][0]:intronRanges[i][0] - PrevInts(intronRanges,i)])
        tmpCDFR.extend(-1*numpy.ones(intronRanges[i][1] - intronRanges[i][0]))
    tmpCDFR.extend(newLoc.cdfR[intronRanges[len(intronRanges)-1][0]-PrevInts(intronRanges,len(intronRanges)-1):])
    
    newLoc.fullScoreVec = numpy.asarray(tmpvec)
    newLoc.cdfL = numpy.asarray(tmpCDFL)
    newLoc.cdfR = numpy.asarray(tmpCDFR)
    
    return newLoc
    
def PrevInts(intronRanges,i):
    summer = 0
    for j in range(i):
        summer += intronRanges[j][1] - intronRanges[j][0]
    return summer
    
def recoverAlignSOI(newLoc,intronRanges):
    if newLoc.SOIrange[1] < intronRanges[0][0]:
        shiftRight = 0
    else:
        shiftRight = intronRanges[0][1] - intronRanges[0][0]
        for i in range(1,len(intronRanges)):
            if newLoc.SOIrange[1] > intronRanges[i][0] - PrevInts(intronRanges,i):
                shiftRight += intronRanges[i][1] - intronRanges[i][0]
            else:
                break
    newLoc.SOIrange = [newLoc.SOIrange[0] + shiftRight,newLoc.SOIrange[1] + shiftRight]
    return newLoc
    
def rangeFinder(vector):
    ranges = []
    start = -1
    for i in range(len(vector)):
        if vector[i] == True:
            if start == -1:
                start = i
        if vector[i] == False:
            if start != -1:
                ranges.append([start,i-1])
                start = -1
    if start != -1:
        ranges.append([start,len(vector)-1])
    return ranges
            
def ratioOverlap(Seq1range,Seq2range):
    if Seq1range[0] < Seq2range[0]:
        if Seq1range[1] <= Seq2range[0]:
            return 0.0
        else:
            
            return (Seq1range[1] - Seq2range[0])/(1.0*Seq1range[1] - Seq1range[0])
    else:
        if Seq2range[1] <= Seq1range[0]:
            return 0.0
        else:
            return (Seq2range[1] - Seq1range[0])/(1.0*Seq1range[1] - Seq1range[0])