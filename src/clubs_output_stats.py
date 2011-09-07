import numpy
from PARAMS import *
from VelvetropeClasses import *
import matplotlib.pylab as plt

def generateStatsOut(setOfAligns,parameters):
    logging.info('Making the stat plots...')
    runLengths = []
    densitys = []
    exacts = []
    for align in setOfAligns.aligns:
        for locA in align.localAligns:
            for i in range(len(locA.SOIrange)):
                if align.seqOfInt.readingFrames == 1:
                    runLengths.append(locA.SOIrange[i][1]-locA.SOIrange[i][0])
                    matchTotal = 0
                    matchExact = 0
                    for j in range(locA.SOIrange[i][1]-locA.SOIrange[i][0]):
                        if align.seqOfInt.dataDepths[0] == 0:
                            if align.seqOfInt.seq[j + locA.SOIrange[i][0]] == align.testSeq.seq[j + locA.TSrange[i][0]]:
                                matchTotal += 1
                        else:
                            if align.seqOfInt.seq[j + locA.SOIrange[i][0]] == align.testSeq.seq[j + locA.TSrange[i][0]]:
                                matchTotal += 1
                                matchExact += 1
                            elif align.seqOfInt.seq[j + locA.SOIrange[i][0]]%2**align.seqOfInt.dataDepths[0] == align.testSeq.seq[j + locA.TSrange[i][0]]%2**align.seqOfInt.dataDepths[0]:
                                matchTotal += 1
                    densitys.append(matchTotal/float(locA.SOIrange[i][1]-locA.SOIrange[i][0]))
                    exacts.append(matchExact/float(locA.SOIrange[i][1]-locA.SOIrange[i][0]))
                else:
                    runLengths.append(locA.SOIrange[i][1]-locA.SOIrange[i][0])
                    matchTotal = 0
                    matchExact = 0
                    for j in range(locA.SOIrange[i][1]-locA.SOIrange[i][0]):
                        if align.seqOfInt.dataDepths[0] == 0:
                            if align.seqOfInt.seq[0][j + locA.SOIrange[i][0]] == align.testSeq.seq[locA.offset][j + locA.TSrange[i][0]]:
                                matchTotal += 1
                        else:
                            if align.seqOfInt.seq[0][j + locA.SOIrange[i][0]] == align.testSeq.seq[locA.offset][j + locA.TSrange[i][0]]:
                                matchTotal += 1
                                matchExact += 1
                            elif align.seqOfInt.seq[0][j + locA.SOIrange[i][0]]%2**align.seqOfInt.dataDepths[0] == align.testSeq.seq[locA.offset][j + locA.TSrange[i][0]]%2**align.seqOfInt.dataDepths[0]:
                                matchTotal += 1
                    densitys.append(matchTotal/float(locA.SOIrange[i][1]-locA.SOIrange[i][0]))
                    exacts.append(matchExact/float(locA.SOIrange[i][1]-locA.SOIrange[i][0]))
                
    BINSIZE = 100
    
    bin = numpy.max(runLengths)/BINSIZE+1
    runLenDist = numpy.zeros(BINSIZE)
    for i in range(len(runLengths)):
        runLenDist[runLengths[i]/bin] += 1
        
    denDist = numpy.zeros(BINSIZE)
    exactDist = numpy.zeros(BINSIZE)
    for i in range(len(densitys)):
        denDist[int(densitys[i]*BINSIZE)] += 1
        exactDist[int(exacts[i]*BINSIZE)] += 1
                
    fig = plt.figure(figsize=(17,7))
    ax = fig.add_subplot(121)
    ax.plot(runLengths,densitys,'b*')
    ax.plot(runLengths,exacts,'r*')
    ax.set_ylim(0,1)
    ax.set_title('Run length vs density (AA and exact)')
    ax.set_xlabel('Run Length')
    ax.set_ylabel('Density (AA = blue, exact = red)')
    
    ax = fig.add_subplot(122)
    ax.bar(numpy.arange(len(runLenDist)),runLenDist,color='blue',align='center')
    ax.set_title('Run length distribution')
    ax.set_xlabel('Run Length')
    labs = []
    for i in range(0,10):
        labs.append(str(len(runLenDist)/10*i*3))
    ax.set_xticklabels(labs)
    ax.set_ylabel('Number of runs at that run length (bin size = ' + str(bin) + ')')
    
    fig.savefig(parameters['OUTFILE'] + '/images/STATS_OUT.png', bbox_inches='tight', pad_inches=0.03)
    fig.savefig(parameters['OUTFILE'] + '/images/STATS_OUT.pdf')