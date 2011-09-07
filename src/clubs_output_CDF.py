import numpy
import matplotlib.pylab as plt
import logging

def GeneCDFplots(setOfAligns, numGenes,parameters):
    logging.info('Making the CDF plots (This may take a while)...')
    lensoi = setOfAligns.seqOfInt.length
    colors = ['b','g','r','c','m']
    
    if numGenes > len(setOfAligns.testSeqs)/setOfAligns.seqOfInt.readingFrames or numGenes < 1:
        numGenes = len(setOfAligns.testSeqs)*setOfAligns.seqOfInt.readingFrames
    geneCount = 0
    
    plots = 0
    
    CDFdata = open(parameters['OUTFILE'] + '/images/CDFdata.txt','wb')
    CDFdata.write('CDFX.pdf:soi!testseq\n')
    
    for pAlign in setOfAligns.aligns:
        if len(pAlign.localAligns) > 0:
            if geneCount < numGenes:
                geneCount += 1
                plotNum = 0
                totalForLocA = 0
                for locA in pAlign.localAligns:
                    totalForLocA += 1
                    # get density ratio
                    seqRatio = []
                    if setOfAligns.seqOfInt.dataDepths != 0:
                        for i in range(len(locA.SOIrange)):
                            scoreSum = (sum(locA.fullScoreVec[locA.SOIrange[i][0]:locA.SOIrange[i][1]]))
                            seqRatio.append(scoreSum/float(locA.SOIrange[i][1]-locA.SOIrange[i][0]))
                    else:
                        for i in range(len(locA.SOIrange)):
                            scoreSum = (sum(locA.fullScoreVec[locA.SOIrange[i][0]:locA.SOIrange[i][1]]))
                            seqRatio.append(scoreSum/float(locA.SOIrange[i][1]-locA.SOIrange[i][0]))
                            
                    # make the empty plot, 4 per page
                    if plotNum%4 == 0:
                        fig = plt.figure(figsize=(22,17))
                        plt.figtext(0.5, 0.965,  'Local Alignments for soi: ' + setOfAligns.seqOfInt.name + ' vs ' + pAlign.testSeq.name + ' rf ' + str(pAlign.localAligns[0].offset), ha='center', color='black', weight='bold', size='large')
                        
                        ax = []
                        ax.append(fig.add_subplot(811))
                        ax.append(fig.add_subplot(812))
                        plotNum = 0
                        plots += 1
                        
                    else:
                        ax.append(fig.add_subplot(811+(plotNum%4)*2))
                        ax.append(fig.add_subplot(811+(plotNum%4)*2+1))
                        
                    # make the upper plot, the match vector
                    temp = numpy.zeros(lensoi)
                    for i in range(len(locA.SOIrange)):
                        temp[locA.SOIrange[i][0]:locA.SOIrange[i][1]] = 2*numpy.ones(locA.SOIrange[i][1]-locA.SOIrange[i][0])
                    ax[plotNum*2].pcolor(numpy.reshape(numpy.array(locA.fullScoreVec) + temp,(1,lensoi)))
                    
                    # plot the cdfs
                    ax[plotNum*2+1].plot(locA.cdfL)
                    ax[plotNum*2+1].plot(locA.cdfR)
                    for i in range(len(locA.SOIrange)):
                        ax[plotNum*2+1].axvspan(locA.SOIrange[i][0],locA.SOIrange[i][1], facecolor = 'r', alpha = 0.35)
                    
                    # mark off the intron area
                    intronArea = numpy.asarray(numpy.asarray(locA.cdfL) == numpy.asarray(locA.cdfR),dtype=numpy.int32)
                    ranges = rangeFinder(intronArea)
                    
                    for i in range(len(ranges)):
                        ax[plotNum*2+1].axvspan(ranges[i][0],ranges[i][1], facecolor = 'm', alpha = 0.35)
                    
                    #ax[plotNum*2+1].text((locA.SOIrange[0]+locA.SOIrange[1])/2.0,(max(max(locA.cdfL),max(locA.cdfR)) + min(min(locA.cdfL),min(locA.cdfR)))/2*1.66,str(seqRatio)[0:4],verticalalignment='center',horizontalalignment='center')
                    
                    
                    # make the tick marks
                    ax[plotNum].set_yticklabels([])
                    ax[plotNum*2].set_yticks([])
                    prevTicks = []
                    prevLabs = []
                    for i in range(len(locA.SOIrange)):
                        prevTicks.extend([locA.SOIrange[i][0],(locA.SOIrange[i][0]+locA.SOIrange[i][1])/2.0,locA.SOIrange[i][1]])
                        prevLabs.extend([str(locA.SOIrange[i][0]),str(seqRatio[i])[0:4],str(locA.SOIrange[i][1])])
                    newXticks, newXtickLabels = intronAdderXticks(prevTicks,prevLabs,ranges)
                    ax[plotNum*2].set_xticks(newXticks)
                    ax[plotNum*2].set_xticklabels(newXtickLabels)
                    #ax[plotNum*2].set_xticks([locA.SOIrange[0],(locA.SOIrange[0]+locA.SOIrange[1])/2.0,locA.SOIrange[1]])
                    #ax[plotNum*2].set_xticklabels([str(locA.SOIrange[0]),str(seqRatio)[0:4],str(locA.SOIrange[1])])
                    ax[plotNum*2+1].set_yticklabels([])
                    ax[plotNum*2+1].set_yticks([])
                    ax[plotNum*2].set_xlim(0,lensoi)
                    ax[plotNum*2+1].set_xlim(0,lensoi)
                    plotNum += 1
                    #plt.show()
                    if ((plotNum)%4 == 0 and plotNum > 1) or totalForLocA == len(pAlign.localAligns):
                        #print plotNum, (plotNum)%4, totalForLocA, len(pAlign.localAligns)
                        CDFdata.write('CDF' + str(plots) + '.pdf:' + setOfAligns.seqOfInt.name + '!' + pAlign.testSeq.name + '\n')
                        fig.savefig(parameters['OUTFILE'] + '/images/CDF' + str(plots) + '.png', bbox_inches='tight', pad_inches=0.03)
                        fig.savefig(parameters['OUTFILE'] + '/images/CDF' + str(plots) + '.pdf')
                        logging.debug('Made CDF plot ' + str(plots) + '...')
                    #return 0
            else:
                #plt.show()
                CDFdata.close()
                return 0
    CDFdata.close()
    #plt.show()
    
def intronAdderXticks(prevTicks,prevLabs,ranges):
    for i in range(len(ranges)):
        prevTicks.extend(ranges[i])
        prevLabs.extend([str(ranges[i][0]),str(ranges[i][1])])
    return prevTicks,prevLabs

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