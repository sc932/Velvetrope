import numpy
import matplotlib.pylab as plt
import logging

def GeneByGeneOutput(setOfAligns,parameters):
    logging.info('Making the GBG plots...')
    h = 2.0 # plot hieght
    l = 4.0 # plot width
    MB = 0.05 # percent buffer from side
    LB = 0.01 # percent buffer for corner
    seq1len = setOfAligns.seqOfInt.length
    colors = ['b','g','r','c','m']
    
    plotNum = 0
    plots = 0
    
    GBGdata = open(parameters['OUTFILE'] + '/images/GBGdata.txt','wb')
    GBGdata.write('GBGX.pdf:ts1!ts2!ts3!ts4!\n')
    
    numToExpect = 0
    for pAlign in setOfAligns.aligns:
        if len(pAlign.localAligns) > 0:
            numToExpect += 1
    
    for pAlign in setOfAligns.aligns:
        if len(pAlign.localAligns) > 0:
            if plotNum%4 == 0:
                fig = plt.figure(figsize=(22,17))
                plots += 1
                plt.figtext(0.5, 0.965,  'Local Alignments for soi: ' + setOfAligns.seqOfInt.name + ', test seqs (' + str(1+(plots-1)*4) + '-' + str(min(4+(plots-1)*4,numToExpect)) + '/' + str(numToExpect), ha='center', color='black', weight='bold', size='large')
                gbgstring = 'GBG' + str(plots) + '.pdf:'
                ax = []
                ax.append(fig.add_subplot(411))
                plotNum = 0
            else:
                ax.append(fig.add_subplot(411+plotNum%4))
            seq2len = pAlign.testSeq.length
            seq1las = []
            seq2las = []
            seqRatio = []
            scoreRolling = 0
            for locA in pAlign.localAligns:
                
                for i in range(len(locA.SOIrange)):
                    seq1las.append(locA.SOIrange[i])
                    seq2las.append(locA.TSrange[i])
                    if setOfAligns.seqOfInt.dataDepths != 0:
                        scoreSum = sum(locA.fullScoreVec[locA.SOIrange[i][0]:locA.SOIrange[i][1]])
                    else:
                        scoreSum = sum(locA.fullScoreVec[locA.SOIrange[i][0]:locA.SOIrange[i][1]])
                    seqRatio.append(scoreSum/float(locA.SOIrange[i][1]-locA.SOIrange[i][0]))
                    scoreRolling += scoreSum
            totRatio = [scoreRolling/float(seq1len),scoreRolling/float(seq2len)]

            if seq1len <= seq2len:
                rL = 0.3 # relitive y start of lower
                rU = 0.7 # relitive y start of higher
                ratio = (seq2len - seq1len)/float(2*seq2len)
                s1xm = MB + ratio
                s1xx = 1 - MB - ratio
            else:
                rL = 0.7 # relitive y start of lower
                rU = 0.3 # relitive y start of higher
                ratio = (seq1len - seq2len)/float(2*seq1len)
                s1xm = MB + ratio
                s1xx = 1 - MB - ratio
            # top, smaller line
            ax[plotNum].axhline(y = rU*h, xmin = s1xm, xmax = s1xx, color='black')
            ax[plotNum].plot([s1xm*l,s1xm*l],[rU*h - 3*MB, rU*h + 3*MB],color='black') # left bound
            ax[plotNum].plot([s1xx*l,s1xx*l],[rU*h - 3*MB, rU*h + 3*MB],color='black') # right bound
            ax[plotNum].axhline(y = rU*h - 3*MB, xmin = s1xm, xmax = s1xm + LB, color='black') # left lower corner
            ax[plotNum].axhline(y = rU*h + 3*MB, xmin = s1xm, xmax = s1xm + LB, color='black') # left upper corner
            ax[plotNum].axhline(y = rU*h - 3*MB, xmin = s1xx, xmax = s1xx - LB, color='black') # right lower corner
            ax[plotNum].axhline(y = rU*h + 3*MB, xmin = s1xx, xmax = s1xx - LB, color='black') # right upper corner
            
            # bottom, larger line
            ax[plotNum].axhline(y = rL*h, xmin = MB, xmax = 1 - MB, color='black') # main line
            ax[plotNum].plot([MB*l,MB*l],[rL*h - 3*MB, rL*h + 3*MB],color='black') # left bound
            ax[plotNum].plot([(1 - MB)*l,(1 - MB)*l],[rL*h - 3*MB, rL*h + 3*MB],color='black') # right bound
            ax[plotNum].axhline(y = rL*h - 3*MB, xmin = MB, xmax = MB + LB, color='black') # left lower corner
            ax[plotNum].axhline(y = rL*h + 3*MB, xmin = MB, xmax = MB + LB, color='black') # left upper corner
            ax[plotNum].axhline(y = rL*h - 3*MB, xmin = 1 - MB, xmax = 1 - (MB + LB), color='black') # right lower corner
            ax[plotNum].axhline(y = rL*h + 3*MB, xmin = 1 - MB, xmax = 1 - (MB + LB), color='black') # right upper corner
            if seq1len <= seq2len:
                ax[plotNum].text(MB*l - MB/2*l, rL*h, '0',verticalalignment='center') # left text, bot
                ax[plotNum].text((1 - MB)*l + MB/4*l, rL*h + 3*MB, str(seq2len),verticalalignment='center') # right text1, bot
                ax[plotNum].text((1 - MB)*l + MB/4*l, rL*h - 3*MB, str(totRatio[1])[0:5],verticalalignment='center') # right text2, bot
                ax[plotNum].text(s1xm*l - MB/2*l, rU*h, '0',verticalalignment='center') # left text, top
                ax[plotNum].text(s1xx*l + MB/4*l, rU*h + 3*MB, str(seq1len),verticalalignment='center') # right text1, top
                ax[plotNum].text(s1xx*l + MB/4*l, rU*h - 3*MB, str(totRatio[0])[0:5],verticalalignment='center') # right text2, top
            else:
                ax[plotNum].text(MB*l - MB/2*l, rL*h, '0',verticalalignment='center') # left text, bot
                ax[plotNum].text((1 - MB)*l + MB/4*l, rL*h + 3*MB, str(seq1len),verticalalignment='center') # right text1, bot
                ax[plotNum].text((1 - MB)*l + MB/4*l, rL*h - 3*MB, str(totRatio[0])[0:5],verticalalignment='center') # right text2, bot
                ax[plotNum].text(s1xm*l - MB/2*l, rU*h, '0',verticalalignment='center') # left text, top
                ax[plotNum].text(s1xx*l + MB/4*l, rU*h + 3*MB, str(seq2len),verticalalignment='center') # right text1, top
                ax[plotNum].text(s1xx*l + MB/4*l, rU*h - 3*MB, str(totRatio[1])[0:5],verticalalignment='center') # right text2, top
            # local alignments
            for i in range(len(seq1las)):
                if seq1len <= seq2len:
                    ax[plotNum].text(l*(MB + (1 - 2*MB)*seq2las[i][0]/float(seq2len) + MB + (1 - 2*MB)*seq2las[i][1]/float(seq2len) + s1xm + (1 - 2*MB)*seq1las[i][0]/float(seq2len) + s1xm + (1 - 2*MB)*seq1las[i][1]/float(seq2len))/4.0,(rL*h - 10*MB + rL*h - 10*MB + rU*h + 5*MB + rU*h + 5*MB)/4.0,str(seqRatio[i])[0:4],horizontalalignment='center')
                    ax[plotNum].axhspan(rU*h - 3*MB, rU*h + 3*MB, xmin = s1xm + (1 - 2*MB - 2*ratio)*seq1las[i][0]/float(seq1len), xmax = s1xm + (1 - 2*MB - 2*ratio)*seq1las[i][1]/float(seq1len), facecolor = colors[i%5], alpha = 0.35) # top seq
                    ax[plotNum].axhspan(rL*h - 3*MB, rL*h + 3*MB, xmin = MB + (1 - 2*MB)*seq2las[i][0]/float(seq2len), xmax = MB + (1 - 2*MB)*seq2las[i][1]/float(seq2len), facecolor = colors[i%5], alpha = 0.35) # bot seq
                    ax[plotNum].plot([l*(s1xm + (1 - 2*MB - 2*ratio)*seq1las[i][0]/float(seq1len)),l*(MB + (1 - 2*MB)*seq2las[i][0]/float(seq2len))],[rU*h - 3*MB,rL*h + 3*MB],color = 'black', alpha = 0.35) # left line
                    ax[plotNum].plot([l*(s1xm + (1 - 2*MB - 2*ratio)*seq1las[i][1]/float(seq1len)),l*(MB + (1 - 2*MB)*seq2las[i][1]/float(seq2len))],[rU*h - 3*MB,rL*h + 3*MB],color = 'black', alpha = 0.35) # right line
                    ax[plotNum].text(l*(MB + .01 + (1 - 2*MB)*seq2las[i][0]/float(seq2len)),rL*h - 10*MB,str(seq2las[i][0]),horizontalalignment='center') # lower left text
                    ax[plotNum].text(l*(MB - .01 + (1 - 2*MB)*seq2las[i][1]/float(seq2len)),rL*h - 10*MB,str(seq2las[i][1]),horizontalalignment='center') # lower right text
                    ax[plotNum].text(l*(s1xm + .01 + (1 - 2*MB - 2*ratio)*seq1las[i][0]/float(seq1len)),rU*h + 5*MB,str(seq1las[i][0]),horizontalalignment='center') # upper left text
                    ax[plotNum].text(l*(s1xm - .01 + (1 - 2*MB - 2*ratio)*seq1las[i][1]/float(seq1len)),rU*h + 5*MB,str(seq1las[i][1]),horizontalalignment='center') # upper right text
                else:
                    ax[plotNum].text(l*(MB + (1 - 2*MB)*seq1las[i][0]/float(seq1len) + MB + (1 - 2*MB)*seq1las[i][1]/float(seq1len) + s1xm + (1 - 2*MB)*seq2las[i][0]/float(seq1len) + s1xm + (1 - 2*MB)*seq2las[i][1]/float(seq1len))/4.0,(rL*h - 10*MB + rL*h - 10*MB + rU*h + 5*MB + rU*h + 5*MB)/4.0,str(seqRatio[i])[0:4],horizontalalignment='center') 
                    ax[plotNum].axhspan(rU*h - 3*MB, rU*h + 3*MB, xmin = s1xm + (1 - 2*MB - 2*ratio)*seq2las[i][0]/float(seq2len), xmax = s1xm + (1 - 2*MB - 2*ratio)*seq2las[i][1]/float(seq2len), facecolor = colors[i%5], alpha = 0.35) # top seq
                    ax[plotNum].axhspan(rL*h - 3*MB, rL*h + 3*MB, xmin = MB + (1 - 2*MB)*seq1las[i][0]/float(seq1len), xmax = MB + (1 - 2*MB)*seq1las[i][1]/float(seq1len), facecolor = colors[i%5], alpha = 0.35) # bot seq
                    ax[plotNum].plot([l*(s1xm + (1 - 2*MB - 2*ratio)*seq2las[i][0]/float(seq2len)),l*(MB + (1 - 2*MB)*seq1las[i][0]/float(seq1len))],[rU*h + 3*MB,rL*h - 3*MB],color = 'black', alpha = 0.35) # left line
                    ax[plotNum].plot([l*(s1xm + (1 - 2*MB - 2*ratio)*seq2las[i][1]/float(seq2len)),l*(MB + (1 - 2*MB)*seq1las[i][1]/float(seq1len))],[rU*h + 3*MB,rL*h - 3*MB],color = 'black', alpha = 0.35) # right line
                    ax[plotNum].text(l*(MB + .01 + (1 - 2*MB)*seq1las[i][0]/float(seq1len)),rL*h + 5*MB,str(seq1las[i][0]),horizontalalignment='center') # lower left text
                    ax[plotNum].text(l*(MB - .01 + (1 - 2*MB)*seq1las[i][1]/float(seq1len)),rL*h + 5*MB,str(seq1las[i][1]),horizontalalignment='center') # lower right text
                    ax[plotNum].text(l*(s1xm + .01 + (1 - 2*MB - 2*ratio)*seq2las[i][0]/float(seq2len)),rU*h - 10*MB,str(seq2las[i][0]),horizontalalignment='center') # upper left text
                    ax[plotNum].text(l*(s1xm - .01 + (1 - 2*MB - 2*ratio)*seq2las[i][1]/float(seq2len)),rU*h - 10*MB,str(seq2las[i][1]),horizontalalignment='center') # upper right text
            ax[plotNum].set_ylim(0,h)
            ax[plotNum].set_xlim(0,l)
            ax[plotNum].set_yticklabels([])
            ax[plotNum].set_xticklabels([])
            ax[plotNum].set_yticks([])
            ax[plotNum].set_xticks([])
            plt.xlabel(pAlign.testSeq.name + ' rf ' + str(pAlign.localAligns[0].offset), multialignment='center')
            gbgstring += pAlign.testSeq.name
            if (plotNum+1)%4 == 0 and plotNum > 1:
                gbgstring += '\n'
                GBGdata.write(gbgstring)
                fig.savefig(parameters['OUTFILE'] + '/images/GBG' + str(plots) + '.png', bbox_inches='tight', pad_inches=0.03)
                fig.savefig(parameters['OUTFILE'] + '/images/GBG' + str(plots) + '.pdf')
                logging.debug('Made GBG plot ' + str(plots))
            else:
                gbgstring += '!'
            plotNum += 1
    gbgstring += '\n'
    GBGdata.write(gbgstring)
    fig.savefig(parameters['OUTFILE'] + '/images/GBG' + str(plots) + '.png', bbox_inches='tight', pad_inches=0.03)
    fig.savefig(parameters['OUTFILE'] + '/images/GBG' + str(plots) + '.pdf')
    logging.debug('Made GBG plot ' + str(plots))
    GBGdata.close()
    #MakeLine(0.4,0.45,rU,ax[0],MB,h)
    #plt.show()
    
def MakeLine(start,stop,y,ax,MB,h):
    ax.axhspan(y*h - 3*MB, y*h + 3*MB, xmin = MB + (1 - 2*MB)*start, xmax = MB + (1 - 2*MB)*stop, facecolor = 'black', alpha = 0.35) # bot seq