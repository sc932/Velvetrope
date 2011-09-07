import numpy
import logging
import sys
from Bio.Seq import Seq
from Bio import SeqIO
from PARAMS import *
from VelvetropeClasses import *
import time
import clubs_output_genes

def generateStdOut(setOfAligns,parameters):
    output = open(parameters['OUTFILE'] + '/VRstdOut.txt', 'wb')
    output.write("Velvetrope standard output\n")
    output.write("Parameters:\n")
    output.write(str(parameters) + '\n\n')
    output.write("Regular expression data from Sequence of Interest (SOI):\n")
    regs = clubs_output_genes.RegExpOut(setOfAligns, parameters)
    #print regs
    for i in range(len(regs)):
        output.write(regs[i] + '\n')
    output.write('\n')
    for align in setOfAligns.aligns:
        if len(align.localAligns) > 0:
            output.write('> ' + align.testSeq.name + '\n')
            output.write(' Length: ' + str(align.testSeq.length) + '\n\n')
        for locA in align.localAligns:
            
            for i in range(len(locA.SOIrange)):
                stringSOI = ''
                stringLink = ''
                stringTS = ''
                matchTotal = 0
                matchExact = 0
                for j in range(locA.SOIrange[i][1]-locA.SOIrange[i][0]):
                    if align.seqOfInt.dataDepths[0] == 0:
                        stringSOI += FIVEBIT_TO_AA[align.seqOfInt.seq[j + locA.SOIrange[i][0]]]
                        stringTS += FIVEBIT_TO_AA[align.testSeq.seq[j + locA.TSrange[i][0]]]
                        if align.seqOfInt.seq[j + locA.SOIrange[i][0]] == align.testSeq.seq[j + locA.TSrange[i][0]]:
                            stringLink += '|'
                            matchTotal += 1
                        else:
                            stringLink += ' '
                    else:
                        if align.seqOfInt.readingFrames == 1:
                            stringSOI += FIVEBIT_TO_AA[align.seqOfInt.seq[j + locA.SOIrange[i][0]]%2**align.seqOfInt.dataDepths[0]]
                            stringTS += FIVEBIT_TO_AA[align.testSeq.seq[j + locA.TSrange[i][0]]%2**align.seqOfInt.dataDepths[0]]
                            if align.seqOfInt.seq[j + locA.SOIrange[i][0]] == align.testSeq.seq[j + locA.TSrange[i][0]]:
                                stringLink += '!'
                                matchTotal += 1
                                matchExact += 1
                            elif align.seqOfInt.seq[j + locA.SOIrange[i][0]]%2**align.seqOfInt.dataDepths[0] == align.testSeq.seq[j + locA.TSrange[i][0]]%2**align.seqOfInt.dataDepths[0]:
                                stringLink += '|'
                                matchTotal += 1
                            else:
                                stringLink += ' '
                        else:
                            stringSOI += FIVEBIT_TO_AA[align.seqOfInt.seq[0][j + locA.SOIrange[i][0]]%2**align.seqOfInt.dataDepths[0]]
                            stringTS += FIVEBIT_TO_AA[align.testSeq.seq[locA.offset][j + locA.TSrange[i][0]]%2**align.seqOfInt.dataDepths[0]]
                            if align.seqOfInt.seq[0][j + locA.SOIrange[i][0]] == align.testSeq.seq[locA.offset][j + locA.TSrange[i][0]]:
                                stringLink += '!'
                                matchTotal += 1
                                matchExact += 1
                            elif align.seqOfInt.seq[0][j + locA.SOIrange[i][0]]%2**align.seqOfInt.dataDepths[0] == align.testSeq.seq[locA.offset][j + locA.TSrange[i][0]]%2**align.seqOfInt.dataDepths[0]:
                                stringLink += '|'
                                matchTotal += 1
                            else:
                                stringLink += ' '
                if align.seqOfInt.dataDepths[0] == 0:
                    output.write('Match of length ' + str(locA.TSrange[i][1]-locA.TSrange[i][0]) + '\nfrom test sequence region: [' + str(locA.TSrange[i][0]) + ',' + str(locA.TSrange[i][1]) + '] to sequence of interest region [' + str(locA.SOIrange[i][0]) + ',' + str(locA.SOIrange[i][1]) + ']\nwith identity density = ' + str(matchTotal/float(locA.SOIrange[i][1]-locA.SOIrange[i][0])) + '\n')
                else:
                    output.write('Match of length ' + str(locA.TSrange[i][1]-locA.TSrange[i][0]) + '\nfrom test sequence region: [' + str(locA.TSrange[i][0]) + ',' + str(locA.TSrange[i][1]) + '] to sequence of interest region [' + str(locA.SOIrange[i][0]) + ',' + str(locA.SOIrange[i][1]) + ']\n with AA identity density = ' + str(matchTotal/float(locA.SOIrange[i][1]-locA.SOIrange[i][0])) + '\nexact (codon) identity density = ' + str(matchExact/float(locA.SOIrange[i][1]-locA.SOIrange[i][0])) + '\n')
                if len(str(locA.SOIrange[i][0])) > len(str(locA.TSrange[i][0])):
                    spacerTS = ' '
                    spacerLink = ' '
                    for j in range(len(str(locA.TSrange[i][0]))):
                        spacerLink += ' '
                    spcDiff = len(str(locA.SOIrange[i][0])) - len(str(locA.TSrange[i][0]))
                    for j in range(spcDiff):
                        spacerTS += ' '
                        spacerLink += ' '
                    output.write(str(locA.SOIrange[i][0]) + ' ' + stringSOI + ' ' + str(locA.SOIrange[i][1]) + '\n')
                    output.write(spacerLink + stringLink + '\n')
                    output.write(str(locA.TSrange[i][0]) + spacerTS + stringTS + ' ' + str(locA.TSrange[i][1]) + '\n')
                else:
                    spacerSOI = ' '
                    spacerLink = ' '
                    for j in range(len(str(locA.SOIrange[i][0]))):
                        spacerLink += ' '
                    spcDiff = len(str(locA.TSrange[i][0])) - len(str(locA.SOIrange[i][0]))
                    for j in range(spcDiff):
                        spacerSOI += ' '
                        spacerLink += ' '
                    output.write(str(locA.SOIrange[i][0]) + spacerSOI + stringSOI + ' ' + str(locA.SOIrange[i][1]) + '\n')
                    output.write(spacerLink + stringLink + '\n')
                    output.write(str(locA.TSrange[i][0]) + ' ' + stringTS + ' ' + str(locA.TSrange[i][1]) + '\n')
                output.write("\n")
    output.close()