import numpy
import logging
import sys
from Bio.Seq import Seq
from Bio import SeqIO
from PARAMS import *
from VelvetropeClasses import *
import time

################
# Clubs Reader #
################
"""
    This program reads in an input file and parses data to be read into clubs_generator.py
    It is designed to read in files in fasta format using biopython
        >gene1
        SATANSSPAWN...
        >gene2
        SATANSPRAWN...
    It returns seqOfInt, testSeqSet, each of the form:
        seqOfInt = [[1,4,5,...],[17,2,5,...],...] where len(seqOfInt) = seqOfInt.readingframes
        and the numbers are the 5/8-bit translations of the AA/DNA data in numpy arrays
"""

def reader(parameters):
    """
        Reads in data from parameters['INPUT_FILE'] defined in Velvetrope.py
    """

    logging.info('Reading data...')
    t0 = time.time()
    
    if parameters['SEQ_TYPE'] == 'AA': # read in amino acid
        data = readInAA(parameters['INPUT_FILE'],parameters['REF_GENE'])
    elif parameters['SEQ_TYPE'][0:3] == 'DNA': # read in DNA 1,3,6 reading frames
        data = readInDNA(parameters)
    else:
        logging.critical('Unrecognized filetype: ' + parameters['SEQ_TYPE'] + '! Expecting AA or DNA')
        return 0
    
    logging.info('Reading took ' + str(time.time() - t0) + ' seconds.')
    
    return data # (seqOfInt, testSeqSet) kick data back to Velvetrope.py
    
def readInAA(infile,refGene):
    """
        Reads in data from infile containing AA data
        Reads the file in using biopython and then converts the sequence using AA_TO_5BIT from PARAMS.py
    """
    testSeqSet = []
    seqOfInt = 0
    handle = open(infile)
    parsed = SeqIO.parse(handle, "fasta")
    for seq_record in parsed:
        convertedSeq = []
        for amino in seq_record.seq:
            if amino in AA_TO_5BIT.keys():
                convertedSeq.append(AA_TO_5BIT[amino])
        if seq_record.description == refGene:
            seqOfInt = Sequence('AA',AA_DATA_DEPTHS,AA_READING_FRAMES,refGene,convertedSeq)
        else:
            testSeqSet.append(Sequence('AA',AA_DATA_DEPTHS,AA_READING_FRAMES,seq_record.description,convertedSeq))
    if seqOfInt == 0:
        logging.critical('Never found sequence of interest: ' + refGene)
        return 0
    return seqOfInt, testSeqSet
    
def readInDNA(parameters):
    """
        Reads in data from infile containing DNA data in either single frame, 3 or 6 frame translations.
        Reads the file in using biopython and then converts the sequence using DNA_TO_5BIT from PARAMS.py
    """
    infile = parameters['INPUT_FILE']
    refGene = parameters['REF_GENE']
    seqT = parameters['SEQ_TYPE']
    
    handle = open(infile)
    
    parsed = SeqIO.parse(handle, "fasta")
    
    testSeqSet = []
    seqOfInt = 0
    if seqT == 'DNA6': # 6 frame translation
        for seq_record in parsed:
            seqLen = len(seq_record.seq)
            convertedSeq = [[],[],[],[],[],[]]
            for j in range(3):
                convertedSeq[j] = numpy.zeros((seqLen-j)/3 - j, dtype = numpy.uint8) # only need a byte per AA
                convertedSeq[j+3] = numpy.zeros((seqLen-j)/3 - j, dtype = numpy.uint8) # only need a byte per AA
            theSeq = seq_record.seq[:]
            theSeqRev = theSeq[::-1] # only need to flip the array once
            for j in range(3):
                # read forward 3 frames
                k = 0
                for i in range(j,(seqLen-j)/3):
                    if str(theSeq[i*3+j:(i+1)*3+j]) in DNA_TO_8BIT.keys():
                        convertedSeq[j][k] = DNA_TO_8BIT[str(theSeq[i*3+j:(i+1)*3+j])]
                        k += 1
                    else:
                        pass # could not find AA, ie AGX
                # read backwards 3 frames
                k = 0
                for i in range(j,(seqLen-j)/3):
                    if str(theSeqRev[i*3+j:(i+1)*3+j]) in DNA_TO_8BIT.keys():
                        convertedSeq[j+3][k] = DNA_TO_8BIT[str(theSeqRev[i*3+j:(i+1)*3+j])]
                        k += 1
                    else:
                        pass # could not find AA, ie AGX
            if seq_record.description == refGene:
                seqOfInt = Sequence('DNA6',DNA_DATA_DEPTHS,6,refGene,convertedSeq)
            else:
                testSeqSet.append(Sequence('DNA6',DNA_DATA_DEPTHS,6,seq_record.description,convertedSeq))
    elif seqT == 'DNA3': # just read in the forward 3 frames
        for seq_record in parsed:
            seqLen = len(seq_record.seq)
            convertedSeq = [[],[],[]]
            for j in range(3):
                convertedSeq[j] = numpy.zeros((seqLen-j)/3 - j, dtype = numpy.uint8) # only need a byte per AA
            theSeq = seq_record.seq[:]
            for j in range(3):
                # read forward 3 frames
                k = 0
                for i in range(j,(seqLen-j)/3):
                    if str(theSeq[i*3+j:(i+1)*3+j]) in DNA_TO_8BIT.keys():
                        convertedSeq[j][k] = DNA_TO_8BIT[str(theSeq[i*3+j:(i+1)*3+j])]
                        k += 1
                    else:
                        pass # could not find AA, ie AGX
            if seq_record.description == refGene:
                seqOfInt = Sequence('DNA3',DNA_DATA_DEPTHS,3,refGene,convertedSeq)
            else:
                testSeqSet.append(Sequence('DNA3',DNA_DATA_DEPTHS,3,seq_record.description,convertedSeq))
    else: # just read in a single frame
        for seq_record in parsed:
            seqLen = len(seq_record.seq)
            convertedSeq = numpy.zeros(seqLen/3, dtype = numpy.uint8) # only need a byte per AA
            # read forward single frame
            k = 0
            for i in range(0,(len(seq_record.seq))/3):
                if str(seq_record.seq[i*3:(i+1)*3]) in DNA_TO_8BIT.keys():
                    convertedSeq[k] = DNA_TO_8BIT[str(seq_record.seq[i*3:(i+1)*3])]
                    k += 1
                else:
                    pass # could not find AA, ie AGX
            if seq_record.description == refGene:
                seqOfInt = Sequence('DNA',DNA_DATA_DEPTHS,1,refGene,convertedSeq)
            else:
                testSeqSet.append(Sequence('DNA',DNA_DATA_DEPTHS,1,seq_record.description,convertedSeq))
        
    if seqOfInt == 0:
        logging.critical('Never found sequence of interest: ' + refGene)
        return 0
    return seqOfInt, testSeqSet