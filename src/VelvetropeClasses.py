import numpy

class Sequence():
    """Sequence Class"""
    def __init__(self, sType, dDepths, rFrames, nm, theSeq, linksto = 0):
        self.sequenceType = sType # sequence type, ie DNA, AA, etc
        self.dataDepths = dDepths # bit-depths where the info is, ie for DNA [5,8]
        self.readingFrames = rFrames # number of reading frames, ie HIV/DNA is usually 6
        self.name = nm # sequence name, ie 'E.Coli K12'
        self.seq = theSeq
        self.linksTo = linksto # if sequence had to be chopped up, link to this sequence next
        if self.readingFrames != 1:
            self.length = len(theSeq[0])
        else:
            self.length = len(theSeq)
            
class SetOfAlignments():
    """Set of Alignments"""
    def __init__(self, soi):
        self.seqOfInt = soi
        self.aligns = []
        self.testSeqs = []
        if self.seqOfInt.dataDepths != [0]:
            if self.seqOfInt.readingFrames > 1:
                self.scoreSum = numpy.zeros((len(self.seqOfInt.dataDepths),len(self.seqOfInt.seq[0])), dtype = numpy.int8)
                self.inClubSum = self.inClub = numpy.zeros(len(self.seqOfInt.seq[0]), dtype = numpy.int16)
            else:
                self.scoreSum = numpy.zeros((len(self.seqOfInt.dataDepths),len(self.seqOfInt.seq)), dtype = numpy.int8)
                self.inClubSum = self.inClub = numpy.zeros(len(self.seqOfInt.seq), dtype = numpy.int16)
        else:
            self.scoreSum = numpy.zeros((1,len(self.seqOfInt.seq)), dtype = numpy.int8)
            self.inClubSum = self.inClub = numpy.zeros(len(self.seqOfInt.seq), dtype = numpy.int16)
    def AddAlign(self, align): # add an alignment
        self.aligns.append(align)
        self.testSeqs.append(align.testSeq)
        self.inClubSum += align.inClub
        self.scoreSum += align.score
            
class Alignment():
    """Alignment Class"""
    def __init__(self, s1, s2, offset = 0):
        self.seqOfInt = s1 # sequence of interest
        self.testSeq = s2 # test sequence
        self.localAligns = [] # initialize local alignments
        self.offset = offset
        # initialize the club and scores
        if self.seqOfInt.dataDepths != [0]:
            if self.seqOfInt.readingFrames > 1:
                self.score = numpy.zeros((len(self.seqOfInt.dataDepths),len(self.seqOfInt.seq[0])), dtype = numpy.bool)
                self.inClub = numpy.zeros(len(self.seqOfInt.seq[0]), dtype = numpy.bool)
            else:
                self.score = numpy.zeros((len(self.seqOfInt.dataDepths),len(self.seqOfInt.seq)), dtype = numpy.bool)
                self.inClub = numpy.zeros(len(self.seqOfInt.seq), dtype = numpy.bool)
        else:
            self.score = numpy.zeros((len(self.seqOfInt.dataDepths),len(self.seqOfInt.seq)), dtype = numpy.bool)
            self.inClub = numpy.zeros(len(self.seqOfInt.seq), dtype = numpy.bool)
    def AddLocalAlign(self,locA): # add a local alignment to the list and club/score
        self.localAligns.append(locA)
        #for i in range(len(locA.SOIrange)):
        #    self.inClub[locA.SOIrange[i][0]:locA.SOIrange[i][1]] += numpy.ones((locA.SOIrange[i][1]-locA.SOIrange[i][0]))
        #self.score += locA.score
        
class LocalAlignment():
    """Local Alignment Class"""
    def __init__(self,s1range=[[0,0]],s2range=[[0,0]],fullRangeSOI=[0,0],fullRangeTS=[0,0],offset=0,score=[0],cdfL=[0],cdfR=[0],fullScoreVec=[0]):
        self.SOIrange = s1range # range in soi where the loc align is
        self.TSrange = s2range # range in the test seq where it is
        self.fSOIrange = fullRangeSOI
        self.fTSrange = fullRangeTS
        self.score = score # score in range
        self.cdfL = cdfL
        self.cdfR = cdfR
        self.fullScoreVec = fullScoreVec
        self.offset = offset