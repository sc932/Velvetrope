import matplotlib.pylab as plt
import numpy
from PARAMS import *

def makeLatexMat(align,w,h,offset):
    locA = align.localAligns[0]
    subsoi = []
    subts = []
    for i in range(-h/2,w+h/2):
        subts.append(FIVEBIT_TO_AA[align.testSeq.seq[0][locA.TSrange[0]+i+offset]%2**5])
    for i in range(w):
        subsoi.append(FIVEBIT_TO_AA[align.seqOfInt.seq[0][locA.SOIrange[0]+i+offset]%2**5])
    print subts
    print subsoi
    string = "\\begin{tabular}{|"
    for i in range(w+2):
        string = string + 'c'
    string = string + '|}'
    print string
    print '\hline'
    string = "\\footnotesize{residue} &"
    for i in range(len(subsoi)):
        string = string + " \\tiny{" + str(locA.SOIrange[0]+i+offset) + "} & "
    string = string + " \\tiny{}$\\cdots$ \\\\"
    print string
    string = "$\\cdots$ &"
    for i in range(len(subsoi)):
        string = string + " \\textcolor{blue}{" + subsoi[i] + "} & "
    string = string + " $\\cdots$ \\\\"
    print string
    print '\hline'
    string = "\\footnotesize{row} &"
    for i in range(len(subsoi)):
        string = string + " $\\vdots$ & "
    string = string + "\\\\"
    print string
    for i in range(h):
        string = "\\tiny{" + str(locA.TSrange[0] - h/2 + i) + "}$\\cdots$ &"
        for j in range(w):
            if subts[j+i] == subsoi[j]:
                string = string + " \\textcolor{red}{" + subts[j+i] + "} "
            else:
                string = string + " " + subts[j+i] + " "
            if j < w - 1:
                string = string + "&"
            else:
                string = string + ' & $\\cdots$ \\\\'
        print string
    string = " &"
    for i in range(len(subsoi)):
        string = string + " $\\vdots$ & "
    string = string + "\\\\"
    print string
    print '\hline'
    print '\end{tabular}'