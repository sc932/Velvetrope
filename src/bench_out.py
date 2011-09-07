import benchmarks
import sys

if __name__ == '__main__':
    aligns = sys.argv[1]
    seqDat = sys.argv[2]
    blaster = sys.argv[3]
    hmmer = sys.argv[4]
    benchmarks.plotBenchMulti(aligns,seqDat,blaster,hmmer)