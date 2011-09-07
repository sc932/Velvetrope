#!/bin/bash
# input soiLen, seqLen, nSub, // add later nCat
echo python src/benchmarks.py  $1 $2 $3
python src/benchmarks.py $1 $2 $3

echo makeblastdb -in benchmarks/benchmark${1}/testSeqs.faa -dbtype 'prot' -title 'bencher'
makeblastdb -in benchmarks/benchmark${1}/testSeqs.faa -dbtype 'prot' -title 'bencher'

let NUMAL=${3}*10
echo time blastp -num_alignments $NUMAL -query benchmarks/benchmark${1}/soiSeq.faa -db benchmarks/benchmark${1}/testSeqs.faa -out benchmarks/benchmark${1}/BLASTout
time blastp -num_alignments $NUMAL -query benchmarks/benchmark${1}/soiSeq.faa -db benchmarks/benchmark${1}/testSeqs.faa -out benchmarks/benchmark${1}/BLASTout

echo time phmmer -o HMMout --cpu 1 benchmarks/benchmark${1}/soiSeq.faa benchmarks/benchmark${1}/testSeqs.faa
time phmmer -o benchmarks/benchmark${1}/HMMout --cpu 1 benchmarks/benchmark${1}/soiSeq.faa benchmarks/benchmark${1}/testSeqs.faa

echo time python Velvetrope.py -odir benchmarks/benchmark${1} -rd rg -o benchmarks/benchmark${1}/VRout benchmarks/benchmark${1}/benchVR.faa
time python Velvetrope.py -odir benchmarks/benchmark${1} -rd rg -o benchmarks/benchmark${1}/VRout benchmarks/benchmark${1}/benchVR.faa

echo python src/bench_out.py benchmarks/benchmark${1}/VRout.vrg.pkl benchmarks/benchmark${1}/tsInfo.pkl benchmarks/benchmark${1}/BLASTout benchmarks/benchmark${1}/HMMout
python src/bench_out.py benchmarks/benchmark${1}/VRout.vrg.pkl benchmarks/benchmark${1}/tsInfo.pkl benchmarks/benchmark${1}/BLASTout benchmarks/benchmark${1}/HMMout