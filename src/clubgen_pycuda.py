# -*- coding: utf-8 -*-
import pycuda.autoinit
import pycuda.driver as cuda
import pycuda.gpuarray as gpuarray
import numpy
import logging
from PARAMS import *
import time
import math

from pycuda.compiler import SourceModule

def GenClubs(seqOfInt, testSeqs, parameters):
    
    logging.info('Starting clubgen on CUDA (from pyCUDA)...')
    logging.info('Device found: ' + str(pycuda.driver.Device(0).name()))
    logging.info('of ' + str(pycuda.driver.Device(0).count()) + ' devices')
    logging.info('Device memory (free,total): ' + str(pycuda.driver.mem_get_info()))
    
    logging.info('Setting up GPU...')
    t0 = time.time()
    
    MAX_MATCH = 20;
    MAX_INLINE = 5;
    
    soiPartTotal = seqOfInt.length/1000+1
    if seqOfInt.readingFrames > 1:
        soiseq = seqOfInt.seq[0]
    else:
        soiseq = seqOfInt.seq
    soi = numpy.zeros(soiPartTotal*1000,dtype = numpy.uint8)
    if seqOfInt.sequenceType == 'AA':
        soi[0:seqOfInt.length] += numpy.asarray(soiseq,dtype = numpy.uint8)
    else:
        soi[0:seqOfInt.length] += numpy.asarray(soiseq,dtype = numpy.uint8)%2**seqOfInt.dataDepths[0]
    
    tsPartTotal = 0
    tsPartInfo = []
    tsPartNum = 0
    if seqOfInt.readingFrames == 1:
        for j in range(len(testSeqs)):
            tmpPart = testSeqs[j].length/1000+1
            tsPartTotal += tmpPart
            for i in range(tmpPart):
                tsPartInfo.append(tsPartNum)
            tsPartNum += 1
    else:
        for k in range(seqOfInt.readingFrames):
            for j in range(len(testSeqs)):
                tmpPart = len(testSeqs[j].seq[k])/1000+1
                tsPartTotal += tmpPart
                for i in range(tmpPart):
                    tsPartInfo.append(tsPartNum)
                tsPartNum += 1
    tsPartInfo = numpy.asarray(tsPartInfo,dtype = numpy.uint16)
    
    thePart = -1
    counter = 0
    tsPartNum = 0
    
    if seqOfInt.readingFrames == 1:
        if seqOfInt.sequenceType == 'AA':
            ts = numpy.zeros(tsPartTotal*1000,dtype = numpy.uint8)
            for i in range(tsPartTotal):
                if thePart == tsPartInfo[i]:
                    counter += 1
                else:
                    thePart = tsPartInfo[i]
                    counter = 0
                for j in range(min(1000,testSeqs[tsPartInfo[i]].length-1000*counter)):
                    ts[j+1000*i] = testSeqs[tsPartInfo[i]].seq[j+1000*counter]
                for j in range(min(1000,testSeqs[tsPartInfo[i]].length-1000*counter),1000):
                    ts[j+1000*i] = 0
        else:
            ts = numpy.zeros(tsPartTotal*1000,dtype = numpy.uint8)
            for i in range(tsPartTotal):
                if thePart == tsPartInfo[i]:
                    counter += 1
                else:
                    thePart = tsPartInfo[i]
                    counter = 0
                for j in range(min(1000,testSeqs[tsPartInfo[i]].length-1000*counter)):
                    ts[j+1000*i] = testSeqs[tsPartInfo[i]].seq[j+1000*counter]%2**seqOfInt.dataDepths[0]
                for j in range(min(1000,testSeqs[tsPartInfo[i]].length-1000*counter),1000):
                    ts[j+1000*i] = 0
    else:
        ts = numpy.zeros(tsPartTotal*1000,dtype = numpy.uint8)
        for k in range(seqOfInt.readingFrames):
            for i in range(tsPartTotal/seqOfInt.readingFrames):
                infoLen = len(testSeqs[tsPartInfo[i]%len(testSeqs)].seq[k])
                if thePart == tsPartInfo[i]:
                    counter += 1
                else:
                    thePart = tsPartInfo[i]
                    counter = 0
                    
                ts[1000*i+k*tsPartTotal/seqOfInt.readingFrames*1000:1000*i+k*tsPartTotal/seqOfInt.readingFrames*1000+min(1000,infoLen-1000*counter)] = testSeqs[tsPartInfo[i]%len(testSeqs)].seq[k][1000*counter:1000*counter+min(1000,infoLen-1000*counter)]%2**seqOfInt.dataDepths[0]
                
                ts[1000*i+k*tsPartTotal/seqOfInt.readingFrames*1000+min(1000,infoLen-1000*counter):1000*i+k*tsPartTotal/seqOfInt.readingFrames*1000+1000] = numpy.zeros(1000-min(1000,infoLen-1000*counter),dtype = numpy.uint8)
                #for j in range(min(1000,infoLen-1000*counter)):
                #    ts[j+1000*i+k*tsPartTotal/seqOfInt.readingFrames*1000] = testSeqs[tsPartInfo[i]%len(testSeqs)].seq[k][j+1000*counter]%2**seqOfInt.dataDepths[0]
                #for j in range(min(1000,infoLen-1000*counter),1000):
                #    ts[j+1000*i+k*tsPartTotal/seqOfInt.readingFrames*1000] = 0
                    
    ts = numpy.asarray(ts,dtype = numpy.uint8)
    
    if seqOfInt.readingFrames == 1:
        tsLens = numpy.zeros(len(testSeqs),dtype= numpy.uint16)
        for i in range(len(testSeqs)):
            tsLens[i] = testSeqs[i].length
    else:
        tsLens = numpy.zeros(len(testSeqs)*seqOfInt.readingFrames,dtype= numpy.uint16)
        for k in range(seqOfInt.readingFrames):
            for i in range(len(testSeqs)):
                tsLens[i+k*len(testSeqs)] = len(testSeqs[i].seq[k])
        
    counter = numpy.zeros(tsPartTotal, dtype = numpy.uint16)
    count = 0
    for i in range(1,tsPartTotal):
        if tsPartInfo[i] == tsPartInfo[i-1]:
            count += 1
            counter[i] += count
        else:
            count = 0
            counter[i] += count
    
    #print list(counter)
    #print list(tsLens)
    #print list(tsPartInfo)
    #print len(ts)
    #print len(soi)
    #print seqOfInt.readingFrames
    #print tsPartTotal
    mod = SourceModule("""
__global__ void GenClubCuda(const unsigned char isEnd, const unsigned char *d_soi, const unsigned char *d_ts, const unsigned short int *tsLens, const unsigned short int *tsPartInfo, const unsigned short int *counter, int *alignmentInfo);
    
__device__ bool xEQy8(const unsigned char x, const unsigned char y){
    //return !(x-y);
    return x==y;
    //if(x == y){return 1;}else{return 0;}
    //return !((x&16)^(y&16))&!((x&8)^(y&8))&!((x&4)^(y&4))&!((x&2)^(y&2))&!((x&1)^(y&1));
}
    
__device__ int LocalFilter(const bool *matchVec, const unsigned short int matLen, const float locSig, const unsigned short int locWin, const float sinExp, const float gapBridge, short int *tempMatch){
    const float MIN_DENSITY = 0.35f;

    // have to retool this a little to fit onto cuda (have like 10KB left on the core and 1000 floats overwhelms that, even 1000 int64s overwhelm it)
    unsigned short int i,j;
    unsigned char matches = 0;
    float significant = locSig*sqrt(sinExp*locWin);
    short int start = -1, stop = 0;
    short int densityMatches = 0;
    //int extension = 0;
    
    for(i = 0; i < 10; i++){
        tempMatch[i] = -1;
    }
    
    // make a cumulative sum of the matchvec, subtracting off the expected match
    //const unsigned char jumper = (int)((1.0)/sinExp);
    unsigned short int cumsum[1000];
    cumsum[0] = matchVec[0];
    for(i = 1; i < matLen; i++){
        cumsum[i] = cumsum[i-1] + matchVec[i];
    }
    //tempMatch[0]=cumsum[0];tempMatch[1]=cumsum[matLen-1];return 1;
    // from the left find a significant jump in the cumsum
    for(i = locWin; i < matLen - locWin; i++){
        if(cumsum[i] - cumsum[i-locWin] - sinExp*locWin > significant && cumsum[i+locWin] - cumsum[i] - sinExp*locWin > significant){
            if(start == -1){
                start = i - locWin;
            }
        }else{
            if(start != -1){
                
                stop = i + locWin;
                
                for(j = 0; j < locWin-1; j++){ // need 2/4 to start/end
                     if(matchVec[start] + matchVec[start+1] + matchVec[start+2] + matchVec[start+3] < 2){start++;}
                     if(matchVec[stop] + matchVec[stop-1] + matchVec[stop-2] + matchVec[stop-3] < 2){stop--;}
                }
                start += !matchVec[start];
                start += !matchVec[start];
                stop -= !matchVec[stop];
                stop -= !matchVec[stop];
                stop++;
                densityMatches = 0;
                for(j = start;j < stop; j++){
                    densityMatches += matchVec[j];
                }
                
                if((densityMatches)/(float)(stop-start) > MIN_DENSITY){ // MINIMUM DENSITY IMPOSED
                    
                    // MINIMUM CUTTOF IMPOSED (to find 11-mers etc)
                    if(stop > start && (stop - start > locWin)){
                        
                        // new local align found
                        if(matches == 0){
                            tempMatch[matches] = start;
                            tempMatch[matches+1] = stop;
                            matches += 2;
                        }else{
                            if(start - tempMatch[matches - 1] < locWin*gapBridge){
                                tempMatch[matches-1] = stop;
                            }else{
                                tempMatch[matches] = start;
                                tempMatch[matches+1] = stop;
                                matches += 2;
                            }
                        }
                        // dont need to search in the area just defined in
                        i = stop;
                    }
                }
                start = -1;
            }
        }
    }
    if(start != -1){
                stop = i + locWin;
                for(j = 0; j < locWin-1; j++){ // need 2/4 to start/end
                     if(matchVec[start] + matchVec[start+1] + matchVec[start+2] + matchVec[start+3] < 2){start++;}
                     if(matchVec[stop] + matchVec[stop-1] + matchVec[stop-2] + matchVec[stop-3] < 2){stop--;}
                }
                start += !matchVec[start];
                start += !matchVec[start];
                stop -= !matchVec[stop];
                stop -= !matchVec[stop];
                densityMatches = 0;
                for(j = start;j < stop; j++){
                    densityMatches += matchVec[j]&1;
                }
                
                if((densityMatches)/(float)(stop-start) > MIN_DENSITY){ // MINIMUM DENSITY IMPOSED
                    // MINIMUM CUTTOF IMPOSED (to find 11-mers etc)
                    
                    if(stop > start && (stop - start > locWin)){
                        
                        if(matches == 0){
                            tempMatch[matches] = start;
                            tempMatch[matches+1] = stop;
                            matches += 2;
                        }else{
                            if(start - tempMatch[matches - 1] < locWin*gapBridge){
                                tempMatch[matches-1] = stop;
                                //extension++;
                            }else{
                                tempMatch[matches] = start;
                                tempMatch[matches+1] = stop;
                                matches += 2;
                            }
                        }
                    }
                }
                start = 0;
            }
    if(matches > 0){
        return 1;
    }else{
        return 0;
    }
}

__global__ void GenClubCuda(const unsigned char isEnd, const unsigned char *d_soi, const unsigned char *d_ts, const unsigned short int *tsLens, const unsigned short int *tsPartInfo, const unsigned short int *counter, int *alignmentInfo){
    const unsigned short int maxMatch = """ + str(MAX_MATCH) + """;
    const unsigned short int MAX_INLINE = """ + str(MAX_INLINE) + """;
    const unsigned int soiLenTotal = """ + str(numpy.uint16(seqOfInt.length)) + """;
    const unsigned int soiLen = """ + str(soiPartTotal) + """;
    const unsigned int tsLen = """ + str(tsPartTotal) + """;
    const float globSig = """ + str(numpy.float32(parameters['GLOBAL_SIG_LEVEL'])) + """;
    const float locSig = """ + str(numpy.float32(parameters['LOCAL_SIG_LEVEL'])) + """;
    const unsigned short int locWin = """ + str(numpy.uint16(parameters['LOCAL_WINDOW'])) + """;
    const float gapBridge = """ + str(numpy.float32(parameters['LOCAL_BRIDGE_WIDTH'])) + """;
    const unsigned char comber = """ + str(numpy.uint8(parameters['COMB_SPEED'])) + """;
    
    unsigned short int i = 0, j = 0, k = 0;

    // find the locations on the alignment grid
    const unsigned int idx = threadIdx.x + blockIdx.x*blockDim.x;
    unsigned int tsID = idx%(tsLen);
    unsigned int soiID = idx/tsLen;
    if(isEnd == 1){
        soiID = soiLen - 1 - (blockDim.x - threadIdx.x)/tsLen;
        tsID = tsLen - (blockDim.x - threadIdx.x)%(tsLen);
    }
    
    // find the relative lengths of the soi and ts we are comparing
    const unsigned short int lensoi = min(1000,soiLenTotal - 1000*soiID);
    const unsigned short int lents = min(1000,tsLens[tsPartInfo[tsID]] - counter[tsID]*1000);
    
    // copy this kernels soi and ts chunks into local/register memory
    unsigned char soi[1000];
    unsigned char ts[1000];
    for(i = 0; i < lensoi; i++){
        soi[i] = d_soi[1000*soiID + i];
    }
    for(i = 0; i < lents; i++){
        ts[i] = d_ts[1000*tsID + i];
    }
    
    
    // find where on the alignmentInfo vector we will be writing
    unsigned int outputLoc = soiID*3*maxMatch + tsID*3*maxMatch*soiLen;
    // this is incremented 3 after each local alignment is found
    
    // set up local variables
    
    int c_start[maxMatch];
    int c_end[maxMatch];
    int c_shift[maxMatch];
    
    bool matchVec[1000];
    
    int soiVecComp[21];
    int tsVecComp[21];
    
    int numMatches = 0;
    float expMatches = 0.0;
    short int actMatches = -1;
    
    short int matchNum = 0;
    short int tempMatch[2*MAX_INLINE];
    
    // initialize local variables
    for(i = 0; i < maxMatch; i++){
        c_start[i] = -1;
        c_end[i] = -1;
        c_shift[i] = -1;
    }
    
    for(i = 0; i < 2*MAX_INLINE; i++){
        tempMatch[i] = -1;
    }
    
    // Builds up to compostion of the first locWin residues
    for(k = 0; k < 21; k++){
        soiVecComp[k] = 0;
        tsVecComp[k] = 0;
    }

    for(k = 0; k < locWin; k++){
        soiVecComp[soi[k]]++;
        tsVecComp[ts[lents-1-k]]++;
    }

    // finds the maximum number of matches over all offsets
    for(k = 0; k < 21; k++){
        numMatches += soiVecComp[k]*tsVecComp[k];
    }
    
    // loops unrolled to get rid of the 4 min,max operations per offset
    if(lents > lensoi){
        // start of the comparison
        for(i = 0; i < lensoi - locWin; i++){
            actMatches = 0;
            // this is the number of matches needed to trigger the global filter
            expMatches = numMatches/(float)(locWin+i) + globSig*sqrt(numMatches/(float)(locWin+i));
            expMatches = expMatches/(float)comber;
            // check soi[0:i+locWin] == ts[lents-locWin-i:lents]
            for(j = 0; j < i + locWin; j+=comber){
                // if there is a match add it
                actMatches += xEQy8(soi[j],ts[lents-locWin-i+j]);
                //if(soi[j] == ts[lents-locWin-i+j]){
                //    actMatches++;
                //}
            }
            // Global test
            if(actMatches > expMatches){
                if(actMatches > locWin/2){ // quick check to make sure there are a reasonable number of actual matches with respect to the local window
                    // build matchvec
                    for(j = 0; j < i + locWin; j++){
                        matchVec[j] = xEQy8(soi[j],ts[lents-locWin-i+j]);
                        // if there is a match add it
                        //if(soi[j] == ts[lents-locWin-i+j]){
                        //    matchVec[j] = 1;
                        //}else{
                        //    matchVec[j] = 0;
                        //}
                    }
                    // local test
                    if(LocalFilter(matchVec,locWin+i,locSig,locWin,numMatches/(float)((locWin+i)*(locWin+i)),gapBridge,tempMatch)){
                        for(j = 0; j < MAX_INLINE; j++){
                            if(tempMatch[j*2+1] != -1){
                                c_start[matchNum] = tempMatch[j*2]; // start place
                                c_end[matchNum] = tempMatch[j*2+1]; // end place
                                c_shift[matchNum] = i; // offset
                                // these values are decomposed in python to fit into the alignment objects
                                matchNum++;
                            }else{
                                j = MAX_INLINE;
                            }
                        }
                    }
                }
            }
            // Change vector compositions and inrement total number of possible matches
            numMatches += soiVecComp[ts[lents-locWin-i-1]];
            numMatches += tsVecComp[soi[i+locWin]];
            soiVecComp[soi[i+locWin]]++;
            tsVecComp[ts[lents-locWin-i-1]]++;
            numMatches += xEQy8(soi[i+locWin],ts[lents-locWin-i-1]);
            //if(soi[i+locWin] == ts[lents-locWin-i-1]){numMatches++;}

        }
        // middle part of the comparison
        for(i = lensoi - locWin; i < lents - locWin; i++){
            actMatches = 0;
            expMatches = numMatches/(float)lensoi + globSig*sqrt(numMatches/(float)lensoi);
            expMatches = expMatches/(float)comber;
            // check soi[0:lensoi] == ts[lents - locWin - i: lents -locWin -i +lensoi]
            for(j = 0; j < lensoi; j+=comber){
                actMatches += xEQy8(soi[j],ts[lents - locWin -i+j]);
                //if(soi[j] == ts[lents - locWin -i+j]){
                //    actMatches++;
                //}
            }
            if(actMatches > expMatches){
                    // build matchvec
                    for(j = 0; j < lensoi; j++){
                        // if there is a match add it
                        matchVec[j] = xEQy8(soi[j],ts[lents - locWin -i+j]);
                        //if(soi[j] == ts[lents - locWin -i+j]){
                        //    matchVec[j] = 1;
                        //}else{
                        //    matchVec[j] = 0;
                        //}
                    }
                    if(LocalFilter(matchVec,lensoi,locSig,locWin,numMatches/(float)(lensoi*lensoi),gapBridge,tempMatch)){
                        for(j = 0; j < MAX_INLINE; j++){
                            if(tempMatch[j*2+1] != -1){
                                c_start[matchNum] = tempMatch[j*2]; // start place
                                c_end[matchNum] = tempMatch[j*2+1]; // end place
                                c_shift[matchNum] = i; // offset
                                // these values are decomposed in python to fit into the alignment objects
                                matchNum++;
                            }else{
                                j = MAX_INLINE;
                            }
                        }
                    }
            }
            // Change vector compositions and inrement total number of possible matches
            numMatches -= soiVecComp[ts[lents - locWin -i + lensoi-1]];
            numMatches += soiVecComp[ts[lents - locWin -i-1]];
            tsVecComp[ts[lents - locWin -i + lensoi-1]]--; //shift out
            tsVecComp[ts[lents - locWin -i-1]]++; //shift in
        }
        // end of the comparison
        for(i = lents - locWin; i < lensoi + lents - locWin - locWin; i++){
            actMatches = 0;
            expMatches = numMatches/(float)(lensoi-i+lents-locWin) + globSig*sqrt(numMatches/(float)(lensoi-i+lents-locWin));
            expMatches = expMatches/(float)comber;
            // check soi[i - lents + locWin: lensoi] == ts[0:lensoi -i+lents-locWin]
            for(j = 0; j < lensoi -i+lents-locWin; j+=comber){
                actMatches += xEQy8(soi[j+i-lents+locWin],ts[j]);
                //if(soi[j+i-lents+locWin] == ts[j]){
                //    actMatches++;
                //}
            }
            // Global test
            if(actMatches > expMatches){
                if(actMatches > locWin/2){
                    // build matchvec
                    for(j = 0; j < lensoi -i+lents-locWin; j++){
                        // if there is a match add it
                        matchVec[j] = xEQy8(soi[j+i-lents+locWin],ts[j]);
                        //if(soi[j+i-lents+locWin] == ts[j]){
                        //    matchVec[j] = 1;
                        //}else{
                        //    matchVec[j] = 0;
                        //}
                    }
                    if(LocalFilter(matchVec,lensoi-i+lents-locWin,locSig,locWin,numMatches/(float)((lensoi-i+lents-locWin)*(lensoi-i+lents-locWin)),gapBridge,tempMatch)){
                        for(j = 0; j < MAX_INLINE; j++){
                            if(tempMatch[j*2+1] != -1){
                                c_start[matchNum] = tempMatch[j*2]; // start place
                                c_end[matchNum] = tempMatch[j*2+1]; // end place
                                c_shift[matchNum] = i; // offset
                                // these values are decomposed in python to fit into the alignment objects
                                matchNum++;
                            }else{
                                j = MAX_INLINE;
                            }
                        }
                    }
                }
            }
            // Change vector compositions and inrement total number of possible matches
            numMatches -= soiVecComp[ts[lents - locWin -i+lensoi-1]];
            numMatches -= tsVecComp[soi[i-lents+locWin]];
            soiVecComp[soi[i-lents+locWin]]--;
            tsVecComp[ts[lents - locWin -i+lensoi-1]]--;
            numMatches += xEQy8(soi[i-lents+locWin],ts[lents - locWin -i+lensoi-1]);
            //if(soi[i-lents+locWin] == ts[lents - locWin -i+lensoi-1]){numMatches++;}
        }
    }else{
        // the other option: lents <= lensoi
        for(i = 0; i < lents - locWin; i++){
            actMatches = 0;
            expMatches = numMatches/(float)(locWin+i) + globSig*sqrt(numMatches/(float)(locWin+i));
            expMatches = expMatches/(float)comber;
            // check soi[0:i+locWin-1] == ts[lents-locWin-i:lents-1]
            for(j = 0; j < i + locWin; j+=comber){
                // if there is a match add it
                actMatches += xEQy8(soi[j],ts[lents-locWin-i+j]);
                //if(soi[j] == ts[lents-locWin-i+j]){
                //    actMatches++;
                //}
            }
            // Global test
            if(actMatches > expMatches){
                if(actMatches > locWin/2){
                    // build matchvec
                    for(j = 0; j < i + locWin; j++){
                        matchVec[j] = xEQy8(soi[j],ts[lents-locWin-i+j]);
                        // if there is a match add it
                        //if(soi[j] == ts[lents-locWin-i+j]){
                        //    matchVec[j] = 1;
                        //}else{
                        //    matchVec[j] = 0;
                        //}
                    }
                    if(LocalFilter(matchVec,locWin+i,locSig,locWin,numMatches/(float)((locWin+i)*(locWin+i)),gapBridge,tempMatch)){
                        for(j = 0; j < MAX_INLINE; j++){
                            if(tempMatch[j*2+1] != -1){
                                c_start[matchNum] = tempMatch[j*2]; // start place
                                c_end[matchNum] = tempMatch[j*2+1]; // end place
                                c_shift[matchNum] = i; // offset
                                // these values are decomposed in python to fit into the alignment objects
                                matchNum++;
                            }else{
                                j = MAX_INLINE;
                            }
                        }
                    }
                }
            }
            // Change vector compositions and inrement total number of possible matches
            numMatches += soiVecComp[ts[lents-locWin-i-1]];
            numMatches += tsVecComp[soi[i+locWin]];
            soiVecComp[soi[i+locWin]]++;
            tsVecComp[ts[lents-locWin-i-1]]++;
            numMatches += xEQy8(soi[i+locWin],ts[lents-locWin-i-1]);
            //if(soi[i+locWin] == ts[lents-locWin-i-1]){numMatches++;}
        }
        // middle
        for(i = lents - locWin; i < lensoi - locWin; i++){
            actMatches = 0;
            expMatches = numMatches/(float)lents + globSig*sqrt(numMatches/(float)lents);
            expMatches = expMatches/(float)comber;
            // check soi[0:lensoi-1] == ts[lents - locWin - i: lents -1 -locWin -i +lensoi]
            for(j = 0; j < lents; j+=comber){
                actMatches += xEQy8(soi[i-lents+locWin+j],ts[j]);
                //if(soi[i-lents+locWin+j] == ts[j]){
                //    actMatches++;
                //}
            }
            // Global test
            if(actMatches > expMatches){
                    // build matchvec
                    for(j = 0; j < lents; j++){
                        matchVec[j] = xEQy8(soi[i-lents+locWin+j],ts[j]);
                        // if there is a match add it
                        //if(soi[i-lents+locWin+j] == ts[j]){
                        //    matchVec[j] = 1;
                        //}else{
                        //    matchVec[j] = 0;
                        //}
                    }
                    if(LocalFilter(matchVec,lents,locSig,locWin,numMatches/(float)(lents*lents),gapBridge,tempMatch)){
                        for(j = 0; j < MAX_INLINE; j++){
                            if(tempMatch[j*2+1] != -1){
                                c_start[matchNum] = tempMatch[j*2]; // start place
                                c_end[matchNum] = tempMatch[j*2+1]; // end place
                                c_shift[matchNum] = i; // offset
                                // these values are decomposed in python to fit into the alignment objects
                                matchNum++;
                            }else{
                                j = MAX_INLINE;
                            }
                        }
                    }
            }
            // Change vector compositions and inrement total number of possible matches
            numMatches += tsVecComp[soi[i-lents+locWin]];
            numMatches -= tsVecComp[soi[i+locWin-1]];
            soiVecComp[soi[i-lents+locWin]]++; //shift in
            soiVecComp[soi[i+locWin-1]]--; //shift out
        }
        //end
        for(i = lensoi - locWin; i < lensoi + lents - locWin - locWin; i++){
            actMatches = 0;
            expMatches = numMatches/(float)(lensoi-i+lents-locWin) + globSig*sqrt(numMatches/(float)(lensoi-i+lents-locWin));
            //expMatches = expMatches/(float)comber;
            // check soi[i - lents + locWin: lensoi -1] == ts[lents-locWin-i:lents-1]
            for(j = 0; j < lensoi -i+lents-locWin; j+=comber){
                actMatches += xEQy8(soi[j+i-lents+locWin],ts[j]);
                //if(soi[j+i-lents+locWin] == ts[j]){
                //    actMatches++;
                //}
            }
            // Global test
            if(actMatches > expMatches){
                if(actMatches > locWin/2){
                    // build matchvec
                    for(j = 0; j < lensoi -i+lents-locWin; j++){
                        matchVec[j] = xEQy8(soi[j+i-lents+locWin],ts[j]);
                        // if there is a match add it
                        //if(soi[j+i-lents+locWin] == ts[j]){
                        //    matchVec[j] = 1;
                        //}else{
                        //    matchVec[j] = 0;
                        //}
                    }
                    if(LocalFilter(matchVec,lensoi -i+lents-locWin,locSig,locWin,numMatches/(float)((lensoi -i+lents-locWin)*(lensoi -i+lents-locWin)),gapBridge,tempMatch)){
                        for(j = 0; j < MAX_INLINE; j++){
                            if(tempMatch[j*2+1] != -1){
                                c_start[matchNum] = tempMatch[j*2]; // start place
                                c_end[matchNum] = tempMatch[j*2+1]; // end place
                                c_shift[matchNum] = i; // offset
                                // these values are decomposed in python to fit into the alignment objects
                                matchNum++;
                            }else{
                                j = MAX_INLINE;
                            }
                        }
                    }
                }
            }
            // Change vector compositions and inrement total number of possible matches
            numMatches -= soiVecComp[ts[lents-(i+ locWin-lensoi)-1]];
            numMatches -= tsVecComp[soi[i-lents+locWin]];
            soiVecComp[soi[i-lents+locWin]]--;
            tsVecComp[ts[lents-(i+ locWin-lensoi)-1]]--;
            numMatches += xEQy8(soi[i-lents+locWin],ts[lents-(i+ locWin-lensoi)-1]);
            //if(soi[i-lents+locWin] == ts[lents-(i+ locWin-lensoi)-1]){numMatches++;}
        }
    }
    for(i = 0; i < maxMatch; i++){
        if(c_start[i] != -1){
            //atomicAdd(int* address, int val);
            atomicAdd(&alignmentInfo[outputLoc], c_start[i]);
            atomicAdd(&alignmentInfo[outputLoc+1], c_end[i]);
            atomicAdd(&alignmentInfo[outputLoc+2], c_shift[i]);
            //alignmentInfo[outputLoc] = c_start[i];
            //alignmentInfo[outputLoc+1] = c_end[i];
            //alignmentInfo[outputLoc+2] = c_shift[i];
            outputLoc += 3;
        }else{
            i = maxMatch;
        }
    }
}
""")
    
    # see if problem will fit on device, if not then break it up
    
    alignmentInfo = numpy.zeros(MAX_MATCH*3*tsPartTotal*soiPartTotal, dtype = numpy.int32)
    
    size = 0
    size += soi.size * soi.dtype.itemsize # soi
    size += ts.size * ts.dtype.itemsize # ts
    size += tsLens.size * tsLens.dtype.itemsize # tsLens
    size += tsPartInfo.size * tsPartInfo.dtype.itemsize # tsPartInfo
    size += counter.size * counter.dtype.itemsize # counter
    alignSize = alignmentInfo.size * alignmentInfo.dtype.itemsize
    
    deviceFree = pycuda.driver.mem_get_info()[0]*.5 # only let it fill 1/2 full
    
    if deviceFree > size + alignSize: # no problem, put it all on the card, otherwise chop it up
    
        # set up optimal thread/block combos
        
        numThreads = soiPartTotal*tsPartTotal
        MP = 24
        MAX_THREADS = 128
        
        MPneeded = math.ceil(numThreads/float(MAX_THREADS))
        
        if MPneeded <= MP:
            # split the threads evenly among MP blocks
            threadsOnBlock = math.floor(numThreads/float(MP))
            threadsOnLast = numThreads - math.floor(numThreads/float(MP))*MP
            grider = (int(MP),1)
        else:
            # split the threads evenly among MPneeded blocks
            threadsOnBlock = math.floor(numThreads/float(MPneeded))
            threadsOnLast = numThreads - math.floor(numThreads/float(MPneeded))*MPneeded
            grider = (int(MPneeded),1)
        blocker = (int(threadsOnBlock),1,1)
        
        logging.info('Main executing block: ' + str(int(threadsOnBlock)) + ' threads on ' + str(grider[0]) + ' blocks.')
    
        # SET UP ARRAYS AND COPY THEM TO GPU
        
        d_soi = cuda.mem_alloc(soi.size * soi.dtype.itemsize)
        cuda.memcpy_htod(d_soi,soi)
        d_soi_dtype = numpy.ctypeslib.ndpointer(dtype=soi.dtype,ndim=soi.size)
        
        d_ts = cuda.mem_alloc(ts.size * ts.dtype.itemsize)
        cuda.memcpy_htod(d_ts,ts)
        d_ts_dtype = numpy.ctypeslib.ndpointer(dtype=ts.dtype,ndim=ts.size)
        
        d_tsLens = cuda.mem_alloc(tsLens.size * tsLens.dtype.itemsize)
        cuda.memcpy_htod(d_tsLens,tsLens)
        d_tsLens_dtype = numpy.ctypeslib.ndpointer(dtype=tsLens.dtype,ndim=tsLens.size)
        
        d_tsPartInfo = cuda.mem_alloc(tsPartInfo.size * tsPartInfo.dtype.itemsize)
        cuda.memcpy_htod(d_tsPartInfo,tsPartInfo)
        d_tsPartInfo_dtype = numpy.ctypeslib.ndpointer(dtype=tsPartInfo.dtype,ndim=tsPartInfo.size)
        
        d_counter = cuda.mem_alloc(counter.size * counter.dtype.itemsize)
        cuda.memcpy_htod(d_counter,counter)
        d_counter_dtype = numpy.ctypeslib.ndpointer(dtype=counter.dtype,ndim=counter.size)
        
        d_alignmentInfo = cuda.mem_alloc(alignmentInfo.size * alignmentInfo.dtype.itemsize)
        cuda.memcpy_htod(d_alignmentInfo,alignmentInfo)
        d_alignmentInfo_dtype = numpy.ctypeslib.ndpointer(dtype=alignmentInfo.dtype,ndim=alignmentInfo.size)
        
        # SET UP FUNCTION MAIN RUN
        
        GenClubCuda = mod.get_function("GenClubCuda")
        
        GenClubCuda.prepare((numpy.uint8,d_soi_dtype,d_ts_dtype,d_tsLens_dtype,d_tsPartInfo_dtype,d_counter_dtype,d_alignmentInfo_dtype),block = blocker)
        
        # TIMING
        
        runner = cuda.Stream()
        runner.synchronize()
        start = cuda.Event()
        stop = cuda.Event()
        start = start.record()
        start.synchronize()
        
        #print GenClubCuda.local_size_bytes
        #print GenClubCuda.shared_size_bytes
        #print GenClubCuda.num_regs
        
        # RUN KERNEL ON GPU
        
        #GenClubCuda.prepared_call((tsPartTotal,soiPartTotal),d_soi, d_ts, numpy.uint16(seqOfInt.length), numpy.uint32(soiPartTotal), numpy.uint32(tsPartTotal), numpy.float32(parameters['GLOBAL_SIG_LEVEL']), numpy.float32(parameters['LOCAL_SIG_LEVEL']), numpy.uint16(parameters['LOCAL_WINDOW']), numpy.uint8(20), numpy.uint8(5), numpy.float32(parameters['LOCAL_BRIDGE_WIDTH']), numpy.uint8(parameters['COMB_SPEED']), d_tsLens, d_tsPartInfo, d_counter, d_alignmentInfo)
        
        GenClubCuda.prepared_async_call(grider,runner,0,d_soi, d_ts, d_tsLens, d_tsPartInfo, d_counter, d_alignmentInfo)
        
        #GenClubCuda(d_soi, d_ts, numpy.uint16(seqOfInt.length), numpy.uint32(soiPartTotal), numpy.uint32(tsPartTotal), numpy.float32(parameters['GLOBAL_SIG_LEVEL']), numpy.float32(parameters['LOCAL_SIG_LEVEL']), numpy.uint16(parameters['LOCAL_WINDOW']), numpy.uint8(20), numpy.uint8(5), numpy.float32(parameters['LOCAL_BRIDGE_WIDTH']), numpy.uint8(parameters['COMB_SPEED']), d_tsLens, d_tsPartInfo, d_counter, d_alignmentInfo, block = (soiPartTotal,1,1), grid = (tsPartTotal,1))
    
        # TIMING P2
    
        runner.synchronize()
        stop = stop.record()
        stop.synchronize()
    
        logging.info('Done executing main on GPU: ' + str(stop.time_since(start)) + '(ms)')
    
        logging.info('Running leftover ' + str(int(threadsOnLast)) + ' threads.')
    
        # SET UP FUNCTION LEFTOVER RUN, this should really be done in serial using the C code..
    
        GenClubCuda.prepare((numpy.uint8,d_soi_dtype,d_ts_dtype,d_tsLens_dtype,d_tsPartInfo_dtype,d_counter_dtype,d_alignmentInfo_dtype),block = (int(threadsOnLast),1,1))
    
        # TIMING LO
    
        runner = cuda.Stream()
        runner.synchronize()
        start = cuda.Event()
        stop = cuda.Event()
        start = start.record()
        start.synchronize()
    
        # RUN LO KERNEL ON GPU
    
        #GenClubCuda.prepared_call((tsPartTotal,soiPartTotal),d_soi, d_ts, numpy.uint16(seqOfInt.length), numpy.uint32(soiPartTotal), numpy.uint32(tsPartTotal), numpy.float32(parameters['GLOBAL_SIG_LEVEL']), numpy.float32(parameters['LOCAL_SIG_LEVEL']), numpy.uint16(parameters['LOCAL_WINDOW']), numpy.uint8(20), numpy.uint8(5), numpy.float32(parameters['LOCAL_BRIDGE_WIDTH']), numpy.uint8(parameters['COMB_SPEED']), d_tsLens, d_tsPartInfo, d_counter, d_alignmentInfo)
    
        GenClubCuda.prepared_async_call((1,1),runner,1,d_soi, d_ts, d_tsLens, d_tsPartInfo, d_counter, d_alignmentInfo)
    
        # COPY INFO FROM GPU TO HOST
    
        # TIMING LO P2
    
        runner.synchronize()
        stop = stop.record()
        stop.synchronize()
    
        logging.info('Done executing leftover on GPU: ' + str(stop.time_since(start)) + '(ms)')
    
        logging.info('Device memory (free,total): ' + str(pycuda.driver.mem_get_info()))
    
        start = cuda.Event()
        stop = cuda.Event()
        start.record()
        start.synchronize()
        getter = cuda.Stream()

        cuda.memcpy_dtoh(alignmentInfo, d_alignmentInfo)
    
        stop = stop.record()
        stop.synchronize()
    
        logging.info('Copying memory from GPU time: ' + str(stop.time_since(start)) + '(ms)')
    else: # too big to fit on card, break up the test sequences, right now the soi has to at least be able to fit...
        alignmentInfo = []
        GenClubCuda = mod.get_function("GenClubCuda")
    
        numberOfBigRunsNeeded = 10000*int(math.floor((size + alignSize)/float(deviceFree)))
        
        logging.info('Need ' + str(numberOfBigRunsNeeded) + ' big runs.')
            
        numThreads = soiPartTotal*tsPartTotal/numberOfBigRunsNeeded
        lastBigBlock = soiPartTotal*tsPartTotal - numberOfBigRunsNeeded*(soiPartTotal*tsPartTotal/numberOfBigRunsNeeded)
        
        numberTSonBigBlock = ts.size/numberOfBigRunsNeeded
        numberTSonLastBlock = int(math.floor(lastBigBlock/float(soiPartTotal))) # conservative
        
        MP = 24
        MAX_THREADS = 128
        
        MPneeded = math.ceil(numThreads/float(MAX_THREADS))
        
        
        if MPneeded <= MP:
            # split the threads evenly among MP blocks
            threadsOnBlock = math.floor(numThreads/float(MP))
            threadsOnLast = numThreads - math.floor(numThreads/float(MP))*MP
            grider = (int(MP),1)
        else:
            # split the threads evenly among MPneeded blocks
            threadsOnBlock = math.floor(numThreads/float(MPneeded))
            threadsOnLast = numThreads - math.floor(numThreads/float(MPneeded))*MPneeded
            grider = (int(MPneeded),1)
        blocker = (int(threadsOnBlock),1,1)
        
        d_soi = cuda.mem_alloc(soi.size * soi.dtype.itemsize)
        cuda.memcpy_htod(d_soi,soi)
        d_soi_dtype = numpy.ctypeslib.ndpointer(dtype=soi.dtype,ndim=soi.size)
        
        d_tsLens = cuda.mem_alloc(tsLens.size * tsLens.dtype.itemsize)
        cuda.memcpy_htod(d_tsLens,tsLens)
        d_tsLens_dtype = numpy.ctypeslib.ndpointer(dtype=tsLens.dtype,ndim=tsLens.size)
        
        d_tsPartInfo = cuda.mem_alloc(tsPartInfo.size * tsPartInfo.dtype.itemsize)
        cuda.memcpy_htod(d_tsPartInfo,tsPartInfo)
        d_tsPartInfo_dtype = numpy.ctypeslib.ndpointer(dtype=tsPartInfo.dtype,ndim=tsPartInfo.size)
        
        d_counter = cuda.mem_alloc(counter.size * counter.dtype.itemsize)
        cuda.memcpy_htod(d_counter,counter)
        d_counter_dtype = numpy.ctypeslib.ndpointer(dtype=counter.dtype,ndim=counter.size)
        
        runner = cuda.Stream()
        
        for i in range(numberOfBigRunsNeeded):
            logging.info('Main executing block ' + str(i) + ': ' + str(int(threadsOnBlock)) + ' threads on ' + str(grider[0]) + ' blocks.')
    
            # SET UP ARRAYS AND COPY THEM TO GPU
        
            tsTemp = ts[i*numberTSonBigBlock:(i+1)*numberTSonBigBlock]
            d_ts = cuda.mem_alloc(ts.size * ts.dtype.itemsize)
            cuda.memcpy_htod(d_ts,ts)
            d_ts_dtype = numpy.ctypeslib.ndpointer(dtype=ts.dtype,ndim=ts.size)
        
            tempAlignment = numpy.zeros(MAX_MATCH*3*numberTSonBigBlock*soiPartTotal, dtype = numpy.int32)
            
            d_alignmentInfo = cuda.mem_alloc(tempAlignment.size * tempAlignment.dtype.itemsize)
            cuda.memcpy_htod(d_alignmentInfo,tempAlignment)
            d_alignmentInfo_dtype = numpy.ctypeslib.ndpointer(dtype=tempAlignment.dtype,ndim=tempAlignment.size)
            
            size = 0
            size += soi.size * soi.dtype.itemsize # soi
            size += tsTemp.size * tsTemp.dtype.itemsize # ts
            size += tsLens.size * tsLens.dtype.itemsize # tsLens
            size += tsPartInfo.size * tsPartInfo.dtype.itemsize # tsPartInfo
            size += counter.size * counter.dtype.itemsize # counter
            alignSize = tempAlignment.size * tempAlignment.dtype.itemsize
            #print size,alignSize
            #return 0
            
            # SET UP FUNCTION MAIN RUN
        
            GenClubCuda.prepare((numpy.uint8,d_soi_dtype,d_ts_dtype,d_tsLens_dtype,d_tsPartInfo_dtype,d_counter_dtype,d_alignmentInfo_dtype),block = blocker)
        
            # TIMING
        
            
            runner.synchronize()
            start = cuda.Event()
            stop = cuda.Event()
            start = start.record()
            start.synchronize()
        
            #print GenClubCuda.local_size_bytes
            #print GenClubCuda.shared_size_bytes
            #print GenClubCuda.num_regs
        
            # RUN KERNEL ON GPU
        
            #GenClubCuda.prepared_call((tsPartTotal,soiPartTotal),d_soi, d_ts, numpy.uint16(seqOfInt.length), numpy.uint32(soiPartTotal), numpy.uint32(tsPartTotal), numpy.float32(parameters['GLOBAL_SIG_LEVEL']), numpy.float32(parameters['LOCAL_SIG_LEVEL']), numpy.uint16(parameters['LOCAL_WINDOW']), numpy.uint8(20), numpy.uint8(5), numpy.float32(parameters['LOCAL_BRIDGE_WIDTH']), numpy.uint8(parameters['COMB_SPEED']), d_tsLens, d_tsPartInfo, d_counter, d_alignmentInfo)
        
            GenClubCuda.prepared_async_call(grider,runner,0,d_soi, d_ts, d_tsLens, d_tsPartInfo, d_counter, d_alignmentInfo)
        
            #GenClubCuda(d_soi, d_ts, numpy.uint16(seqOfInt.length), numpy.uint32(soiPartTotal), numpy.uint32(tsPartTotal), numpy.float32(parameters['GLOBAL_SIG_LEVEL']), numpy.float32(parameters['LOCAL_SIG_LEVEL']), numpy.uint16(parameters['LOCAL_WINDOW']), numpy.uint8(20), numpy.uint8(5), numpy.float32(parameters['LOCAL_BRIDGE_WIDTH']), numpy.uint8(parameters['COMB_SPEED']), d_tsLens, d_tsPartInfo, d_counter, d_alignmentInfo, block = (soiPartTotal,1,1), grid = (tsPartTotal,1))
    
            # TIMING P2
    
            runner.synchronize()
            stop = stop.record()
            stop.synchronize()
    
            logging.info('Done executing main ' + str(i) + ' on GPU: ' + str(stop.time_since(start)) + '(ms)')
    
            logging.info('Running leftover ' + str(i) + ', ' + str(int(threadsOnLast)) + ' threads.')
    
            # SET UP FUNCTION LEFTOVER RUN, this should really be done in serial using the C code..
    
            GenClubCuda.prepare((numpy.uint8,d_soi_dtype,d_ts_dtype,d_tsLens_dtype,d_tsPartInfo_dtype,d_counter_dtype,d_alignmentInfo_dtype),block = (int(threadsOnLast),1,1))
    
            # TIMING LO
    
            runner = cuda.Stream()
            runner.synchronize()
            start = cuda.Event()
            stop = cuda.Event()
            start = start.record()
            start.synchronize()
    
            # RUN LO KERNEL ON GPU
    
            #GenClubCuda.prepared_call((tsPartTotal,soiPartTotal),d_soi, d_ts, numpy.uint16(seqOfInt.length), numpy.uint32(soiPartTotal), numpy.uint32(tsPartTotal), numpy.float32(parameters['GLOBAL_SIG_LEVEL']), numpy.float32(parameters['LOCAL_SIG_LEVEL']), numpy.uint16(parameters['LOCAL_WINDOW']), numpy.uint8(20), numpy.uint8(5), numpy.float32(parameters['LOCAL_BRIDGE_WIDTH']), numpy.uint8(parameters['COMB_SPEED']), d_tsLens, d_tsPartInfo, d_counter, d_alignmentInfo)
    
            GenClubCuda.prepared_async_call((1,1),runner,1,d_soi, d_ts, d_tsLens, d_tsPartInfo, d_counter, d_alignmentInfo)
    
            # COPY INFO FROM GPU TO HOST
    
            # TIMING LO P2
    
            runner.synchronize()
            stop = stop.record()
            stop.synchronize()
    
            logging.info('Done executing leftover ' + str(i) + ' on GPU: ' + str(stop.time_since(start)) + '(ms)')
    
            logging.info('Device memory (free,total): ' + str(pycuda.driver.mem_get_info()))
    
            start = cuda.Event()
            stop = cuda.Event()
            start.record()
            start.synchronize()
            getter = cuda.Stream()

            cuda.memcpy_dtoh(tempAlignment, d_alignmentInfo)
            
    
            stop = stop.record()
            stop.synchronize()
            
            d_ts.free()
            d_alignmentInfo.free()
    
            logging.info('Copying memory from GPU time: ' + str(stop.time_since(start)) + '(ms)')
            
            alignmentInfo.extend(tempAlignment)
        # now run the leftover
            
        numThreads = soi.size*numberTSonLastBlock
        
        MP = 24
        MAX_THREADS = 128
        
        MPneeded = math.ceil(numThreads/float(MAX_THREADS))
        
        
        if MPneeded <= MP:
            # split the threads evenly among MP blocks
            threadsOnBlock = math.floor(numThreads/float(MP))
            threadsOnLast = numThreads - math.floor(numThreads/float(MP))*MP
            grider = (int(MP),1)
        else:
            # split the threads evenly among MPneeded blocks
            threadsOnBlock = math.floor(numThreads/float(MPneeded))
            threadsOnLast = numThreads - math.floor(numThreads/float(MPneeded))*MPneeded
            grider = (int(MPneeded),1)
        blocker = (int(threadsOnBlock),1,1)
        
        tsTemp = ts[-numberTSonLastBlock]
        d_ts = cuda.mem_alloc(ts.size * ts.dtype.itemsize)
        cuda.memcpy_htod(d_ts,ts)
        d_ts_dtype = numpy.ctypeslib.ndpointer(dtype=ts.dtype,ndim=ts.size)
        
        tempAlignment = numpy.zeros(MAX_MATCH*3*numberTSonLastBlock*soiPartTotal, dtype = numpy.int32)
        d_alignmentInfo = cuda.mem_alloc(tempAlignment.size * tempAlignment.dtype.itemsize)
        cuda.memcpy_htod(d_alignmentInfo,tempAlignment)
        d_alignmentInfo_dtype = numpy.ctypeslib.ndpointer(dtype=tempAlignment.dtype,ndim=tempAlignment.size)
        
        GenClubCuda.prepare((numpy.uint8,d_soi_dtype,d_ts_dtype,d_tsLens_dtype,d_tsPartInfo_dtype,d_counter_dtype,d_alignmentInfo_dtype),block = blocker)
        
        # TIMING
        
        runner = cuda.Stream()
        runner.synchronize()
        start = cuda.Event()
        stop = cuda.Event()
        start = start.record()
        start.synchronize()
        
        #print GenClubCuda.local_size_bytes
        #print GenClubCuda.shared_size_bytes
        #print GenClubCuda.num_regs
        
        # RUN KERNEL ON GPU
        
        #GenClubCuda.prepared_call((tsPartTotal,soiPartTotal),d_soi, d_ts, numpy.uint16(seqOfInt.length), numpy.uint32(soiPartTotal), numpy.uint32(tsPartTotal), numpy.float32(parameters['GLOBAL_SIG_LEVEL']), numpy.float32(parameters['LOCAL_SIG_LEVEL']), numpy.uint16(parameters['LOCAL_WINDOW']), numpy.uint8(20), numpy.uint8(5), numpy.float32(parameters['LOCAL_BRIDGE_WIDTH']), numpy.uint8(parameters['COMB_SPEED']), d_tsLens, d_tsPartInfo, d_counter, d_alignmentInfo)
        
        GenClubCuda.prepared_async_call(grider,runner,0,d_soi, d_ts, d_tsLens, d_tsPartInfo, d_counter, d_alignmentInfo)
        
        #GenClubCuda(d_soi, d_ts, numpy.uint16(seqOfInt.length), numpy.uint32(soiPartTotal), numpy.uint32(tsPartTotal), numpy.float32(parameters['GLOBAL_SIG_LEVEL']), numpy.float32(parameters['LOCAL_SIG_LEVEL']), numpy.uint16(parameters['LOCAL_WINDOW']), numpy.uint8(20), numpy.uint8(5), numpy.float32(parameters['LOCAL_BRIDGE_WIDTH']), numpy.uint8(parameters['COMB_SPEED']), d_tsLens, d_tsPartInfo, d_counter, d_alignmentInfo, block = (soiPartTotal,1,1), grid = (tsPartTotal,1))
    
        # TIMING P2
    
        runner.synchronize()
        stop = stop.record()
        stop.synchronize()
    
        logging.info('Done executing main leftover on GPU: ' + str(stop.time_since(start)) + '(ms)')
    
        logging.info('Running leftover^2 ' + str(int(threadsOnLast)) + ' threads.')
    
        # SET UP FUNCTION LEFTOVER RUN, this should really be done in serial using the C code..
    
        GenClubCuda.prepare((numpy.uint8,d_soi_dtype,d_ts_dtype,d_tsLens_dtype,d_tsPartInfo_dtype,d_counter_dtype,d_alignmentInfo_dtype),block = (int(threadsOnLast),1,1))
    
        # TIMING LO
    
        runner = cuda.Stream()
        runner.synchronize()
        start = cuda.Event()
        stop = cuda.Event()
        start = start.record()
        start.synchronize()
    
        # RUN LO KERNEL ON GPU
    
        #GenClubCuda.prepared_call((tsPartTotal,soiPartTotal),d_soi, d_ts, numpy.uint16(seqOfInt.length), numpy.uint32(soiPartTotal), numpy.uint32(tsPartTotal), numpy.float32(parameters['GLOBAL_SIG_LEVEL']), numpy.float32(parameters['LOCAL_SIG_LEVEL']), numpy.uint16(parameters['LOCAL_WINDOW']), numpy.uint8(20), numpy.uint8(5), numpy.float32(parameters['LOCAL_BRIDGE_WIDTH']), numpy.uint8(parameters['COMB_SPEED']), d_tsLens, d_tsPartInfo, d_counter, d_alignmentInfo)
    
        GenClubCuda.prepared_async_call((1,1),runner,1,d_soi, d_ts, d_tsLens, d_tsPartInfo, d_counter, d_alignmentInfo)
    
        # COPY INFO FROM GPU TO HOST
    
        # TIMING LO P2
    
        runner.synchronize()
        stop = stop.record()
        stop.synchronize()
    
        logging.info('Done executing leftover^2 on GPU: ' + str(stop.time_since(start)) + '(ms)')
    
        logging.info('Device memory (free,total): ' + str(pycuda.driver.mem_get_info()))
    
        start = cuda.Event()
        stop = cuda.Event()
        start.record()
        start.synchronize()
        getter = cuda.Stream()

        cuda.memcpy_dtoh(tempAlignment, d_alignmentInfo)
    
        stop = stop.record()
        stop.synchronize()
    
        logging.info('Copying memory from GPU time: ' + str(stop.time_since(start)) + '(ms)')
        
        alignmentInfo.extend(tempAlignment)
        alignmentInfo = numpy.asarray(alignmentInfo, dtype = numpy.int32)
    return alignmentInfo, soiPartTotal, tsPartTotal, counter
        
def finishOnC():
    import clubgen_c
    
