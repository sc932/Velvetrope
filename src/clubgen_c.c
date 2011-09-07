/* File: clubgen_c.c */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "clubgen_c.h"
// int mm_Aligner_cu(const int soiLen, int soi[][1000], const int tsLen, int ts[][1000], const int tsPartInfo[], const int soiLenTotal, const int tsLens[], const float globSig, const float locSig, const int locWin, const int maxMatch, const int maxInline, const float gapBridge, const int comber, int alignmentInfo[]);

// int CUDA_Align_c(int soiLen, int soi[][1000], int tsLen, int ts[][1000], int tsPartInfo[], int alignmentInfo[], int soiLenTotal, int tsLens[], double globSig, double locSig, int locWin, int maxMatch, int maxInline, double gapBridge, int comber){
//     mm_Aligner_cu(soiLen, soi, tsLen, ts, tsPartInfo, soiLenTotal, tsLens, (float)globSig, (float)locSig, locWin, maxMatch, maxInline, (float)gapBridge, comber, alignmentInfo);
//     return 1;
// }

int GenAlignment_c(int *soi, int *ts, int lensoi, int lents, double globSig, double locSig, int locWin, int maxMatch, int maxInline, double gapBridge, int comber, int *c_start, int *c_end, int *c_shift){
    //printf("in gen align 1\n");
    int MAX_INLINE = maxInline;
    
    int matchNum = 0;
    int tempMatch[2*MAX_INLINE];
    
    
    double expMatches = 0.0;
    
    int i,j,k;
    for(i = 0; i < 2*MAX_INLINE; i++){tempMatch[i] = -1;}
    for(i = 0; i < maxMatch; i++){
        c_start[i] = 0;
        c_end[i] = 0;
        c_shift[i] = 0;
    }
    
//     for(i = 0; i < 40; i++){
//         printf("%d ",ts[i]);
//     }
//     printf("\n");
//     for(i = 0; i < 40; i++){
//         printf("%d ",soi[i]);
//     }
//     printf("\n");
    
    int soiVecComp[21];
    int tsVecComp[21];
    
    int numMatches = 0;

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

    int actMatches = -1;
    //printf("in gen align 2\n");
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
                if(soi[j] == ts[lents-locWin-i+j]){
                    actMatches++;
                }
            }
            // Global test
            if(actMatches > expMatches){
                if(actMatches > locWin/2){ // quick check to make sure there are a reasonable number of actual matches with respect to the local window
                    // build matchvec
                    int matchVec[i + locWin];
                    for(j = 0; j < i + locWin; j++){
                        // if there is a match add it
                        if(soi[j] == ts[lents-locWin-i+j]){
                            matchVec[j] = 1;
                        }else{
                            matchVec[j] = 0;
                        }
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
            if(soi[i+locWin] == ts[lents-locWin-i-1]){numMatches++;}

        }
        // middle part of the comparison
        for(i = lensoi - locWin; i < lents - locWin; i++){
            actMatches = 0;
            expMatches = numMatches/(float)lensoi + globSig*sqrt(numMatches/(float)lensoi);
            expMatches = expMatches/(float)comber;
            // check soi[0:lensoi] == ts[lents - locWin - i: lents -locWin -i +lensoi]
            for(j = 0; j < lensoi; j+=comber){
                if(soi[j] == ts[lents - locWin -i+j]){
                    actMatches++;
                }
            }
            if(actMatches > expMatches){
                    // build matchvec
                    int matchVec[lensoi];
                    for(j = 0; j < lensoi; j++){
                        // if there is a match add it
                        if(soi[j] == ts[lents - locWin -i+j]){
                            matchVec[j] = 1;
                        }else{
                            matchVec[j] = 0;
                        }
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
                if(soi[j+i-lents+locWin] == ts[j]){
                    actMatches++;
                }
            }
            // Global test
            if(actMatches > expMatches){
                if(actMatches > locWin/2){
                    // build matchvec
                    int matchVec[lensoi -i+lents-locWin];
                    for(j = 0; j < lensoi -i+lents-locWin; j++){
                        // if there is a match add it
                        if(soi[j+i-lents+locWin] == ts[j]){
                            matchVec[j] = 1;
                        }else{
                            matchVec[j] = 0;
                        }
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
            if(soi[i-lents+locWin] == ts[lents - locWin -i+lensoi-1]){numMatches++;}
        }
    }else{
        // the other option: lents < lensoi
        for(i = 0; i < lents - locWin; i++){
            actMatches = 0;
            expMatches = numMatches/(float)(locWin+i) + globSig*sqrt(numMatches/(float)(locWin+i));
            expMatches = expMatches/(float)comber;
            // check soi[0:i+locWin-1] == ts[lents-locWin-i:lents-1]
            for(j = 0; j < i + locWin; j+=comber){
                // if there is a match add it
                if(soi[j] == ts[lents-locWin-i+j]){
                    actMatches++;
                }
            }
            // Global test
            if(actMatches > expMatches){
                if(actMatches > locWin/2){
                    // build matchvec
                    int matchVec[i + locWin];
                    for(j = 0; j < i + locWin; j++){
                        // if there is a match add it
                        if(soi[j] == ts[lents-locWin-i+j]){
                            matchVec[j] = 1;
                        }else{
                            matchVec[j] = 0;
                        }
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
            if(soi[i+locWin] == ts[lents-locWin-i-1]){numMatches++;}
        }
        // middle
        for(i = lents - locWin; i < lensoi - locWin; i++){
            actMatches = 0;
            expMatches = numMatches/(float)lents + globSig*sqrt(numMatches/(float)lents);
            expMatches = expMatches/(float)comber;
            // check soi[0:lensoi-1] == ts[lents - locWin - i: lents -1 -locWin -i +lensoi]
            for(j = 0; j < lents; j+=comber){
                if(soi[i-lents+locWin+j] == ts[j]){
                    actMatches++;
                }
            }
            // Global test
            if(actMatches > expMatches){
                    // build matchvec
                    int matchVec[lents];
                    for(j = 0; j < lents; j++){
                        // if there is a match add it
                        if(soi[i-lents+locWin+j] == ts[j]){
                            matchVec[j] = 1;
                        }else{
                            matchVec[j] = 0;
                        }
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
            expMatches = expMatches/(float)comber;
            // check soi[i - lents + locWin: lensoi -1] == ts[lents-locWin-i:lents-1]
            for(j = 0; j < lensoi -i+lents-locWin; j+=comber){
                if(soi[j+i-lents+locWin] == ts[j]){
                    actMatches++;
                }
            }
            // Global test
            if(actMatches > expMatches){
                if(actMatches > locWin/2){
                    // build matchvec
                    int matchVec[lensoi -i+lents-locWin];
                    for(j = 0; j < lensoi -i+lents-locWin; j++){
                        // if there is a match add it
                        if(soi[j+i-lents+locWin] == ts[j]){
                            matchVec[j] = 1;
                        }else{
                            matchVec[j] = 0;
                        }
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
            if(soi[i-lents+locWin] == ts[lents-(i+ locWin-lensoi)-1]){numMatches++;}
        }
    }
    //printf("in gen align 3\n");
    return matchNum;
}

int LocalFilter(int *matchVec,int matLen,double locSig,int locWin,double sinExp, double gapBridge, int *tempMatch){
    double cumSum[matLen];
    int i,j;
    int matches = 0;
    double significant = locSig*sqrt(sinExp*locWin);
    int start = 0, stop = 0;
    int densityMatches = 0;
    
    for(i = 0; i < 10; i++){
        tempMatch[i] = -1;
    }
    
//     printf("Testing!\n");
//     for(i = 0; i < matLen; i++){
//         printf("%d ",matchVec[i]);
//     }
//     printf("\n");

    // make a cumulative sum of the matchvec, subtracting off the expected match
    
    cumSum[0] = matchVec[0] - sinExp;
    for(i = 1; i < matLen;i++){
        cumSum[i] = cumSum[i-1] + (double)matchVec[i] - sinExp;
    }

    // from the left find a significant jump in the cumsum
    for(i = locWin; i < matLen - locWin; i++){
        if(cumSum[i] - cumSum[i-locWin] > significant && cumSum[i+locWin] - cumSum[i] > significant){
            if(start == 0){
                start = i - locWin;
            }
        }else{
            if(start != 0){
                stop = i + locWin;
                for(j = 0; j < locWin-1; j++){ // need 2/4 to start/end
                     if(matchVec[start] + matchVec[start+1] + matchVec[start+2] + matchVec[start+3] < 2){start++;}
                     if(matchVec[stop] + matchVec[stop-1] + matchVec[stop-2] + matchVec[stop-3] < 2){stop--;}
                }
                if(matchVec[start] == 0){start++;}
                if(matchVec[stop] == 0){stop--;}
                if(matchVec[start] == 0){start++;}
                if(matchVec[stop] == 0){stop--;}
//                 for(j = 0; j < locWin-1; j++){ // need 1 to start
//                      if(matchVec[start] == 0){start++;}
//                      if(matchVec[stop] == 0){stop--;}
//                 }
                stop++;
                densityMatches = 0;
                for(j = start;j < stop; j++){
                    if(matchVec[j] == 1){densityMatches++;}
                }
                if((densityMatches)/(float)(stop-start) > 0.35){ // MINIMUM DENSITY IMPOSED
                    // MINIMUM CUTTOF IMPOSED (to find 11-mers etc)
                    if(stop > start && (stop - start > locWin)){
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
                        i = stop;
                    }
                }
                start = 0;
            }
        }
    }
    
    // Catches them if they are at the end
    if(start != 0){
        stop = i + locWin;
//         printf("End Catch! %d %d\n", start, stop);
//         printf("Testing!\n");
//         for(i = start; i < stop; i++){
//             printf("%d ",matchVec[i]);
//         }
//         printf("\n");
        for(j = 0; j < locWin-1; j++){ // need 2/4 to start/end
             if(matchVec[start] + matchVec[start+1] + matchVec[start+2] + matchVec[start+3] < 2){start++;}
             if(matchVec[stop] + matchVec[stop-1] + matchVec[stop-2] + matchVec[stop-3] < 2){stop--;}
        }
        if(matchVec[start] == 0){start++;}
        if(matchVec[stop] == 0){stop--;}
        if(matchVec[start] == 0){start++;}
        if(matchVec[stop] == 0){stop--;}
//         for(j = 0; j < locWin-1; j++){ // need 1 to start
//              if(matchVec[start] == 0){start++;}
//              if(matchVec[stop] == 0){stop--;}
//         }
        stop++;
//         printf("End Catch! %d %d\n", start, stop);
        densityMatches = 0;
        for(j = start;j < stop; j++){
             if(matchVec[j] == 1){densityMatches++;}
        }
//         printf("End Catch! %d\n", densityMatches);
        if((densityMatches)/(float)(stop-start) > 0.35){ // MINIMUM DENSITY IMPOSED

            if(stop > start && (stop - start > locWin)){
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
                i = stop;
            }
        }
    }
    //for(i = 0; i < 10; i++){printf("%d ",tempMatch[i]);}
    //printf("\n");
    if(matches > 0){
        return 1;
    }else{
        return 0;
    }
}

int TesterAlign_c(int soiLen, int soi[][1000], int tsLen, int ts[][1000], int tsPartInfo[], int alignmentInfo[], int soiLenTotal, int tsLens[], double globSig, double locSig, int locWin, int maxMatch, int maxInline, double gapBridge, int comber){
    int i,j,k;

    int localStart[maxMatch];
    int localEnd[maxMatch];
    int localShift[maxMatch];
    
    
//     for(j = 0; j < tsLen; j++){
//         printf("\n%d\n",j);
//         for(i = 0; i < 1000; i++){
//             printf("%d ",ts[j][i]);
//         }
//         printf("\n");
//     }

    int totalNumTs = tsPartInfo[tsLen-1]+1;
    int tsOutPos[totalNumTs];
    for(i = 0; i < totalNumTs; i++){
        tsOutPos[i] = 5*i*maxMatch;
    }
    
    int thePart = -1, modify;
    int counter = 0;
    for(i = 0; i < soiLen; i++){
        thePart = -1;
        for(j = 0; j < tsLen; j++){
            // printf("Comparing %d to %d\n",i,j);
            if(thePart == tsPartInfo[j]){
                counter++;
            }else{
                thePart = tsPartInfo[j];
                counter = 0;
            }
//             for(k = 0; k < 40; k++){
//                 printf("%d ",ts[j][k]);
//             }
//             printf("\n");
//             for(k = 0; k < 40; k++){
//                 printf("%d ",soi[i][k]);
//             }
//             printf("\n");
            GenAlignment_c(soi[i], ts[j], min(1000,soiLenTotal-i*1000), min(1000,tsLens[tsPartInfo[j]]-1000*counter), globSig, locSig, locWin, maxMatch, maxInline, gapBridge, comber, localStart, localEnd, localShift);
            for(k = 0; k < maxMatch; k++){
                //printf("Start: %d  Stop: %d\n",localStart[k],localEnd[k]);
                if(localStart[k] != localEnd[k]){
                    // links across gaps
                    modify = 0;
                    if(tsOutPos[tsPartInfo[j]] != 5*tsPartInfo[j]*maxMatch){
                        if(tsOutPos[tsPartInfo[j]-3] == localShift[k] + i*1000 + max(tsLens[tsPartInfo[j]] -1000*(counter+1),0)){
                            if(localStart[k] - tsOutPos[tsPartInfo[j]-4] < locWin*gapBridge){
                                //printf("EXTENSION!\n");
                                alignmentInfo[tsOutPos[tsPartInfo[j]]-4] = localEnd[k];
                                modify = 1;
                            }
                        }
                    }
                    if(modify == 0){
                        alignmentInfo[tsOutPos[tsPartInfo[j]]] = localStart[k];
                        alignmentInfo[tsOutPos[tsPartInfo[j]]+1] = localEnd[k];
                        alignmentInfo[tsOutPos[tsPartInfo[j]]+2] = localShift[k];
                        alignmentInfo[tsOutPos[tsPartInfo[j]]+3] = i;
                        alignmentInfo[tsOutPos[tsPartInfo[j]]+4] = counter;
                        //printf("%d %d %d %d %d\n",localStart[k],localEnd[k],localShift[k],i,j);
                        //printf("%d %d %d %d %d\n",alignmentInfo[tsOutPos[tsPartInfo[j]]],alignmentInfo[tsOutPos[tsPartInfo[j]]+1],alignmentInfo[tsOutPos[tsPartInfo[j]]+2],i,j);
                        tsOutPos[tsPartInfo[j]] += 5;
                    }
                }else{
                    k = maxMatch;
                }
            }
        }
    }
    
//     outputPos = 0;
//     for(i = 0; i < soiLen; i++){
//         for(j = 0; j < tsLen; j++){
//             for(k = 0; k < maxMatch; k++){
//                 printf("%d %d %d\n",alignmentInfo[outputPos],alignmentInfo[outputPos+1],alignmentInfo[outputPos+2]);
//                 outputPos += 3;
//             }
//         }
//     }
    return 1;
}

int min(int a, int b){
    //return (1/2)*(a + b - abs(a-b));
    if(a < b){
        return a;
    }else{
        return b;
    }
}

int max(int a, int b){
    //return (1/2)*(a + b - abs(a-b));
    if(a < b){
        return b;
    }else{
        return a;
    }
}

/* Create any sort of [size] array */
int *int_array(int size) {
    return (int *) malloc(size*sizeof(int));
}
/* Create a two-dimension array [size][1000] */
int (*int_array_1000(int size))[1000] {
    return (int (*)[1000]) malloc(size*1000*sizeof(int));
}

void a_set(int i, int j, int val, int a[][1000]) {
   a[i][j] = val;
}

int a_get(int i, int j, int a[][1000]) {
   return a[i][j];
}