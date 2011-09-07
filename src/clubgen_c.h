/* File: clubgen_c.h */

int GenAlignment_c(int soi[], int ts[], int lensoi, int lents, double globSig, double locSig, int locWin, int maxMatch, int maxInline, double gapBridge, int comber, int c_start[], int c_end[], int c_shift[]);
int LocalFilter(int matchVec[],int matLen,double locSig,int locWin,double sinExp,double gapBridge,int tempMatch[]);
int TesterAlign_c(int soiLen, int soi[][1000], int tsLen, int ts[][1000], int tsPartInfo[], int alignmentInfo[], int soiLenTotal, int tsLens[], double globSig, double locSig, int locWin, int maxMatch, int maxInline, double gapBridge, int comber);
// int CUDA_Align_c(int soiLen, int soi[][1000], int tsLen, int ts[][1000], int tsPartInfo[], int alignmentInfo[], int soiLenTotal, int tsLens[], double globSig, double locSig, int locWin, int maxMatch, int maxInline, double gapBridge, int comber);
int *int_array(int size);
int (*int_array_1000(int size))[1000];
void a_set(int i, int j, int val, int a[][1000]);
int a_get(int i, int j, int a[][1000]);
int min(int a, int b);
int max(int a, int b);