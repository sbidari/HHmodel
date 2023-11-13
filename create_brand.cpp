/*********************************************************************
 * Keep in mind:
 * <> Use 0-based indexing as always in C or C++
 * <> Indexing is column-based as in Matlab (not row-based as in C)
 * <> Use linear indexing.  [i-1+(j-1)*num_row] in C++ instead of [i][j] in Matlab (note starting index)

 ********************************************************************/
#include <matrix.h>
#include <mex.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <random>
#include <vector>
#include <algorithm>
using namespace std;

void mexFunction(int nlmxhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
//declare variables
    mxArray *mxX, *mxY;//input
    mxArray *mxZ;//output
    const mwSize *dims;
    double *H_prev,*P_any;//input
    double *N_any;//output
    int L;//length of input

//associate inputs
    mxX = mxDuplicateArray(prhs[0]);
    mxY = mxDuplicateArray(prhs[1]);
    
//figure out dimensions
    dims = mxGetDimensions(prhs[0]);//dimensional of the first input
    L = (int)dims[0];
    
//associate outputs
    mxZ = plhs[0] = mxCreateDoubleMatrix(L,1,mxREAL);
    
//associate pointers
    H_prev = mxGetPr(mxX);
    P_any = mxGetPr(mxY);
    N_any = mxGetPr(mxZ);
    
    //Generate random binomial using input param
    int i;
    for (i=0; i<L; i++) {
        //unsigned seed = chrono::system_clock::now().time_since_epoch().count();
//        default_random_engine generator (42);
        random_device rd;
        mt19937 gen(rd());
        binomial_distribution<int> distribution(H_prev[i],P_any[i]);
        N_any[i] = distribution(gen);
        distribution.reset();
    }
    return;
}
