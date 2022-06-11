#include "mex.h"
#include "matrix.h"
#include <string.h>
#include <omp.h>
#include <bits/stdc++.h> 

using namespace std;

double** mat;
int** mot;

int* QDXY;

double max(double a, double b)
{
	if(a> b)
		return a;
	else
		return b;
}

void transMotifs(int beginData, int l, int c)
{
	double last, blast, val;

	for (int i = beginData + 2; i < l; ++i)
	{
		for (int j = 0; j < c; ++j)
		{
			val = mat[i][j];
			blast = mat[i-1][j];
			last = mat[i-2][j];

			if (last > blast && blast > val)
	        	mot[i-2][j] = 1;
	        else if (last > blast && blast < val)
	        	mot[i-2][j] = 2;
	        else if (last < blast && blast > val)
	        	mot[i-2][j] = 3;
	        else if (last < blast && blast < val)
	        	mot[i-2][j] = 4;
	        else
	        	mot[i-2][j] = 5;
		}	
	}
}

void motifSync(int begin, int end, int slidWindow, int c, int taoMin, int taoMax, int dim)
{
	// Initialize variables
	int soma_XY, soma_YX, max_CXY, max_CYX;

	// Third dimension index of the 3d array
	int dim3 = dim * (c * c);

	for (int i = 0; i < c; ++i)
	{
		for (int j = i+1; j < c; ++j)
		{
			max_CXY = 0;
			max_CYX = 0;

			for (int l = taoMin; l <= taoMax ; ++l)
			{
				soma_XY = 0;
				soma_YX = 0;

				for (int k = begin; k < begin + slidWindow; ++k)
				{
					soma_XY += mot[k][i] == mot[k+l][j];
					soma_YX += mot[k][j] == mot[k+l][i];
				}
				if(soma_XY > max_CXY )
				 	max_CXY = soma_XY;

				if(soma_YX > max_CYX )
					max_CYX = soma_YX;

			}

			if (max_CYX > max_CXY)
				QDXY[j + i*c + dim3] = QDXY[i + j*c + dim3] = max_CYX;
			else
				QDXY[j + i*c + dim3] = QDXY[i + j*c + dim3] = max_CXY;

			
		}
	}
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	double *Data, threshold;
	int beginData, endData, numNodes, slidWindow, taoMin,taoMax;    /* pointers to input matrices */
	size_t mrows,ncols;
	mxArray *Awork, *mxPivot;

  
  	// Check the number of inputs
	if(nrhs!=7) {
		mexErrMsgIdAndTxt( "MATLAB:Motif_Sync_mex:invalidNumInputs",
	        "Seven input required.");
	} 

	Data = mxGetPr(prhs[0]); /* pointer to first input matrix */
	beginData = mxGetPr(prhs[1])[0];  /* pointer to second input matrix */
	endData = mxGetPr(prhs[2])[0];  /* pointer to third input matrix */
	numNodes = mxGetPr(prhs[3])[0];  /* pointer to fourth input matrix */
	slidWindow = mxGetPr(prhs[4])[0];  /* pointer to fifth input matrix */
	taoMin = mxGetPr(prhs[5])[0];  /* pointer to sixth input matrix */
    taoMax = mxGetPr(prhs[6])[0];  /* pointer to seventh input matrix */

	mrows = mxGetM(prhs[0]);   /* get number of rows of data */
    ncols = mxGetN(prhs[0]);   /* get number of columns of data */

	/* Initialize matrix pointer which will receive the data */
	mat = new double*[mrows];
	for(size_t i = 0 ; i < (mrows) ; ++i)
	    mat[i] = new double[ncols];

	/* Initialize matrix pointer which will receive the motifs */
	mot = new int*[mrows];
	for(size_t i = 0 ; i < (mrows) ; ++i)
	    mot[i] = new int[ncols];

	/* Receive the data */
    for (int col=0; col < ncols; col++)
    	for (int row=0; row < mrows; row++)
            mat[row][col] = Data[row+col*mrows];

    /* Calculate the slide window and the size of the TVG */
	double window = slidWindow + taoMax;
	int sizeTVG = endData - window - beginData ;

    /* Transform the data into motifs */
	transMotifs(beginData, endData, numNodes);

    /* Define the dimension of the 3d array (tensor) which will receive the outputs */
	mwSize dims[3] = {numNodes,numNodes,sizeTVG};


    /* Atach the pointer to the array which will receive the TVG */
    plhs[0] = mxCreateNumericArray(3, dims, mxINT32_CLASS, mxREAL);
    QDXY = (int *)mxGetData(plhs[0]);

    int i = 0;

    /* Call the functon to calculate the TVG by the motifs synchronization */
    #pragma omp parallel for
	for (i = 0; i < sizeTVG; ++i)
		motifSync(i, i+ window, slidWindow, numNodes, taoMin,taoMax, i);

	int n = sizeof(QDXY)/sizeof(QDXY[0]);

	sort(QDXY, QDXY+n);

	/* Destroy the matrix pointer */
    for(size_t i = 0 ; i < endData+1 ; ++i)
	    delete mat[i];
	delete[] mat;

	/* Destroy the matrix pointer */
	for(size_t i = 0 ; i < endData+1 ; ++i)
	    delete mot[i];
	delete[] mot;
}
