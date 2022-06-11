#include "mex.h"
#include "matrix.h"
#include <string.h>
#include <omp.h>

using namespace std;

double** mat;
bool m;
int** mot;

bool* QXY;

double max(double a, double b)
{
	if(a > b)
		return a;
	else
		return b;
}


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	double *data;
	double threshold;
	size_t mrows,ncols;
	const mwSize *pDims;
	mxArray *Awork, *mxPivot;

  
  	// Check the number of inputs
	if(nrhs!=2) {
		mexErrMsgIdAndTxt( "MATLAB:Motif_Sync_mex:invalidNumInputs",
	        "Two input required.");
	}

	data = mxGetPr(prhs[0]); /* pointer to first input matrix */
	threshold = mxGetPr(prhs[1])[0];  /* pointer to second input matrix */

	mrows = mxGetM(prhs[0]);   /* get number of rows of data */
    ncols = mxGetN(prhs[0]);   /* get number of columns of data */

	/* Initialize matrix pointer which will receive the data */
	mat = new double*[mrows];
	for(size_t i = 0 ; i < (mrows) ; ++i)
	    mat[i] = new double[ncols];

	/* Receive the data */
    // for (int col=0; col < ncols; col++)
    // 	for (int row=0; row < mrows; row++)
    //         mat[row][col] = data[row+col*mrows];

    m = mat[0][0];

    plhs[0] = mxCreateLogicalArray(3, dims);
    QXY = (bool *)mxGetData(plhs[0]);

    //QXY[0] = (double) mat[0][0];
}
