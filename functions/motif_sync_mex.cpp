#include "mex.h"
#include "matrix.h"
#include <string.h>
#include <omp.h>
#include <stdio.h>

using namespace std;

double** mat;
int** mot;

bool* QXY;
bool* QDH;
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

			if (last > blast && blast > val)        /* \  */
	        	mot[i-2][j] = 1;
	        else if (last > blast && blast < val)   /* V  */
	        	mot[i-2][j] = 2;
            else if (last < blast && blast > val)   /* V invertido */
	        	mot[i-2][j] = 3;
            
	        else if (last < blast && blast < val)   /* /   */
	        	mot[i-2][j] = 4;
            else if (last==blast && blast == val && val==0)   /* Null motif */
	        	mot[i-2][j] = 0;
	        else
	        	mot[i-2][j] = 5;
		}	
	}
}

void motifSync(int begin, int end, int slidWindow, int c, int taoMin, int taoMax, double threshold, int dim)
{
	// Initialize variables
	int soma_XY, soma_YX, max_CXY, max_CYX, tao_CXY, tao_CYX, qXY, max;

	// Third dimension index of the 3d array
	int dim3 = dim * (c * c);

	// Calculate the real threshold
	double rthreshold = threshold * slidWindow;
    
    double sumI,sumJ;

	for (int i = 0; i < c; ++i)
	{
		for (int j = i+1; j < c; ++j)
		{
			tao_CXY = 0;
			tao_CYX = 0;
			max_CXY = 0;
			max_CYX = 0;
			max = 0;

			for (int l = taoMin; l <= taoMax ; ++l)
			{
				soma_XY = 0;
				soma_YX = 0;
                sumI=sumJ=0;
				for (int k = begin; k < begin + slidWindow; ++k)
				{
                    sumI+=(mat[k+l][i]!=0);
                    sumJ+=(mat[k+l][j]!=0);
					soma_XY += mot[k][i] == mot[k+l][j];
					soma_YX += mot[k][j] == mot[k+l][i];
				}
                if (sumI<slidWindow || sumJ<slidWindow) // disregard flat data
                    soma_XY = soma_YX = 0;
				if(soma_XY > max_CXY )
				{
				 	max_CXY = soma_XY;
				 	tao_CXY = l;
				}
				if(soma_YX > max_CYX )
				{
					max_CYX = soma_YX;
					tao_CYX = l;
				}
			}

			// Return max logical value if the synchronization value is bigger than the trheshold
            max = max_CYX >= rthreshold || max_CXY >= rthreshold;
//             if(i==103 && j==254 && begin==0)
//             {
//                 FILE * pFile;
//                 pFile = fopen ("debug.txt","a");
//                 fprintf(pFile,"motifs i %d %d %d %d\n",mot[0][i],mot[1][i],mot[2][i],mot[3][i]);
//                 fprintf(pFile,"motifs j %d %d %d %d\n",mot[0][j],mot[1][j],mot[2][j],mot[3][j]);
//                 fprintf(pFile,"%d %d %d %g %d %d\n",max,max_CYX,max_CXY,rthreshold,tao_CXY,tao_CYX);
//                 fclose(pFile);
//             }
			// Binary TVG
			QXY[i + j*c + dim3] = QXY[j + i*c + dim3] = max;

			if(max)
			{
				qXY = max_CXY - max_CYX;
				if(qXY > 0 || (qXY==0 && tao_CXY<tao_CYX)) // caso a sincronizaçao de X para Y for maior ou (seja a mesma mas tenha um tempo de atraso menor)
				{
                    if(tao_CXY==0)      // non directed edge
                        QDXY[j + i*c + dim3] = QDXY[i + j*c + dim3] = 1;
                    else
                    {
					// Weigthed oriented TVG
					QDXY[i + j*c + dim3] =  tao_CXY + 1;
					// Oriented Binary TVG
					QDH[i + j*c + dim3] = max;
                    }
					
				}
				else if (qXY < 0 || (qXY==0 && tao_CXY>=tao_CYX))
				{
                    if(tao_CYX==0)      // non directed edge
                        QDXY[j + i*c + dim3] = QDXY[i + j*c + dim3] = 1;
                    else
                    {
                        // Weigthed oriented TVG
                        QDXY[j + i*c + dim3] =  tao_CYX + 1;
                        // Oriented Binary TVG
                        QDH[j + i*c + dim3] = max;
                    }
				}
                else 
                {    
					QDXY[j + i*c + dim3] = QDXY[i + j*c + dim3] = 1;
				}
			}

			
		}
	}
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	double *Data, threshold;
	int beginData, endData, numNodes, slidWindow, taoMin, taoMax;    /* pointers to input matrices */
	size_t mrows,ncols;
	mxArray *Awork, *mxPivot;

  
  	// Check the number of inputs
	if(nrhs!=8) {
		mexErrMsgIdAndTxt( "MATLAB:Motif_Sync_mex:invalidNumInputs",
	        "Eigth input required.");
	} 

	Data = mxGetPr(prhs[0]); /* pointer to first input matrix */
	beginData = mxGetPr(prhs[1])[0];  /* pointer to second input matrix */
	endData = mxGetPr(prhs[2])[0];  /* pointer to third input matrix */
	numNodes = mxGetPr(prhs[3])[0];  /* pointer to fourth input matrix */
	slidWindow = mxGetPr(prhs[4])[0];  /* pointer to fifth input matrix */
	taoMin = mxGetPr(prhs[5])[0];  /* pointer to sixth input matrix */
    taoMax = mxGetPr(prhs[6])[0];  /* pointer to seventh input matrix */    
	threshold = mxGetPr(prhs[7])[0];  /* pointer to eighth input matrix */

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
	int sizeTVG = endData - window - beginData;

    /* Transform the data into motifs */
	transMotifs(beginData, endData, numNodes);

    /* Define the dimension of the 3d array (tensor) which will receive the outputs */
	mwSize dims[3] = {numNodes,numNodes,sizeTVG};

	/* Atach the pointer to the array which will receive the TVG */
	plhs[0] = mxCreateLogicalArray(3, dims);
    QXY = (bool *)mxGetData(plhs[0]);

    /* Atach the pointer to the array which will receive the oriented TVG */
    plhs[1] =mxCreateLogicalArray(3, dims);
    QDH = (bool *)mxGetData(plhs[1]);

    /* Atach the pointer to the array which will receive the weighted oriented TVG */
    plhs[2] = mxCreateNumericArray(3, dims, mxINT32_CLASS, mxREAL);
    QDXY = (int *)mxGetData(plhs[2]);

    int i = 0;

    /* Call the functon to calculate the TVG by the motifs synchronization */
    #pragma omp parallel for
	for (i = 0; i < sizeTVG; ++i)
		motifSync(i, i+ window, slidWindow, numNodes, taoMin, taoMax, threshold, i);


	plhs[3] = mxCreateNumericMatrix(endData-2, numNodes, mxINT32_CLASS, mxREAL);
    int* motif = (int *)mxGetData(plhs[3]);
    
	for (int col=0; col < numNodes; col++) {
    	for (int row=0; row < endData-2; row++) {
            motif[row + col*(endData-2)] =  mot[row][col];
        }
    }
    
	/* Destroy the matrix pointer */
    for(size_t i = 0 ; i < endData+1 ; ++i)
	    delete mat[i];
	delete[] mat;

	/* Destroy the matrix pointer */
	for(size_t i = 0 ; i < endData+1 ; ++i)
	    delete mot[i];
	delete[] mot;
}
