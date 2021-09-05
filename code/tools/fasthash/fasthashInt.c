/*******************************************************************************
 ** Fast C++ based has function.
 ** 
 ** Input:
 ** keyArray
 ** inputVector
 ** 
 ** Output:
 ** idxVector
 **
 *******************************************************************************/

#include "mex.h"
#include "uthash.h"


//#define IS_REAL_2D_FULL_DOUBLE(P) (!mxIsComplex(P) && mxGetNumberOfDimensions(P) == 2 && !mxIsSparse(P) && mxIsDouble(P))
//#define IS_REAL_SCALAR(P) (IS REAL 2D FULL DOUBLE(P) && mxGetNumberOfElements(P) == 1)
//#define MM(m,n,D) (m + (D*n))
//#define xMM(A,m,n,D) (A[m+(D*n)])
//#define keySize 16
//#define MIN(a,b) (((a)<(b))?(a):(b))


struct keyID {
    uint64_T name;        /* key (string is WITHIN the structure) */
	uint64_T id;
    UT_hash_handle hh;         /* makes this structure hashable */
};


/*******************************************************************************
 **entry point for the mex call
 **nlhs - number of outputs
 **plhs - pointer to array of outputs
 **nrhs - number of inputs
 **prhs - pointer to array of inputs
 *******************************************************************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/*****************************************************************************
	 ** this function takes a flag and a variable number of arguments
	 ** depending on the value of the flag and returns either a construct 
	 ** containing probability estimates, a merged vector or a double value 
	 ** representing an entropy or mutual information
	 *****************************************************************************/


	const uint64_T *keyVector = (uint64_T *)mxGetData(prhs[0]);
	const uint64_T *inputVector = (uint64_T *)mxGetData(prhs[1]);

	mwSize lenKeys = mxGetNumberOfElements(prhs[0]);
	mwSize lenInputVector = mxGetNumberOfElements(prhs[1]);


	struct keyID *keyMap = NULL;
	struct keyID *curr = NULL;

	struct keyID *mapMemPool = (struct keyID*)mxCalloc(lenKeys,sizeof(struct keyID));	
//	struct keyID **memPtr = &mapMemPool;
//	uint64_T currKey = 0;
//	const mxArray* cellP;
//	mwSize keyLen = 0;
	
	mwSize cidx = 0;
	
	//mexPrintf("Start READ -- %d\n",lenKeys);
	for (mwSize i = 0; i < lenKeys;++i) 
	{
		
        uint64_T currKey = keyVector[i];
//    	mexPrintf("String %d\n",currKey);
		if (currKey)
		{
		    HASH_FIND(hh, keyMap, &currKey, sizeof(uint64_T), curr);  /* id already in the hash? */
//            mexPrintf("After find %d\n",currKey);
		    if (curr==NULL) 
		    {
		      curr = mapMemPool+cidx;
//		      //mexPrintf("Debug: %d %x\n",cidx,curr);
//			  strncpy(curr->name, currKey,keySize);
//		      // curr->name = currKey;
              curr->name = currKey;
		      curr->id = cidx++;
              HASH_ADD_KEYPTR(hh, keyMap, &(curr->name), sizeof(uint64_T), curr);  /* id: name of key field */
		    }
		    else
		    {
		    	mexPrintf("Warn: duplicate key %s\n",currKey);
		    }
		}			
	}

	plhs[0] = mxCreateDoubleMatrix(lenInputVector,1,mxREAL);
	double* idxVector = (double *)mxGetPr(plhs[0]);
	double nan = mxGetNaN();
	
	uint64_T currOut;
	
	// mexPrintf("Start data\n");
	for (mwSize i = 0;i<lenInputVector;++i)
	{
        currOut = inputVector[i];
		HASH_FIND(hh, keyMap, &currOut, sizeof(uint64_T), curr);

		if (curr)
		{
			idxVector[i] = curr->id+1;
		}
		else
		{
			 // mexPrintf("Warn: key %s not found\n",currOut);
			idxVector[i] = nan;
		}
	}

	return;
}/*mexFunction()*/

