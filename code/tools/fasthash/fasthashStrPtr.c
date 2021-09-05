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
#include "string.h"

#define IS_REAL_2D_FULL_DOUBLE(P) (!mxIsComplex(P) && mxGetNumberOfDimensions(P) == 2 && !mxIsSparse(P) && mxIsDouble(P))
#define IS_REAL_SCALAR(P) (IS REAL 2D FULL DOUBLE(P) && mxGetNumberOfElements(P) == 1)
#define MM(m,n,D) (m + (D*n))
#define xMM(A,m,n,D) (A[m+(D*n)])
#define keySize 16
#define MIN(a,b) (((a)<(b))?(a):(b))


struct keyID {
    //char name[keySize];        /* key (string is WITHIN the structure) */
    char* name;
    unsigned int id;
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


	const mxArray *keyVector = prhs[0];
	const mxArray *inputVector = prhs[1];

	mwSize lenKeys = mxGetNumberOfElements(keyVector);
	mwSize lenInputVector = mxGetNumberOfElements(inputVector);
	

	struct keyID *keyMap = NULL;	
	
	struct keyID *curr = NULL;
	struct keyID *mapMemPool = (struct keyID*)mxCalloc(lenKeys,sizeof(struct keyID));	
	struct keyID **memPtr = &mapMemPool;
	char * currKey = NULL;
	const mxArray* cellP;
	mwSize keyLen = 0;
	
	unsigned int cidx = 0;
	
	//mexPrintf("Start READ\n");
	for (mwSize i = 0; i < lenKeys;++i) 
	{
		
		cellP = mxGetCell(keyVector,i);
		if (cellP != NULL)
		{
			currKey = mxArrayToString(cellP);
		}
		//mexPrintf("String %s\n",currKey);
		if (currKey)
		{
			keyLen = strlen(currKey);
			
		    HASH_FIND_STR(keyMap, currKey, curr);  /* id already in the hash? */

		    if (curr==NULL) 
		    {
		      curr = mapMemPool+cidx;
		      //mexPrintf("Debug: %d %x\n",cidx,curr);
			  curr->name = currKey;
		      curr->id = cidx++;
		      // HASH_ADD_STR( keyMap, name, curr );  /* id: name of key field */
		      HASH_ADD_KEYPTR( hh, keyMap, curr->name, strlen(curr->name), curr );
		    }
		    else
		    {
		    	mexPrintf("Warn: duplicate key %s\n",currKey);
		    }
		}			
	}

	//int ddims[] = { lenInputVector, 1 };
	//plhs[0] = mxCreateNumericArray(2,ddims,mxUINT32_CLASS,0);
	plhs[0] = mxCreateDoubleMatrix(lenInputVector,1,mxREAL);
	double* idxVector = (double *)mxGetPr(plhs[0]);
	double nan = mxGetNaN();
	
	char* currOut;
	
	//mexPrintf("Start data\n");
	for (mwSize i = 0;i<lenInputVector;++i)
	{
		cellP = mxGetCell(inputVector,i);
		if (cellP != NULL)
		{
			currOut = mxArrayToString(cellP);
		}
		if (currOut)
		{
			HASH_FIND_STR( keyMap, currOut, curr);			
		}		
		if (curr) 
		{
			idxVector[i] = curr->id+1;
		}
		else
		{
			mexPrintf("Warn: key %s not found\n",currOut);
			idxVector[i] = nan;
		}
	}

	return;
}/*mexFunction()*/

