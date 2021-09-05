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


    // const uint64_T *keyVector = (uint64_T *)mxGetData(prhs[0]);
    const uint64_T *inputVector = (uint64_T *)mxGetData(prhs[0]);
    mwSize lenInputVector = mxGetNumberOfElements(prhs[0]);

    mwSize lenKeys = mxGetScalar(prhs[1]);



    struct keyID *keyMap = NULL;
    struct keyID *curr = NULL;

    struct keyID *mapMemPool = (struct keyID*)mxCalloc(lenKeys,sizeof(struct keyID));
    mwSize cidx = 0;

    uint64_T* occCounter = (uint64_T *)mxCalloc(lenKeys,sizeof(uint64_T));
    memset(occCounter, 0, lenKeys*sizeof(uint64_T));

    plhs[0] = mxCreateNumericMatrix(lenKeys, 1, mxUINT64_CLASS, mxREAL );
    plhs[1] = mxCreateNumericMatrix(lenKeys, 1, mxUINT64_CLASS, mxREAL );
    plhs[2] = mxCreateNumericMatrix(lenInputVector, 1, mxUINT64_CLASS, mxREAL );
    plhs[3] = mxCreateNumericMatrix(lenKeys, 1, mxUINT64_CLASS, mxREAL );
    plhs[4] = mxCreateNumericMatrix(lenInputVector, 1, mxUINT64_CLASS, mxREAL );

    uint64_T* keyVect = (uint64_T*)mxGetData(plhs[0]);
    uint64_T* orderA = (uint64_T*)mxGetData(plhs[1]);
    uint64_T* orderVect = (uint64_T*)mxGetData(plhs[2]);
    uint64_T* countVect = (uint64_T*)mxGetData(plhs[3]);
    uint64_T* dupVect = (uint64_T*)mxGetData(plhs[4]);

    //mexPrintf("Start READ -- %d %d\n",lenKeys,lenInputVector);
    for (mwSize i = 0; i < lenInputVector;++i) {

        uint64_T currKey = inputVector[i];

        HASH_FIND(hh, keyMap, &currKey, sizeof(uint64_T), curr);  /* id already in the hash? */
        if (curr == NULL) {
            keyVect[cidx] = currKey;

            occCounter[cidx] = 1;
            orderA[cidx] = i + 1;

            curr = mapMemPool + cidx;
            curr->name = currKey;
            curr->id = cidx;
            HASH_ADD_KEYPTR(hh, keyMap, &(curr->name), sizeof(uint64_T), curr);  /* id: name of key field */
            ++cidx;
        }
        else {
            ++occCounter[curr->id];
        }
    }

    uint64_T* occCumSum = (uint64_T *)mxCalloc(lenKeys+1,sizeof(uint64_T));
    occCumSum[0] = 0;

    for (mwSize zi = 1; zi<=cidx;++zi)
        occCumSum[zi] = occCumSum[zi-1] + occCounter[zi-1];

    uint64_T currOut;

    //mexPrintf("Start data -- cidx %d \n",cidx);
    for (mwSize i = 0;i<lenInputVector;++i)
    {
        currOut = inputVector[i];
        HASH_FIND(hh, keyMap, &currOut, sizeof(uint64_T), curr);

        if (curr)
        {
            uint64_T cPos = curr->id;
            // mexPrintf("%d) Pos: %d\n",i,cPos);
            orderVect[i] = cPos+1;

            mwSize zP = occCumSum[cPos]+countVect[cPos];

            //mexPrintf("%d) Assign dupvect %d\n",i,zP);

            if (zP < lenInputVector) {
                dupVect[zP] = i + 1;
            }
            else
            { mexPrintf("Error!\n"); }

            //mexPrintf("%d) Dup\n",i);
            ++countVect[cPos];
            //mexPrintf("%d) count\n",i);
        }
        else
        {
            mexPrintf("Warn: key %s not found\n",currOut);
        }
    }

    return;
}/*mexFunction()*/


