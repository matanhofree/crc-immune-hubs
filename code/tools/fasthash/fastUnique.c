/*******************************************************************************
 ** Fast C++ based unique function.
 **
 ** Input:
 ** inputVector A
 **
 ** Output: [ keyVect, orderA, orderKey, countVect, dupList
 **  keyVect - a unique key vector in order of occurrence
 **  keyVect = A[orderA]
 **  A = keyVect[orderKey]
 **  countVect = count occurances of keyVect in A
 **  dupList[i] = positions of keyVect[i] in A
 *******************************************************************************/

#include "mex.h"
#include "utarray.h"
#include "uthash.h"
#include "string.h"

#define KEYSIZE sizeof(uint64_T)

struct keyID {
    //char name[keySize];        /* key (string is WITHIN the structure) */
    uint64_T key;
    mwSize idx;
    UT_hash_handle hh;         /* makes this structure hashable */
};


struct keyID *keyMap = NULL;

UT_icd uint64_icd = {sizeof(uint64_T), NULL, NULL, NULL };
// UT_icd uint64_icd = {sizeof(uint64_T), NULL, NULL, NULL };


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

    uint64_T *A = (uint64_T *) mxGetData(prhs[0]);
    mwSize lenA = mxGetNumberOfElements(prhs[0]);


//    struct keyID *keyMap = NULL;
    struct keyID *curr = NULL;
    struct keyID *mapMemPool = (struct keyID*)mxCalloc(lenA,sizeof(struct keyID));
//    uint64_T* currKey = NULL;
//    const mxArray* cellP;
//    mwSize keyLen = 0;

    mwSize cidx = 0;
//   std::vector<uint64_T> occCounter;
//    std::vector<std::vector<mwSize> > posDupList;
    UT_array *occCounter;
    utarray_new(occCounter,&uint64_icd);
    uint64_T ul, *currOcc;

    mexPrintf("Start index\n");
    for (mwSize i = 0; i < lenA;++i)
    {

            HASH_FIND(hh,keyMap,A[i],KEYSIZE,curr);
            if (curr==NULL)
            {
                // Id Is not in hash
                curr = mapMemPool+cidx;

                curr->key = A[i];
                curr->id = cidx++;
                HASH_ADD( hh, keyMap, curr->key, KEYSIZE, curr );

                ul = 0;

                utarray_push_back(occCounter,&ul);
//                posDupList[i].push_back(i);
            }
            else
            {
                mwSize cPos = curr->idx;
                p=(long*)utarray_next(occCounter,);
  //              posDupList[cPos].push_back(i);
            }
        }
    }

    plhs[0] = mxCreateNumericMatrix(cidx, 1, mxUINT64_CLASS , mxREAL );
    plhs[1] = mxCreateNumericMatrix(cidx, 1, mxUINT64_CLASS , mxREAL );
    plhs[2] = mxCreateNumericMatrix(lenA, 1, mxUINT64_CLASS , mxREAL );
    plhs[3] = mxCreateNumericMatrix(cidx, 1, mxUINT64_CLASS , mxREAL );
    plhs[4] = mxCreateCellArray(cidx, 1);


    uint64_T* keyVect = (uint64_T*)mxGetPr(plhs[0]);
    uint64_T* orderA = (uint64_T*)mxGetPr(plhs[0]);
    uint64_T* orderKey = (uint64_T*)mxGetPr(plhs[0]);
    uint64_T* countVect = (uint64_T*)mxGetPr(plhs[0]);
    mxArray* dupList = (uint64_T*)mxGetPr(plhs[0]);

    mexPrintf("Start output \n");
    struct my_struct *tmp;

    mwSize zi = 0;
    uint64_T pIdx = 0;
    HASH_ITER(hh, keyMap, curr, tmp) {
        keyVect[zi] = curr->key;
        countVect[zi] = mxCreateNumericMatrix(cidx, 1, mxUINT64_CLASS , mxREAL );
//        std::copy(posDupList[zi].begin(), posDupList[zi].end(), mxGetPr(countVect[zi]));
    }

    return;
}/*mexFunction()*/


