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
#include <mex.h>
#include <vector>
#include <unordered_map>
#include <mexplus.h>
#include <stdlib.h> // size_t, malloc, free
#include <new> // bad_alloc, bad_array_new_length
#include <string>

template <class T> struct Mallocator {
    typedef T value_type;
    Mallocator() noexcept { } // default ctor not required
    template <class U> Mallocator(const Mallocator<U>&) noexcept { }
    template <class U> bool operator==(
            const Mallocator<U>&) const noexcept { return true; }
    template <class U> bool operator!=(
            const Mallocator<U>&) const noexcept { return false; }

    T * allocate(const size_t n) const {
        if (n == 0) { return nullptr; }
        if (n > static_cast<size_t>(-1) / sizeof(T)) {
            throw std::bad_array_new_length();
        }
        void * const pv = mxMalloc(n * sizeof(T));
        if (!pv) { throw std::bad_alloc(); }
        return static_cast<T *>(pv);
    }
    void deallocate(T * const p, size_t) const noexcept {
        mxFree(p);
    }
};

using namespace mexplus;
using namespace std;

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

    //string *A = (string *) mxGetData(prhs[0]);
    vector<string> A = MxArray::to<vector<string> >(prhs[0]);
    mwSize lenA = mxGetNumberOfElements(prhs[0]);
    std::unordered_map<string, mwSize,std::hash<string>,std::equal_to<string>,Mallocator<std::pair<const string,mwSize> > > keyMap;
    mwSize cidx = 0;
    std::vector<mwSize> occCounter;
    std::vector<std::vector<mwSize> > posDupList;
    std::vector<mwSize> orderA;
    std::vector<mwSize> orderVect = std::vector<mwSize>(lenA,0);
    string keyVal;

    //mexPrintf("Start index %d\n",lenA);
    for (mwSize i = 0; i < lenA;++i)
    {
        //auto curr = keyMap.find(A[i]);
        keyVal = A[i];
        std::unordered_map<string,mwSize>::const_iterator curr = keyMap.find(A[i]);
        //mexPrintf("After find %s\n",keyVal.c_str());


        if (curr==keyMap.end())
        {
            keyMap[keyVal] = cidx;
            occCounter.push_back(1);
            posDupList.push_back(std::vector<mwSize> (1,i+1));
            //exPrintf("After new index\n");

            orderA.push_back(i+1);
            orderVect[i] = cidx+1;
            ++cidx;
        }
        else
        {
            //mexPrintf("Before old index\n");

            mwSize cPos = curr->second;
            occCounter[cPos]++;
            posDupList[cPos].push_back(i+1);
            orderVect[i] = cPos+1;
        }

    }


    //mexPrintf("End index\n");

    plhs[0] = mxCreateCellMatrix(cidx, 1);
    //plhs[1] = mxCreateNumericMatrix(cidx, 1, mxUINT64_CLASS, mxREAL );
    //plhs[2] = mxCreateNumericMatrix(lenA, 1, mxUINT64_CLASS, mxREAL );
    plhs[3] = mxCreateNumericMatrix(cidx, 1, mxUINT64_CLASS, mxREAL );
    plhs[4] = mxCreateCellMatrix(cidx, 1);


    mxArray *keyVect = plhs[0];
    mwSize *countVect = (mwSize*)mxGetData(plhs[3]);
    mxArray *dupList = plhs[4];
    //mexPrintf("Start output \n");

    plhs[1] = mexplus::MxArray::from(orderA);
    plhs[2] = mexplus::MxArray::from(orderVect);

    mwSize zi = 0;
    //mexPrintf("Start output loop \n");
    for (auto curr=keyMap.begin(); curr!=keyMap.end(); ++curr) {
        zi = curr->second;
        // keyVect[zi] = curr->first;
        //string keyVal = curr->first;
        mxSetCell(keyVect,zi,mxCreateString((curr->first).c_str()));
        countVect[zi] = occCounter[zi];
        // dupList[zi] = MxArray::from<<std::vector<mwSize> >(posDupList[zi])

        //countVect[zi] = mxCreateNumericMatrix(cidx, 1, mxUINT64_CLASS , mxREAL );
        //std::copy(posDupList[zi].begin(), posDupList[zi].end(), mxGetPr(countVect[zi]));
        //std::copy(posDupList[zi].begin(), posDupList[zi].end(), mxGetCell(plhs[4],zi));
        mxArray *outArray = mxCreateNumericMatrix(cidx, 1, mxUINT64_CLASS , mxREAL );
        outArray = mexplus::MxArray::from(posDupList[zi]);
        mxSetCell(dupList,zi,outArray);
    }
    mexPrintf("Done - return\n");
    return;
}/*mexFunction()*/


