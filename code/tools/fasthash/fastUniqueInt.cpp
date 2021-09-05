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
    std::unordered_map<uint64_T, mwSize,std::hash<uint64_T>,std::equal_to<uint64_T>,Mallocator<std::pair<const uint64_T,mwSize> > > keyMap;
    //std::unordered_map<uint64_T, mwSize> keyMap;
    mwSize cidx = 0;
    std::vector<uint64_T, Mallocator<uint64_t> > occCounter;
    std::vector<std::vector<mwSize, Mallocator<mwSize> >, Mallocator<std::vector<mwSize, Mallocator<mwSize> > > > posDupList;
    std::vector<mwSize, Mallocator<mwSize> > orderA;
    std::vector<mwSize, Mallocator<mwSize> > orderVect = std::vector<mwSize, Mallocator<mwSize> >(lenA,0);
    uint64_T keyVal;

    mexPrintf("Start index %d %d\n",lenA,sizeof(mwSize));
    for (mwSize i = 0; i < lenA;++i)
    {
        // auto curr = keyMap.find(keyVal);
        keyVal = A[i];
        if (i%100 == 0)
            mexPrintf("Line %d -- %u\n",i,keyVal);

        std::unordered_map<uint64_T,mwSize>::const_iterator curr = keyMap.find(keyVal);

        if (curr==keyMap.end())
        {
            keyMap[keyVal] = cidx;
            occCounter.push_back(1);
            posDupList.push_back(std::vector<mwSize, Mallocator<mwSize> > (1,i+1));
            //mexPrintf("After new index\n");

            orderA.push_back(i+1);
            orderVect[i] = cidx+1;
            ++cidx;
            //mexPrintf("After creating key %d\t%u\n",cidx,keyMap.max_size() );
        }
        else
        {
            //mexPrintf("Before old index\n");

            mwSize cPos = curr->second;
            occCounter[cPos]++;
            posDupList[cPos].push_back(i+1);
            orderVect[i] = cPos+1;
        }
        //mexPrintf("Line %d\n",i);
    }


    mexPrintf("End index\n");

    plhs[0] = mxCreateNumericMatrix(cidx, 1, mxUINT64_CLASS, mxREAL );
    //plhs[1] = mxCreateNumericMatrix(cidx, 1, mxUINT64_CLASS, mxREAL );
    //plhs[2] = mxCreateNumericMatrix(lenA, 1, mxUINT64_CLASS, mxREAL );
    plhs[3] = mxCreateNumericMatrix(cidx, 1, mxUINT64_CLASS, mxREAL );
    plhs[4] = mxCreateCellMatrix(cidx, 1);


    uint64_T* keyVect = (uint64_T*)mxGetData(plhs[0]);
    //uint64_T* orderA = (uint64_T*)mxGetData(plhs[1]);
    //uint64_T* orderKey = (uint64_T*)mxGetData(plhs[2]);
    uint64_T* countVect = (uint64_T*)mxGetData(plhs[3]);
    mxArray* dupList = plhs[4];
    //mexPrintf("Start output \n");

    plhs[1] = mexplus::MxArray::from(orderA);
    plhs[2] = mexplus::MxArray::from(orderVect);

    mwSize zi = 0;
    uint64_T pIdx = 0;

    for (auto curr=keyMap.begin(); curr!=keyMap.end(); ++curr) {
        zi = curr->second;
        keyVect[zi] = curr->first;
        countVect[zi] = occCounter[zi];
        // dupList[zi] = MxArray::from<<std::vector<mwSize> >(posDupList[zi])

        //countVect[zi] = mxCreateNumericMatrix(cidx, 1, mxUINT64_CLASS , mxREAL );
        //std::copy(posDupList[zi].begin(), posDupList[zi].end(), mxGetPr(countVect[zi]));
        //std::copy(posDupList[zi].begin(), posDupList[zi].end(), mxGetCell(plhs[4],zi));
        mxArray *outArray = mxCreateNumericMatrix(cidx, 1, mxUINT64_CLASS , mxREAL );
        outArray = mexplus::MxArray::from(posDupList[zi]);
        mxSetCell(dupList,zi,outArray);
    }

    return;
}/*mexFunction()*/


