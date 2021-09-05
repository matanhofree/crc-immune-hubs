/*******************************************************************************
 ** Fast C++ based unique function.
 **
 ** Input:
 ** inputVector A
 **
 ** Output: [ keyVect, orderA, orderKey, countVect, dupListIdx, dupVect ] = fastUniqueIntFlat
 **  keyVect - a unique key vector in order of occurrence
 **  orderA -  keyVect = A[orderA]
 **  orederVect - A = keyVect[orderKey]
 **  countVect -= count occurrences of keyVect in A
 **
 **  Alternative for improved perormance (?):
 **  dupListIdx -- (length of keyVect) Index in dupVec of stat occuranc
 **  dupVect -- occurances of keyVect in A
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

    // uint64_T *keyVectInput = (uint64_T *) mxGetData(prhs[0]);
    mwSize lenKey = mxGetNumberOfElements(prhs[0]);
    std::vector<uint64_T > keyVect;
    mexplus::MxArray::to<std::vector<uint64_T> >(prhs[0], &keyVect);

    mwSize lenA = mxGetNumberOfElements(prhs[1]);
    std::vector<uint64_T > A;
    mexplus::MxArray::to<std::vector<uint64_T> >(prhs[1], &A);

    std::unordered_map<uint64_T, mwSize,std::hash<uint64_T>,std::equal_to<uint64_T>,Mallocator<std::pair<const uint64_T,mwSize> > > keyMap;
    mwSize cidx = 1;
    for (auto const & x : keyVect) { keyMap[x] = cidx++; }


    plhs[0] = mxCreateDoubleMatrix(lenA,1,mxREAL);
    double* idxVector = (double *)mxGetPr(plhs[0]);
    double nanVal = mxGetNaN();

    for (mwSize i = 0; i < lenA;++i)
    {
        std::unordered_map<uint64_T,mwSize>::const_iterator curr = keyMap.find(A[i]);

        if (curr==keyMap.end())
        {
            idxVector[i] = nanVal;
        }
        else
        {
            idxVector[i] = curr->second;
        }

    }

    return;


}/*mexFunction()*/


