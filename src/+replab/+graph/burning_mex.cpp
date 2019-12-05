#include "mex.h"
#include <iostream>
#include <vector>
#include <set>
#include <algorithm>

/*
  This files is part of the RepLAB library. It provides a faster 
  implementation of the algorithm contained in burning.m.
*/

using namespace std;


/*// Embed an array into a c++ vector without copying
// From https://stackoverflow.com/questions/7278347/c-pointer-array-to-vector
template <class T>
class vectorWrapper : public std::vector<T>
{   
public:
  vectorWrapper() {
    this->_M_impl _M_start = this->_M_impl _M_finish = this->_M_impl _M_end_of_storage = NULL;
  }   

  vectorWrapper(T* sourceArray, int arraySize)
  {   
    this->_M_impl _M_start = sourceArray;
    this->_M_impl _M_finish = this->_M_impl _M_end_of_storage = sourceArray + arraySize;
  }   

  ~vectorWrapper() {
    this->_M_impl _M_start = this->_M_impl _M_finish = this->_M_impl _M_end_of_storage = NULL;
  }   

  void wrapArray(T* sourceArray, int arraySize)
  {   
    this->_M_impl _M_start = sourceArray;
    this->_M_impl _M_finish = this->_M_impl _M_end_of_storage = sourceArray + arraySize;
  }   
};*/


// The following function is supposed to deal with all the memory allocation
// by itself.
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  // We check that the parameters are correct
  if (nrhs != 1)
    mexErrMsgTxt("burning_mex: Unexpected number of arguments.");
  if (nlhs != 1)
    mexErrMsgTxt("burning_mex: Unexpected number of outputs.");

  // The matlab object is supposed to be an array
  if (!mxIsDouble(prhs[0]))
    mexErrMsgTxt("burning_mex: The argument should be an array of double.");
  
  // We extract the input data
  /* Get the size and pointers to input data */
  mwSize m(mxGetM(prhs[0]));
  mwSize n(mxGetN(prhs[0]));
  double* pr(mxGetPr(prhs[0]));
  double* pi(mxGetPi(prhs[0])); // Not used since we don't care about imaginary part
  bool isComplex = (pi==NULL ? 0 : 1);

  // The input should be real
  if (isComplex != 0)
    mexErrMsgTxt("burning_mex: The argument should not be complex.");

  // Second dimension should be 2
  if (n != 2)
    mexErrMsgTxt("burning_mex: The input should be of dimension m x 2.");
  
  
  // We initialize the graph data structure
  vector < set < long int > > graphData (0);
  
  // We first iterate on the rows
  for (mwIndex i = 0; i < m; ++i) {
    long int a(*pr - 1);
    long int b(*(pr+m) - 1);
    long int maxab(max(a,b));
    
    // If required, we make the data structure big enough to contain the elements we need
    if (graphData.size() <= maxab)
      graphData.resize(maxab+1, set < long int > ());
    
    graphData[a].insert(b);
    graphData[b].insert(a);
  }


  // TODO: Perform the actual burning algorithm


  // We allocate space for the result
  plhs[0] = mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL);

  // We check where the output data should be places
  int* outputMatrix = (int*)mxGetData(plhs[0]);

  // And stock it at the right place
  //outputMatrix[0] = mpfr::mpreal::get_default_prec();

  return;
}
