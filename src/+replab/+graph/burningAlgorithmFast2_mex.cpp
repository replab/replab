#include "mex.h"
#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <functional>

#undef DEBUG
#ifdef DEBUG
  #include <chrono>
#endif

/*
  This files is part of the RepLAB library. It provides a faster
  implementation of the algorithm contained in burningAlgorithm.m.

  For optimal efficiency, the edges are not checked for redundency.
  Hence, it is best to provide each edge only once (in particular, the
  edge "1-2" already implies "2-1" so the latter should not be given).
  If we wish to change the code to eliminate redundancy, the data
  structure could be changed from vector < vector < long int > > to
  vector < set < long int > >

  This implementation is optimized for dense inputs: inputs like
  [1 2; 2 n] will trigger the creation of n vertices.

  Moreover, this implementation directly encodes the result into a
  matlab array. This is possible only for the list of orbits. The
  corresponding cell array still needs to be copied at the end (matlab
  does not support dynamic arrays).
*/

using namespace std;


// The data type we use for indices of vertices
typedef uint64_t Index;


/* This is the function that is called from matlab. It has just one
   possible calling pattern:
    - prhs should point to:
      - a double scalar containing the number of vertices
      - a double matrix of size n x 2 of undirected edges
    - nlhs will be returned as:
      - a vector of double defining the orbit that each vertex belongs to
      - a cell array of groups of connex vertices in double (optional)
   Note that the second output is optional.

   Remember that the following function is supposed to deal with all the memory allocation
   by itself.
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  //-//-// Argument checking //-//-//

  // We check that the parameters are correct
  if (nrhs != 2)
    mexErrMsgTxt("burningAlgorithmFast2_mex: Unexpected number of arguments.");
  if ((nlhs != 1) && (nlhs != 2))
    mexErrMsgTxt("burningAlgorithmFast2_mex: Unexpected number of outputs.");

  // First, we get the number of vertices
  // It should be a scalar
  if ((mxGetM(prhs[0]) != 1) || (mxGetN(prhs[0]) != 1))
    mexErrMsgTxt("burningAlgorithmFast2_mex: Number of vertices should be a scalar.");
  double* pr0(mxGetPr(prhs[0]));
  double* pi0(mxGetPi(prhs[0]));
  bool isComplex0 = (pi0==NULL ? 0 : 1);

  // The input should be real
  if (isComplex0 != 0)
    mexErrMsgTxt("burningAlgorithmFast2_mex: The second argument should not be complex.");

  // The number of vertices
  Index nbVertices(*pr0);


  // Second, the list of edges
  // The matlab object is supposed to be an array
  if (!mxIsDouble(prhs[1]))
    mexErrMsgTxt("burningAlgorithmFast2_mex: The argument should be an array of double.");

  // Get the size and pointers to input data
  mwSize m(mxGetM(prhs[1]));
  mwSize n(mxGetN(prhs[1]));
  double* pr1(mxGetPr(prhs[1]));
  double* pi1(mxGetPi(prhs[1]));
  bool isComplex = (pi1==NULL ? 0 : 1);

  // The input should be real
  if (isComplex != 0)
    mexErrMsgTxt("burningAlgorithmFast2_mex: The argument should not be complex.");

  // Second dimension should be 2
  if (n != 2)
    mexErrMsgTxt("burningAlgorithmFast2_mex: The number of input should be of dimension m x 2.");



  //-//-// Data initialization //-//-//
#ifdef DEBUG
  auto t0 = std::chrono::system_clock::now();

  cout << "Number of vertices : " << nbVertices << endl << flush;
  cout << "Number of edges : " << m << endl << flush;
#endif

  // We initialize the graph data structure
  vector < vector < Index > > graphData(nbVertices);
  for (mwIndex i = 0; i < m; ++i) {
    Index a(*(pr1+i));
    Index b(*(pr1+i+m));

    // We map these numbers to the new compact indices
    Index newA(a-1);
    Index newB(b-1);

    // We save the link in both directions
    graphData[newA].push_back(newB);
    graphData[newB].push_back(newA);
  }

#ifdef DEBUG
  auto t1 = std::chrono::system_clock::now();
  std::chrono::duration<double> delta01 = t1 - t0;
  cout << "Initialization finished (" << delta01.count() << " s)" << endl << flush;
#endif

  //-//-// Algorithm //-//-//

  // Now we perform the actual burning algorithm
  plhs[0] = mxCreateNumericMatrix(1, nbVertices, mxDOUBLE_CLASS, mxREAL); // We directly save this info in matlab format
  double* reached(mxGetPr(plhs[0]));
//  vector < Index > reached(nbVertices, 0);
  vector < Index > neighbors [2];
  short int ptr(0);
  Index lastStart(0);    // We monitor the last starting point and begin the algorithm by reaching the first site 0.
  neighbors[ptr].push_back(lastStart);
  Index nbSets(0);       // The same set number is assigned to each vertices belonging to a connex group

  vector < vector < Index > > allSets(0); // We keep track of which vertex ends up in which set, with a numbering of vertices starting at 1

#ifdef DEBUG
  // For debugging purpose
  Index nbTouchedVertices(0);
  Index lastPercentage(0);
#endif

  // Let's "burn" all the sites that touch a reached site recursively until there none is left.
  do {
    ++nbSets;
    allSets.push_back(vector < Index >(0));
    do {
      neighbors[1-ptr].clear();
      for (Index i(0); i < neighbors[ptr].size(); ++i) {
        if (reached[neighbors[ptr][i]] == 0) {
          reached[neighbors[ptr][i]] = nbSets;
          allSets[nbSets-1].push_back(1+neighbors[ptr][i]);
          for (Index j(0); j < graphData[neighbors[ptr][i]].size(); ++j)
            if (reached[graphData[neighbors[ptr][i]][j]] == 0)
              neighbors[1-ptr].push_back(graphData[neighbors[ptr][i]][j]);
        }
      }
      ptr = 1-ptr;
    } while (neighbors[ptr].size() > 0);

#ifdef DEBUG
    // Update on advancement, for debugging purpose
    nbTouchedVertices += allSets[nbSets-1].size();
    if (nbTouchedVertices*100/nbVertices > lastPercentage) {
      lastPercentage = nbTouchedVertices*100/nbVertices;
      cout << lastPercentage << "% : " << nbTouchedVertices << "/" << nbVertices << endl << flush;
    }
#endif

    // We look for the next un-attained vertex
    for (Index i(lastStart+1); i < nbVertices; ++i) {
      if (reached[i] == 0) {
        lastStart = i;
        neighbors[ptr].push_back(i);
        break;
      }
    }
  } while (neighbors[ptr].size() > 0);

#ifdef DEBUG
  auto t2 = std::chrono::system_clock::now();
  std::chrono::duration<double> delta12 = t2 - t1;
  cout << "Orbits identified (" << delta12.count() << " s)" << endl << flush;
#endif

  if (nlhs == 2)
  {
    //-//-// Preparing additional output field //-//-//
    plhs[1] = mxCreateCellMatrix(1, nbSets);

    // Now we iterate on all the elements of the cell array
    for (Index i = 0; i < nbSets; ++i) {
      mxArray* oneFamily(mxCreateNumericMatrix(1, allSets[i].size(), mxDOUBLE_CLASS, mxREAL));
      copy(allSets[i].begin(), allSets[i].end(), mxGetPr(oneFamily)); // copy the data
      mxSetCell(plhs[1], i, oneFamily); // Assign the data to the cell element
    }
  }

  return;
}
