#include "mex.h"
#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <functional>

/*
  This files is part of the RepLAB library. It provides a faster 
  implementation of the algorithm contained in burning.m.
  
  For optimal efficiency, the edges are not checked for redundency. 
  Hence, it is best to provide each edge only once (in particular, the 
  edge "1-2" already implies "2-1" so the latter should not be given).
  If we wish to change the code to eliminate redundancy, the data 
  structure should be changed from vector < vector < long int > > to
  vector < set < long int > >.  
  
  This implementation supports sparse inputs, such as [1 2; 1 n] with n
  very large; with an overhead in memory and time much smaller than n.
*/

using namespace std;


// The data type we use for indices of vertices
typedef unsigned long int Index;


// This class comparator compares the first elements of pairs
struct compareFirstPart {
  bool operator() (const pair < Index, Index >& lhs, const pair < Index, Index >& rhs) const
  {
    return lhs.first < rhs.first;
  }
};


/* This is the function that is called from matlab. It has just one 
   possible calling pattern:
    - prhs should point to a matrix of size n x 2 of undirected edges
    - nlhs will be returned as a cell array of groups of connex vertices

   Remember that the following function is supposed to deal with all the memory allocation
   by itself.
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  //-//-// Argument checking //-//-//
  
  // We check that the parameters are correct
  if (nrhs != 1)
    mexErrMsgTxt("burningAlgorithm_mex: Unexpected number of arguments.");
  if (nlhs != 1)
    mexErrMsgTxt("burningAlgorithm_mex: Unexpected number of outputs.");

  // The matlab object is supposed to be an array
  if (!mxIsDouble(prhs[0]))
    mexErrMsgTxt("burningAlgorithm_mex: The argument should be an array of double.");
  
  // Get the size and pointers to input data
  mwSize m(mxGetM(prhs[0]));
  mwSize n(mxGetN(prhs[0]));
  double* pr(mxGetPr(prhs[0]));
  double* pi(mxGetPi(prhs[0]));
  bool isComplex = (pi==NULL ? 0 : 1);

  // The input should be real
  if (isComplex != 0)
    mexErrMsgTxt("burningAlgorithm_mex: The argument should not be complex.");

  // Second dimension should be 2
  if (n != 2)
    mexErrMsgTxt("burningAlgorithm_mex: The input should be of dimension m x 2.");
  
  
  
  //-//-// Data initialization //-//-//
  
  // First, we quickly list the vertices numbers in a compact vector
  set < Index > uniqueVertices(pr, pr + m*n); // This efficiently removes duplicates
  vector < Index > vertices(uniqueVertices.begin(), uniqueVertices.end()); // We extract the unique vertices into a vector
  Index nbVertices(vertices.size());

  // And create a fast lookup for vertices' numbers from their index in 'vector'
  set < pair < Index, Index >, compareFirstPart > initialIndex;
  for (unsigned int i = 0; i < nbVertices; ++i)
    initialIndex.insert(std::make_pair(vertices[i], i));

  // This function is such that verticesInverse(vertices[i]) gives back i
  function < Index (Index) > verticesInverse = [=](Index i){ return (*initialIndex.find(pair < Index, Index >(i,0))).second; };
  /*for (unsigned int i = 0; i < nbVertices; ++i)
  {
    cout << i << " == " << verticesInverse(vertices[i]) << endl;
  }*/
  
  // We initialize the graph data structure
  vector < vector < Index > > graphData(nbVertices);
  for (mwIndex i = 0; i < m; ++i) {
    Index a(*(pr+i));
    Index b(*(pr+i+m));
    
    // We map these numbers to the new compact indices
    Index newA(verticesInverse(a));
    Index newB(verticesInverse(b));

    // We save the link in both directions
    graphData[newA].push_back(newB);
    graphData[newB].push_back(newA);
  }



  //-//-// Algorithm //-//-//
  
  // Now we perform the actual burning algorithm
  Index reached [nbVertices];
  fill_n(reached, nbVertices, 0);

  vector < Index > neighbors [2];
  short int ptr(0);
  neighbors[ptr].push_back(0); // We begin the algorithm by reaching the first site 0.
  Index nbSets(0);       // The same set number is assigned to each vertices belonging to a connex group

  vector < vector < Index > > allSets(0); // We keep track of which vertex ends up in which set, here with the original numbering

  // Let's "burn" all the sites that touch a reached site recursively until there is no left.
  do {
    ++nbSets;
    allSets.push_back(vector < Index >(0));
    do {
      neighbors[1-ptr].clear();
      for (Index i(0); i < neighbors[ptr].size(); ++i) {
        if (reached[neighbors[ptr][i]] == 0) {
          reached[neighbors[ptr][i]] = nbSets;
          allSets[nbSets-1].push_back(vertices[neighbors[ptr][i]]);
          for (Index j(0); j < graphData[neighbors[ptr][i]].size(); ++j)
            if (reached[graphData[neighbors[ptr][i]][j]] == 0)
              neighbors[1-ptr].push_back(graphData[neighbors[ptr][i]][j]);
        }
      }
      ptr = 1-ptr;
    } while (neighbors[ptr].size() > 0);
    
    // We look for the next un-attained vertex
    for (Index i(0); i < nbVertices; ++i) {
      if (reached[i] == 0) {
        neighbors[ptr].push_back(i);
        break;
      }
    }
  } while (neighbors[ptr].size() > 0);
  


  //-//-// Preparing output fields //-//-//
  plhs[0] = mxCreateCellMatrix(1, nbSets);

  // Now we iterate on all the elements of the cell array
  for (Index i = 0; i < nbSets; ++i) {
    mxArray* oneFamily(mxCreateNumericMatrix(1, allSets[i].size(), mxDOUBLE_CLASS, mxREAL));
    copy(allSets[i].begin(), allSets[i].end(), mxGetPr(oneFamily)); // copy the data
    mxSetCell(plhs[0], i, oneFamily); // Assign the data to the cell element
  }
  
  return;
}
