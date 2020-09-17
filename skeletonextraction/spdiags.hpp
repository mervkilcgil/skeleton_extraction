#pragma once
#include <vector>
#include "include\Eigen\Sparse"
#include <algorithm>
using namespace Eigen;
//#pragma begin<0>
//! similar to matlab function A = spdiags(B,d,m,n)
//! spdiags(B,d,m,n) creates an m-by-n sparse matrix by taking the 
//! columns of B and placing them along the diagonals specified by d.
template <class numeric_t> 
SparseMatrix<numeric_t> spdiags(const Matrix<numeric_t,-1,-1> &B, 
					const VectorXi &d, const int m, const int n) {					
	typedef Triplet<numeric_t> triplet_t;
	std::vector<triplet_t> triplets;
	triplets.reserve(std::min(m,n)*d.size());
	for (int k = 0; k < d.size(); ++k) {
		int diag = d(k);	// get diagonal
		int i_start = std::max(-diag, 0); // get row of 1st element
		int i_end = std::min(m, m-diag-(m-n)); // get row of last element
		int j = -std::min(0, -diag); // get col of 1st element
		int B_i; // start index i in matrix B
		if(m < n)
			B_i = std::max(-diag,0); // m < n
		else
			B_i = std::max(0,diag); // m >= n
		for(int i = i_start; i < i_end; ++i, ++j, ++B_i){
			triplets.push_back( {i, j,  B(B_i,k)} );
		}
	}
	SparseMatrix<numeric_t> A(m,n);
	A.setFromTriplets(triplets.begin(), triplets.end());
	return A;
}
//#pragma end<0>
