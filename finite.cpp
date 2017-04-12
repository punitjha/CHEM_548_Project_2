#include <armadillo>
#include <iostream>

int main() 
{
	int N = 4;
	arma::vec diags = 1.2 * arma::ones(N);
	arma::vec offdiags = 2.3 * arma::ones(N-1);
	arma::mat M(N, N, arma::fill::zeros);
	M.diag() = diags;
	M.diag(-1) = offdiags;
	M.diag(1) = offdiags;
	M.print("M = ");

	arma::vec eigvals; 
	arma::mat eigvecs;
	arma::eig_sym(eigvals, eigvecs, M);
	eigvals.print("eigvals = ");
	eigvecs.print("eigvecs = ");

	arma::uvec idx = arma::sort_index(eigvals);
	eigvals = eigvals(idx);
	eigvecs = eigvecs(idx);

	eigvals.print("sorted vals");
	eigvecs.print("sorted vecs");

 return 0;
}

 // if undefined symbols at compilation:
 // g++ prog.cpp -o prog -O2 -larmadillo -llapack -lblas
 // normally just
// g++ -larmadillo *.cpp 
