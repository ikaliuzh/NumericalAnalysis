#pragma once

#include "Matrix.hpp"
#include <optional>

template <size_t N>
class LinearSystemSolver {
public:
	LinearSystemSolver(SquareMatrix<N> mat, Vector<N> vec)
		: A(mat.transpose() * mat), b(mat.transpose() * vec){}

	// Iterative Solver
	std::optional< Vector<N> > solve(long double eps) {
		Vector<N> solution[2];
		for (size_t i = 0; i < N; ++i){
			solution[0][i][0] = 1;
		}
		long double tau = 0.016;
		int64_t k = 0;
		while (true){
			if (k > 1000) {return {};}
			solution[(k + 1) % 2] = solution[k % 2] + tau * (b - A * solution[k % 2]);
			if ((solution[0] - solution[1]).norm() < eps) {
				std::cout << k << '\n';
				return solution[(k+1)%2];
			}
			++k;
		}
	}
private:
	SquareMatrix<N> A;
	Vector<N> b;
};
