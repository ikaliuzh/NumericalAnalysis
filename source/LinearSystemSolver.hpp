#pragma once

#include "Matrix.hpp"
#include <optional>

template <size_t N>
class LinearSystemSolver {
public:
	LinearSystemSolver(SquareMatrix<N> mat, Vector<N> vec)
		: A(mat), b(vec){}

	// Iterative Solver
	std::optional< Vector<N> > solve(long double eps) {

		SquareMatrix<N> B = A - SquareMatrix<N>().identity();
		if (B.norm() >= 1){
			return {};
		}
		Vector<N> solution[2];
		for (size_t i = 0; i < N; ++i){
			solution[0][i][0] = 1;
		}
		int64_t k = 0;
		while (true){
			if (k > 100) return {};
			solution[(k + 1) % 2] = B * solution[k % 2] + b;
			if ((solution[0] - solution[1]).norm() < eps) {
				return solution[(k+1)%2];
			}
			++k;
		}
	}
private:
	SquareMatrix<N> A;
	Vector<N> b;
};
