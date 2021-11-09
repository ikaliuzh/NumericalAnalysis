#include "PolynomialInterpolator.hpp"
#include "CubicSplines.hpp"
#include "LinearSystemSolver.hpp"

int main(int argc, char *argv[]){
	Nodes::Uniform nodes{40, -4, 4};
	CubicSpline cs(&nodes);
	cs.report();

//	SquareMatrix<4> A{{2, 4, -1, -2}, {6, 1, 3, 4}, {-7, 1, 1, -1}, {0, 2, -3, 5}};
//	Vector<4> b{{1}, {1}, {3}, {5}};
//
//	LinearSystemSolver lss(A, b);
//	auto solution = lss.solve(1e-7);
//	if (solution.has_value()) {
//		std::cout << "Решение:\n" << solution.value() << "Невязка:\n" << b - A * solution.value();
//	} else {
//		std::cout << "Failed to solve!";
//	}
//	return 0;
}
