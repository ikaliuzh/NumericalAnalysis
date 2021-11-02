#include "PolynomialInterpolator.hpp"
#include "CubicSplines.hpp"
#include "LinearSystemSolver.hpp"

int main(int argc, char *argv[]){
	/* Nodes::Uniform nodes{20, -4, 4};
	CubicSpline cs(&nodes);
	cs.report(); */

	SquareMatrix<2> A{{1, 0}, {0, 0.5}};
	Vector<2> b{{1}, {2}};

	LinearSystemSolver lss(A, b);
	std::cout << *lss.solve(0.001);
	return 0;
}

