#include "PolynomialInterpolator.hpp"
#include "CubicSplines.hpp"
#include "Matrix.hpp"

int main(int argc, char *argv[]){
	/* Nodes::Uniform nodes{20, -4, 4};
	CubicSpline cs(&nodes);
	cs.report(); */
	Matrix<double, 2, 2> A{{1, 2}, {3, 4}}, B{{1, 2}, {0 , 0}};
	std::cout << A;
	return 0;
}

