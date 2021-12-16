#include "PolynomialInterpolator.hpp"
#include "CubicSplines.hpp"
#include "LinearSystemSolver.hpp"
#include "NumericIntegration.hpp"

long double g (long double x) {
	return (x - sin(x)) / x;
}

int main(int argc, char *argv[]){
	std::cout << integrate(g, 0, pi/2, 1e-5, 't').value() << std::endl;
}

