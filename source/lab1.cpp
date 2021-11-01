#include "PolynomialInterpolator.hpp"
#include "CubicSplines.hpp"

int main(){
	Nodes::Uniform nodes{20, -4, 4};
	CubicSpline cs(&nodes);
	cs.report();
	return 0;
}

