#include "PolynomialInterpolator.hpp"


int main(){
	LagrangeInterpolator lp2;
	Nodes::Chebyshev nodes2{30, -4, 4};
	lp2.setNodes(&nodes2);
	lp2.fit();
	lp2.report();
	return 0;
}

