#pragma once

#include <iostream>
#include <iomanip>

#include "Nodes.hpp"
#include "Polynoms.hpp"

template<typename T> T abs(T x){
	if (x < T(0)) return x * T(-1);
	return x;
}


class LagrangeInterpolator{
public:
	static long double f(long double x);

	LagrangeInterpolator()
		: P(), nodes() {}

	void setNodes(Nodes::Root* v){
		delete nodes;
		nodes = v;
	}

	void fit() {
		for (size_t i = 0; i < nodes->size(); ++i){
			Polynomial prod{1};
			for (size_t j = 0; j < nodes->size(); ++j){
				if (j == i)
					continue;
				prod = prod * Polynomial{-nodes->at(j), 1} / (nodes->at(i) - nodes->at(j));
			}
			P = P + prod * Polynomial{f(nodes->at(i))};
		}
	}

	void addNode(long double x){
		long double a = f(x) - P(x);
		Polynomial Ptmp{1};
		for (size_t i = 0; i < nodes->size(); ++i){
			a /= (x - nodes->at(i));
			Ptmp = Ptmp * Polynomial{-nodes->at(i), 1};
		}
		P = P + a * Ptmp;

		if (nodes->size() > 0){
			Nodes::Root* ndsnew = nodes->addNode(x);
			if (nodes != ndsnew){
				// delete old node; memory leak
				nodes = ndsnew;
			}
		} else {
			nodes = new Nodes::Ordinary({x});
		}
	}

	long double precision() const{
		srand(time(NULL));
		long double eps = 0;
		for (size_t i = 0; i < nodes->size(); ++i){
			long double x = ((long double) rand() / (RAND_MAX)) *
					(nodes->interval().second - nodes->interval().first)
					+ nodes->interval().first;
			eps = std::max<long double>(eps, abs(f(x)-P(x)));
		}
		return eps;
	}

	Polynomial getPolynom() const{
		return P;
	}

	void report(std::ostream& stream = std::cout) const{
		stream << "\t x_i\t|\tf(x_i)\t|\tP(x_i)\t|\t|P(x_i) - f(x_i)|\n"
				<< std::string(75, '-') << '\n';
		for (size_t i = 0; i < nodes->size(); ++i){
			stream << '\t' << std::setprecision(3) << nodes->at(i) <<
					"\t|\t" << f(nodes->at(i)) << "\t|\t" <<
					P(nodes->at(i)) << "\t|\t" <<
					abs(P(nodes->at(i)) - f(nodes->at(i))) << '\n';
		}
		stream << std::string(75, '=') << "\n\n\n";

		srand(time(NULL));
		stream << "\t y_i\t|\t f(y_i)\t|\tP(y_i)\t|\t|P(y_i) - f(y_i)|\n" <<
				std::string(75, '-') << '\n';
		for (size_t i = 0; (i + 1) < nodes->size(); ++i){
			// long double alpha = (long double) rand() / (RAND_MAX);
			long double alpha = 0.5;
			long double xi = nodes->at(i) + alpha *
					(nodes->at(i + 1) - nodes->at(i));
			stream << '\t' << std::setprecision(3) << xi << "\t|\t" << f(xi)
					<< "\t|\t" << P(xi) << "\t|\t" <<
					abs(P(xi) - f(xi)) << '\n';
		}
		stream << std::string(75, '=') << "\n";
	}

	virtual ~LagrangeInterpolator() = default;
private:
	Polynomial P;
	Nodes::Root *nodes;
};

long double LagrangeInterpolator::f(long double x){
	return exp(-x*x) * sin(x) + x*x;
}
