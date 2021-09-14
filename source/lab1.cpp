#include <vector>
#include <initializer_list>
#include <sstream>
#include <utility>
#include <memory>
#include <cmath>
#include <limits>
#include <cassert>
#include <iostream>
#include <iomanip>


template<typename T> T abs(T x){
	if (x < T(0)) return x * T(-1);
	return x;
}

class Polynomial{
public:
	Polynomial()
		: coeffs({0}) {}

	Polynomial(size_t N)
		: coeffs(N){}

	Polynomial(std::initializer_list<double> initl)
		: coeffs(initl) {}

	size_t size() const{
		if (coeffs.size() > 0 && coeffs[coeffs.size() - 1] == 0){
			size_t size = 0;
			for (size_t i = 0; i < coeffs.size(); ++i){
				if (coeffs[i] != 0){
					size = i;
				}
			}
			return size + 1;
		}
		return coeffs.size();
	}

	void shrink(){
		coeffs.resize(size());
	}

	double operator[](size_t i) const{
		if (i < coeffs.size())
			return coeffs[i];
		return 0;
	}

	double& operator[](size_t i){
		if (i < coeffs.size())
			return coeffs[i];
		coeffs.resize(i + 1);
		return coeffs[i];
	}

	double operator()(double x) const{
		double res = 0;
		for (size_t i = 0; i < size(); ++i){
			res += pow(x, i) * coeffs[i];
		}
		return res;
	}


private:
	std::vector<double> coeffs;
};


Polynomial operator+(const Polynomial& lhs, const Polynomial& rhs){
	Polynomial res;
	for (size_t i = 0; i < std::max(lhs.size(), rhs.size()); ++i){
		res[i] = lhs[i] + rhs[i];
	}
	return res;
}

Polynomial operator*(const Polynomial& lhs, const Polynomial& rhs){
	Polynomial res;
	for (size_t i = 0; i < lhs.size() + rhs.size(); ++i){
		res[i] = 0;
		for (size_t j = 0; j <= i; ++j)
			res[i] += lhs[j] * rhs[i-j];
	}
	return res;
}


Polynomial operator*(Polynomial lhs, const double& rhs){
	for (size_t i = 0; i < lhs.size(); ++i){
		lhs[i] *= rhs;
	}
	return lhs;
}

Polynomial operator/(Polynomial lhs, const double& rhs){
	for (size_t i = 0; i < lhs.size(); ++i){
		lhs[i] /= rhs;
	}
	return lhs;
}

std::ostream& operator<<(std::ostream& stream, const Polynomial& p){
	for (size_t i = 0; i < p.size(); ++i){
		int sign = 1;
		if (i != 0 && p[i] > 0)
			stream << " + ";
		if (i != 0 && p[i] <= 0){
			stream << " - ";
			sign = -1;
		}
		stream << p[i] * sign;
		if (i == 1)
			stream << "x";
		else if (i > 1)
			stream << "x^" << i;
	}
	return stream;
}


namespace Nodes{
	class Root{
	public:
		Root(size_t N, double A, double B)
			: n(N), a(A), b(B){};
		Root(std::initializer_list<double> ilist){
			size_t i = 0;
			for (auto v : ilist){
				if (i == 0){
					a = v;
				}
				b = v;
				n = i;
				++i;
			}
			++n;
		}
		virtual ~Root() = default;
		virtual double operator[](size_t) const =0;
		double at(size_t i) const {return this->operator[](i);};
		void transform(size_t m, double c, double d) {
			n = m; a = c; b = d;
		}

		size_t size() const{
			return n;
		}

		std::pair<double, double> interval() const{
			return {a, b};
		}
	protected:
		size_t n;
		double a, b;
	};


	class Uniform : public Nodes::Root{
	public:
		Uniform(size_t N, double A, double B)
			: Root(N, A, B) {}

		double operator[](size_t i) const{
			return a + i * (b - a) / n;
		}
	};

	class Ordinary : public Nodes::Root{
	public:
		Ordinary(std::initializer_list<double> ilist)
			:Root(ilist),
			 vals(ilist)
			 {}

		double operator[](size_t i) const{
				return vals[i];
			}
	private:
		std::vector<double> vals;
	};
}


class LagrangeInterpolator{
public:
	static double f(double x);

	LagrangeInterpolator()
		: P(), nodes() {}

	void setNodes(Nodes::Root* v){
		delete nodes;
		nodes = v;
	}

	void fitUniform() {
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

	double precision() const{
		srand(time(NULL));
		double eps = 0;
		for (size_t i = 0; i < nodes->size(); ++i){
			double x = ((double) rand() / (RAND_MAX)) *
					(nodes->interval().second - nodes->interval().first)
					+ nodes->interval().first;
			eps = std::max<double>(eps, abs(f(x)-P(x)));
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
			double alpha = (double) rand() / (RAND_MAX);
			double xi = nodes->at(i) + alpha *
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

double LagrangeInterpolator::f(double x){
	return exp(-x*x) * sin(x) + x*x;
}

int main(){
	LagrangeInterpolator lp;
	Nodes::Uniform nodes{20, -1, 1};
	// Nodes::Ordinary nodes{0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
	lp.setNodes(&nodes);
	lp.fitUniform();
	lp.report();
	return 0;
}
