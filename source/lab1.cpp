#include <vector>
#include <initializer_list>
#include <algorithm>
#include <cmath>
#include <limits>
#include <cassert>
#include <iostream>
#include <iomanip>


template<typename T> T abs(T x){
	if (x < T(0)) return x * T(-1);
	return x;
}

class Pynomial{
public:
	Pynomial()
		: coeffs({0}) {}

	Pynomial(size_t N)
		: coeffs(N){}

	Pynomial(std::initializer_list<double> initl)
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


Pynomial operator+(const Pynomial& lhs, const Pynomial& rhs){
	Pynomial res;
	for (size_t i = 0; i < std::max(lhs.size(), rhs.size()); ++i){
		res[i] = lhs[i] + rhs[i];
	}
	return res;
}

Pynomial operator*(const Pynomial& lhs, const Pynomial& rhs){
	Pynomial res;
	for (size_t i = 0; i < lhs.size() + rhs.size(); ++i){
		res[i] = 0;
		for (size_t j = 0; j <= i; ++j)
			res[i] += lhs[j] * rhs[i-j];
	}
	return res;
}


Pynomial operator*(Pynomial lhs, const double& rhs){
	for (size_t i = 0; i < lhs.size(); ++i){
		lhs[i] *= rhs;
	}
	return lhs;
}

Pynomial operator/(Pynomial lhs, const double& rhs){
	for (size_t i = 0; i < lhs.size(); ++i){
		lhs[i] /= rhs;
	}
	return lhs;
}

std::ostream& operator<<(std::ostream& stream, const Pynomial& p){
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

class LagrangeInterpolator{
public:
	static double f(double x);

	LagrangeInterpolator()
		: P(), a(0), b(1), niters(10), npoints(0) {}

	void setInterval(double A, double B){
		assert(a < b);
		a = A; b = B;
	}
	void fitUniform(size_t N) {
		npoints = N;
		for (size_t i = 0; i < npoints; ++i){
			double xi = a + i * (b - a) / N;
			Pynomial prod{1};
			for (size_t j = 0; j < N; ++j){
				if (j == i)
					continue;
				double xj = a + j * (b - a) / N;
				prod = prod * Pynomial{-xj, 1} / (xi - xj);
			}
			P = P + prod * Pynomial{f(xi)};
		}
	}

	double precision() const{
		srand(time(NULL));
		double eps = 0;
		for (size_t i = 0; i < niters; ++i){
			double x = ((double) rand() / (RAND_MAX)) * (b - a) + a;
			eps = std::max<double>(eps, abs(f(x)-P(x)));
		}
		return eps;
	}

	bool fitPrecision(double eps){
		// not stable
		for (size_t k = 1; k < 100; ++k){
			fitUniform(10 * k);
			if (precision() < eps){
				return true;
			}
		}
		return false;
	}

	Pynomial getPolynom() const{
		return P;
	}

	void report(std::ostream& stream = std::cout) const{
		stream << "x_i\t|\tf(x_i)\t|\tP(x_i)\t|\t|P(x_i) - f(x_i)|\n" << std::string(64, '-') << '\n';
		for (size_t i = 0; i < npoints; ++i){
			double xi = a + i * (b - a) / npoints;
			stream << xi << "\t|\t" << f(xi) << "\t|\t" << P(xi) << "\t|\t" << abs(P(xi) - f(xi)) << '\n';
		}
		stream << std::string(64, '=') << "\n\n\n";

		srand(time(NULL));
		stream << "y_i\t|\tf(y_i)\t|\tP(y_i)\t|\t|P(y_i) - f(y_i)|\n" << std::string(64, '-') << '\n';
		for (size_t i = 0; i < npoints; ++i){
			double alpha = (double) rand() / (RAND_MAX);
			double xi = a + (i + alpha) * (b - a) / npoints;
			stream << std::setprecision(5) << xi << "\t|\t" << f(xi) << "\t|\t" << P(xi) << "\t|\t" << abs(P(xi) - f(xi)) << '\n';
		}
		stream << std::string(64, '=') << "\n";
	}

	virtual ~LagrangeInterpolator() = default;
private:
	Pynomial P;
	double a, b;
	size_t niters, npoints;
};

double LagrangeInterpolator::f(double x){
	return pow(2, x);
}

int main(){
	LagrangeInterpolator lp;
	lp.setInterval(0, 10);
	lp.fitUniform(20);
	lp.report();
	return 0;
}
