#pragma once

#include "Nodes.hpp"


class CubicSpline {
public:
	static long double f(long double x);

	CubicSpline(Nodes::Root *X)
		: x(X->size()),
		  a(X->size()),
		  b(X->size()),
		  c(X->size()),
		  d(X->size()) {
		for (size_t i = 0; i < x.size(); ++i){
			x[i] = X->at(i);
			a[i] = f(x[i]);
		}

		int n = x.size()-1;
	    std::vector<long double> h;

	    for(int i = 0; i < n; ++i)
	        h.push_back(x[i+1]-x[i]);

	    std::vector<long double> alpha;
	    alpha.push_back(0);
	    for(int i = 1; i < n; ++i)
	        alpha.push_back( 3*(a[i+1]-a[i])/h[i] - 3*(a[i]-a[i-1])/h[i-1]  );

	    std::vector<long double> c(n+1), l(n+1), mu(n+1), z(n+1);
	    l[0] = 1; mu[0] = 0; z[0] = 0;

	    for(int i = 1; i < n; ++i)
	    {
	        l[i] = 2 *(x[i+1]-x[i-1])-h[i-1]*mu[i-1];
	        mu[i] = h[i]/l[i];
	        z[i] = (alpha[i]-h[i-1]*z[i-1])/l[i];
	    }

	    l[n] = 1;
	    z[n] = 0;
	    c[n] = 0;

	    for(int j = n-1; j >= 0; --j)
	    {
	        c[j] = z [j] - mu[j] * c[j+1];
	        b[j] = (a[j+1]-a[j])/h[j]-h[j]*(c[j+1]+2*c[j])/3;
	        d[j] = (c[j+1]-c[j])/3/h[j];
	    }
	}

	long double operator()(long double t) const{
		size_t j = std::lower_bound(x.begin(), x.end(), t) - x.begin();
		return a[j] + b[j]*(t - x[j])
				+ c[j]*(t - x[j])*(t - x[j])
				+ d[j]*(t - x[j])*(t - x[j])*(t - x[j]);
	}

	void report(std::ostream& stream = std::cout) const{
		stream << "\t x_i\t|\tf(x_i)\t|\t S(x_i)\t|\t|S(x_i) - f(x_i)|\n"
				<< std::string(75, '-') << '\n';
		for (size_t i = 0; i < x.size(); ++i){
			stream << '\t' << std::setprecision(3) << x[i] <<
					"\t|\t" << f(x[i]) << "\t|\t" <<
					operator()(x[i]) << "\t|\t" <<
					abs(operator()(x[i]) - f(x[i])) << '\n';
		}
		stream << std::string(75, '=') << "\n\n\n";

		srand(time(NULL));
		stream << "\t y_i\t|\t f(y_i)\t|\t S(y_i)\t|\t|S(y_i) - f(y_i)|\n" <<
				std::string(75, '-') << '\n';
		for (size_t i = 0; (i + 1) < x.size(); ++i){
			long double alpha = 0.5;
			long double xi = x[i] + alpha *
					(x.at(i + 1) - x[i]);
			stream << '\t' << std::setprecision(3) << xi << "\t|\t" << f(xi)
					<< "\t|\t" << operator()(xi) << "\t|\t" <<
					abs(operator()(xi) - f(xi)) << '\n';
		}
		stream << std::string(75, '=') << "\n";
	}
private:
	std::vector<long double> x, a, b, c, d;
};


long double CubicSpline::f(long double x){
	return exp(-x*x) * sin(x) + x*x;
}
