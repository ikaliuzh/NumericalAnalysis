#pragma once

#include <optional>
#include <cmath>


long double __integrate_trapezoidal(
		long double (*f)(long double),
		long double a,
		long double b,
		size_t n){
	long double h = (b - a) / n;
	long double c = a, d = b;

	while (std::isnan(f(c)) && c < a + h) {
		c += h / 100.;
	}
	while (std::isnan(f(d)) && d > b - h) {
		d -= h / 100.;
	}
	long double result = (f(c) + f(d))/2.;
	// h = (d - c) / n;

	for (size_t i = 1; i < n; ++i) {
		long double xi = a + h * i;
		// while (std::isnan(f(xi)) && xi < a + h * i + h) {
		// 	xi += h / 100.;
		// }
		result += f(xi);
	}
	return result * h;
}

long double __integrate_composite_simpson(
		long double (*f)(long double),
		long double a,
		long double b,
		size_t n){
	long double h = (b - a) / n;
	long double c = a, d = b;

	while (std::isnan(f(c)) && c < a + h) {
		c += h / 100.;
	}
	while (std::isnan(f(d)) && d > b - h) {
		c -= h / 100.;
	}
	// h = (d - c) / n;


	long double sum_odds = 0, sum_evens = 0; 	

	for (size_t i = 1; i < n; ++i) {
		long double xi = a + h * i;
		// while (std::isnan(f(xi)) && xi < c + h * i + h) {
		// 	xi += h / 100.;
		// }
		if (i % 2)
			sum_odds += f(xi);
		else
			sum_evens += f(xi);
	}

	return (f(c) + f(d) + 2 * sum_evens + 4 * sum_odds) * h / 3;
}


std::optional<long double> integrate(
		long double (*f)(long double),
		long double a,
		long double b,
		long double eps,
		char method){
	int p;
	long double (*integrator)(
		long double (*f)(long double),
		long double a,
		long double b,
		size_t n);

	if (method == 's'){
		p = 2;
		integrator = __integrate_composite_simpson;
	} else if (method == 't'){
		p = 1;
		integrator = __integrate_trapezoidal;
	}
	
	long double k = 1. / (pow(2, p) - 1);
	int n = 1;
	while (true) {
		long double In = integrator(f, a, b, n);
		long double I2n = integrator(f, a, b, 2 * n);
		if (k * abs(In - I2n) < eps) {
			std::cout << n << std::endl;
			return I2n;
		}
		if (n > 100'000){
			return {};
		}
		n *= 2;
	}
}
