#pragma once

#include <optional>
#include <cmath>


long double __integrate(
		long double (*f)(long double),
		long double a,
		long double b,
		size_t n){
	long double h = (b - a) / n;
	long double result = 0;
	long double c = a, d = b;

	while (std::isnan(f(c)) && c < a + h) {
		c += h / 100.;
	}
	while (std::isnan(f(d)) && d > b - h) {
		c -= h / 100.;
	}
	h = (d - c) / n;

	for (size_t i = 1; i < n; ++i) {
		long double xi = c + h * i;
		while (std::isnan(f(xi)) && xi < c + h * i + h) {
			xi += h / 20.;
		}
		result += f(xi);
	}
	return result * h;
}


std::optional<long double> integrate(
		long double (*f)(long double),
		long double a,
		long double b,
		long double eps){

	int p = 2;
	long double k = 1. / (pow(2, p) - 1);
	int n = 1;
	while (true) {
		long double In = __integrate(f, a, b, n);
		long double I2n = __integrate(f, a, b, 2 * n);
		if (k * abs(In - I2n) < eps) {
			return I2n;
		}
		if (n > 100'000){
			return {};
		}
		n += 200;
	}
}
