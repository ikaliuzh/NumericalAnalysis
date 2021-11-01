#pragma once

#include <vector>

class Polynomial{
public:
	Polynomial()
		: coeffs({0}) {}

	Polynomial(size_t N)
		: coeffs(N){}

	Polynomial(std::initializer_list<long double> initl)
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

	long double operator[](size_t i) const{
		if (i < coeffs.size())
			return coeffs[i];
		return 0;
	}

	long double& operator[](size_t i){
		if (i < coeffs.size())
			return coeffs[i];
		coeffs.resize(i + 1);
		return coeffs[i];
	}

	long double operator()(long double x) const{
		long double res = 0;
		for (size_t i = 0; i < size(); ++i){
			res += pow(x, i) * coeffs[i];
		}
		return res;
	}


private:
	std::vector<long double> coeffs;
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


Polynomial operator*(const long double& lhs, Polynomial rhs){
	for (size_t i = 0; i < rhs.size(); ++i){
		rhs[i] *= lhs;
	}
	return rhs;
}

Polynomial operator*(Polynomial lhs, const long double& rhs){
	for (size_t i = 0; i < lhs.size(); ++i){
		lhs[i] *= rhs;
	}
	return lhs;
}

Polynomial operator/(Polynomial lhs, const long double& rhs){
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
