#pragma once

#include <vector>
#include <iostream>
#include <cassert>
#include <utility>

template <typename T, size_t M, size_t N>
class Matrix {
public:
	Matrix() : data(M, std::vector<T>(N)) {};
	Matrix(std::initializer_list<std::initializer_list<T>> _data) {
		assert(M == _data.size());
		size_t i = 0;
		for (const auto& v : _data){
			if (i == 0){
				assert(N == v.size());
			}
			data[i] = std::vector<T>(v);
			++i;
		}
	};

	std::vector<T>& operator[](size_t i) { return data[i]; }
	const std::vector<T>& operator[](size_t i) const { return data[i]; }

	std::pair<size_t, size_t> size() const { return {data.size(), data[0].size()}; }
private:
	std::vector<std::vector<T>> data;
};

template <typename T, size_t M, size_t N>
std::ostream& operator<<(std::ostream& stream, const Matrix<T, M, N>& matrix) {
	for (size_t i = 0; i < M; ++i){
		for (size_t j = 0; j < N; ++j){
			stream << matrix[i][j] << ' ';
		}
		stream << '\n';
	}
	return stream;
}

template <typename T, size_t M, size_t N>
Matrix<T, M, N> operator+(Matrix<T, M, N> lhs, const Matrix<T, M, N>& rhs) {
	for (size_t i = 0; i < M; ++i) {
		for (size_t j = 0; j < N; ++i) {
			lhs[i][j] += rhs[i][j];
		}
	}
	return lhs;
}
