#pragma once

#include <vector>
#include <iostream>
#include <cassert>
#include <utility>

template <typename T, size_t M, size_t N>
class Matrix {
public:
	Matrix() : data(M, std::vector<T>(N, 0)) {}

	Matrix(std::initializer_list<std::initializer_list<T>> _data)
	 : data(M, std::vector<T>(N)) {
		assert(M == _data.size());
		size_t i = 0;
		for (const auto& row : _data){
			assert(N == row.size());
			size_t j = 0;
			for (const auto& element : row){
				data[i][j] = element;
				++j;
			}
			++i;
		}
	};

	std::vector<T>& operator[](const size_t i) { return data[i]; }
	const std::vector<T>& operator[](const size_t i) const { return data[i]; }

	std::pair<size_t, size_t> size() const { return {data.size(), data[0].size()}; }

	long double norm() const {
		long double result = 0;
		for (size_t i = 0; i < M; ++i){
			long double row_result = 0;
			for (size_t j = 0; j < N; ++j){
				row_result += std::abs(data[i][j]);
			}
			result = std::max(row_result, result);
		}
		return result;
	}

	Matrix<T, N, M> transpose() {
		Matrix<T, N, M> result;
		for(size_t i = 0; i < M; ++i){
			for (size_t j = 0; j < N; ++j){
				result[j][i] = data[i][j];
			}
		}
		return result;
	}

	static Matrix<T, N, N> identity() {
		Matrix<T, N, N> result;
		for (size_t i = 0; i < std::min(M, N); ++i){
			result[i][i] = 1;
		}
		return result;
	}
private:
	std::vector<std::vector<T>> data;
};

template <size_t N>
using SquareMatrix = Matrix<long double, N, N>;

template <size_t N>
using Vector = Matrix<long double, N, 1>;


template <typename T, size_t M, size_t N>
std::ostream& operator<<(std::ostream& stream, const Matrix<T, M, N>& matrix) {
	for (size_t i = 0; i < M; ++i){
		for (size_t j = 0; j < N; ++j){
			stream << matrix[i][j] << '\t';
		}
		stream << '\n';
	}
	return stream;
}

template <typename T, size_t M, size_t N>
Matrix<T, M, N> operator+(Matrix<T, M, N> lhs, const Matrix<T, M, N>& rhs) {
	for (size_t i = 0; i < M; ++i) {
		for (size_t j = 0; j < N; ++j) {
			lhs[i][j] += rhs[i][j];
		}
	}
	return lhs;
}

template <typename T, size_t M, size_t N>
Matrix<T, M, N> operator-(Matrix<T, M, N> lhs, const Matrix<T, M, N>& rhs) {
	for (size_t i = 0; i < M; ++i) {
		for (size_t j = 0; j < N; ++j) {
			lhs[i][j] -= rhs[i][j];
		}
	}
	return lhs;
}

template <typename T, size_t M1, size_t M2, size_t M3>
Matrix<T, M1, M3> operator*(Matrix<T, M1, M2> lhs, const Matrix<T, M2, M3>& rhs) {
	Matrix<T, M1, M3> result;
	for (size_t i = 0; i < M1; ++i) {
		for (size_t j = 0; j < M3; ++j) {
			for (size_t k = 0; k < M2; ++k){
				result[i][j] += lhs[i][k] * rhs[k][j];
			}
		}
	}
	return result;
}

template <typename T, size_t M, size_t N>
Matrix<T, M, N> operator* (const T& alpha, Matrix<T, M, N> rhs) {
	for (size_t i = 0; i < M; ++i) {
		for (size_t j = 0; j < N; ++j) {
			rhs[i][j] *= alpha;
		}
	}
	return rhs;
}
