#pragma once

#include <vector>
#include <algorithm>
#include <cmath>

const long double pi = 3.14159265358979323846;

namespace Nodes{
	class Root{
	public:
		Root(size_t N, long double A, long double B)
			: n(N), a(A), b(B){};
		Root(std::initializer_list<long double> ilist){
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

		virtual long double operator[](size_t) const =0;
		long double at(size_t i) const {return this->operator[](i);};
		void transform(size_t m, long double c, long double d) {
			n = m; a = c; b = d;
		}

		size_t size() const{
			return n;
		}
		std::pair<long double, long double> interval() const{
			return {a, b};
		}

		virtual Root* addNode(long double) =0;

	protected:
		size_t n;
		long double a, b;
	};

	class Ordinary : public Nodes::Root{
	public:
		Ordinary(std::initializer_list<long double> ilist)
			:Root(ilist),
			 vals(ilist)
			 {}

		Ordinary(Root *nodes)
			:Root(*nodes)
			{
			vals.resize(nodes->size());
			for (size_t i = 0; i < nodes->size(); ++i){
				vals[i] = nodes->at(i);
			}
		}

		long double operator[](size_t i) const{
			return vals[i];
		}

		Ordinary* addNode(long double x){
			vals.insert(std::upper_bound(vals.begin(), vals.end(), x), x);
			return this;
		}
	private:
		std::vector<long double> vals;
	};


	class Chebyshev : public Nodes::Root{
	public:
		Chebyshev(size_t N, long double A, long double B)
			: Root(N, A, B), vals(N), a(A), b(B){
			for (size_t i = 1; i <= N; ++i){
				vals[i - 1] = cos(pi * (1.0 * i - 0.5) / N);
			}
		}

		long double operator[](size_t i) const{
			return a + (vals[i] + 1.0) * (b - a) / 2.0;
		}

		Ordinary* addNode(long double x){
			Ordinary* nodes = new Ordinary(this);
			nodes->addNode(x);
			return nodes;
		}
	private:
		std::vector<long double> vals;
		long double a, b;
	};

	class Uniform : public Nodes::Root{
	public:
		Uniform(size_t N, long double A, long double B)
			: Root(N, A, B) {}

		long double operator[](size_t i) const{
			return a + i * (b - a) / n;
		}

		Ordinary* addNode(long double x){
			Ordinary* nodes = new Ordinary(this);
			nodes->addNode(x);
			return nodes;
		}
	};
}
