#include <vector>
#include <cmath>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/numeric/ublas/vector.hpp>
using namespace boost::numeric::ublas;
using namespace boost::math;
using namespace std;

#pragma once
class basis
{
public:
	basis(void);
	~basis(void);
	static matrix<double> polynomial(std::vector<double>& x,double& alpha,double& beta,int N);
	static matrix<double> gradient(std::vector<double>& r,double& alpha, double& beta, int N);
	static void gauss(double alpha,double beta,int N, std::vector<double>&x, std::vector<double>& w);
};

