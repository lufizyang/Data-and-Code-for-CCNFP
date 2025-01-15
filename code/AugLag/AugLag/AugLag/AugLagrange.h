#ifndef _AUGLAGRANGE_H_
#define _AUGLAGRANGE_H_

#include<functional>
#include<map>
#include<tuple>
#include<iostream>

#include "utils.h"


class AugmentedLagrange
{
public:

	void GenerateAugLagrange();
	std::vector<double> solve();
	AugmentedLagrange(vec_l tail, vec_l head, vec_l capa, vec_l supd, std::function<double(const std::vector<double>&)> func,
		std::function<std::vector<double>(const std::vector<double>&)> jacf, int max_iter = 100, double tol = 1e-8,
		double init_mu = 0.1, double init_lamb = 0.1, double init_rho = 0.1, double rho_up = 1.1);
	~AugmentedLagrange() {}

	long NODENUM, EDGENUM;
	std::function<std::vector<double>(const vec_d& x)> fcons, bnds;
	mat_d track;

private:

	std::map<std::string, long> edge2index_;
	vec_l tail_, head_, capa_, supd_;
	std::function<double(const std::vector<double>&)> func_;
	std::function<std::vector<double>(const std::vector<double>&)> jacf_;
	double tol_, init_mu_, init_lamb_, init_rho_, rho_up_;
	int max_iter_;
	std::vector<double> x0;

	// Define functions
	std::function<mat_d(const vec_d& x)> jcons;
	std::function<mat_d(const vec_d& x)> jbnds;
	std::function<double(const vec_d& x, const vec_d& lamb, const vec_d& mu, const double& rho)> Lag;
	std::function<vec_d(const vec_d& x, const vec_d& lamb, const vec_d& mu, const double& rho)> dLag;
};



#endif // !_AUGLAGRANGE_H_

