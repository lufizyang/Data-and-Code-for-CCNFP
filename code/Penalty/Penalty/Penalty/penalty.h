#ifndef _PENALTY_H_
#define _PENALTY_H_

#include<functional>
#include<map>
#include<tuple>
#include<iostream>

#include"utils.h"


class Penalty
{
public:
	
	void GeneratePenalty();
	vec_d solve();
	Penalty(vec_l tail, vec_l head, vec_l capa, vec_l supd, std::function<double(const vec_d&)> func, 
		std::function<vec_d(const vec_d&)> jacf, int max_iter = 100, double tol = 1e-6, 
		double init_sigma = 1.0, double rho = 3.0);
	~Penalty() {}

	long NODENUM, EDGENUM;
	std::function<vec_d(const vec_d& x)> fcons, bnds;
	mat_d track;

private:

	std::map<std::string, long> edge2index_;
	vec_l tail_, head_, capa_, supd_;
	std::function<double(const vec_d&)> func_;
	std::function<vec_d(const vec_d&)> jacf_;
	double tol_, init_sigma_, rho_;
	int max_iter_;
	vec_d x0;

	// function definition
	std::function<mat_d(const vec_d& x)> jcons, jbnds;
	std::function<double(const vec_d& x, const double& sigma)> P;
	std::function<vec_d(const vec_d& x, const double& sigma)> dP;
};







#endif // !_PENALTY_H_

