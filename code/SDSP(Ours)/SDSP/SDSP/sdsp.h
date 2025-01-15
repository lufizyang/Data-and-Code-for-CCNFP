#ifndef _SDSP_H_
#define _SDSP_H_

#include<functional>
#include"dijkstra.h"


class SDSP
{
public:
	SDSP(const std::vector<long>& tail, const std::vector<long>& head, const std::vector<double>& capa, 
		const std::vector<double>& supd, const std::function<std::vector<double>(const std::vector<double>& x)>& func_in, 
		const std::function<std::vector<double>(const std::vector<double>& x)>& jac_in);
	~SDSP() {}

	long edgenum, nodenum;
	std::vector<std::vector<double>> track;

	std::vector<double> SSP(const std::vector<double>& weight, std::vector<double>& ek, std::vector<long>& Sk, std::vector<long>& Dk);
	std::vector<double> solve(const int& sample_number = 5,const std::string& ini_med = "lpu", const double& reduce_coef = 0.5, const vec_d& wkini = vec_d(), const bool& verbose = false);

	std::function<std::vector<double>(const std::vector<double>&)> func, jacf;
	std::vector<double> capa_, e_, x_, p_;
	std::vector<long> tail_, head_, S_, D_;
private:
	
	
};

#endif // !_SDSP_H_

