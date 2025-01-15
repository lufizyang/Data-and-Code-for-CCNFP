#ifndef DYNAMIC_SLOPE_SCALING_PROCEDURE_H_
#define DYNAMIC_SLOPE_SCALING_PROCEDURE_H_

#include"utils.h"

#include<string>
#include<vector>
#include<functional>



using std::string;

class DSSP
{
public:
	// Solve function
	std::vector<double> Solve();

	DSSP(const std::vector<long>& tail, const std::vector<long>& head, const std::vector<double>& capa, 
		const std::vector<double>& supd, const std::function<std::vector<double>(const std::vector<double>& x)>& func_in, 
		const std::function<std::vector<double>(const std::vector<double>& x)>& jacl_in, 
		string cost_ini = "maxflow", int max_iter = 30);
	~DSSP() {};

private:
	int MAXITER, EDGENUM, NODENUM;
	std::function<std::vector<double>(const std::vector<double>&)> func, jacl;
	std::vector<long> tail_, head_;
	std::vector<double> capa_, supd_, ecost;
};



#endif // !DYNAMIC_SLOPE_SCALING_PROCEDURE_H_
