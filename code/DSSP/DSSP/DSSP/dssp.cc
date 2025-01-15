#include"dssp.h"
#include"sdsp.h"

#include<algorithm>
#include<iostream>

using namespace std;

DSSP::DSSP(const vector<long>& tail, const vector<long>& head, const vector<double>& capa,
	const std::vector<double>& supd, const std::function<std::vector<double>(const std::vector<double>& x)>& func_in, 
	const std::function<std::vector<double>(const std::vector<double>& x)>& jacl_in, string cost_ini, int max_iter)
{
	tail_ = tail;
	head_ = head;
	capa_ = capa;
	supd_ = supd;
	func = func_in;
	jacl = jacl_in;
	MAXITER = max_iter;
	NODENUM = supd.size();
	EDGENUM = tail.size();
	if (cost_ini == "unitflow")
	{
		ecost = func(vector<double>(tail.size(), 1.0));
	}
	else if (cost_ini == "unitcost")
	{
		ecost = jacl(capa);
	}
	else
	{
		ecost = func(capa) / capa;
	}
}

vector<double> DSSP::Solve()
{
	size_t k = 0;
	vector<vector<double>> Solus;
	vector<double> MaxhisEdgeCost; // save the max edgecost in history of a edge.

	for (size_t i = 0; i < EDGENUM; i++)
	{
		MaxhisEdgeCost.push_back(0.0);
	}

	while (k <= MAXITER)
	{
		k++;
		SDSP mcf(tail_, head_, capa_, supd_, func, jacl);
		vector<double> x = mcf.SSP(ecost, mcf.e_, mcf.S_, mcf.D_);
		double tcost = vsum(func(x));
		// stop criterion
		if (find(Solus.begin(),Solus.end(), x) != Solus.end())
		{
			break;
		}
		Solus.push_back(x);
		// Update the ecost
		for (size_t i = 0; i < EDGENUM; i++)
		{
			if (x[i] > 1e-5)
			{
				ecost[i] = func(x)[i] / x[i];
				if (MaxhisEdgeCost[i] < ecost[i])
				{
					MaxhisEdgeCost[i] = ecost[i];
				}
			}
			else
			{
				if (MaxhisEdgeCost[i] > 1e-5)
				{
					ecost[i] = MaxhisEdgeCost[i];
				}
			}
		}
	}
	cout << "The iteration number is " << k << ". " << endl;
	return Solus.back();
}