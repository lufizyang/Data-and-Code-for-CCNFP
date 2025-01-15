#include"sdsp.h"

#include<map>
#include<string>
#include<iostream>
#include<numeric>
#include<algorithm>



using namespace std;

vector<long> sortbb(const vector<long>& v1, const vector<double>& v2, const bool& accending)
{
	vector<pair<long, double>> temp;
	for (size_t i = 0; i < v1.size(); i++)
	{
		temp.push_back(make_pair(v1[i], v2[i]));
	}
	sort(temp.begin(), temp.end(), [accending](pair<long, double>a, pair<long, double>b) {
		if (accending)
		{
			return a.second < b.second;
		}
		else {
			return a.second > b.second;
		}
		});
	vector<long> newv;
	for (size_t i = 0; i < temp.size(); i++)
	{
		newv.push_back(temp[i].first);
	}
	return newv;
}


SDSP::SDSP(const std::vector<long>& tail, const std::vector<long>& head, const std::vector<double>& capa,
	const std::vector<double>& supd, const std::function<std::vector<double>(const std::vector<double>& x)>& func_in,
	const std::function<std::vector<double>(const std::vector<double>& x)>& jac_in)
{
	tail_ = tail;
	head_ = head;
	capa_ = capa;
	e_ = supd;
	func = func_in;
	jacf = jac_in;
	edgenum = tail.size();
	nodenum = supd.size();
	vector<double> Se, De;
	for (size_t i = 0; i < nodenum; i++)
	{
		if (e_[i] > 0)
		{
			S_.push_back(i);
			Se.push_back(e_[i]);
		}
		if (e_[i] < 0)
		{
			D_.push_back(i);
			De.push_back(e_[i]);
		}
	}
	// sorting S_,D_ by the value of e_
	S_ = sortbb(S_, Se, true);
	D_ = sortbb(D_, De, false);
}

vector<double> SDSP::SSP(const vector<double>& weight_, vector<double>& ek, vector<long>& Sk, vector<long>& Dk) {
	double MaxWeight = accumulate(weight_.begin(), weight_.end(), 0.0) + 1;
	long edgenum_solve = 2 * edgenum;
	vector<long> tail = tail_, head = head_;
	vector<double> weight = weight_, r = capa_, x = vector<double>(edgenum_solve, 0), p = vector<double>(nodenum, 0);
	vector<double> c_pi = vector<double>(edgenum_solve, 0);
	tail.insert(tail.end(), head_.begin(), head_.end());
	head.insert(head.end(), tail_.begin(), tail_.end());
	for (size_t i = 0; i < weight_.size(); i++)
	{
		weight.push_back(-weight_[i]);
	}
	for (size_t i = 0; i < edgenum; i++)
	{
		r.push_back(0);
	}
	// solve process
	for (size_t i = 0; i < edgenum_solve; i++)
	{
		if (r[i] < 1e-6)
		{
			c_pi[i] = MaxWeight;
		}
		else {
			c_pi[i] = weight[i] - p[tail[i]] + p[head[i]];
		}
	}
	map<string, vector<long>> edge2index;
	for (size_t i = 0; i < edgenum_solve; i++)
	{
		edge2index[to_string(tail[i]) + "_" + to_string(head[i])].push_back(i);
	}
	while (Sk.size() > 0)
	{
		long nodek = Sk.back(), nodel = Dk.back();
		Dijkstra dijk(tail, head, c_pi, nodek);
		vector<long> fix_edge_index_vec = dijk.solve(nodel);
		vector<double> dis = dijk.dis;
		for (size_t j = 0; j < nodenum; j++)
		{
			p[j] -= dis[j];
		}
		vector<double> compare_vector = { ek[nodek],-ek[nodel] };
		for (long i : fix_edge_index_vec)
		{
			compare_vector.push_back(r[i]);
		}
		double delta = *min_element(compare_vector.begin(), compare_vector.end());
		for (long i : fix_edge_index_vec)
		{
			x[i] += delta;
			r[i] -= delta;
			if (i < edgenum)
			{
				r[i + edgenum] += delta;
			}
			else {
				r[i - edgenum] += delta;
			}
		}
		ek[nodek] -= delta;
		ek[nodel] += delta;
		if (ek[nodek] < 1e-8)
		{
			Sk.pop_back();
		}
		if (ek[nodel] > -1e-8)
		{
			Dk.pop_back();
		}
		// update the c_pi
		for (size_t i = 0; i < edgenum_solve; i++)
		{
			if (r[i] < 1e-6)
			{
				c_pi[i] = MaxWeight;
			}
			else {
				c_pi[i] = weight[i] - p[tail[i]] + p[head[i]];
			}
		}
	}
	vector<double> resx;
	for (size_t i = 0; i < edgenum; i++)
	{
		resx.push_back(x[i] - x[i + edgenum]);
	}
	return resx;
}
