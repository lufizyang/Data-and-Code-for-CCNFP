#ifndef _SDSP_DIJKSTRA_H_
#define _SDSP_DIJKSTRA_H_

#include"utils.h"

#include<vector>
#include<tuple>

class ArgHeap
{
public:
	ArgHeap(const std::vector<double>& list);
	~ArgHeap() {}

	void start_heaplify();
	void heaplify(const size_t& i);
	long pop();
	long length();
	void update(const size_t& i, const double& val);

private:
	std::vector<double> list_;
	std::vector<size_t> ordered;
	void swap(const size_t&, const size_t&);

};

class Dijkstra
{
public:
	Dijkstra(const std::vector<long>& tail, const std::vector<long>& head, const std::vector<double>& weight, const long& start_node);
	~Dijkstra() {}

	std::vector<long> solve(const long& end_node);

	std::vector<double> dis;
	std::vector<long> parent;

private:
	std::vector<std::vector<long>> nbh;
	std::vector<long> tail_, head_, visited;
	std::vector<double> weight_;
	long start_node_, nodenum_;
	double Max_Weight;

	void dijkstra();
	void dijkstrasparse();
	void dijkstranotsparse();
};

class MultiDijkstra
{
public:
	MultiDijkstra(const std::vector<long>& tail, const std::vector<long>& head, const std::vector<double>& weight, const long& start_node, const long& sp_num);
	~MultiDijkstra() {}

	std::vector <std::vector<long>> solve(const long& end_node);
	std::vector<std::vector<double>> Dis;
	std::vector<std::vector<long>> Parent;

private:
	std::vector<std::vector<double>> _Dis;
	std::vector<std::vector<long>> _Parent, _Path;
	std::vector<long> tail_, head_;
	std::vector<double> weight_;
	long start_node_, nodenum_, sp_num_;
	double Max_Weight;
};


#endif // !_SDSP_DIJKSTRA_H_
