#include"dijkstra.h"

#include<algorithm>
#include<numeric>
#include<math.h>
#include<iostream>

using namespace std;

ArgHeap::ArgHeap(const vector<double>& list) {
	list_ = list;
	for (size_t i = 0; i < list_.size(); i++)
	{
		ordered.push_back(i);
	}
	start_heaplify();
}

void ArgHeap::start_heaplify() {
	for (int i = list_.size()/2; i >= 0; i--)
	{
		heaplify(i);
	}
}

void ArgHeap::swap(const size_t& i, const size_t& j) {
	size_t temp = ordered[j];
	ordered[j] = ordered[i];
	ordered[i] = temp;
}

void ArgHeap::heaplify(const size_t& i) {
	if (i >= ordered.size()) { return; }
	int left = 2 * i + 1, right = 2 * i + 2;
	if ((left < ordered.size()) && (list_[ordered[left]] < list_[ordered[i]]))
	{
		swap(i, left);
		heaplify(left);
	}
	if ((right < ordered.size()) && (list_[ordered[right]] < list_[ordered[i]]))
	{
		swap(i, right);
		heaplify(right);
	}
}

long ArgHeap::pop() {
	if (ordered.size() == 0) { return -1; }
	swap(0, ordered.size() - 1);
	long i = ordered.back();
	ordered.pop_back();
	heaplify(0);
	return i;
}

long ArgHeap::length() {
	return ordered.size();
}

void ArgHeap::update(const size_t& i, const double& val) {
	list_[i] = val;
	size_t k = find(ordered.begin(), ordered.end(), i) - ordered.begin();
	while ((k > 0) && (list_[ordered[(k-1)/2]] > list_[ordered[k]]))
	{
		swap(k, (k - 1) / 2);
		k = (k - 1) / 2;
		heaplify(k);
	}
}


Dijkstra::Dijkstra(const vector<long>& tail, const vector<long>& head, const vector<double>& weight, const long& start_node) {
	tail_ = tail;
	head_ = head;
	weight_ = weight;
	start_node_ = start_node;
	vector<long> nodeset;
	for (size_t i = 0; i < tail.size(); i++)
	{
		if (find(nodeset.begin(),nodeset.end(),tail[i]) == nodeset.end())
		{
			nodeset.push_back(tail[i]);
		}
		if (find(nodeset.begin(),nodeset.end(),head[i]) == nodeset.end())
		{
			nodeset.push_back(head[i]);
		}
	}
	nodenum_ = nodeset.size();
	// initialize the vectors dis and parent;
	Max_Weight = accumulate(weight.begin(), weight.end(), 0.0);
	dis = vector<double>(nodenum_, Max_Weight);
	parent = vector<long>(nodenum_, tail_.size() + 1);
	visited = vector<long>(nodenum_, 0);
	dis[start_node_] = 0;
	parent[start_node_] = -1;
	nbh = vector<vector<long>>(nodenum_, vector<long>());
	for (size_t i = 0; i < tail_.size(); i++)
	{
		nbh[tail_[i]].push_back(i);
	}
}

void Dijkstra::dijkstrasparse() {
	ArgHeap h(dis);
	while (h.length() > 0)
	{
		long i = h.pop();
		visited[i] = 1;
		for (auto k : nbh[i])
		{
			if ((visited[head_[k]]) || (weight_[k] == Max_Weight)) { continue; }
			if (dis[head_[k]] > dis[i] + weight_[k])
			{
				dis[head_[k]] = dis[i] + weight_[k];
				parent[head_[k]] = k;
				h.update(head_[k], dis[head_[k]]);
			}
		}
	}
}

void Dijkstra::dijkstranotsparse() {
	visited[start_node_] = 1;
	long last_point = start_node_;
	for (size_t i = 0; i < nodenum_ - 1; i++)
	{
		double min_dis = Max_Weight;
		vector<long> nvl;
		for (size_t j = 0; j < nodenum_; j++)
		{
			if (visited[j] == 0)
			{
				nvl.push_back(j);
			}
		}
		for (size_t j : nvl)
		{
			if (dis[j] < min_dis)
			{
				min_dis = dis[j];
				last_point = j;
			}
		}
		visited[last_point] = 1;
		for (size_t k : nbh[last_point])
		{
			if ((visited[head_[k]]) || (weight_[k] == Max_Weight)) { continue; }
			if (dis[head_[k]] > dis[last_point] + weight_[k])
			{
				dis[head_[k]] = dis[last_point] + weight_[k];
				parent[head_[k]] = k;
			}
		}
	}
}

void Dijkstra::dijkstra() {
	long edgenum = tail_.size();
	if (edgenum >= nodenum_ * log(nodenum_))
	{
		dijkstranotsparse();
	}
	else {
		dijkstrasparse();
	}
}

vector<long> Dijkstra::solve(const long& end_node) {
	dijkstra();
	vector<long> path = vector<long>(1, parent[end_node]);
	while (parent[tail_[path[0]]] != -1)
	{
		//cout << "Judge " << end_node << " " << endl;
		//vector<long> add_node = vector<long>(1, parent[tail_[path[0]]]);
		path.insert(path.begin(), parent[tail_[path[0]]]);
	}
	/*cout << "Current path is: " << endl;
	for (size_t i = 0; i < path.size(); i++)
	{
		cout << "The edge " << path[i] << ": " << tail_[path[i]] << " -> " << head_[path[i]] << endl;
	}
	cout << endl;*/
	return path;
}


MultiDijkstra::MultiDijkstra(const std::vector<long>& tail, const std::vector<long>& head, const std::vector<double>& weight, const long& start_node, const long& sp_num)
{
	tail_ = tail;
	head_ = head;
	weight_ = weight;
	start_node_ = start_node;
	sp_num_ = sp_num;
	vector<long> nodeset;
	for (size_t i = 0; i < tail.size(); i++)
	{
		if (find(nodeset.begin(), nodeset.end(), tail[i]) == nodeset.end())
		{
			nodeset.push_back(tail[i]);
		}
		if (find(nodeset.begin(), nodeset.end(), head[i]) == nodeset.end())
		{
			nodeset.push_back(head[i]);
		}
	}
	nodenum_ = nodeset.size();
	Max_Weight = accumulate(weight.begin(), weight.end(), 0.0);
}

mat_l MultiDijkstra::solve(const long& end_node) {
	mat_l _Path, Path;
	Dijkstra sp(tail_, head_, weight_, start_node_);
	vec_l path = sp.solve(end_node);
	mat_l Del_edge;
	vec_l visited;
	_Path.push_back(path);
	_Dis.push_back(sp.dis);
	_Parent.push_back(sp.parent);
	Del_edge.push_back(vec_l());
	visited.push_back(0);
	while (Path.size() < sp_num_)
	{
		// 找到当前未拆解的最短路的下标
		size_t min_index = 0;
		double min_dist = Max_Weight;
		for (size_t i = 0; i < _Path.size(); i++)
		{
			if ((visited[i] == 0) && _Dis[i][end_node] < min_dist)
			{
				min_index = i;
				min_dist = _Dis[i][end_node];
			}
		}
		//cout << "the minimal distance is: " << min_dist << ", the max weight is: " << Max_Weight << endl;
		/*if (min_dist == Max_Weight)
		{
			cout << "There are " << Path.size() << " paths between node " << start_node_ << " and " << end_node << endl;
			break;
		}*/
		Path.push_back(_Path[min_index]);
		Dis.push_back(_Dis[min_index]);
		Parent.push_back(_Parent[min_index]);
		visited[min_index] = 1;
		vec_l del_edge = Del_edge[min_index];
		vec_l temp_path = _Path[min_index];
		for (size_t j = 0; j < temp_path.size(); j++)
		{
			// 构造去边的新图求解最短路
			vec_d wt = weight_;
			for (size_t m = 0; m < del_edge.size(); m++)
			{
				wt[m] = Max_Weight;
			}
			wt[temp_path[j]] = Max_Weight;
			Dijkstra spp(tail_, head_, wt, start_node_);
			vec_l add_path = spp.solve(end_node);
			bool Is_already_adding = false;
			for (size_t m = 0; m < _Path.size(); m++)
			{
				if (add_path == _Path[m])
				{
					Is_already_adding = true;
				}
			}
			if (Is_already_adding) { continue; }
			_Path.push_back(add_path);
			_Dis.push_back(spp.dis);
			_Parent.push_back(spp.parent);
			Del_edge.push_back(del_edge);
			Del_edge.back().push_back(temp_path[j]);
			visited.push_back(0);
		}
	}
	//// 找出Path中的前sp_num_个最短路
	//vector<pair<size_t,double>> dis_pair;
	//for (size_t i = 0; i < _Path.size(); i++)
	//{
	//	dis_pair.push_back(make_pair(i,_Dis[i][end_node]));
	//}
	//// 排序
	//sort(dis_pair.begin(), dis_pair.end(),
	//	[&](const auto& a, const auto& b)
	//	{
	//		return a.second < b.second;
	//	});
	//// 提取出前sp_num_个最短路输出
	//for (size_t i = 0; i < sp_num_; i++)
	//{
	//	size_t index_i = dis_pair[i].first;
	//	Path.push_back(_Path[index_i]);
	//	Dis.push_back(_Dis[index_i]);
	//	Parent.push_back(_Parent[index_i]);
	//}
	return Path;
}