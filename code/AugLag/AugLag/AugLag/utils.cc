#include"utils.h"
#include<iostream>
#include<fstream>
#include<sstream>

#include<tuple>

using namespace std;

void readMatFromFile(string filename, vector<vector<double>>& output) {
	ifstream fs(filename);
	fs.unsetf(ios::skipws);
	while (fs.good())
	{
		string line_str;
		getline(fs, line_str);
		if (line_str.size() < 2) { continue; }
		istringstream ss(line_str);
		vector<double> line_data;
		while (ss.good()) {
			float v;
			ss >> v;
			line_data.push_back(v);
		}
		output.push_back(line_data);
	}
}

readflatreturn ReadFlatFile(string edta_path, string ndta_path) {
	mat_d edta, ndta;
	vec_l tail, head, cap, bvec;
	vec_d ucost, fcost;

	readMatFromFile(edta_path, edta);
	for (size_t i = 0; i < edta[0].size(); i++)
	{
		tail.push_back((long)edta[0][i]);
		head.push_back((long)edta[1][i]);
		cap.push_back((long)edta[2][i]);
		fcost.push_back(edta[3][i]);
		ucost.push_back(edta[4][i]);
	}
	readMatFromFile(ndta_path, ndta);
	for (size_t i = 0; i < ndta[0].size(); i++)
	{
		bvec.push_back((long)ndta[0][i]);
	}
	return make_tuple(tail, head, bvec, cap, fcost, ucost);
}

vec_d operator+(const vec_d& v1, const vec_d& v2) {
	if (v1.size() != v2.size())
	{
		std::cout << "These two vectors can not add itemwise!" << std::endl;
		exit;
	}
	vec_d newv;
	for (size_t i = 0; i < v1.size(); i++)
	{
		newv.push_back(v1[i] + v2[i]);
	}
	return newv;
}

vec_d operator+(const vec_d& v, const double& a) {
	vec_d newv;
	for (size_t i = 0; i < v.size(); i++)
	{
		newv.push_back(v[i] + a);
	}
	return newv;
}

vec_d operator*(const double& a, const vec_d& v) {
	vec_d newv;
	for (size_t i = 0; i < v.size(); i++)
	{
		newv.push_back(a * v[i]);
	}
	return newv;
}

vec_d operator*(const vec_d& al, const vec_d& v) {
	if (al.size() != v.size())
	{
		std::cout << "These two vectors can not add itemwise!" << std::endl;
		exit;
	}
	vec_d newv;
	for (size_t i = 0; i < al.size(); i++)
	{
		newv.push_back(al[i] * v[i]);
	}
	return newv;
}

vec_d operator-(const vec_d& v1) {
	vec_d newv;
	for (size_t i = 0; i < v1.size(); i++)
	{
		newv.push_back(-v1[i]);
	}
	return newv;
}

vec_d operator-(const vec_d& v1, const vec_d& v2) {
	if (v1.size() != v2.size())
	{
		std::cout << "These two vectors can not add itemwise!" << std::endl;
		exit;
	}
	vec_d newv;
	for (size_t i = 0; i < v1.size(); i++)
	{
		newv.push_back(v1[i] - v2[i]);
	}
	return newv;
}

vec_d& operator+=(vec_d& v1, const vec_d& v2) {
	if (v1.size() != v2.size())
	{
		std::cout << "These two vectors can not add itemwise!" << std::endl;
		exit;
	}
	for (size_t i = 0; i < v1.size(); i++)
	{
		v1[i] += v2[i];
	}
	return v1;
}

vec_d& operator+=(vec_d& v, double& a) {
	for (size_t i = 0; i < v.size(); i++)
	{
		v[i] += a;
	}
	return v;
}

vec_d& operator*=(vec_d& v1, const vec_d& v2) {
	if (v1.size() != v2.size())
	{
		std::cout << "These two vectors can not add itemwise!" << std::endl;
		exit;
	}
	for (size_t i = 0; i < v1.size(); i++)
	{
		v1[i] *= v2[i];
	}
	return v1;
}

vec_d& operator*=(vec_d& v, const double& a) {
	for (size_t i = 0; i < v.size(); i++)
	{
		v[i] *= a;
	}
	return v;
}

double VecMultiple(const vec_d& v1, const vec_d& v2) {
	if (v1.size() != v2.size())
	{
		cout << "The entering vectors have different length!" << endl;
		exit;
	}
	double s = 0.0;
	for (size_t i = 0; i < v1.size(); i++)
	{
		s += v1[i] * v2[i];
	}
	return s;
}

vec_d VecMultiMat(const vec_d& v1, const mat_d& m2) {
	if (v1.size() != m2.size() or m2.size() == 0)
	{
		cout << "The length of v1 and the row number of m2 are different! " << endl;
		exit;
	}
	vec_d res;
	for (size_t i = 0; i < m2[0].size(); i++)
	{
		double s = 0.0;
		for (size_t j = 0; j < v1.size(); j++)
		{
			s += v1[j] * m2[j][i];
		}
		res.push_back(s);
	}
	return res;
}