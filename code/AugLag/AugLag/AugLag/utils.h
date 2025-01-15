#pragma once
#ifndef _UTILS_H_
#define _UTILS_H_

#include<string>
#include<vector>

using mat_l = std::vector<std::vector<long>>;
using mat_d = std::vector<std::vector<double>>;
using vec_d = std::vector<double>;
using vec_l = std::vector<long>;
using readflatreturn = std::tuple<vec_l, vec_l, vec_l, vec_l, vec_d, vec_d>;


void readMatFromFile(std::string filename, std::vector<std::vector<double>>& output);

readflatreturn ReadFlatFile(std::string edta_path, std::string ndta_path);

double VecMultiple(const vec_d& v1, const vec_d& v2);
vec_d VecMultiMat(const vec_d& v1, const mat_d& m2);

vec_d operator+(const vec_d& v1, const vec_d& v2);
vec_d operator+(const vec_d& v, const double& a);
vec_d operator*(const double& a, const vec_d& v);
vec_d operator*(const vec_d& al, const vec_d& v);
vec_d operator-(const vec_d& v1);
vec_d operator-(const vec_d& v1, const vec_d& v2);
vec_d& operator+=(vec_d& v1, const vec_d& v2);
vec_d& operator+=(vec_d& v, double& a);
vec_d& operator*=(vec_d& v1, const vec_d& v2);
vec_d& operator*=(vec_d& v, const double& a);

#endif // !_UTILS_H_
