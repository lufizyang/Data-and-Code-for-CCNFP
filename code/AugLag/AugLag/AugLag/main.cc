#include "AugLagrange.h"

#include<ctime>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<numeric>
#include<Windows.h>
#include<Psapi.h>


using namespace std;
PROCESS_MEMORY_COUNTERS pms;

int main()
{
	// data load path setting
	string ndta_path = "../../../test_case/ndta.txt";
	string edta_path = "../../../test_case/edta.txt";

	// obj-type setting
	string obj_type = "log"; // log, sqrt, exp, mix

	// Data loading
	readflatreturn dat = ReadFlatFile(edta_path, ndta_path);
	vec_l tail = get<0>(dat), head = get<1>(dat), supd = get<2>(dat), capa = get<3>(dat);
	vec_d ucost = get<5>(dat), fcost = get<4>(dat);
	double mf = accumulate(fcost.begin(), fcost.end(), 0.0) / fcost.size();
	// Output the memory usage of data and objective functions
	double memo1;
	if (GetProcessMemoryInfo(GetCurrentProcess(), &pms, sizeof(pms)))
	{
		cout << "The Memory Usage of Data and objectives is " << pms.WorkingSetSize / 1024 << "." << endl;
	}
	

	// log(x)
	auto cons_log_a = [](const vec_d& a) -> function<double(const vec_d&)>
		{
			auto log_a = [a](const vec_d& x) -> double
				{
					double s = 0.0;
					for (size_t i = 0; i < x.size(); i++)
					{
						s += log(x[i] + 1) / log(a[i]);
					}
					return s;
				};
			return log_a;
		};
	auto jac_log_a = [](const vec_d& a) -> function<vec_d(const vec_d&)>
		{
			auto jac_log = [a](const vec_d& x) -> vec_d
				{
					vec_d newv;
					for (size_t i = 0; i < a.size(); i++)
					{
						newv.push_back(1 / (log(a[i]) * (x[i] + 1)));
					}
					return newv;
				};
			return jac_log;
		};

	// sqrt(x)
	auto cons_sqrt = [](const vec_d& a) -> function<double(const vec_d&)>
		{
			auto sqrt_x = [a](const vec_d& x) -> double
				{
					double s = 0.0;
					for (size_t i = 0; i < x.size(); i++)
					{
						s += pow(x[i] + 1, 1 / a[i]) - 1;
					}
					return s;
				};
			return sqrt_x;
		};
	auto jac_sqrt = [](const vec_d& a) -> function<vec_d(const vec_d&)>
		{
			auto jac_sqrt_x = [a](const vec_d& x) -> vec_d
				{
					vec_d newv;
					for (size_t i = 0; i < x.size(); i++)
					{
						newv.push_back(pow(x[i] + 1, 1.0 / a[i] - 1) / a[i]);
					}
					return newv;
				};
			return jac_sqrt_x;
		};

	// exp(x)
	auto cons_exp = [](const vec_d& a) -> function<double(const vec_d&)>
		{
			auto a_exp_x = [a](const vec_d& x) -> double
				{
					vec_d newv;
					for (size_t i = 0; i < x.size(); i++)
					{
						newv.push_back(1.0 / (1.0 + pow(a[i], -x[i])) - 1.0 / 2);
					}
					return accumulate(newv.begin(), newv.end(), 0.0);
				};
			return a_exp_x;
		};
	auto jac_exp = [](const vec_d& a) -> function<vec_d(const vec_d&)>
		{
			auto jac_exp_x = [a](const vec_d& x) -> vec_d
				{
					vec_d newv;
					for (size_t i = 0; i < x.size(); i++)
					{
						newv.push_back(1.0 / pow(1.0 + pow(a[i], -x[i]), 2.0) * log(a[i]) * pow(a[i], -x[i]));
					}
					return newv;
				};
			return jac_exp_x;
		};

	// mix(x)
	vector<unsigned> lu;
	int n0 = 0, n1 = 0, n2 = 0;
	for (size_t i = 0; i < tail.size(); i++)
	{
		if (fcost[i] < 4)
		{
			lu.push_back(0);
			n0++;
		}
		else if (fcost[i] > 7.5)
		{
			lu.push_back(2);
			n2++;
		}
		else
		{
			lu.push_back(1);
			n1++;
		}
	}
	//cout << "...nx are " << n0 << ", " << n1 << ", " << n2 << endl;
	auto cons_mix_ = [lu](const vec_d& a) -> function<double(const vec_d&)>
		{
			auto mix_x = [lu, a](const vec_d& x) -> double
				{
					vec_d newv;
					for (size_t i = 0; i < x.size(); i++)
					{
						if (lu[i] == 0)
						{
							newv.push_back(pow(x[i] + 1, 1.0 / a[i]) - 1.0);
						}
						else if (lu[i] == 1)
						{
							newv.push_back(1.0 / (1.0 + pow(a[i], -x[i])) - 1.0 / 2.0);
						}
						else
						{
							newv.push_back(log(x[i] + 1) / log(a[i]));
						}
					}
					return accumulate(newv.begin(), newv.end(), 0.0);
				};
			return mix_x;
		};
	auto cons_mix_jac = [lu](const vec_d& a) -> function<vec_d(const vec_d&)>
		{
			auto jac_x = [lu, a](const vec_d& x) -> vec_d
				{
					vec_d newv;
					for (size_t i = 0; i < x.size(); i++)
					{
						if (lu[i] == 0)
						{
							newv.push_back(pow(x[i] + 1, 1.0 / a[i] - 1.0) / a[i]);
						}
						else if (lu[i] == 1)
						{
							newv.push_back(1.0 / pow(1.0 + pow(a[i], -x[i]), 2) * log(a[i]) * pow(a[i], -x[i]));
						}
						else
						{
							newv.push_back(1.0 / (log(a[i]) * (x[i] + 1)));
						}
					}
					return newv;
				};
			return jac_x;
		};
	
	// obj. setting
	function<double(const vec_d&)> fucl;
	function<vec_d(const vec_d&)> jacl;
	if (obj_type == "log") { fucl = cons_log_a(ucost + 1.0); jacl = jac_log_a(ucost + 1.0); }
	else if (obj_type == "sqrt") { fucl = cons_sqrt(ucost + 1.0); jacl = jac_sqrt(ucost + 1.0); }
	else if (obj_type == "exp") { fucl = cons_exp(ucost + 1.0); jacl = jac_exp(ucost + 1.0); }
	else if (obj_type == "mix") { fucl = cons_mix_(ucost + 1.0); jacl = cons_mix_jac(ucost + 1.0); }
	else { cout << "The obj. setting gets something wrong!" << endl; }

	// AugLag
	memo1 = pms.WorkingSetSize / 1024;
	AugmentedLagrange alg(tail, head, capa, supd, fucl, jacl, 500, 1e-6);
	const clock_t begin_time = clock();
	vec_d sol = alg.solve();
	double seconds = double(clock() - begin_time) / 1000;
	cout << "The objective value of ALg is: " << std::setprecision(12) << fucl(sol) << "." << endl;
	cout << "The solving time of ALg is: " << seconds << "s." << endl;
	cout << "The number of iterations is: " << alg.track.size() << "." << endl;
	// Output the memory usage of AugmentedLagrange
	if (GetProcessMemoryInfo(GetCurrentProcess(), &pms, sizeof(pms)))
	{
		cout << "The Memory Usage of AugLag is " << pms.WorkingSetSize / 1024 - memo1 << "." << endl;
	}
	return 0;
}