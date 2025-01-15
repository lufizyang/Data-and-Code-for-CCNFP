#include"dssp.h"

#include<random>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<functional>
#include<Windows.h>
#include<Psapi.h>

using namespace std;
using vec_d = vector<double>;
using vec_l = vector<long>;
PROCESS_MEMORY_COUNTERS pms;

int main() {

	// Data load path setting
	string ndta_path = "../../../test_case/ndta.txt";
	string edta_path = "../../../test_case/edta.txt";

	// obj.-type setting
	string obj_type = "log"; // log, sqrt, exp, mix

	// data loading
	readflatreturn dat = ReadFlatFile(edta_path, ndta_path);
	vec_l tail = get<0>(dat), head = get<1>(dat);
	vec_d ucost = get<5>(dat), supd = get<2>(dat), capa = get<3>(dat);
	vec_d fcost = get<4>(dat);
	double mf = vsum(fcost) / fcost.size();
	cout << "The number of node is " << supd.size() << ", and the number of edge is " << tail.size() << "." << endl;
	int count_source = std::count_if(supd.begin(), supd.end(), [](int x) { return x > 0; });
	int count_sink = std::count_if(supd.begin(), supd.end(), [](int x) { return x < 0; });
	cout << "The number of source is " << count_source << ", and the number of sink is " << count_sink << "." << endl;
	//// Output the memory usage of data and objective functions
	//if (GetProcessMemoryInfo(GetCurrentProcess(), &pms, sizeof(pms)))
	//{
	//	cout << "The Memory Usage of Data and objectives is " << pms.WorkingSetSize / 1024 << "." << endl;
	//}
	double memo1;

	// log-type obj. construction
	auto cons_log_a = [](const vec_d& a) -> function<vec_d(const vec_d&)>
	{
		auto log_a = [a](const vec_d& x) -> vec_d
		{
			vec_d newv;
			for (size_t i = 0; i < x.size(); i++)
			{
				newv.push_back(log(x[i] + 1) / log(a[i]));
			}
			return newv;
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

	// sqrt-type obj. construction
	auto cons_sqrt = [](const vec_d& a) -> function<vec_d(const vec_d&)>
	{
		auto sqrt_x = [a](const vec_d& x) -> vec_d
		{
			vec_d newv;
			for (size_t i = 0; i < x.size(); i++)
			{
				newv.push_back(pow(x[i] + 1,1 / a[i]) - 1);
			}
			return newv;
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
				newv.push_back(pow(x[i] + 1, 1.0 / a[i] - 1.0) / a[i]);
			}
			return newv;
		};
		return jac_sqrt_x;
	};

	// exp-type obj.construction
	auto cons_exp = [](const vec_d& a) -> function<vec_d(const vec_d&)>
	{
		auto a_exp_x = [a](const vec_d& x) -> vec_d
		{
			vec_d newv;
			for (size_t i = 0; i < x.size(); i++)
			{
				newv.push_back(1.0 / (1 + pow(a[i], -x[i])) - 1.0 / 2.0);
			}
			return newv;
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
				newv.push_back(1.0 / pow(1 + pow(a[i], -x[i]), 2.0) * log(a[i]) * pow(a[i], -x[i]));
			}
			return newv;
		};
		return jac_exp_x;
	};

	// mix-type obj. construction
	vector<unsigned> lu; // Auxiliary variable: Use to define the mix-type objective function.
	int n0 = 0, n1 = 0, n2 = 0; // Counting the numbers of above three types obj. that used in the mix-type
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
	auto cons_mix_ = [lu](const vec_d& a) -> function<vec_d(const vec_d&)>
	{
		auto mix_x = [lu,a](const vec_d& x) -> vec_d
		{
			vec_d newv;
			for (size_t i = 0; i < x.size(); i++)
			{
				if (lu[i] == 0)
				{
					newv.push_back(pow(x[i] + 1, 1.0 / a[i]) - 1);
				}
				else if (lu[i] == 1)
				{
					newv.push_back(1.0 / (1.0 + pow(a[i], -x[i])) - 1.0 / 2);
				}
				else
				{
					newv.push_back(log(x[i] + 1) / log(a[i]));
				}
			}
			return newv;
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
	function<vec_d(const vec_d&)> fucl, jacl;
	if (obj_type == "log")
	{
		fucl = cons_log_a(ucost + 1.0);
		jacl = jac_log_a(ucost + 1.0);
	}
	else if (obj_type == "exp")
	{
		fucl = cons_exp(ucost + 1.0);
		jacl = jac_exp(ucost + 1.0);
	}
	else if (obj_type == "sqrt")
	{
		fucl = cons_sqrt(ucost + 1.0);
		jacl = jac_sqrt(ucost + 1.0);
	}
	else if (obj_type == "mix")
	{
		fucl = cons_mix_(ucost + 1.0); 
		jacl = cons_mix_jac(ucost + 1.0);
	}
	else
	{
		cout << "Please input the correct parameter. " << endl;
	}

	// DSSP
	memo1 = pms.WorkingSetSize / 1024;
	DSSP mdl(tail, head, capa, supd, fucl, jacl, "unitcost");
	const clock_t btime_dssp = clock();
	vec_d sol_dssp = mdl.Solve();
	double sec_dssp = double(clock() - btime_dssp) / 1000;
	cout << "The objectibe value of DSSP is: " << vsum(fucl(sol_dssp)) << "." << endl;
	cout << "The solving time of DSSP is: " << sec_dssp << "s." << endl;
	cout << endl;
	// Output the memory usage of DSSP
	if (GetProcessMemoryInfo(GetCurrentProcess(), &pms, sizeof(pms)))
	{
		cout << "The Memory Usage of DSSP is " << pms.WorkingSetSize / 1024 - memo1 << "." << endl;
	}
}