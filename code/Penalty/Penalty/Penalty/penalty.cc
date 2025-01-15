#include"penalty.h"

#include<iostream>

using namespace std;

Penalty::Penalty(vec_l tail, vec_l head, vec_l capa, vec_l supd, function<double(const vec_d&)> func,
	function<vec_d(const vec_d&)> jacf, int max_iter, double tol, double init_sigma, double rho)
{
	tail_ = tail;
	head_ = head;
	capa_ = capa;
	supd_ = supd;
	func_ = func;
	jacf_ = jacf;
	max_iter_ = max_iter;
	tol_ = tol;
	init_sigma_ = init_sigma;
	rho_ = rho;
	NODENUM = supd_.size();
	EDGENUM = tail_.size();
	for (size_t i = 0; i < EDGENUM; i++)
	{
		edge2index_[to_string(tail_[i]) + "_" + to_string(head_[i])] = i;
	}
	x0 = vec_d(EDGENUM, 1.0);
}

void Penalty::GeneratePenalty()
{
	map<long, vector<long>> arci, arco;
	for (size_t i = 0; i < NODENUM; i++)
	{
		arci[i] = vector<long>();
		arco[i] = vector<long>();
	}
	for (size_t i = 0; i < EDGENUM; i++)
	{
		arci[head_[i]].push_back(tail_[i]);
		arco[tail_[i]].push_back(head_[i]);
	}
	// generate the constraint function
	auto generate_cons = [](const vector<double>& ai, const long& bi) -> function<double(const vector<double>&)>
	{
		auto eq_cons = [ai, bi](const vector<double>& x) -> double
		{
			double s = 0.0;
			for (size_t i = 0; i < x.size(); i++)
			{
				s += ai[i] * x[i];
			}
			return s - bi;
		};
		return eq_cons;
	};

	vector<function<double(const vector<double>&)>> ais;
	mat_d jacs;
	for (size_t i = 0; i < NODENUM; i++)
	{
		vec_d ai = vector<double>(EDGENUM, 0.0);
		for (long j : arci[i])
		{
			ai[edge2index_[to_string(j) + "_" + to_string(i)]] = -1;
		}
		for (long j : arco[i])
		{
			ai[edge2index_[to_string(i) + "_" + to_string(j)]] = 1;
		}
		ais.push_back(generate_cons(ai, supd_[i]));
		jacs.push_back(ai);
	}
	// Ax - b
	fcons = [ais](const vec_d& x) -> vec_d
	{
		vec_d res;
		for (size_t i = 0; i < ais.size(); i++)
		{
			res.push_back(ais[i](x));
		}
		return res;
	};
	jcons = [jacs](const vec_d& x) -> mat_d
	{
		return jacs;
	};
	// bounds
	vec_l cap = capa_;
	bnds = [cap](const vec_d& x) -> vec_d
	{
		vector<double> bnd;
		for (size_t i = 0; i < x.size(); i++)
		{
			bnd.push_back(-x[i]);
		}
		for (size_t i = 0; i < x.size(); i++)
		{
			bnd.push_back(x[i] - cap[i]);
		}
		return bnd;
	};
	// jac of bounds
	long edgenum = EDGENUM;
	jbnds = [edgenum](const vec_d& x) -> mat_d
	{
		mat_d jac_bnds = mat_d(2 * edgenum, vec_d(edgenum, 0));
		for (size_t i = 0; i < edgenum; i++)
		{
			jac_bnds[i][i] = -1;
			jac_bnds[i + edgenum][i] = 1;
		}
		return jac_bnds;
	};
	// Penalty function
	P = [this](const vec_d& x, const double& sigma) -> double
	{	
		vec_d fcx = fcons(x), bdx = bnds(x);
		double bds = 0.0;
		for (size_t i = 0; i < bdx.size(); i++)
		{
			double maxv = max(bdx[i], 0.0);
			bds += maxv * maxv;
		}
		return func_(x) + 0.5 * sigma * (VecMultiple(fcx, fcx) + bds);
	};
	// Jacobean function of Penalty
	dP = [this](const vec_d& x, const double& sigma) -> vec_d
	{
		mat_d jcx = jcons(x), jbd = jbnds(x);
		vec_d rcx = fcons(x), bdx = bnds(x), bds;
		for (size_t i = 0; i < bdx.size(); i++)
		{
			bds.push_back(max(bdx[i], 0.0));
		}
		return jacf_(x) + sigma * (VecMultiMat(rcx, jcx) + VecMultiMat(bds, jbd));
	};
}

vec_d Penalty::solve()
{
	GeneratePenalty();
	vec_d x = x0;
	track.push_back(x);
	double sigma = init_sigma_;
	for (size_t _ = 0; _ < max_iter_; _++)
	{
		size_t k = 0;
		mat_d track_inner = mat_d(1, x0);
		mat_d trackd_inner = mat_d(1, dP(x0, sigma));
		while (true)
		{
			k += 1;
			vec_d d = -dP(x, sigma);
			vec_d alpha = vec_d(d.size(), 1.0);
			for (size_t i = 0; i < d.size(); i++)
			{
				if (x[i] + d[i] <= 0.0)
				{
					alpha[i] = -x[i] / d[i];
				}
				if (x[i] + d[i] >= capa_[i])
				{
					alpha[i] = (capa_[i] - x[i]) / d[i];
				}
			}
			// Line search
			if (true)
			{
				while (P(x+alpha*d,sigma) > P(x,sigma) + 0.1 * VecMultiple(-alpha * d,d))
				{
					alpha *= 0.5;
				}
				x += alpha * d;
			}
			else {
				vec_d yk_1 = trackd_inner.back() - trackd_inner[trackd_inner.size() - 2];
				vec_d sk_1 = track_inner.back() - track_inner[track_inner.size() - 2];
				double alpha1 = VecMultiple(sk_1, yk_1) / VecMultiple(yk_1, yk_1);
				double alpha2 = VecMultiple(sk_1, sk_1) / VecMultiple(sk_1, yk_1);
				double alpha_base = (alpha1 + alpha2) / 2.0;
				for (size_t i = 0; i < EDGENUM; i++)
				{
					if (alpha[i] > alpha_base)
					{
						alpha[i] = alpha_base;
					}
					if (alpha[i] < 0)
					{
						alpha[i] = 0;
					}
				}
				x += alpha * d;
			}
			double error_iter = func_(x) - func_(track_inner.back());
			if (max(error_iter, -error_iter) < tol_ / 10.0)
			{
				break;
			}
			track_inner.push_back(x);
			trackd_inner.push_back(-d);
			if (k > 500) { break; }
		}
		track.push_back(x);
		double dcx = VecMultiple(fcons(x), fcons(x));
		if (sqrt(dcx) < tol_) { break; }
		//cout << "Iteraion: " << _ << ", the error is " << dcx << ", the sqrt of the sum of above two is " << sqrt(dcx) << ", the tol is " << tol_ << ". " << endl;
		sigma *= rho_;
	}
	return track.back();
}