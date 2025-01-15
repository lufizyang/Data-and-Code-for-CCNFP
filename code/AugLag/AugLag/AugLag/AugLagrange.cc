#include "AugLagrange.h"

#include<iostream>
#include<math.h>

using namespace std;

AugmentedLagrange::AugmentedLagrange(vec_l tail, vec_l head, vec_l capa, vec_l supd, std::function<double(const std::vector<double>&)> func,
	std::function<std::vector<double>(const std::vector<double>&)> jacf, int max_iter, double tol,
	double init_mu, double init_lamb, double init_rho, double rho_up) {
	
	tail_ = tail;
	head_ = head;
	capa_ = capa;
	supd_ = supd;
	func_ = func;
	jacf_ = jacf;
	max_iter_ = max_iter;
	NODENUM = supd_.size();
	EDGENUM = tail_.size();
	x0 = vector<double>(EDGENUM, 1.0);
	init_mu_ = init_mu;
	init_lamb_ = init_lamb;
	init_rho_ = init_rho;
	rho_up_ = rho_up;
	tol_ = tol;
	for (size_t i = 0; i < EDGENUM; i++)
	{
		edge2index_[to_string(tail_[i]) + "_" + to_string(head_[i])] = i;
	}
}


void AugmentedLagrange::GenerateAugLagrange() {

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
		auto eq_cons = [ai,bi](const vector<double>& x) -> double
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
		ais.push_back(generate_cons(ai,supd_[i]));
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
	// Lagrange function
	Lag = [this](const vec_d& x, const vec_d& lamb, const vec_d& mu, const double& rho) -> double
	{
		vec_d fcx = fcons(x);
		vec_d bdx = bnds(x);
		double thrits = 0.0;
		if (mu.size() != bdx.size())
		{
			cout << "The length of 'mu' and 'bounds' are not equal! " << endl;
			exit;
		}
		for (size_t i = 0; i < bdx.size(); i++)
		{
			double maxtempvalue = max(mu[i] / rho + bdx[i], 0.0);
			thrits += maxtempvalue * maxtempvalue - mu[i] * mu[i] / rho / rho;
		}
		return func_(x) + VecMultiple(lamb, fcx) + rho * VecMultiple(fcx, fcx) / 2.0 + rho * thrits / 2.0;
	};
	// Jacobean function of Lagrange
	dLag = [this](const vec_d& x, const vec_d& lamb, const vec_d& mu, const double& rho) -> vec_d
	{
		mat_d jcx = jcons(x), jbd = jbnds(x);
		vec_d rcx = fcons(x);
		vec_d bdx = bnds(x), thrits;
		for (size_t i = 0; i < rcx.size(); i++)
		{
			rcx[i] *= rho;
		}
		for (size_t i = 0; i < bdx.size(); i++)
		{
			thrits.push_back(rho * max(mu[i] / rho + bdx[i], 0.0));
		}
		return jacf_(x) + VecMultiMat(lamb, jcx) + VecMultiMat(rcx, jcx) + VecMultiMat(thrits, jbd);
	};
}


vec_d AugmentedLagrange::solve() {
	// Generate the Lagrange function
	GenerateAugLagrange();
	vec_d x = x0;
	track.push_back(x);
	vec_d lamb = vec_d(NODENUM, init_lamb_);
	vec_d mu = vec_d(2 * EDGENUM, init_mu_);
	double rho = init_rho_;
	for (size_t _ = 0; _ < max_iter_; _++)
	{
		size_t k = 0;
		mat_d track_inner = mat_d(1, x0);
		mat_d trackd_inner = mat_d(1, dLag(x0, lamb, mu, rho));
		while (true)
		{
			vec_d d = -dLag(x, lamb, mu, rho);
			vec_d alpha = vec_d(d.size(), 1);
			for (size_t i = 0; i < d.size(); i++)
			{
				if (x[i] + d[i] <= 0.0)	     { alpha[i] = -x[i] / d[i]; }
				if (x[i] + d[i] >= capa_[i]) { alpha[i] = (capa_[i] - x[i]) / d[i]; }
			}
			// One-dim line search
			if (true) // track_inner.size() < 2
			{
				while (Lag(x + alpha * d, lamb, mu, rho) > Lag(x, lamb, mu, rho) + 0.1 * VecMultiple(alpha * dLag(x, lamb, mu, rho), d))
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
			if (max(error_iter, -error_iter) < tol_ / 10.0) { break; }
			track_inner.push_back(x);
			trackd_inner.push_back(-d);
			k += 1;
			if (k > 500) { break; }
		}
		track.push_back(x);
		vec_d fcx = fcons(x), bdx = bnds(x);
		double dcx = VecMultiple(fcx, fcx), s_ = 0.0;
		for (size_t i = 0; i < bdx.size(); i++)
		{
			double maxit = max(bdx[i], -mu[i] / rho);
			s_ += maxit * maxit;
		}
		if (sqrt(dcx + s_) <= tol_) { break; } // Stopping criterion
		//cout << "Iteraion: " << _ << ", the first item is " << dcx << ", and the second item is " << s_ << ", the sqrt of the sum of above two is " << sqrt(dcx + s_) << ", the tol is " << tol_ << ". " << endl;
		lamb += rho * fcx;
		for (size_t i = 0; i < mu.size(); i++)
		{
			mu[i] = max(mu[i] + rho * bdx[i], 0.0);
		}
		rho *= rho_up_;
	}
	return track.back();
}