#pragma once

#include <string>
#include <array>

class Thelen2003MuscleGazebo
{
	public:
		struct ReturnValue
		{
			double q_dot;
			double lm_dot;
			double f;
		};

		Thelen2003MuscleGazebo(const std::string name,
		                 double fmax,
		                 double l_0,
		                 double l_t0,
		                 double alpha,
		                 double l_init,
		                 double q_init);

		// compute active length dependent gain
		double comp_gamma_a(double l_n);

		// compute passive gain
		double comp_gamma_p(double l_n);

		// compute tendon (serial) gain
		double comp_gamma_s(double eps_t);

		// compute q_dot, the activation dynamics
		double comp_q_dot(double q,
		                  double u);

		// compute ln_dot, the normalized muscle velocity
		double comp_ln_dot(double gamma_a,
		                   double gamma_c,
		                   double q);

		// compute alpha, the pennation angle
		double comp_alpha(double l);

		ReturnValue compute_muscle_dynamics(double q,
		                                    double u,
		                                    double l,
		                                    double l_bar);

		// compute muscle equilibrium (e.g. after initialization)
		double compute_muscle_equilibrium(double l,
		                                  double q,
		                                  double l_bar);

		std::array<double,2> equalibrate_muscle(double l,
		                          double q,
		                          double l_bar);

	private:
		// params for ln_dot
		double p1;
		double p2;
		double p3;
		double p4;
		double p5;
		double p6;
		double p7;
		double p8;

		// params for gamma_s
		double eps_th;
		double a1;
		double a2;
		double a3;

		// params for gamma_a and gamma_p
		double a0;
		double k;
		double e_0m;
		double g;

		// params for q_dot
		double tau_a;
		double tau_d;

		std::string name;
		double fmax;
		double l_0;
		double l_t0;
		double alpha;

		// constraints
		double q_min;
		double q_max;

	public:
		// muscle state
		double q;
		double l;
};
