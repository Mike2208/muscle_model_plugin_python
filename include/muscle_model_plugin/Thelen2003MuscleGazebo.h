#pragma once

#include <math.h>
#include <fstream>

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

		/**
		 * @brief Thelen2003MuscleGazebo
		 * @param name
		 * @param fmax
		 * @param l_0
		 * @param l_t0
		 * @param alpha
		 * @param l_init
		 * @param q_init
		 */
		Thelen2003MuscleGazebo(const std::string name,
		                 double fmax,
		                 double l_0,
		                 double l_t0,
		                 double alpha,
		                 double l_init,
		                 double q_init);

		/**
		 * @brief comp_gamma_a, length dependent gain
		 * @param l_n
		 * @return gamma_a
		 */
		double comp_gamma_a(double l_n);

		/**
		 * @brief Compute passive elastic gain
		 * @param l_n normalized muscle length
		 * @return gamma_p passive elastic gain
		 */
		double comp_gamma_p(double l_n);

		/**
		 * @brief Compute serial elastic (tendon) gain
		 * @param eps_t the tendon strain
		 * @return gamma_s series element gain
		 */
		double comp_gamma_s(double eps_t);

		/**
		 * @brief comp_q_dot, activation dynamics
		 * @param q, activation
		 * @param u, stimulus
		 * @return q_dot
		 */
		double comp_q_dot(double q,
		                  double u);

		/**
		 * @brief comp_ln_dot, the muscle dynamics
		 * @param gamma_a
		 * @param gamma_c
		 * @param q
		 * @return ln_dot
		 */
		double comp_ln_dot(double gamma_a,
		                   double gamma_c,
		                   double q);

		/**
		 * @brief comp_alpha, the pennation angle between muscle and tendon
		 * @param l
		 * @return alpha
		 */
		double comp_alpha(double l);

		ReturnValue compute_muscle_dynamics(double q,
		                                    double u,
		                                    double l,
		                                    double l_bar);

		/**
		 * @brief compute_muscle_equilibrium
		 * @param l
		 * @param q
		 * @param l_bar
		 * @return
		 */
		double compute_muscle_equilibrium(double l,
		                                  double q,
		                                  double l_bar);

		/**
		 * @brief equalibrate_muscle
		 * @param l
		 * @param q
		 * @param l_bar
		 * @return
		 */
		std::array<double,2> equalibrate_muscle(double l,
		                                        double q,
		                                        double l_bar);

		/**
		 * @brief compute_rs
		 * @param theta_s
		 * @return r
		 */
		double compute_rs(double theta_s); // TODO: improve accuracy

		/**
		 * @brief compute_re
		 * @param theta_e
		 * @return r
		 */
		double compute_re(double theta_e);

		/**
		 * @brief compute_l_bar
		 * @param theta_s
		 * @param theta_e
		 * @return l_bar
		 */
		double compute_l_bar(double theta_s,
		                     double theta_e);

	private:
		// params for ln_dot
		double p1 = 1.0/4.0;
		double p2 = 3.0/4.0;
		double p3 = 1.8;
		double p4 = 0.95;
		double p5 = 3.0/10.0;
		double p6 = this->p3 - 1.0;
		double p7 = 2 + 2 / this->p5;
		double p8 = 10.0;

		// params for gamma_s
		double eps_0 = 0.033;

		double eps_th = 99.0*this->eps_0*std::exp(3.0) / (166.0*std::exp(3) - 67);;
		double a1 = 0.33;
		double a2 = 67.0 / (100.0 * (eps_0 - (99.0*eps_0*exp(3.0))/(166.0*exp(3.0)-67.0)));
		double a3 = 3.0;

		// params for gamma_a and gamma_p
		double a0 = 0.5;
		double k = 4.0;
		double e_0m = 0.6;
		double g = 0.5;

		// params for q_dot
		double tau_a = 0.01;
		double tau_d = 0.04;

		std::string name;
		double fmax;
		double l_0;
		double l_t0;
		double alpha;

		// constraints
		double q_min = 0.01;
		double q_max = 1.0;

		// polynomial coefficients for moment arm (r) computation
		std::fstream rs_p_file;
		std::array<double, 13> rs_p_vec;
		std::fstream re_p_file;
		std::array<double, 13> re_p_vec;

		// polynomial coefficients for l_bar computation
		std::fstream l_bar_p_file;
		std::array<double, 64> l_bar_p_vec;

	public:
		// muscle state
		double q;
		double l;
};

/**
 * @brief rotate_p_urh
 * @param p1
 * @param p2
 * @param theta
 * @return p*g
 */
std::array<double, 2> rotate_p_urh(double p1,
                                   double p2,
                                   double theta);
/**
 * @brief comp_d
 * @param p1
 * @param p2
 * @return d
 */
double comp_d(double p1,
              double p2);

// TODO: extend to 3D
double comp_l_bar_BRA(double theta);
double comp_l_bar_TRIlong(double theta);
double comp_l_bar_TRIlat(double theta);
double comp_l_bar_TRImed(double theta);
double comp_l_bar_BIClong(double theta);
double comp_l_bar_BICshort(double theta);
