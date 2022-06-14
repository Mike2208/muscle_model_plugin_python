#include "muscle_model_plugin/Millard2012Equilibrium.h"
#include <math.h>
#include <iostream>



Thelen2003MuscleGazebo::Thelen2003MuscleGazebo(const std::string name,
                                   double fmax,
                                   double l_0,
                                   double l_t0,
                                   double alpha,
                                   double l_init,
                                   double q_init)
{
	this->p1 = 1.0/4.0;
	this->p2 = 3.0/4.0;
	this->p3 = 1.8;
	this->p4 = 0.95;
	this->p5 = 3.0/10.0;
	this->p6 = this->p3 - 1;
	this->p7 = 2.0 + 2.0 / this->p5;
	this->p8 = 1e1;

	double eps_0 = 0.033;

	this->eps_th = 99.0*eps_0*exp(3.0) / (166.0*exp(3) - 67);
	this->a1 = 0.33;
	this->a2 = 67.0 / (100.0 * (eps_0 - (99.0*eps_0*exp(3.0))/(166.0*exp(3.0)-67.0)));
	this->a3 = 3.0;

	this->a0 = 0.5;
	this->k = 4.0;
	this->e_0m = 0.6;
	this->g = 0.5;

	this->tau_a = 0.01;
	this->tau_d = 0.04;

	this->name = name;
	this->fmax = fmax;
	this->l_0 = l_0;
	this->l_t0 = l_t0;
	this->alpha = alpha;

	// constraints
	this->q_min = 0.01;
	this->q_max = 1.0;

	// muscle state
	this->q = q_init;
	this->l = l_init;

};

double Thelen2003MuscleGazebo::comp_gamma_a(double l_n)
{
	return exp ( -( pow((l_n - 1), 2) ) / this->g );
}

double Thelen2003MuscleGazebo::comp_gamma_p(double l_n)
{
	double gamma_p = 0;

	if(l_n <= 1){
		gamma_p = 0;
	}else {
		gamma_p = ( exp( this->k*(l_n - 1) / this->e_0m) - 1) / (exp(this->k) - 1);
	}

	return gamma_p;
}

double Thelen2003MuscleGazebo::comp_gamma_s(double eps_t)
{
	double gamma_s = 0;

	if(eps_t <= 0) {
		gamma_s = 0;
	}else if(eps_t <= this->eps_th) {
		gamma_s = this->a1 * (exp(this->a3*eps_t/this->eps_th) - 1) / (exp(this->a3) - 1);
	} else {
		gamma_s = this->a2*(eps_t - this->eps_th) + this->a1;
	}

	return gamma_s;
}

double Thelen2003MuscleGazebo::comp_alpha(double l)
{
	double alpha = acos(sqrt(1 - pow((this->l_0*sin(this->alpha)/l), 2)));

	return alpha;
}

Thelen2003MuscleGazebo::ReturnValue Thelen2003MuscleGazebo::compute_muscle_dynamics(double q,
                                                                                    double u,
                                                                                    double l,
                                                                                    double l_bar)
{
	// inputs: q (activation), u (excitation), l (muscle length) and l_bar(muscle tendon length)
	// outputs: q_dot, ln_dot, f (force)
	ReturnValue res;

	double q_dot = Thelen2003MuscleGazebo::comp_q_dot(q, u);
	res.q_dot = q_dot; // directly assign q_dot

	double l_n = l / l_0; // normalized by optimal length l_0

	// passive and active force gains
	double gamma_p = Thelen2003MuscleGazebo::comp_gamma_p(l_n);
	double gamma_a = Thelen2003MuscleGazebo::comp_gamma_a(l_n);

	//series force gain
	double c_alpha = cos(comp_alpha(l));
	double eps_t = (l_bar - l*c_alpha - l_t0) / l_t0;
	double gamma_s = Thelen2003MuscleGazebo::comp_gamma_s(eps_t);

	// compute contractile element force gain
	double gamma_c = gamma_s / c_alpha - gamma_p;

	// constrain length dynamics
	if(q < q_min) {
		q = q_min;
	}else if(q > q_max) {
		q = 1;
	}

	// compute normalized length dynamics
	double ln_dot = Thelen2003MuscleGazebo::comp_ln_dot(gamma_a,
	                                              gamma_c,
	                                              q);

	double lm_dot = p8 * ln_dot * l_0;

	res.f = gamma_s * fmax;
	res.lm_dot = lm_dot;

	return res;
}

double Thelen2003MuscleGazebo::comp_q_dot(double q, double u)
{
	double q_dot = 0;
	double tau = 0;

	if(u>q) {
		tau = this->tau_a * (0.5 + 1.5*q);
	}else {
		tau = this->tau_d / (0.5 + 1.5*q);
	}

	q_dot = (u - q)/tau;

	return q_dot;
}

double Thelen2003MuscleGazebo::comp_ln_dot(double gamma_a, double gamma_c, double q)
{
	double ln_dot = 0;

	if( (gamma_c > 0) && (gamma_c < q*gamma_a*this->p3*this->p4)) {
		if(gamma_c <= q*gamma_a) {
			ln_dot = (this->p1 + this->p2*q)*(gamma_c - q*gamma_a)
			        / (q*gamma_a + gamma_c / this->p5);
		} else {
			ln_dot = (this->p1 + this->p2 * q)*( gamma_c - q*gamma_a)
			        * this->p6 / (this->p7*(q*gamma_a*this->p3 - gamma_c));
		}
	}else {
		if(gamma_c <= 0) {
			ln_dot = (this->p1 + this->p2*q)*( (1+1/this->p5) * gamma_c - q*gamma_a)
			        / (q*gamma_a);
		} else {
			double c1 = p1*p6/ (p7*p3*(1-p4));
			double c2 = p2*p6/(p7*p3*(1-p4));
			double c3 = (p3*p4 - 1)*(1-2*p4)/(1-p4);
			double c4 = p3*p4;
			double c5 = (p3*p4 - 1)/(p3*(1 - p4));

			ln_dot = (c1 + c2*q) * (c3 - gamma_c + c4*q*gamma_a + c5 *gamma_c/(q*gamma_a) );
		}
	}

	return ln_dot;
}

double Thelen2003MuscleGazebo::compute_muscle_equilibrium(double l,
                                                    double q,
                                                    double l_bar)
{
	double l_n = l / l_0;

	double c_alpha = cos(Thelen2003MuscleGazebo::comp_alpha(l));
	double eps_t = ( l_bar - l*c_alpha - l_t0) / l_t0;

	double gamma_s = Thelen2003MuscleGazebo::comp_gamma_s(eps_t);
	double gamma_p = Thelen2003MuscleGazebo::comp_gamma_p(l_n);
	double gamma_a = Thelen2003MuscleGazebo::comp_gamma_a(l_n);
	double gamma_fv = (gamma_s/c_alpha - gamma_p) / (q * gamma_a);

	double val = (q * gamma_a * gamma_fv + gamma_p) * c_alpha - gamma_s;

	return val;
}


