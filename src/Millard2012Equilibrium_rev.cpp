#include "muscle_model_plugin/Millard2012Equilibrium.h"
#include <math.h>
#include <iostream>
#include <fstream>


Thelen2003MuscleGazebo::Thelen2003MuscleGazebo(const std::string name,
											   double fmax,
											   double l_0,
											   double l_t0,
											   double alpha,
											   double l_init,
											   double q_init)
{
	//this->p1 = 1.0/4.0;
	this->p2 = 3.0/4.0;
	this->p3 = 1.8;
	this->p4 = 0.95;
	this->p5 = 3.0/10.0;
	this->p6 = this->p3 - 1.0;
	this->p7 = 2.0 + 2.0 / this->p5;
	this->p8 = 10.0;

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
	return exp ( -( pow((l_n - 1.0), 2) ) / this->g );
}

double Thelen2003MuscleGazebo::comp_gamma_p(double l_n)
{
	double gamma_p = 0;

	if(l_n <= 1.0)
	{
		gamma_p = 0.0;
	}
	else
	{
		gamma_p = ( exp( this->k*(l_n - 1.0) / this->e_0m) - 1.0) / (exp(this->k) - 1.0);
	}

	return gamma_p;
}

double Thelen2003MuscleGazebo::comp_gamma_s(double eps_t)
{

	double gamma_s = 0;

	if(eps_t <= 0.0)
	{
		gamma_s = 0.0;
	}
	else if(eps_t <= this->eps_th)
	{
		gamma_s = this->a1 * (exp(this->a3*eps_t/this->eps_th) - 1.0) / (exp(this->a3) - 1.0);
	}
	else
	{
		gamma_s = this->a2*(eps_t - this->eps_th) + this->a1;
	}

	return gamma_s;
}

double Thelen2003MuscleGazebo::comp_alpha(double l)
{
	double alpha = acos(sqrt(1 - pow((this->l_0*sin(this->alpha)/l), 2)));

	return alpha;
}

Thelen2003MuscleGazebo::ReturnValue Thelen2003MuscleGazebo:: compute_muscle_dynamics(double q,
                                                                                     double u,
                                                                                     double l,
                                                                                     double l_bar)
{
	// inputs: q (activation), u (excitation), l (muscle length) and l_bar(muscle tendon length)
	// outputs: q_dot, l_dot, f (force)
	ReturnValue res;

	res.q_dot = this->comp_q_dot(q, u); // directly assign q_dot

	const double l_n = l / l_0; // normalized by optimal length l_0

	// passive and active force gains
	// gamma_p: Definition
	const double gamma_p = this->comp_gamma_p(l_n);
	const double gamma_a = this->comp_gamma_a(l_n);

	//series force gain
	double c_alpha = cos(comp_alpha(l));
	double eps_t = (l_bar - l*c_alpha - l_t0) / l_t0;
	double gamma_s = this->comp_gamma_s(eps_t);

	// compute contractile element force gain
	double gamma_c = gamma_s / c_alpha - gamma_p;

	// constrain length dynamics
	if(q < q_min) {
		q = q_min;
	}else if(q > q_max) {
		q = 1.0;
	}

	// compute normalized length dynamics
	double ln_dot = this->comp_ln_dot(gamma_a,
	                                  gamma_c,
	                                  q);

	double l_dot = p8 * ln_dot * l_0;

	res.f = gamma_s * fmax;
	res.lm_dot = l_dot;

	return res;
}

double Thelen2003MuscleGazebo::comp_q_dot(double q, double u)
{
	double q_dot = 0.0;
	double tau = 0.0;

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
	double ln_dot = 0.0;

	if( (gamma_c > 0.0) && (gamma_c < q*gamma_a*this->p3*this->p4)) {
		if(gamma_c <= q*gamma_a) {
			ln_dot = (this->p1 + this->p2*q)*(gamma_c - q*gamma_a)
			        / (q*gamma_a + gamma_c / this->p5);
		} else {
			ln_dot = (this->p1 + this->p2 * q)*( gamma_c - q*gamma_a)
			        * this->p6 / (this->p7*(q*gamma_a*this->p3 - gamma_c));
		}
	}else {
		if(gamma_c <= 0.0) {
			ln_dot = (this->p1 + this->p2*q)*( (1.0 + 1.0/this->p5) * gamma_c - q*gamma_a)
			        / (q*gamma_a);
		} else {
			double c1 = p1*p6/ (p7*p3*(1.0-p4));
			double c2 = p2*p6/(p7*p3*(1.0-p4));
			double c3 = (p3*p4 - 1.0)*(1.0-2.0*p4)/(1.0-p4);
			double c4 = p3*p4;
			double c5 = (p3*p4 - 1.0)/(p3*(1.0 - p4));

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

double Thelen2003MuscleGazebo::compute_moment_arm(double theta)
{
	// not properly debugged yet
	std::fstream p_file;
	std::string filename = "/home/devel/Downloads/r_models/" + this->name + ".dat";
	p_file.open(filename, std::ios::in); // load polynomial coefficients

	std::array<double, 16> p_vec;

	for(size_t i=0; i<16; ++i)
	{
		p_file >> p_vec[i];
	}

	double r = p_vec[15];

	for(size_t i=0; i<15; ++i)
	{
		r += p_vec[i] * std::pow(theta, (15 - i));
	}

	p_file.close();

	std::cout << r << std::endl;

	return r;
}


double comp_l_bar_TRIlong(double theta)
{
	double p11 = -0.0537 + 0.0; // adding ground offset
	double p12 = -0.0137 + 0.0;

	double p21 = -0.0271 -0.017545; //adding humerus offset
	double p22 = -0.1144 -0.007;

	double p31 = -0.0318 -0.017545;
	double p32 = -0.2264 -0.007;

	double p41 = -0.0174 -0.017545;
	double p42 = -0.2676 -0.007;

	//double p51 = -0.0219 -0.011445; // adding ulna_radius_hand offset
	//double p52 = 0.0105 -0.2974;

	double p51l = -0.0219;
	double p52l = 0.0105;
	// rotate point and add offset to get global point
	double p51g = std::cos(theta) * p51l - std::sin(theta) * p52l -0.011445;
	double p52g = std::sin(theta) * p51l +std::cos(theta) * p52l -0.2974;

	// compute distance between each pair of points that form a direct connection globally
	double d1 = std::sqrt( std::pow((p21 - p11), 2) + std::pow((p22 - p12), 2) );
	double d2 = std::sqrt( std::pow((p31 - p21), 2) + std::pow((p32 - p22), 2) );
	double d3 = std::sqrt( std::pow((p41 - p31), 2) + std::pow((p42 - p32), 2) );
	double d4 = std::sqrt( std::pow((p51g - p41), 2) + std::pow((p52g - p42), 2) );

	double total_d = d1 + d2 + d3 + d4;

	return total_d;
}

double comp_r_TRIlong(double theta)
{

//	double p41 = -0.0174 -0.017545;
//	double p42 = -0.2676 -0.007;

//	double p51l = -0.0219;
//	double p52l = 0.0105;

//	// rotate point and add offset to get global point
//	double p51g = std::cos(theta) * p51l - std::sin(theta) * p52l -0.011445;
//	double p52g = std::sin(theta) * p51l +std::cos(theta) * p52l -0.2974;

//	double b1 = std::cos(theta) * p51l - std::sin(theta) * p52l;
//	double b2 = std::cos(theta) * p51l - std::sin(theta) * p52l;

//	double m1 = p51g - p41;
//	double m2 = p52g - p42;

//	double nb = std::sqrt( std::pow((b1), 2) + std::pow((b2), 2) );
//	double nm = std::sqrt( std::pow((m1), 2) + std::pow((m2), 2) );

//	double dot_product = m1 * b1 + m2 * b2;
//	double cos_phi = dot_product / (nb + nm);
//	double sin_phi = std::sqrt(std::pow((1 - cos_phi), 2));

//	double r = - nb * sin_phi;

	std::fstream p_file;
	p_file.open("/home/devel/Downloads/r_models/TRIlong.dat", std::ios::in); // load polynomial coefficients

	std::array<double, 16> p_vec;

	for(size_t i=0; i<16; ++i)
	{
		p_file >> p_vec[i];
	}

	double r = p_vec[15];

	for(size_t i=0; i<15; ++i)
	{
		r += p_vec[i] * std::pow(theta, (15 - i));
	}

	p_file.close();

	return r;
}

double comp_l_bar_TRIlat(double theta)
{
	double p11 = -0.0060 -0.017545;
	double p12 = -0.1265 -0.007;

	double p21 = -0.0234 -0.017545; //adding humerus offset
	double p22 = -0.1453 -0.007;

	double p31 = -0.0318 -0.017545;
	double p32 = -0.2264 -0.007;

	double p41 = -0.0174 -0.017545;
	double p42 = -0.2676 -0.007;

	//double p51 = -0.0219 -0.011445; // adding ulna_radius_hand offset
	//double p52 = 0.0105 -0.2974;

	double p51l = -0.0219;
	double p52l = 0.0105;
	// rotate point and add offset to get global point
	double p51g = std::cos(theta) * p51l - std::sin(theta) * p52l -0.011445;
	double p52g = std::sin(theta) * p51l + std::cos(theta) * p52l -0.2974;

	// compute distance between each pair of points that form a direct connection globally
	double d1 = std::sqrt( std::pow((p21 - p11), 2) + std::pow((p22 - p12), 2) );
	double d2 = std::sqrt( std::pow((p31 - p21), 2) + std::pow((p32 - p22), 2) );
	double d3 = std::sqrt( std::pow((p41 - p31), 2) + std::pow((p42 - p32), 2) );
	double d4 = std::sqrt( std::pow((p51g - p41), 2) + std::pow((p52g - p42), 2) );

	double total_d = d1 + d2 + d3 + d4;

	return total_d;
}

double comp_r_TRIlat(double theta)
{

//	double p41 = -0.0174 -0.017545;
//	double p42 = -0.2676 -0.007;

//	double p51l = -0.0219;
//	double p52l = 0.0105;

//	// rotate point and add offset to get global point
//	double p51g = std::cos(theta) * p51l - std::sin(theta) * p52l -0.011445;
//	double p52g = std::sin(theta) * p51l +std::cos(theta) * p52l -0.2974;

//	double b1 = std::cos(theta) * p51l - std::sin(theta) * p52l;
//	double b2 = std::cos(theta) * p51l - std::sin(theta) * p52l;

//	double m1 = p51g - p41;
//	double m2 = p52g - p42;

//	double nb = std::sqrt( std::pow((b1), 2) + std::pow((b2), 2) );
//	double nm = std::sqrt( std::pow((m1), 2) + std::pow((m2), 2) );

//	double dot_product = m1 * b1 + m2 * b2;
//	double cos_phi = dot_product / (nb + nm);
//	double sin_phi = std::sqrt(std::pow((1 - cos_phi), 2));

//	double r = - nb * sin_phi;

	std::fstream p_file;
	p_file.open("/home/devel/Downloads/r_models/TRIlat.dat", std::ios::in); // load polynomial coefficients

	std::array<double, 16> p_vec;

	for(size_t i=0; i<16; ++i)
	{
		p_file >> p_vec[i];
	}

	double r = p_vec[15];

	for(size_t i=0; i<15; ++i)
	{
		r += p_vec[i] * std::pow(theta, (15 - i));
	}

	p_file.close();

	return r;
}

double comp_l_bar_TRImed(double theta)
{
	double p11 = -0.0084 -0.017545;
	double p12 = -0.1370 -0.007;

	double p21 = -0.0260 -0.017545; //adding humerus offset
	double p22 = -0.1514 -0.007;

	double p31 = -0.0318 -0.017545;
	double p32 = -0.2264 -0.007;

	double p41 = -0.0174 -0.017545;
	double p42 = -0.2676 -0.007;

	//double p51 = -0.0219 -0.011445; // adding ulna_radius_hand offset
	//double p52 = 0.0105 -0.2974;

	double p51l = -0.0219;
	double p52l = 0.0105;
	// rotate point and add offset to get global point
	double p51g = std::cos(theta) * p51l - std::sin(theta) * p52l -0.011445;
	double p52g = std::sin(theta) * p51l + std::cos(theta) * p52l -0.2974;

	// compute distance between each pair of points that form a direct connection globally
	double d1 = std::sqrt( std::pow((p21 - p11), 2) + std::pow((p22 - p12), 2) );
	double d2 = std::sqrt( std::pow((p31 - p21), 2) + std::pow((p32 - p22), 2) );
	double d3 = std::sqrt( std::pow((p41 - p31), 2) + std::pow((p42 - p32), 2) );
	double d4 = std::sqrt( std::pow((p51g - p41), 2) + std::pow((p52g - p42), 2) );

	double total_d = d1 + d2 + d3 + d4;

	return total_d;
}

double comp_r_TRImed(double theta)
{

//	double p41 = -0.0174 -0.017545;
//	double p42 = -0.2676 -0.007;

//	double p51l = -0.0219;
//	double p52l = 0.0105;

//	// rotate point and add offset to get global point
//	double p51g = std::cos(theta) * p51l - std::sin(theta) * p52l -0.011445;
//	double p52g = std::sin(theta) * p51l +std::cos(theta) * p52l -0.2974;

//	double b1 = std::cos(theta) * p51l - std::sin(theta) * p52l;
//	double b2 = std::cos(theta) * p51l - std::sin(theta) * p52l;

//	double m1 = p51g - p41;
//	double m2 = p52g - p42;

//	double nb = std::sqrt( std::pow((b1), 2) + std::pow((b2), 2) );
//	double nm = std::sqrt( std::pow((m1), 2) + std::pow((m2), 2) );

//	double dot_product = m1 * b1 + m2 * b2;
//	double cos_phi = dot_product / (nb + nm);
//	double sin_phi = std::sqrt(std::pow((1 - cos_phi), 2));

//	double r = - nb * sin_phi;

	std::fstream p_file;
	p_file.open("/home/devel/Downloads/r_models/TRImed.dat", std::ios::in); // load polynomial coefficients

	std::array<double, 16> p_vec;

	for(size_t i=0; i<16; ++i)
	{
		p_file >> p_vec[i];
	}

	double r = p_vec[15];

	for(size_t i=0; i<15; ++i)
	{
		r += p_vec[i] * std::pow(theta, (15 - i));
	}

	p_file.close();

	return r;
}

double comp_l_bar_BIClong(double theta)
{
	double p11 = -0.0392 + 0.0; // adding ground offset
	double p12 = 0.0035 + 0.0;

	double p21 = -0.0289 + 0.0; //adding humerus offset
	double p22 = 0.0139 + 0.0;

	double p31 = 0.0213 -0.017545;
	double p32 = 0.0179 -0.007;

	double p41 = 0.0238 -0.017545;
	double p42 = -0.0051 -0.007;

	double p51 = 0.0135 - 0.017545;
	double p52 = -0.0283 - 0.007;

	double p61 = 0.0107 - 0.017545;
	double p62 = -0.0774 - 0.007;

	double p71 = 0.0170 - 0.017545;
	double p72 = -0.1213 - 0.007;

	double p81 = 0.0228 - 0.017545;
	double p82 = -0.1754 - 0.007;

	double p91l = 0.0075;
	double p92l = -0.0484;

	// rotate point and add offset to get global point
	double p91g = std::cos(theta) * p91l - std::sin(theta) * p92l -0.011445;
	double p92g = std::sin(theta) * p91l + std::cos(theta) * p92l -0.2974;

	// compute distance between each pair of points that form a direct connection globally
	double d1 = std::sqrt( std::pow((p21 - p11), 2) + std::pow((p22 - p12), 2) );
	double d2 = std::sqrt( std::pow((p31 - p21), 2) + std::pow((p32 - p22), 2) );
	double d3 = std::sqrt( std::pow((p41 - p31), 2) + std::pow((p42 - p32), 2) );
	double d4 = std::sqrt( std::pow((p51 - p41), 2) + std::pow((p52 - p42), 2) );
	double d5 = std::sqrt( std::pow((p61 - p51), 2) + std::pow((p62 - p52), 2) );
	double d6 = std::sqrt( std::pow((p71 - p61), 2) + std::pow((p72 - p62), 2) );
	double d7 = std::sqrt( std::pow((p81 - p71), 2) + std::pow((p82 - p72), 2) );
	double d8 = std::sqrt( std::pow((p91g - p81), 2) + std::pow((p92g - p82), 2) );


	double total_d = d1 + d2 + d3 + d4 + d5 + d6 + d7 + d8;

	return total_d;
}

double comp_r_BIClong(double theta)
{

//	double p81 = 0.0228 - 0.017545;
//	double p82 = -0.1754 - 0.007;

//	double p91l = 0.0075;
//	double p92l = -0.0484;

//	// rotate point and add offset to get global point
//	double p91g = std::cos(theta) * p91l - std::sin(theta) * p92l -0.011445;
//	double p92g = std::sin(theta) * p91l +std::cos(theta) * p92l -0.2974;

//	double b1 = std::cos(theta) * p91l - std::sin(theta) * p92l;
//	double b2 = std::cos(theta) * p91l - std::sin(theta) * p92l;

//	double m1 = p91g - p81;
//	double m2 = p92g - p82;

//	double nb = std::sqrt( std::pow((b1), 2) + std::pow((b2), 2) );
//	double nm = std::sqrt( std::pow((m1), 2) + std::pow((m2), 2) );

//	double dot_product = m1 * b1 + m2 * b2;
//	double cos_phi = dot_product / (nb + nm);
//	double sin_phi = std::sqrt(std::pow((1 - cos_phi), 2));

//	double r = nb * sin_phi;

	std::fstream p_file;
	p_file.open("/home/devel/Downloads/r_models/BIClong.dat", std::ios::in); // load polynomial coefficients

	std::array<double, 16> p_vec;

	for(int i=0; i<16; ++i)
	{
		p_file >> p_vec[i];
	}

	double r = p_vec[15];

	for(int i=0; i<15; ++i)
	{
		r += p_vec[i] * std::pow(theta, (15 - i));
	}

	p_file.close();

	return r;
}

double comp_l_bar_BICshort(double theta)
{
	double p11 = 0.0047 + 0.0; // adding ground offset
	double p12 = -0.0123 + 0.0;

	double p21 = -0.0071 + 0.0; //adding humerus offset
	double p22 = -0.0400 + 0.0;

	double p31 = 0.0112 -0.017545;
	double p32 = -0.0758 -0.007;

	double p41 = 0.0170 -0.017545;
	double p42 = -0.1213 -0.007;

	double p51 = 0.0228 - 0.017545;
	double p52 = -0.1754 - 0.007;

	double p61l = 0.0075;
	double p62l = -0.0484;

	// rotate point and add offset to get global point
	double p61g = std::cos(theta) * p61l - std::sin(theta) * p62l -0.011445;
	double p62g = std::sin(theta) * p61l + std::cos(theta) * p62l -0.2974;

	// compute distance between each pair of points that form a direct connection globally
	double d1 = std::sqrt( std::pow((p21 - p11), 2) + std::pow((p22 - p12), 2) );
	double d2 = std::sqrt( std::pow((p31 - p21), 2) + std::pow((p32 - p22), 2) );
	double d3 = std::sqrt( std::pow((p41 - p31), 2) + std::pow((p42 - p32), 2) );
	double d4 = std::sqrt( std::pow((p51 - p41), 2) + std::pow((p52 - p42), 2) );
	double d5 = std::sqrt( std::pow((p61g - p51), 2) + std::pow((p62g - p52), 2) );



	double total_d = d1 + d2 + d3 + d4 + d5;

	return total_d;
}

double comp_r_BICshort(double theta)
{

//	double p51 = 0.0228 - 0.017545;
//	double p52 = -0.1754 - 0.007;

//	double p61l = 0.0075;
//	double p62l = -0.0484;

//	// rotate point and add offset to get global point
//	double p61g = std::cos(theta) * p61l - std::sin(theta) * p62l -0.011445;
//	double p62g = std::sin(theta) * p61l +std::cos(theta) * p62l -0.2974;

//	double b1 = std::cos(theta) * p61l - std::sin(theta) * p62l;
//	double b2 = std::cos(theta) * p61l - std::sin(theta) * p62l;

//	double m1 = p61g - p51;
//	double m2 = p62g - p52;

//	double nb = std::sqrt( std::pow((b1), 2) + std::pow((b2), 2) );
//	double nm = std::sqrt( std::pow((m1), 2) + std::pow((m2), 2) );

//	double dot_product = m1 * b1 + m2 * b2;
//	double cos_phi = dot_product / (nb + nm);
//	double sin_phi = std::sqrt(std::pow((1 - cos_phi), 2));

//	double r = nb * sin_phi;

	std::fstream p_file;
	p_file.open("/home/devel/Downloads/r_models/BICshort.dat", std::ios::in); // load polynomial coefficients

	std::array<double, 16> p_vec;

	for(size_t i=0; i<16; ++i)
	{
		p_file >> p_vec[i];
	}

	double r = p_vec[15];

	for(size_t i=0; i<15; ++i)
	{
		r += p_vec[i] * std::pow(theta, (15 - i));
	}

	p_file.close();

	//std::cout << r << std::endl;

	return r;
}

double comp_l_bar_BRA(double theta)
{
	double p11 = 0.0068 -0.017545; // adding ground offset
	double p12 = -0.1739 -0.007;

	double p21l = -0.0032;
	double p22l = -0.0239;

	// rotate point and add offset to get global point
	double p21g = std::cos(theta) * p21l - std::sin(theta) * p22l -0.011445;
	double p22g = std::sin(theta) * p21l + std::cos(theta) * p22l -0.2974;

	// compute distance between each pair of points that form a direct connection globally
	double d1 = std::sqrt( std::pow((p21g - p11), 2) + std::pow((p22g - p12), 2) );

	double total_d = d1;

	return total_d;
}

double comp_r_BRA(double theta)
{

//	double p11 = 0.0068 -0.017545; // adding ground offset
//	double p12 = -0.1739 -0.007;

//	double p21l = -0.0032;
//	double p22l = -0.0239;

//	// rotate point and add offset to get global point
//	double p21g = std::cos(theta) * p21l - std::sin(theta) * p22l -0.011445;
//	double p22g = std::sin(theta) * p21l +std::cos(theta) * p22l -0.2974;

//	double b1 = std::cos(theta) * p21l - std::sin(theta) * p22l;
//	double b2 = std::cos(theta) * p21l - std::sin(theta) * p22l;

//	double m1 = p21g - p11;
//	double m2 = p22g - p12;

//	double nb = std::sqrt( std::pow((b1), 2) + std::pow((b2), 2) );
//	double nm = std::sqrt( std::pow((m1), 2) + std::pow((m2), 2) );

//	double dot_product = m1 * b1 + m2 * b2;
//	double cos_phi = dot_product / (nb + nm);
//	double sin_phi = std::sqrt(std::pow((1 - cos_phi), 2));

//	double r = nb * sin_phi;

	std::fstream p_file;
	p_file.open("/home/devel/Downloads/r_models/BRA.dat", std::ios::in); // load polynomial coefficients

	std::array<double, 16> p_vec;

	for(int i=0; i<16; ++i)
	{
		p_file >> p_vec[i];
	}

	double r = 0.0;

	for(int i=0; i<15; ++i)
	{
		r += p_vec[i] * std::pow(theta, (15 - i));
	}
	r += p_vec[15];
	p_file.close();

	//std::cout << r << std::endl;

	return r;
}

