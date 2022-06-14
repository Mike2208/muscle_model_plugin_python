#include "muscle_model_plugin/Thelen2003MuscleGazebo.h"
#include "muscle_model_plugin/muscle_model_python_config.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <filesystem>


Thelen2003MuscleGazebo::Thelen2003MuscleGazebo(const std::string name,
                                               double fmax,
                                               double l_0,
                                               double l_t0,
                                               double alpha,
                                               double l_init,
                                               double q_init)
{
	// muscle specific parameters
	this->name = name;
	this->fmax = fmax;
	this->l_0 = l_0;
	this->l_t0 = l_t0;
	this->alpha = alpha;

	// muscle state
	this->q = q_init;
	this->l = l_init;

//	std::string filename_rs = "/home/devel/Downloads/data_files/W_s/" + this->name + ".dat";
//	std::string filename_re = "/home/devel/Downloads/data_files/W_e/" + this->name + ".dat";
//	std::string filename_l_bar = "/home/devel/Downloads/data_files/W_l_bar/" + this->name + ".dat";

	const std::filesystem::path module_data_dir = std::filesystem::path(MUSCLE_MODEL_PYTHON_MODULE_DIR) / "data_files";
	std::filesystem::path filename_rs = module_data_dir / "W_s" / (this->name + ".dat");
	std::filesystem::path filename_re = module_data_dir / "W_e" / (this->name + ".dat");
	std::filesystem::path filename_l_bar = module_data_dir / "W_l_bar" / (this->name + ".dat");
	this->rs_p_file.open(filename_rs, std::ios::in); // load polynomial coefficients
	this->re_p_file.open(filename_re, std::ios::in); // load polynomial coefficients
	this->l_bar_p_file.open(filename_l_bar, std::ios::in); // load polynomial coefficients
//	std::cout << filename_rs << std::endl;
//	std::cout << filename_re << std::endl;
//	std::cout << filename_l_bar << std::endl;
	for(size_t i=0; i<13; ++i)
	{
		rs_p_file >> this->rs_p_vec[i];
		re_p_file >> this->re_p_vec[i];
//		std::cout << this->p_vec[i] << std::endl;
	}

	rs_p_file.close();
	re_p_file.close();

	for(size_t i=0; i<64; ++i)
	{
		l_bar_p_file >> this->l_bar_p_vec[i];
		// std::cout << l_bar_p_vec[i] << std::endl; works fine
	}

	l_bar_p_file.close();
};

double Thelen2003MuscleGazebo::comp_gamma_a(double l_n)
{
	return std::exp ( -( pow((l_n - 1.0), 2) ) / this->g ); // gaussian shape centered around 1
}

double Thelen2003MuscleGazebo::comp_gamma_p(double l_n)
{
	double gamma_p = 0;

	if(l_n <= 1.0)
	{
		gamma_p = 0.0; // if muscle shorter than l_0 no elastic force
	}
	else
	{
		gamma_p = ( std::exp( this->k*(l_n - 1.0) / this->e_0m) - 1.0) / (std::exp(this->k) - 1.0);
		// else muscle generates elastic force modeled by exponential
	}

	return gamma_p;
}

double Thelen2003MuscleGazebo::comp_gamma_s(double eps_t)
{

	double gamma_s = 0;

	if(eps_t <= 0.0)
	{
		gamma_s = 0.0; // if tendon strain below zero (no stretch) no force generated
	}
	else if(eps_t <= this->eps_th)
	{
		gamma_s = this->a1 * (std::exp(this->a3*eps_t/this->eps_th) - 1.0) / (std::exp(this->a3) - 1.0);
		// slightly stretched causes exponential increase in tendon gain generated
	}
	else
	{
		gamma_s = this->a2*(eps_t - this->eps_th) + this->a1;
		// above exponential phase tendon gain increases linearly with stretch
	}

	return gamma_s;
}

double Thelen2003MuscleGazebo::comp_alpha(double l)
{
	double alpha = acos(sqrt(1 - pow((this->l_0*sin(this->alpha)/l), 2)));
	// pennation angle of muscle and tendon elements
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

	res.q_dot = Thelen2003MuscleGazebo::comp_q_dot(q, u); // directly assign q_dot

	const double l_n = l / l_0; // normalized by optimal length l_0

	// passive and active force gains
	const double gamma_p = Thelen2003MuscleGazebo::comp_gamma_p(l_n);
	const double gamma_a = Thelen2003MuscleGazebo::comp_gamma_a(l_n);

	//series force gain
	double c_alpha = cos(comp_alpha(l));
	double eps_t = (l_bar - l*c_alpha - l_t0) / l_t0;
	double gamma_s = Thelen2003MuscleGazebo::comp_gamma_s(eps_t);

	// compute contractile element force gain
	double gamma_c = gamma_s / c_alpha - gamma_p;

	// constrain length dynamics
	if(q < q_min)
	{
		q = q_min;
	}
	else if(q > q_max)
	{
		q = 1.0;
	}

	// compute normalized length dynamics
	double ln_dot = Thelen2003MuscleGazebo::comp_ln_dot(gamma_a,
	                                                    gamma_c,
	                                                    q);

	double l_dot = p8 * ln_dot * l_0; // unnormalize

	res.f = gamma_s * fmax;
	res.lm_dot = l_dot;

	return res;
}

double Thelen2003MuscleGazebo::comp_q_dot(double q, double u)
{
	double q_dot = 0.0;
	double tau = 0.0;

	if(u>q)
	{
		tau = this->tau_a * (0.5 + 1.5*q);
	}
	else
	{
		tau = this->tau_d / (0.5 + 1.5*q);
	}

	q_dot = (u - q)/tau;

	return q_dot;
}

double Thelen2003MuscleGazebo::comp_ln_dot(double gamma_a, double gamma_c, double q)
{
	double ln_dot = 0.0;

	if( (gamma_c > 0.0) && (gamma_c < q*gamma_a*this->p3*this->p4)) // tendon stretches elastic element of muscle
	{
		if(gamma_c <= q*gamma_a) // muscle contracts since active length dependent part are greater
		{
			ln_dot = (this->p1 + this->p2*q)*(gamma_c - q*gamma_a)
			        / (q*gamma_a + gamma_c / this->p5);
		}
		else // muscle extends
		{
			ln_dot = (this->p1 + this->p2 * q)*( gamma_c - q*gamma_a)
			        * this->p6 / (this->p7*(q*gamma_a*this->p3 - gamma_c));
		}
	}
	else
	{
		if(gamma_c <= 0.0) // muscle contracts
		{
			ln_dot = (this->p1 + this->p2*q)*( (1.0 + 1.0/this->p5) * gamma_c - q*gamma_a) / (q*gamma_a);
		}
		else // muscle extends fast
		{		// c1-c5: simplification parameters
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
{	// mainly used to check if the muscle is balanced within itself
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

double Thelen2003MuscleGazebo::compute_rs(double theta_s)
{

	double r = 0.0;

	for(int i=0; i<13; ++i)
	{
		r += rs_p_vec[i] * std::pow(theta_s, i);

	}

	return r;
}

double Thelen2003MuscleGazebo::compute_re(double theta_e)
{

	double r = 0.0;

	for(int i=0; i<13; ++i)
	{
		r += re_p_vec[i] * std::pow(theta_e, i);

	}

	return r;
}

double Thelen2003MuscleGazebo::compute_l_bar(double theta_s, double theta_e)
{

	std::array<double, 64> phi;
	double n = 0;
	double a = 0.0;
	double b = 0.0;

	//std::cout << theta_s << std::endl;
	for(int i=0; i<8; ++i)
	{
		a = std::pow(theta_s, i);

		for(int j=0; j<8; ++j)
		{
			b = std::pow(theta_e, j);
			phi[n] = a*b;
			n += 1;
		}
	}

	double l_bar = 0.0;

	for(int i=0; i<64; ++i)
	{
		l_bar += l_bar_p_vec[i] * phi[i];

	}

	return l_bar;

}


std::array<double, 2> rotate_p_urh(double p1, double p2, double theta)
{
	// rotate point and add offset to get global point
	double p1g = std::cos(theta) * p1 - std::sin(theta) * p2 -0.011445;
	double p2g = std::sin(theta) * p1 + std::cos(theta) * p2 -0.2974;

	std::array<double, 2> p_rot;

	p_rot.at(0) = p1g;
	p_rot.at(1) = p2g;

	return p_rot;
}

double comp_d(double p11, double p12, double p21, double p22)
{
	return std::sqrt( std::pow((p21 - p11), 2) + std::pow((p22 - p12), 2) );
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
	std::array<double, 2> p5g = rotate_p_urh(p51l, p52l, theta);

	// compute distance between each pair of points that form a direct connection globally
	double d1 = comp_d(p11, p12, p21, p22);
	double d2 = comp_d(p21, p22, p31, p32);
	double d3 = comp_d(p31, p32, p41, p42);
	double d4 = comp_d(p41, p42, p5g.at(0), p5g.at(1));

	double total_d = d1 + d2 + d3 + d4;

	return total_d;
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
	std::array<double, 2> p5g = rotate_p_urh(p51l, p52l, theta);

	// compute distance between each pair of points that form a direct connection globally
	double d1 = comp_d(p11, p12, p21, p22);
	double d2 = comp_d(p21, p22, p31, p32);
	double d3 = comp_d(p31, p32, p41, p42);
	double d4 = comp_d(p41, p42, p5g.at(0), p5g.at(1));

	double total_d = d1 + d2 + d3 + d4;

	return total_d;
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
	std::array<double, 2> p5g = rotate_p_urh(p51l, p52l, theta);

	// compute distance between each pair of points that form a direct connection globally
	double d1 = comp_d(p11, p12, p21, p22);
	double d2 = comp_d(p21, p22, p31, p32);
	double d3 = comp_d(p31, p32, p41, p42);
	double d4 = comp_d(p41, p42, p5g.at(0), p5g.at(1));

	double total_d = d1 + d2 + d3 + d4;

	return total_d;
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
	std::array<double, 2> p9g = rotate_p_urh(p91l, p92l, theta);

	// compute distance between each pair of points that form a direct connection globally
	double d1 = comp_d(p11, p12, p21, p22);
	double d2 = comp_d(p21, p22, p31, p32);
	double d3 = comp_d(p31, p32, p41, p42);
	double d4 = comp_d(p41, p42, p51, p52);

	double d5 = comp_d(p51, p52, p61, p62);
	double d6 = comp_d(p61, p62, p71, p72);
	double d7 = comp_d(p71, p72, p81, p82);
	double d8 = comp_d(p81, p82, p9g.at(0), p9g.at(1));


	double total_d = d1 + d2 + d3 + d4 + d5 + d6 + d7 + d8;

	return total_d;
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
	std::array<double, 2> p6g = rotate_p_urh(p61l, p62l, theta);

	// compute distance between each pair of points that form a direct connection globally
	double d1 = comp_d(p11, p12, p21, p22);
	double d2 = comp_d(p21, p22, p31, p32);
	double d3 = comp_d(p31, p32, p41, p42);
	double d4 = comp_d(p41, p42, p51, p52);
	double d5 = comp_d(p51, p52, p6g.at(0), p6g.at(1));

	double total_d = d1 + d2 + d3 + d4 + d5;

	return total_d;
}

double comp_l_bar_BRA(double theta)
{
	double p11 = 0.0068 -0.017545; // adding ground offset
	double p12 = -0.1739 -0.007;

	double p21l = -0.0032;
	double p22l = -0.0239;

	// rotate point and add offset to get global point
	std::array<double, 2> p2g = rotate_p_urh(p21l, p22l, theta);

	// compute distance between each pair of points that form a direct connection globally
	double d1 = comp_d(p11, p12, p2g.at(0), p2g.at(1));

	double total_d = d1;

	return total_d;
}


