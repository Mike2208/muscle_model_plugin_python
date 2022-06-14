#include "muscle_model_plugin/Thelen2003MuscleGazebo.h"
//#include "muscle_model_plugin/muscle_model_plugin.h"

#include <boost/numeric/odeint.hpp>

#include <iostream>
using std::endl;
#include <fstream>
#include <cstdlib> // for exit function

/**
 * @brief comp_skeletal_dynamics
 * @param theta
 * @return
 */
double comp_skeletal_dynamics(double theta)
{

	double gr = 9.81;
	double l2 = 0.24;

	double g = - (gr / l2) * std::sin(theta);

	return g;
}


class test_node_two_muscles_class
{
	public:
		test_node_two_muscles_class();

		std::array<Thelen2003MuscleGazebo,7> _muscles = {
		    Thelen2003MuscleGazebo("DELT1", 1218.9, 0.0976, 0.0930, 0.3840, 0.1149, 0.1496),
		    Thelen2003MuscleGazebo("TRIlong", 771.8, 0.1340, 0.1430, 0.2094, 0.1336, 0.848),
		    Thelen2003MuscleGazebo("TRIlat", 717.5, 0.1138, 0.0980, 0.1571, 0.0712, 0.0223),
		    Thelen2003MuscleGazebo("TRImed", 717.5, 0.1138, 0.0908, 0.1571, 0.0662, 0.0239),
		    Thelen2003MuscleGazebo("BIClong", 525.1, 0.1157, 0.2723, 0, 0.1556, 0.0617),
		    Thelen2003MuscleGazebo("BICshort", 316.8, 0.1321, 0.1923, 0, 0.1565, 0.1383),
		    Thelen2003MuscleGazebo("BRA", 1177.4, 0.0858, 0.0535, 0, 0.0867, 0.9767)
		};

		std::array<double, 18> _states = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

		std::fstream thetasf;
		std::fstream tspanf;
		std::fstream stimulusf;

		std::array<double, 100> thetas;
		std::array<double, 100> tspan;
		std::array<double, 7> _stimulus = {0.15,0.01,0.01,0.01,0.15,0.15,0.15};

		using stepper_t = decltype(boost::numeric::odeint::make_controlled(1.E-12, 1.E-12, boost::numeric::odeint::runge_kutta_dopri5< std::array<double,18> >()));
		stepper_t _stepper = boost::numeric::odeint::make_controlled(1.E-12, 1.E-12, boost::numeric::odeint::runge_kutta_dopri5< std::array<double,18> >());


		void IntegrationFcn( const std::array<double,18> &x , std::array<double,18> &dxdt , const double t );

		struct ReturnValue
		{
			std::vector<std::array<double,18> > x_vec;
			std::vector<double> times;
		};

		ReturnValue simulate_episode(double t0, double t1);
		int inspect_muscle_model();

};

test_node_two_muscles_class::test_node_two_muscles_class()
{
	int i = 0;
	for(const auto &muscle : this->_muscles)
	{
		this->_states[2*i + 0] = muscle.q;
		this->_states[2*i + 1] = muscle.l;

		++i;
	}

	this->thetasf.open("/home/devel/Downloads/sim_data/thetas.dat", std::ios::in);
	this->tspanf.open("/home/devel/Downloads/sim_data/tspan.dat", std::ios::in);
	this->stimulusf.open("/home/devel/Downloads/sim_data/stimulus.dat", std::ios::in);

	for(size_t i=0; i<100; ++i)
	{
		std::string line;
		line.resize(100);

		this->thetasf.getline(line.data(), 100);
		this->thetas[i] = std::atof(line.data());
		//this->thetasf >> this->thetas[i];

		this->tspanf.getline(line.data(), 100);
		this->tspan[i] = std::atof(line.data());
		//this->tspanf >> this->tspan[i];
		//std::cout << this->tspan[i] << std::endl;
	}

//	for(size_t i=0; i<100; ++i)
//	{
//		for(size_t j=0; j<6; ++j)
//		{
//			this->stimulusf >> this->stimulus.at(i).at(j);
//			std::cout << this->stimulus.at(i).at(j) << std::endl;
//		}
//	}


}

void test_node_two_muscles_class::IntegrationFcn(const std::array<double, 18> &x, std::array<double, 18> &dxdt, const double t)
{

	double theta_s = x[14];
	//std::cout << theta_s << std::endl;
	double theta_s_dot = x[15];
	double theta_e = x[16];
	double theta_e_dot = x[17];


	double tau_s = 0.0;
	double tau_e = 0.0;

//	if(theta < 1.57)
//		this->_stimulus = {0.01, 0.01, 0.01, 0.15, 0.15, 0.15};
//	else
//		this->_stimulus = {0.15, 0.15, 0.15, 0.01, 0.01, 0.01};
	if(((0.05 < t) && (t < 0.14)) || ((0.5 < t) && (t < 0.53)))
		this->_stimulus = {0.15, 0.01, 0.01, 0.01, 0.15, 0.15, 0.15};
	else
		this->_stimulus = {0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};

	int i = 0;
	for(auto &muscle : this->_muscles)
	{
		double l_bar = muscle.compute_l_bar(theta_s, theta_e);
		std::cout << l_bar << std::endl;
		Thelen2003MuscleGazebo::ReturnValue ret = muscle.compute_muscle_dynamics(x[2*i + 0], this->_stimulus[i], x[2*i + 1], l_bar);
		dxdt[2*i + 0] = ret.q_dot;
		dxdt[2*i + 1] = ret.lm_dot;

//		double r = fcnr[i](theta);
		double rs = muscle.compute_rs(theta_s);
		double re = muscle.compute_re(theta_e);

		tau_s = tau_s + rs * ret.f;
		tau_e = tau_e + re * ret.f;
		++i;
	}

	double theta_s_ddot = comp_skeletal_dynamics(theta_s) + tau_s / (1.86 * std::pow(0.24, 2));
	double theta_e_ddot = comp_skeletal_dynamics(theta_e) + tau_e / (1.53 * std::pow(0.24, 2)); // incorrect inertia but sufficient

	dxdt[14] = theta_s_dot;
	dxdt[15] = theta_s_ddot;
	dxdt[16] = theta_e_dot;
	dxdt[17] = theta_e_ddot;

}

//[ integrate_observer
struct push_back_state_and_time
{
		std::vector< std::array<double,18> >& m_states;
	std::vector< double >& m_times;

	push_back_state_and_time( std::vector< std::array<double,18> > &states , std::vector< double > &times )
	: m_states( states ) , m_times( times ) { }

	void operator()( const std::array<double,18> &x , double t )
	{
		m_states.push_back( x );
		m_times.push_back( t );
	}
};
//]

test_node_two_muscles_class::ReturnValue test_node_two_muscles_class::simulate_episode(double t0, double t1)
{
	ReturnValue ret;

	std::vector<std::array<double,18> > x_vec;
	std::vector<double> times;

	boost::numeric::odeint::integrate_adaptive(this->_stepper,
	                                           std::bind(&test_node_two_muscles_class::IntegrationFcn, this,
	                                                     std::placeholders::_1, std::placeholders::_2,
	                                                     std::placeholders::_3),
	                                           this->_states,
	                                           t0, t1, 0.0001,
	                                           push_back_state_and_time(x_vec, times));


	ret.times = times;
	ret.x_vec = x_vec;

	return ret;
}

int test_node_two_muscles_class::inspect_muscle_model()
{

	Thelen2003MuscleGazebo muscle("BRA", 987.0, 0.0858, 0.0535, 0.0, 0.0979, 0.0399); // use parameters for BRA to yield an example muscle
	std::cout << muscle.compute_muscle_equilibrium(0.0979, 0.0399, 0.1405) << std::endl;

	using fcn_t = double(*)(double);
	fcn_t fcnl[6] = {&comp_l_bar_TRIlong, &comp_l_bar_TRIlat, &comp_l_bar_TRImed, &comp_l_bar_BIClong, &comp_l_bar_BICshort, &comp_l_bar_BRA};

	std::string muscle_names[6] = {"TRIlong", "TRIlat", "TRImed", "BIClong", "BICshort", "BRA"};

	// map out moment arm and l_bar computation for each muscle
	for(int m=0; m<6; ++m)
	{
		std::fstream  lbardata("/home/devel/Thelen2003DebugData/l_bar/"+muscle_names[m]+".dat", std::ios_base::out | std::ios_base::trunc);
		std::fstream  rdata("/home/devel/Thelen2003DebugData/r/"+muscle_names[m]+".dat", std::ios_base::out | std::ios_base::trunc);

		for(int i=0; i<100; ++i)
		{
			lbardata << fcnl[m](0.0314*i) << std::endl;
			rdata << _muscles[m].compute_re(i*0.0314/1.3) << std::endl;

		}
		lbardata.close();
		rdata.close();
	}

	std::fstream  gamma_a("/home/devel/Thelen2003DebugData/BRA/gamma_a.dat", std::ios_base::out | std::ios_base::trunc);
	std::fstream  gamma_p("/home/devel/Thelen2003DebugData/BRA/gamma_p.dat", std::ios_base::out | std::ios_base::trunc);
	std::fstream  gamma_s("/home/devel/Thelen2003DebugData/BRA/gamma_s.dat", std::ios_base::out | std::ios_base::trunc);
	std::fstream  ln_dot_data("/home/devel/Thelen2003DebugData/BRA/ln_dot.dat", std::ios_base::out | std::ios_base::trunc);
	std::fstream  q_dot("/home/devel/Thelen2003DebugData/BRA/q_dot.dat", std::ios_base::out | std::ios_base::trunc);
	std::fstream  q_dot_input("/home/devel/Thelen2003DebugData/BRA/q_dot_input.dat", std::ios_base::out | std::ios_base::trunc);
	std::fstream  l_input("/home/devel/Thelen2003DebugData/BRA/l_input.dat", std::ios_base::out | std::ios_base::trunc);
	std::fstream  l_bar_input("/home/devel/Thelen2003DebugData/BRA/l_bar_input.dat", std::ios_base::out | std::ios_base::trunc);

	double l_0 = 0.0858;
	double l_t0 = 0.0535;
	int i;	int j;
	for(i=1; i<=200; ++i)
	{
		gamma_a << muscle.comp_gamma_a(0.01*i) << endl;
		gamma_p << muscle.comp_gamma_p(0.01*i) << endl;
		gamma_s << muscle.comp_gamma_s(0.001*i) << endl;

		double l_n = 1.0 - (100.0 - i) / 5000.0;
		double l = l_n * l_0;	double q = 1.0;
//		std::cout << l_n << endl;
//		std::cout << l << endl;
		l_input << l << endl;
		for(j=1; j<=200; ++j)
		{
			double l_bar = (1.0 - (100.0 - j) / 5000.0) * (l_0 + l_t0);
			double c_alpha = std::cos(muscle.comp_alpha(l));
			double eps_t = (l_bar - l*c_alpha - l_t0) / l_t0;

			double gamma_c = muscle.comp_gamma_s(eps_t) / c_alpha - muscle.comp_gamma_p(l_n);

			double ln_dot = muscle.comp_ln_dot(muscle.comp_gamma_a(l_n),
			                                   gamma_c,
			                                   q);
			if(i==200)
				l_bar_input << l_bar << endl;

			ln_dot_data << ln_dot << endl;
			//map out q_dot
			q_dot << muscle.comp_q_dot(0.005*i, 0.005*j) << endl;
			q_dot_input << 0.005*i * 0.005*j << endl;
//			std::cout << 0.005*i * 0.005*j << endl;
		}
	}
	gamma_a.close();
	gamma_p.close();
	gamma_s.close();
	ln_dot_data.close();
	q_dot.close();
	q_dot_input.close();
	l_input.close();
	l_bar_input.close();

	return 0;
}


int main(int argc, char **argv)
{


	test_node_two_muscles_class sim_tester = test_node_two_muscles_class();
	//sim_tester.inspect_muscle_model();

	std::fstream  q_sim("/home/devel/Thelen2003DebugData/sim_debug/q.dat", std::ios_base::out | std::ios_base::trunc);
	std::fstream  l_sim("/home/devel/Thelen2003DebugData/sim_debug/l.dat", std::ios_base::out | std::ios_base::trunc);
	std::fstream  theta_sim("/home/devel/Thelen2003DebugData/sim_debug/theta.dat", std::ios_base::out | std::ios_base::trunc);
	std::fstream  theta_dot_sim("/home/devel/Thelen2003DebugData/sim_debug/theta_dot.dat", std::ios_base::out | std::ios_base::trunc);
	std::fstream  tspan("/home/devel/Thelen2003DebugData/sim_debug/tspan.dat", std::ios_base::out | std::ios_base::trunc);

	test_node_two_muscles_class::ReturnValue res = sim_tester.simulate_episode(0.0, 1.0);
	int vec_size = res.x_vec.size();

	for(int i=0; i<vec_size; ++i)
	{
		for(int m=0; m<7; ++m)
		{

			q_sim << res.x_vec.at(i).at(2*m + 0) << endl;
			l_sim << res.x_vec.at(i).at(2*m + 1) << endl;
		}
		theta_sim << res.x_vec.at(i).at(16) << endl;	// elbow joint data
		theta_dot_sim << res.x_vec.at(i).at(17) << endl;
		tspan << res.times.at(i) << endl;
	}

	q_sim.close();
	//f_sim.close();
	l_sim.close();
	theta_sim.close();
	theta_dot_sim.close();
	tspan.close();




}
