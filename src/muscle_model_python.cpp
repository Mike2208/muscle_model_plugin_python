#include "muscle_model_plugin/Thelen2003MuscleGazebo.h"

#include "muscle_model_plugin/muscle_model_python_config.h"

#include <boost/python.hpp>
#include <boost/numeric/odeint.hpp>

namespace bp = boost::python;


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


class TwoMuscleModel
{
	public:
		enum MUSCLE_INDEX
		{
			MINDEX_DELTI     = 0,
			MINDEX_TRILONG   = 1,
			MINDEX_TRILAT    = 2,
			MINDEX_TRIMED    = 3,
			MINDEX_BICLONG   = 4,
			MINDEX_BICSHORT  = 5,
			MINDEX_BRA       = 6,

		};

		enum JOINT_INDEX
		{
			JINDEX_SHOULDER_0 = 0,
			JINDEX_ELBOW = 1,
		};

		enum LINK_INDEX
		{
			LINDEX_UPPER_ARM = 1,
			LINDEX_FOREARM = 2,
		};

		TwoMuscleModel()
		{
			int i = 0;
			for(const auto &muscle : this->_muscles)
			{
				this->_states[2*i + 0] = muscle.q;
				this->_states[2*i + 1] = muscle.l;

				++i;
			}
		}

		void Integrate(double dt)
		{
			boost::numeric::odeint::integrate_adaptive(this->_stepper,
			                                           std::bind(&TwoMuscleModel::IntegrationFcn, this,
			                                                     std::placeholders::_1, std::placeholders::_2,
			                                                     std::placeholders::_3),
			                                           this->_states,
			                                           this->_t, this->_t+dt, 0.0001
			                                           /*push_back_state_and_time(x_vec, times)*/);

			this->_t += dt;

			const auto theta_s = this->_states[14];
			const auto theta_e = this->_states[16];

			// Clear torques
			memset(this->_torques.data(), 0, sizeof(this->_torques));

			int i = 0;

			const std::array<Thelen2003MuscleGazebo*,2> target_muscles{&this->_muscles[3], &this->_muscles[6]};
			for(auto &muscle : target_muscles)
			{
				double l_bar = muscle->compute_l_bar(theta_s, theta_e);
				Thelen2003MuscleGazebo::ReturnValue ret = muscle->compute_muscle_dynamics(this->_states[i * 2 + 0], this->_stimulus[i], this->_states[i*2 + 1], l_bar);
				double rs = muscle->compute_rs(theta_s);
				double re = muscle->compute_re(theta_e);

				this->_torques[JINDEX_SHOULDER_0] = this->_torques[JINDEX_SHOULDER_0] + rs * ret.f;
				this->_torques[JINDEX_ELBOW] = this->_torques[JINDEX_ELBOW] + re * ret.f;

				++i;
			}
		}

		double GetTorque(JOINT_INDEX joint_index) const
		{
			return this->_torques.at(joint_index);
		}

		void SetJointPos(JOINT_INDEX joint_index, double rad)
		{
			switch(joint_index)
			{
				case JINDEX_SHOULDER_0: this->_states[14] = rad; break;
				case JINDEX_ELBOW:      this->_states[16] = rad; break;
				default:                break;
			}
		}

		void SetJointVel(JOINT_INDEX joint_index, double rad)
		{
			switch(joint_index)
			{
				case JINDEX_SHOULDER_0: this->_states[15] = rad; break;
				case JINDEX_ELBOW:      this->_states[17] = rad; break;
				default:                break;
			}
		}

		void SetStimulus(MUSCLE_INDEX muscle_index, double stimulus)
		{
			this->_stimulus.at(muscle_index) = stimulus;
		}

	private:
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

		std::array<double, 2> _torques = {0,0};
		std::array<double, 7> _stimulus = {0.15,0.01,0.01,0.01,0.15,0.15,0.15};

		using stepper_t = decltype(boost::numeric::odeint::make_controlled(1.E-12, 1.E-12, boost::numeric::odeint::runge_kutta_dopri5< std::array<double,18> >()));
		stepper_t _stepper = boost::numeric::odeint::make_controlled(1.E-12, 1.E-12, boost::numeric::odeint::runge_kutta_dopri5< std::array<double,18> >());

		double _t = 0.0;

		void IntegrationFcn(const std::array<double,18> &x , std::array<double,18> &dxdt , const double /*t*/)
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
		//	if(((0.05 < t) && (t < 0.14)) || ((0.5 < t) && (t < 0.53)))
		//		this->_stimulus = {0.15, 0.01, 0.01, 0.01, 0.15, 0.15, 0.15};
		//	else
		//		this->_stimulus = {0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};

			const std::array<Thelen2003MuscleGazebo*,2> target_muscles{&this->_muscles[3], &this->_muscles[6]};

			int i = 0;
			for(auto &muscle : target_muscles)
			{
				double l_bar = muscle->compute_l_bar(theta_s, theta_e);
				//std::cout << l_bar << std::endl;
				Thelen2003MuscleGazebo::ReturnValue ret = muscle->compute_muscle_dynamics(x[2*i + 0], this->_stimulus[i], x[2*i + 1], l_bar);
				dxdt[2*i + 0] = ret.q_dot;
				dxdt[2*i + 1] = ret.lm_dot;

		//		double r = fcnr[i](theta);
				double rs = muscle->compute_rs(theta_s);
				double re = muscle->compute_re(theta_e);

				tau_s = tau_s + rs * ret.f;
				tau_e = tau_e + re * ret.f;
				++i;
			}

//			double theta_s_ddot = comp_skeletal_dynamics(theta_s) + tau_s / (1.86 * std::pow(0.24, 2));
//			double theta_e_ddot = comp_skeletal_dynamics(theta_e) + tau_e / (1.53 * std::pow(0.24, 2)); // incorrect inertia but sufficient
			double theta_s_ddot = tau_s / (1.86 * std::pow(0.24, 2));
			double theta_e_ddot = tau_e / (1.53 * std::pow(0.24, 2)); // incorrect inertia but sufficient

			dxdt[14] = theta_s_dot;
			dxdt[15] = theta_s_ddot;
			dxdt[16] = theta_e_dot;
			dxdt[17] = theta_e_ddot;
		}
};

BOOST_PYTHON_MODULE(MuscleModelPython)
{
//	enum MUSCLE_INDEX
//	{
//		MINDEX_DELTI     = 0,
//		MINDEX_TRILONG   = 1,
//		MINDEX_TRILAT    = 2,
//		MINDEX_TRIMED    = 3,
//		MINDEX_BICLONG   = 4,
//		MINDEX_BICSHORT  = 5,
//		MINDEX_BRA       = 6,

//	};

//	enum JOINT_INDEX
//	{
//		JINDEX_SHOULDER_0 = 0,
//		JINDEX_ELBOW = 1,
//	};

	bp::enum_<TwoMuscleModel::MUSCLE_INDEX>("MUSCLE_INDEX")
	        .value("MINDEX_DELTI", TwoMuscleModel::MINDEX_DELTI)
	        .value("MINDEX_TRILONG", TwoMuscleModel::MINDEX_TRILONG)
	        .value("MINDEX_TRILAT", TwoMuscleModel::MINDEX_TRILAT)
	        .value("MINDEX_TRIMED", TwoMuscleModel::MINDEX_TRIMED)
	        .value("MINDEX_BICLONG", TwoMuscleModel::MINDEX_BICLONG)
	        .value("MINDEX_BICSHORT", TwoMuscleModel::MINDEX_BICSHORT)
	        .value("MINDEX_BRA", TwoMuscleModel::MINDEX_BRA);

	bp::enum_<TwoMuscleModel::JOINT_INDEX>("JOINT_INDEX")
	        .value("JINDEX_SHOULDER_0", TwoMuscleModel::JINDEX_SHOULDER_0)
	        .value("JINDEX_ELBOW", TwoMuscleModel::JINDEX_ELBOW);

	bp::class_<TwoMuscleModel, boost::noncopyable>("TwoMuscleModel")
	        .def("Integrate", &TwoMuscleModel::Integrate)
	        .def("GetTorque", &TwoMuscleModel::GetTorque)
	        .def("SetJointPos", &TwoMuscleModel::SetJointPos)
	        .def("SetJointVel", &TwoMuscleModel::SetJointVel)
	        .def("SetStimulus", &TwoMuscleModel::SetStimulus);
}
