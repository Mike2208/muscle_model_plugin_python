#include "muscle_model_plugin/muscle_model_plugin.h"

#include <gazebo/physics/World.hh>
#include <gazebo/physics/Model.hh>
#include <gazebo/physics/Link.hh>
#include <gazebo/physics/Joint.hh>
#include <sdf/sdf.hh>
#include <math.h>

#include "muscle_model_plugin/Thelen2003MuscleGazebo.h"

void MuscleModelPlugin::Load(gazebo::physics::ModelPtr model, sdf::ElementPtr sdf)
{
	// Create stepper

	this->_model = model;
	this->_onWorldUpdate = gazebo::event::Events::ConnectWorldUpdateBegin(std::bind(&MuscleModelPlugin::OnWorldUpdateBegin, this));

	this->_muscles.clear();
	this->_muscles.emplace_back(Thelen2003MuscleGazebo("DELT1", 1218.9, 0.0976, 0.0930, 0.3840, 0.1149, 0.1496));
	this->_muscles.emplace_back(Thelen2003MuscleGazebo("TRIlong", 771.8, 0.1340, 0.1430, 0.2094, 0.1336, 0.848));
	this->_muscles.emplace_back(Thelen2003MuscleGazebo("TRIlat", 717.5, 0.1138, 0.0980, 0.1571, 0.0712, 0.0223));
	this->_muscles.emplace_back(Thelen2003MuscleGazebo("TRImed", 717.5, 0.1138, 0.0908, 0.1571, 0.0662, 0.0239));
	this->_muscles.emplace_back(Thelen2003MuscleGazebo("BIClong", 525.1, 0.1157, 0.2723, 0, 0.1556, 0.0617));
	this->_muscles.emplace_back(Thelen2003MuscleGazebo("BICshort", 316.8, 0.1321, 0.1923, 0, 0.1565, 0.1383));
	this->_muscles.emplace_back(Thelen2003MuscleGazebo("BRA", 1177.4, 0.0858, 0.0535, 0, 0.0867, 0.9767));

//	this->_muscles = {
//	    Thelen2003MuscleGazebo("TRIlong", 798.52, 0.1340, 0.1430, 0.2094, 0.1362, 0.0399),
//	    Thelen2003MuscleGazebo("TRIlat", 624.3, 0.1138, 0.0980, 0.1571, 0.0702, 0.0399),
//	    Thelen2003MuscleGazebo("TRImed", 624.3, 0.1138, 0.0908, 0.1571, 0.0653, 0.0400),
//	    Thelen2003MuscleGazebo("BIClong", 624.3, 0.1157, 0.2723, 0, 0.1501, 0.0409),
//	    Thelen2003MuscleGazebo("BICshort", 435.56, 0.1321, 0.1923, 0, 0.1473, 0.0404),
//	    Thelen2003MuscleGazebo("BRA", 987.26, 0.0858, 0.0535, 0, 0.0979, 0.0399)
//	};

	int i = 0;
	for(const auto &muscle : this->_muscles)
	{
		this->_states[2*i + 0] = muscle.q;
		this->_states[2*i + 1] = muscle.l;

		++i;
	}

	this->_curTime = this->_model->GetWorld()->SimTime();

	// Integration function
//	boost::numeric::odeint::integrate_adaptive(this->_stepper,
//	                                           std::bind(&MuscleModelPlugin::IntegrationFcn, this,
//	                                                     std::placeholders::_1, std::placeholders::_2,
//	                                                     std::placeholders::_3),
//	                                           this->_states,
//	                                           0.0, 0.001, 0.0001);
}


void MuscleModelPlugin::IntegrationFcn(const std::array<double, 14> &x, std::array<double, 14> &dxdt, const double t)
{

	auto joint_s = this->_model->GetJoint("arm_r_shoulder_joint");
	auto joint_e = this->_model->GetJoint("arm_r_elbow_joint");
	auto theta_s = joint_s->Position(0);
	auto theta_e = joint_e->Position(0);

	int i = 0;
	for(auto &muscle : this->_muscles)
	{
		double l_bar = muscle.compute_l_bar(theta_s, theta_e);
		Thelen2003MuscleGazebo::ReturnValue ret = muscle.compute_muscle_dynamics(x[2*i + 0], this->_stimulus[i], x[2*i + 1], l_bar);
		dxdt[2*i + 0] = ret.q_dot;
		dxdt[2*i + 1] = ret.lm_dot;

		++i;
	}

}

void MuscleModelPlugin::OnWorldUpdateBegin()
{
	//ignition::math::Pose3d pose = this->_model->WorldPose();
	//this->_model->WorldLinearVel();
	//this->_model->WorldAngularVel();
//	auto joint = this->_model->GetJoint("mixamorig_RightForeArm");
//	auto rad = joint->Position(0);

	const auto newTime = this->_model->GetWorld()->SimTime();

	// Integration function
	boost::numeric::odeint::integrate_adaptive(this->_stepper,
	                                           std::bind(&MuscleModelPlugin::IntegrationFcn, this,
	                                                     std::placeholders::_1, std::placeholders::_2,
	                                                     std::placeholders::_3),
	                                           this->_states,
	                                           this->_curTime.Double(), newTime.Double(), 0.0001);

	this->_curTime = newTime;



	auto joint_s = this->_model->GetJoint("arm_r_shoulder_joint");
	auto joint_e = this->_model->GetJoint("arm_r_elbow_joint");
	auto theta_s = joint_s->Position(0);
	auto theta_e = joint_e->Position(0);

	if(theta_e < 1.57)
		this->_stimulus = {0.15, 0.01, 0.01, 0.01, 0.15, 0.15, 0.15};
	else
		this->_stimulus = {0.01, 0.15, 0.15, 0.15, 0.01, 0.01, 0.01};

	int i = 0;
	double tau_s = 0;
	double tau_e = 0;
	for(auto &muscle : this->_muscles)
	{
		double l_bar = muscle.compute_l_bar(theta_s, theta_e);
		Thelen2003MuscleGazebo::ReturnValue ret = muscle.compute_muscle_dynamics(this->_states[i * 2 + 0], this->_stimulus[i], this->_states[i*2 + 1], l_bar);
		double rs = muscle.compute_rs(theta_s);
		double re = muscle.compute_re(theta_e);

		tau_s = tau_s + rs * ret.f;
		tau_e = tau_e + re * ret.f;

		++i;
	}


	joint_s->SetForce(0, tau_s);
	joint_e->SetForce(0, tau_e);
}

GZ_REGISTER_MODEL_PLUGIN(MuscleModelPlugin);
