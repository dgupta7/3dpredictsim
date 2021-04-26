/*  This code describes the OpenSim model and the skeleton dynamics
Author: Antoine Falisse
Contributor: Joris Gillis, Gil Serrancoli, Chris Dembia
*/
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/SimbodyEngine/PlanarJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/PinJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/WeldJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/Joint.h>
#include <OpenSim/Simulation/SimbodyEngine/SpatialTransform.h>
#include <OpenSim/Simulation/SimbodyEngine/CustomJoint.h>
#include <OpenSim/Common/LinearFunction.h>
#include <OpenSim/Common/Constant.h>
#include <OpenSim/Common/SimmSpline.h>
#include <OpenSim/Simulation/Model/ConditionalPathPoint.h>
#include <OpenSim/Simulation/Model/MovingPathPoint.h>
#include <OpenSim/Simulation/Model/HuntCrossleyForce_smooth.h>
#include "SimTKcommon/internal/recorder.h"

#include <iostream>
#include <iterator>
#include <random>
#include <cassert>
#include <algorithm>
#include <vector>
#include <fstream>

using namespace SimTK;
using namespace OpenSim;

/*  The function F describes the OpenSim model and, implicitly, the skeleton
dynamics. F takes as inputs joint positions and velocities (states x),
joint accelerations (controls u), and returns the joint torques as well as
several variables for use in the optimal control problems. F is templatized
using type T. F(x,u)->(r).
*/

// Inputs/outputs of function F
/// number of vectors in inputs/outputs of function F
constexpr int n_in = 2;
constexpr int n_out = 1;
/// number of elements in input/output vectors of function F
constexpr int ndof = 14;           // # degrees of freedom (excluding locked)
constexpr int ndofr = ndof;    // # degrees of freedom (including locked)
constexpr int NX = ndof * 2;       // # states
constexpr int NU = ndof;           // # controls
constexpr int NR = ndof + 6 * 3 + 6 *3;   // 6 joint origins + 6 contact spheres

									   // Helper function value
template<typename T>
T value(const Recorder& e) { return e; }
template<>
double value(const Recorder& e) { return e.getValue(); }


// Function F
template<typename T>
int F_generic(const T** arg, T** res) {

	// OpenSim model: create components
	/// Model
	OpenSim::Model* model;
	/// Bodies
	OpenSim::Body* pelvis;
	OpenSim::Body* femur_r;
	OpenSim::Body* tibia_r;
	OpenSim::Body* talus_r;
	OpenSim::Body* calcn_r;
	OpenSim::Body* central_foot_r; // This is the fusion of metatarsi, cuboid, navicular and cuneiform bones and their surrounding tissue.
	OpenSim::Body* toes_r;
	/// Joints
	OpenSim::CustomJoint* ground_pelvis;
	OpenSim::CustomJoint* hip_r;
	OpenSim::CustomJoint* knee_r;
	OpenSim::CustomJoint* ankle_r;
	OpenSim::CustomJoint* subtalar_r;
	OpenSim::PinJoint* midtarsal_r;
	OpenSim::PinJoint* mtp_r;
	/// Contact elements
	OpenSim::HuntCrossleyForce_smooth* HC_1_r;
	OpenSim::HuntCrossleyForce_smooth* HC_2_r;
	OpenSim::HuntCrossleyForce_smooth* HC_3_r;
	OpenSim::HuntCrossleyForce_smooth* HC_4_r;
	OpenSim::HuntCrossleyForce_smooth* HC_5_r;
	OpenSim::HuntCrossleyForce_smooth* HC_6_r;

	// Inertia of the exoskeleton for right side
	Vec3 COMTibia = Vec3(0, -0.178389, 0);
	Vec3 COMTalus = Vec3(0, 0, 0);
	osim_double_adouble mTibia = 2.8735207629071704;
	osim_double_adouble mTalus = 0.07750561734071933;
	Inertia ITibia(0.0356624, 0.00360871, 0.0361578, 0, 0, 0);
	Inertia ITalus(0.000647371, 0.000647371, 0.000647371, 0, 0, 0);

	// using different way to scale foot model

	// OpenSim model: initialize components
	/// Model
	model = new OpenSim::Model();
	/// Body specifications
	pelvis = new OpenSim::Body("pelvis", 8.84259166189724, Vec3(-0.0682778, 0, 0), Inertia(0.0741799006400181, 0.0741799006400181, 0.0405455944309864, 0, 0, 0));
	femur_r = new OpenSim::Body("femur_r", 6.98382288222561, Vec3(0, -0.170467, 0), Inertia(0.101089610270247, 0.0264992182261812, 0.106600843690507, 0, 0, 0));
	tibia_r = new OpenSim::Body("tibia_r", mTibia, COMTibia, ITibia);
	talus_r = new OpenSim::Body("talus_r", mTalus, COMTalus, ITalus);
	calcn_r = new OpenSim::Body("calcn_r", 0.40857, Vec3(0.036283, 0.017101, -0.004057), Inertia(8.2986e-05, 0.0001696, 0.00019494, 0, 0, 0));
	central_foot_r = new OpenSim::Body("central_foot_r", 0.52998, Vec3(0.03937, -0.01465, 0.0018528), Inertia(0.0013620, 0.002204, 0.0019264, 0, 0, 0));
	toes_r = new OpenSim::Body("toes_r", 0.167877, Vec3(0.0316218, 0.00548355, -0.0159937), Inertia(6.2714-005, 0.0001254, 6.2714e-005, 0, 0, 0));
	/// Joints
	/// Ground-Pelvis transform
	SpatialTransform st_ground_pelvis;
	st_ground_pelvis[0].setCoordinateNames(OpenSim::Array<std::string>("pelvis_tilt", 1, 1));
	st_ground_pelvis[0].setFunction(new LinearFunction());
	st_ground_pelvis[0].setAxis(Vec3(0, 0, 1));
	st_ground_pelvis[1].setCoordinateNames(OpenSim::Array<std::string>("pelvis_list", 1, 1));
	st_ground_pelvis[1].setFunction(new LinearFunction());
	st_ground_pelvis[1].setAxis(Vec3(1, 0, 0));
	st_ground_pelvis[2].setCoordinateNames(OpenSim::Array<std::string>("pelvis_rotation", 1, 1));
	st_ground_pelvis[2].setFunction(new LinearFunction());
	st_ground_pelvis[2].setAxis(Vec3(0, 1, 0));
	st_ground_pelvis[3].setCoordinateNames(OpenSim::Array<std::string>("pelvis_tx", 1, 1));
	st_ground_pelvis[3].setFunction(new LinearFunction());
	st_ground_pelvis[3].setAxis(Vec3(1, 0, 0));
	st_ground_pelvis[4].setCoordinateNames(OpenSim::Array<std::string>("pelvis_ty", 1, 1));
	st_ground_pelvis[4].setFunction(new LinearFunction());
	st_ground_pelvis[4].setAxis(Vec3(0, 1, 0));
	st_ground_pelvis[5].setCoordinateNames(OpenSim::Array<std::string>("pelvis_tz", 1, 1));
	st_ground_pelvis[5].setFunction(new LinearFunction());
	st_ground_pelvis[5].setAxis(Vec3(0, 0, 1));
	/// Hip_r transform
	SpatialTransform st_hip_r;
	st_hip_r[0].setCoordinateNames(OpenSim::Array<std::string>("hip_flexion_r", 1, 1));
	st_hip_r[0].setFunction(new LinearFunction());
	st_hip_r[0].setAxis(Vec3(0, 0, 1));
	st_hip_r[1].setCoordinateNames(OpenSim::Array<std::string>("hip_adduction_r", 1, 1));
	st_hip_r[1].setFunction(new LinearFunction());
	st_hip_r[1].setAxis(Vec3(1, 0, 0));
	st_hip_r[2].setCoordinateNames(OpenSim::Array<std::string>("hip_rotation_r", 1, 1));
	st_hip_r[2].setFunction(new LinearFunction());
	st_hip_r[2].setAxis(Vec3(0, 1, 0));
	/// Knee_r transform
	SpatialTransform st_knee_r;
	st_knee_r[2].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_r", 1, 1));
	st_knee_r[2].setFunction(new LinearFunction());
	st_knee_r[2].setAxis(Vec3(0, 0, 1));
	/// Ankle_r transform
	SpatialTransform st_ankle_r;
	st_ankle_r[0].setCoordinateNames(OpenSim::Array<std::string>("ankle_angle_r", 1, 1));
	st_ankle_r[0].setFunction(new LinearFunction());
	st_ankle_r[0].setAxis(Vec3(-0.10501355, -0.17402245, 0.97912632));
	/// Subtalar_r transform
	SpatialTransform st_subtalar_r;
	st_subtalar_r[0].setCoordinateNames(OpenSim::Array<std::string>("subtalar_angle_r", 1, 1));
	st_subtalar_r[0].setFunction(new LinearFunction());
	st_subtalar_r[0].setAxis(Vec3(0.78717961, 0.60474746, -0.12094949));
	
	/// Joint specifications
	ground_pelvis = new CustomJoint("ground_pelvis", model->getGround(), Vec3(0), Vec3(0), *pelvis, Vec3(0), Vec3(0), st_ground_pelvis); // 6
	hip_r = new CustomJoint("hip_r", *pelvis, Vec3(-0.0682778001711179, -0.0638353973311301, 0.0823306940058688), Vec3(0), *femur_r, Vec3(0), Vec3(0), st_hip_r); // 3
	knee_r = new CustomJoint("knee_r", *femur_r, Vec3(-0.00451221232146798, -0.396907245921447, 0), Vec3(0), *tibia_r, Vec3(0), Vec3(0), st_knee_r); // 1
	ankle_r = new CustomJoint("ankle_r", *tibia_r, Vec3(0, -0.41085882914747662, 0), Vec3(0), *talus_r, Vec3(0), Vec3(0), st_ankle_r); // 1
	subtalar_r = new CustomJoint("subtalar_r", *talus_r, Vec3(-0.044572, -0.051452, 0.0072383), Vec3(0), *calcn_r, Vec3(0), Vec3(0), st_subtalar_r); // 1
	midtarsal_r = new PinJoint("midtarsal_r", *calcn_r, Vec3(0.078024, 0.025452, -0.0068103), Vec3(0), *central_foot_r, Vec3(0), Vec3(0)); // 1
	mtp_r = new PinJoint("mtp_r", *central_foot_r, Vec3(0.085386, -0.02728, 0.0077974), Vec3(0), *toes_r, Vec3(0), Vec3(0)); // 1

	/// Add bodies and joints to model
	model->addBody(pelvis);		    model->addJoint(ground_pelvis);
	model->addBody(femur_r);		model->addJoint(hip_r);
	model->addBody(tibia_r);		model->addJoint(knee_r);
	model->addBody(talus_r);		model->addJoint(ankle_r);
	model->addBody(calcn_r);		model->addJoint(subtalar_r);
	model->addBody(central_foot_r);	model->addJoint(midtarsal_r);
	model->addBody(toes_r);		    model->addJoint(mtp_r);

	/// Contact elements
	/// Parameters
	osim_double_adouble radiusSphere_s1 = 0.03232;
	osim_double_adouble radiusSphere_s2 = 0.03232;
	osim_double_adouble radiusSphere_s3 = 0.023374;
	osim_double_adouble radiusSphere_s4 = 0.020508;
	osim_double_adouble radiusSphere_s5 = 0.016244;
	osim_double_adouble radiusSphere_s6 = 0.018414;
	osim_double_adouble stiffness = 1000000;
	osim_double_adouble dissipation = 2.0;
	osim_double_adouble staticFriction = 0.8;
	osim_double_adouble dynamicFriction = 0.8;
	osim_double_adouble viscousFriction = 0.5;
	osim_double_adouble transitionVelocity = 0.2;
	Vec3 normal = Vec3(0, 1, 0);
	osim_double_adouble offset = 0;
	Vec3 locSphere_1_r(-0.00042152, -0.01, -0.0049972);
	Vec3 locSphere_2_r(0.06, -0.01, 0.020001);
	Vec3 locSphere_3_r(0.086976, -0.035452,  0.027993);
	Vec3 locSphere_4_r(0.086976, -0.035452,  -0.00318967);
	Vec3 locSphere_5_r(0.053154, -0.01, -0.0034173);
	Vec3 locSphere_6_r(1.7381e-06, -0.01, 0.022294);
	/// Right foot contact shere specifications
	HC_1_r = new HuntCrossleyForce_smooth("sphere_1_r", "calcn_r", locSphere_1_r, radiusSphere_s1,
		stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction, transitionVelocity, normal, offset);
	HC_2_r = new HuntCrossleyForce_smooth("sphere_2_r", "calcn_r", locSphere_2_r, radiusSphere_s2,
		stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction, transitionVelocity, normal, offset);
	HC_3_r = new HuntCrossleyForce_smooth("sphere_3_r", "central_foot_r", locSphere_3_r, radiusSphere_s3,
		stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction, transitionVelocity, normal, offset);
	HC_4_r = new HuntCrossleyForce_smooth("sphere_4_r", "central_foot_r", locSphere_4_r, radiusSphere_s4,
		stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction, transitionVelocity, normal, offset);
	HC_5_r = new HuntCrossleyForce_smooth("sphere_5_r", "toes_r", locSphere_5_r, radiusSphere_s5,
		stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction, transitionVelocity, normal, offset);
	HC_6_r = new HuntCrossleyForce_smooth("sphere_6_r", "toes_r", locSphere_6_r, radiusSphere_s6,
		stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction, transitionVelocity, normal, offset);
	/// Add right foot contact spheres to model
	model->addComponent(HC_1_r);
	HC_1_r->connectSocket_body_sphere(*calcn_r);
	model->addComponent(HC_2_r);
	HC_2_r->connectSocket_body_sphere(*calcn_r);
	model->addComponent(HC_3_r);
	HC_3_r->connectSocket_body_sphere(*central_foot_r);
	model->addComponent(HC_4_r);
	HC_4_r->connectSocket_body_sphere(*central_foot_r);
	model->addComponent(HC_5_r);
	HC_5_r->connectSocket_body_sphere(*toes_r);
	model->addComponent(HC_6_r);
	HC_6_r->connectSocket_body_sphere(*toes_r);

	// Initialize system and state
	SimTK::State* state;
	state = new State(model->initSystem());

	// Read inputs
	std::vector<T> x(arg[0], arg[0] + NX);
	std::vector<T> u(arg[1], arg[1] + NU);

	// States and controls
	T ua[ndof]; /// joint accelerations (Qdotdots) - controls (+2 because of locked joints)
	Vector QsUs(NX); /// joint positions (Qs) and velocities (Us) - states

						 // Assign inputs to model variables
						 /// States
	for (int i = 0; i < NX; ++i) QsUs[i] = x[i];
	/// Controls
	for (int i = 0; i < ndof; ++i) ua[i] = u[i];

	// Set state variables and realize
	model->setStateVariableValues(*state, QsUs);
	model->realizeVelocity(*state);

	// Compute residual forces
	/// appliedMobilityForces (# mobilities)
	Vector appliedMobilityForces(ndofr);
	appliedMobilityForces.setToZero();
	/// appliedBodyForces (# bodies + ground)
	Vector_<SpatialVec> appliedBodyForces;
	int nbodies = model->getBodySet().getSize() + 1;
	appliedBodyForces.resize(nbodies);
	appliedBodyForces.setToZero();
	/// Set gravity
	Vec3 gravity(0);
	gravity[1] = -9.81;
	/// Add weights to appliedBodyForces
	for (int i = 0; i < model->getBodySet().getSize(); ++i) {
		model->getMatterSubsystem().addInStationForce(*state,
			model->getBodySet().get(i).getMobilizedBodyIndex(),
			model->getBodySet().get(i).getMassCenter(),
			model->getBodySet().get(i).getMass()*gravity, appliedBodyForces);
	}
	/// Add contact forces to appliedBodyForces
	/// Right foot
	Array<osim_double_adouble> Force_values_1_r = HC_1_r->getRecordValues(*state);
	Array<osim_double_adouble> Force_values_2_r = HC_2_r->getRecordValues(*state);
	Array<osim_double_adouble> Force_values_3_r = HC_3_r->getRecordValues(*state);
	Array<osim_double_adouble> Force_values_4_r = HC_4_r->getRecordValues(*state);
	Array<osim_double_adouble> Force_values_5_r = HC_5_r->getRecordValues(*state);
	Array<osim_double_adouble> Force_values_6_r = HC_6_r->getRecordValues(*state);
	SpatialVec GRF_1_r;
	GRF_1_r[0] = Vec3(Force_values_1_r[9], Force_values_1_r[10], Force_values_1_r[11]);
	GRF_1_r[1] = Vec3(Force_values_1_r[6], Force_values_1_r[7], Force_values_1_r[8]);
	SpatialVec GRF_2_r;
	GRF_2_r[0] = Vec3(Force_values_2_r[9], Force_values_2_r[10], Force_values_2_r[11]);
	GRF_2_r[1] = Vec3(Force_values_2_r[6], Force_values_2_r[7], Force_values_2_r[8]);
	SpatialVec GRF_3_r;
	GRF_3_r[0] = Vec3(Force_values_3_r[9], Force_values_3_r[10], Force_values_3_r[11]);
	GRF_3_r[1] = Vec3(Force_values_3_r[6], Force_values_3_r[7], Force_values_3_r[8]);
	SpatialVec GRF_4_r;
	GRF_4_r[0] = Vec3(Force_values_4_r[9], Force_values_4_r[10], Force_values_4_r[11]);
	GRF_4_r[1] = Vec3(Force_values_4_r[6], Force_values_4_r[7], Force_values_4_r[8]);
	SpatialVec GRF_5_r;
	GRF_5_r[0] = Vec3(Force_values_5_r[9], Force_values_5_r[10], Force_values_5_r[11]);
	GRF_5_r[1] = Vec3(Force_values_5_r[6], Force_values_5_r[7], Force_values_5_r[8]);
	SpatialVec GRF_6_r;
	GRF_6_r[0] = Vec3(Force_values_6_r[9], Force_values_6_r[10], Force_values_6_r[11]);
	GRF_6_r[1] = Vec3(Force_values_6_r[6], Force_values_6_r[7], Force_values_6_r[8]);
	int ncalcn_r = model->getBodySet().get("calcn_r").getMobilizedBodyIndex();
	int ncentral_foot_r = model->getBodySet().get("central_foot_r").getMobilizedBodyIndex();
	int ntoes_r = model->getBodySet().get("toes_r").getMobilizedBodyIndex();
	appliedBodyForces[ncalcn_r] = appliedBodyForces[ncalcn_r] + GRF_1_r + GRF_2_r;
	appliedBodyForces[ncentral_foot_r] = appliedBodyForces[ncentral_foot_r] + GRF_3_r + GRF_4_r;
	appliedBodyForces[ntoes_r] = appliedBodyForces[ntoes_r] + GRF_5_r + GRF_6_r;

	/// knownUdot
	Vector knownUdot(ndofr);
	knownUdot.setToZero();
	for (int i = 0; i < ndofr; ++i) knownUdot[i] = ua[i];

	/// Calculate residual forces
	Vector residualMobilityForces(ndofr);
	residualMobilityForces.setToZero();
	model->getMatterSubsystem().calcResidualForceIgnoringConstraints(*state,
		appliedMobilityForces, appliedBodyForces, knownUdot,
		residualMobilityForces);

	// Extract several joint origins to set constraints in problem
	Vec3 femur_or_r = femur_r->getPositionInGround(*state);
	Vec3 tibia_or_r = tibia_r->getPositionInGround(*state);
	Vec3 calcn_or_r = calcn_r->getPositionInGround(*state);
	Vec3 talus_or_r = talus_r->getPositionInGround(*state);
	Vec3 central_foot_or_r = central_foot_r->getPositionInGround(*state);
	Vec3 toes_or_r = toes_r->getPositionInGround(*state);

	// Extract contact sphere origin positions
	Vec3 locSphere_1_r_GND = calcn_r->expressVectorInGround(*state,locSphere_1_r);
	Vec3 locSphere_2_r_GND = calcn_r->expressVectorInGround(*state,locSphere_2_r);
	Vec3 locSphere_3_r_GND = central_foot_r->expressVectorInGround(*state,locSphere_3_r);
	Vec3 locSphere_4_r_GND = central_foot_r->expressVectorInGround(*state,locSphere_4_r);
	Vec3 locSphere_5_r_GND = toes_r->expressVectorInGround(*state,locSphere_5_r);
	Vec3 locSphere_6_r_GND = toes_r->expressVectorInGround(*state,locSphere_6_r);

	// Extract results
	int nc = 3;

	/// Residual forces (0-13)
	for (int i = 0; i < ndof; ++i) res[0][i] = value<T>(residualMobilityForces[i]);

	/// Joint origins femur
	for (int i = 0; i < nc; ++i) {
		res[0][i + ndof] = value<T>(femur_or_r[i]);
	}

	/// Joint origins tibia
	for (int i = 0; i < nc; ++i) {
		res[0][i + ndof + nc] = value<T>(tibia_or_r[i]);
	}

	/// Joint origins talus
	for (int i = 0; i < nc; ++i) {
		res[0][i + ndof + nc * 2] = value<T>(talus_or_r[i]);
	}

	/// Joint origins calcaneus
	for (int i = 0; i < nc; ++i) {
		res[0][i + ndof + nc * 3] = value<T>(calcn_or_r[i]);
	}

	/// Joint origins central foot
	for (int i = 0; i < nc; ++i) {
		res[0][i + ndof + nc * 4] = value<T>(central_foot_or_r[i]);
	}

	/// Joint origins toes
	for (int i = 0; i < nc; ++i) {
		res[0][i + ndof + nc * 5] = value<T>(toes_or_r[i]);
	}

	/// contact sphere 1
	for (int i = 0; i < nc; ++i) {
		res[0][i + ndof + nc * 6] = value<T>(locSphere_1_r_GND[i]);
	}
	/// contact sphere 2
	for (int i = 0; i < nc; ++i) {
		res[0][i + ndof + nc * 7] = value<T>(locSphere_2_r_GND[i]);
	}

	/// contact sphere 3
	for (int i = 0; i < nc; ++i) {
		res[0][i + ndof + nc * 8] = value<T>(locSphere_3_r_GND[i]);
	}
	/// contact sphere 4
	for (int i = 0; i < nc; ++i) {
		res[0][i + ndof + nc * 9] = value<T>(locSphere_4_r_GND[i]);
	}

	/// contact sphere 5
	for (int i = 0; i < nc; ++i) {
		res[0][i + ndof + nc * 10] = value<T>(locSphere_5_r_GND[i]);
	}
	/// contact sphere 6
	for (int i = 0; i < nc; ++i) {
		res[0][i + ndof + nc * 11] = value<T>(locSphere_6_r_GND[i]);
	}

	
	return 0;

}

/* In main(), the Recorder is used to save the expression graph of function F.
This expression graph is saved as a MATLAB function named foo.m. From this
function, a c-code can be generated via CasADi and then compiled as a dll. This
dll is then imported in MATLAB as an external function. With this workflow,
CasADi can use algorithmic differentiation to differentiate the function F.
*/
int main() {

	Recorder x[NX];
	Recorder u[NU];
	Recorder tau[NR];

	for (int i = 0; i < NX; ++i) x[i] <<= 0;
	for (int i = 0; i < NU; ++i) u[i] <<= 0;

	const Recorder* Recorder_arg[n_in] = { x,u };
	Recorder* Recorder_res[n_out] = { tau };

	F_generic<Recorder>(Recorder_arg, Recorder_res);

	double res[NR];
	for (int i = 0; i < NR; ++i) Recorder_res[0][i] >>= res[i];

	Recorder::stop_recording();

	return 0;

}