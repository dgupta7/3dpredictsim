/*  This code describes the OpenSim model and the skeleton dynamics
Author: Antoine Falisse
Implementation exoskeleton: Maarten Afschrift
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
constexpr int ndof = 31+2;        // # degrees of freedom (excluding locked)
constexpr int ndofr = ndof + 2;   // # degrees of freedom (including locked)
constexpr int NX = ndof * 2;      // # states
constexpr int NU = ndof + 2;        // # controls
constexpr int NR = ndof + 5 * 4 + 2;    // # residual torques + # joint origins + 2 exoskeleton velocities

									// Helper function value
template<typename T>
T value(const Recorder& e) { return e; }
template<>
double value(const Recorder& e) { return e.getValue(); }

// OpenSim and Simbody use different indices for the states/controls when the
// kinematic chain has joints up and down the origin (e.g., lumbar joint/arms
// and legs with pelvis as origin).
// The two following functions allow getting the indices from one reference
// system to the other. These functions are inspired from
// createSystemYIndexMap() in Moco.
// getIndicesOSInSimbody() returns the indices of the OpenSim Qs in the Simbody
// reference system. Note that we only care about the order here so we divide
// by 2 because the states include both Qs and Qdots.
SimTK::Array_<int> getIndicesOSInSimbody(const Model& model) {
	auto s = model.getWorkingState();
	const auto svNames = model.getStateVariableNames();
	SimTK::Array_<int> idxOSInSimbody(s.getNQ());
	s.updQ() = 0;
	for (int iy = 0; iy < s.getNQ(); ++iy) {
		s.updQ()[iy] = SimTK::NaN;
		const auto svValues = model.getStateVariableValues(s);
		for (int isv = 0; isv < svNames.size(); ++isv) {
			if (SimTK::isNaN(svValues[isv])) {
				s.updQ()[iy] = 0;
				idxOSInSimbody[iy] = isv / 2;
				break;
			}
		}
	}
	return idxOSInSimbody;
}
// getIndicesSimbodyInOS() returns the indices of the Simbody Qs in the OpenSim
// reference system.
SimTK::Array_<int> getIndicesSimbodyInOS(const Model& model) {
	auto idxOSInSimbody = getIndicesOSInSimbody(model);
	auto s = model.getWorkingState();
	SimTK::Array_<int> idxSimbodyInOS(s.getNQ());
	for (int iy = 0; iy < s.getNQ(); ++iy) {
		for (int iyy = 0; iyy < s.getNQ(); ++iyy) {
			if (idxOSInSimbody[iyy] == iy) {
				idxSimbodyInOS[iy] = iyy;
				break;
			}
		}
	}
	return idxSimbodyInOS;
}

// Torque actuator. Adds torque to the appliedBodyForces
void ExoActuation(Model& model, const OpenSim::Body& bodyA, const OpenSim::Body& bodyB, State& s, const Vec3& ExoTorque, Vector_<SpatialVec>& appliedBodyForces)
{

	/// Get the mobilized bodies
	SimTK::MobilizedBodyIndex IndexBodyA = bodyA.getMobilizedBodyIndex();
	SimTK::MobilizedBodyIndex IndexBodyB = bodyB.getMobilizedBodyIndex();
	SimTK::MobilizedBody Mobil_BodyA = bodyA.getMobilizedBody();
	SimTK::MobilizedBody Mobil_BodyB = bodyB.getMobilizedBody();

	/// Express the torque in ground frame.
	/// Note that we assume that the ExoTorque is always expressed in  the coordinate system of body A
	Vec3 Exo_Ground = Mobil_BodyA.expressVectorInGroundFrame(s, ExoTorque);

	/// apply the torque on body A
	model.getMatterSubsystem().addInBodyTorque(s, IndexBodyA, Exo_Ground, appliedBodyForces);

	/// and torque in opposite direction on body B
	model.getMatterSubsystem().addInBodyTorque(s, IndexBodyB, -Exo_Ground, appliedBodyForces);
}

// Find angular velocity of the exoskeleton motor
osim_double_adouble ExoAngularVelocity(Model& model, const OpenSim::Body& bodyA, const OpenSim::Body& bodyB, State& s, const Vec3& ExoTorque)
{

	/// Get the mobilized bodies in the coordinate system of body A
	SimTK::MobilizedBodyIndex IndexBodyA = bodyA.getMobilizedBodyIndex();
	SimTK::MobilizedBodyIndex IndexBodyB = bodyB.getMobilizedBodyIndex();
	SimTK::MobilizedBody Mobil_BodyA = bodyA.getMobilizedBody();
	SimTK::MobilizedBody Mobil_BodyB = bodyB.getMobilizedBody();

	/// Get the relative angular velocity vector
	Vec3 Exo_Vel = Mobil_BodyB.findBodyAngularVelocityInAnotherBody(s, Mobil_BodyA);

	/// Note that we assume that the ExoTorque is always expressed in the coordinate system of body A

	/// Project velocity onto torque to get scalar velocity of motor
	return dot(Exo_Vel, ExoTorque.normalize());
}

// Function F
template<typename T>
int F_generic(const T** arg, T** res) {

	// OpenSim model: create components
	/// Model
	OpenSim::Model* model;
	/// Bodies
	OpenSim::Body* pelvis;
	OpenSim::Body* femur_r;
	OpenSim::Body* femur_l;
	OpenSim::Body* tibia_r;
	OpenSim::Body* tibia_l;
	OpenSim::Body* talus_r;
	OpenSim::Body* talus_l;
	OpenSim::Body* calcn_r;
	OpenSim::Body* calcn_l;
	OpenSim::Body* metatarsi_r;
	OpenSim::Body* metatarsi_l;
	OpenSim::Body* toes_r;
	OpenSim::Body* toes_l;
	OpenSim::Body* torso;
	OpenSim::Body* humerus_r;
	OpenSim::Body* humerus_l;
	OpenSim::Body* ulna_r;
	OpenSim::Body* ulna_l;
	OpenSim::Body* radius_r;
	OpenSim::Body* radius_l;
	OpenSim::Body* hand_r;
	OpenSim::Body* hand_l;
	/// Joints
	OpenSim::CustomJoint* ground_pelvis;
	OpenSim::CustomJoint* hip_r;
	OpenSim::CustomJoint* hip_l;
	OpenSim::CustomJoint* knee_r;
	OpenSim::CustomJoint* knee_l;
	OpenSim::CustomJoint* ankle_r;
	OpenSim::CustomJoint* ankle_l;
	OpenSim::CustomJoint* subtalar_r;
	OpenSim::CustomJoint* subtalar_l;
	OpenSim::PinJoint* tmt_r;
	OpenSim::PinJoint* tmt_l;
	OpenSim::PinJoint* mtp_r;
	OpenSim::PinJoint* mtp_l;
	OpenSim::CustomJoint* back;
	OpenSim::CustomJoint* shoulder_r;
	OpenSim::CustomJoint* shoulder_l;
	OpenSim::CustomJoint* elbow_r;
	OpenSim::CustomJoint* elbow_l;
	OpenSim::CustomJoint* radioulnar_r;
	OpenSim::CustomJoint* radioulnar_l;
	OpenSim::WeldJoint* radius_hand_r;
	OpenSim::WeldJoint* radius_hand_l;
	/// Contact elements
	OpenSim::HuntCrossleyForce_smooth* HC_1_r;
	OpenSim::HuntCrossleyForce_smooth* HC_2_r;
	OpenSim::HuntCrossleyForce_smooth* HC_3_r;
	OpenSim::HuntCrossleyForce_smooth* HC_4_r;
	OpenSim::HuntCrossleyForce_smooth* HC_5_r;
	OpenSim::HuntCrossleyForce_smooth* HC_6_r;
	OpenSim::HuntCrossleyForce_smooth* HC_1_l;
	OpenSim::HuntCrossleyForce_smooth* HC_2_l;
	OpenSim::HuntCrossleyForce_smooth* HC_3_l;
	OpenSim::HuntCrossleyForce_smooth* HC_4_l;
	OpenSim::HuntCrossleyForce_smooth* HC_5_l;
	OpenSim::HuntCrossleyForce_smooth* HC_6_l;


	// Inertia of the exoskeleton
	/// Human segment properties
	Vec3 COMTibia = Vec3(0, -0.178389, 0);
	Vec3 COMTalus = Vec3(0, 0, 0);
	osim_double_adouble lTibia = 0.39;
	osim_double_adouble mTibia = 2.8735207629071704;
	osim_double_adouble mTalus = 0.07750561734071933;
	Vec3 ITibia = Vec3(0.0356624, 0.00360871, 0.0361578);
	Vec3 ITalus = Vec3(0.000647371, 0.000647371, 0.000647371);

	/// Ankle foot exoskeleton properties
	osim_double_adouble mShank_Exo = 0.88*0.5;
	osim_double_adouble mFoot_Exo = 0.88*0.5;
	Vec3 COMShank_Exo = Vec3(-0.0381, -lTibia + 0.1016, 0);
	Vec3 COMFoot_Exo = Vec3(0, -0.0381, 0);
	Vec3 IFootExo = Vec3(0.0021, 0.0068, 0.0050);
	Vec3 IShankExo = Vec3(0.0073, 0.0027, 0.0066);

	/// Get combined COM location
	osim_double_adouble mTib_tot = mTibia + mShank_Exo;
	osim_double_adouble mFoot_tot = mTalus + mFoot_Exo;
	Vec3 COM_TibNew, COM_FootNew;
	osim_double_adouble COMs;
	for (int i = 0; i < 3; ++i) {
		COMs = (COMShank_Exo.get(i)*mShank_Exo + COMTibia.get(i)*mTibia) / mTib_tot;
		COM_TibNew.set(i, COMs);
		COMs = (COMFoot_Exo.get(i)*mFoot_Exo + COMTalus.get(i)*mTalus) / mTib_tot;
		COM_FootNew.set(i, COMs);
	}

	// Get the new inertia: x -axis	
	osim_double_adouble dTibia, dExoShank, dTalus, dExoFoot;
	dTibia = sqrt(pow(COM_TibNew.get(1) - COMTibia.get(1), 2) + pow(COM_TibNew.get(2) - COMTibia.get(2), 2));
	dExoShank = sqrt(pow(COM_TibNew.get(1) - COMShank_Exo.get(1), 2) + pow(COM_TibNew.get(2) - COMShank_Exo.get(2), 2));
	dTalus = sqrt(pow(COM_FootNew.get(1) - COMTalus.get(1), 2) + pow(COM_FootNew.get(2) - COMTalus.get(2), 2));
	dExoFoot = sqrt(pow(COM_FootNew.get(1) - COMTalus.get(1), 2) + pow(COM_FootNew.get(2) - COMTalus.get(2), 2));
	osim_double_adouble IFootTotalx = IFootExo.get(0) + mFoot_Exo * pow(dExoFoot, 2) + ITalus.get(0) + mTalus * pow(dTalus, 2);
	osim_double_adouble IShankTotalx = IShankExo.get(0) + mShank_Exo * pow(dExoShank, 2) + ITibia.get(0) + mTibia * pow(dTibia, 2);

	// Get the new inertia: y -axis	
	dTibia = sqrt(pow(COM_TibNew.get(0) - COMTibia.get(0), 2) + pow(COM_TibNew.get(2) - COMTibia.get(2), 2));
	dExoShank = sqrt(pow(COM_TibNew.get(0) - COMShank_Exo.get(0), 2) + pow(COM_TibNew.get(2) - COMShank_Exo.get(2), 2));
	dTalus = sqrt(pow(COM_FootNew.get(0) - COMTalus.get(0), 2) + pow(COM_FootNew.get(2) - COMTalus.get(2), 2));
	dExoFoot = sqrt(pow(COM_FootNew.get(0) - COMFoot_Exo.get(0), 2) + pow(COM_FootNew.get(2) - COMTalus.get(2), 2));
	osim_double_adouble IFootTotaly = IFootExo.get(1) + mFoot_Exo * pow(dExoFoot, 2) + ITalus.get(1) + mTalus * pow(dTalus, 2);
	osim_double_adouble IShankTotaly = IShankExo.get(1) + mShank_Exo * pow(dExoShank, 2) + ITibia.get(1) + mTibia * pow(dTibia, 2);

	// Get the new inertia: z -axis	
	dTibia = sqrt(pow(COM_TibNew.get(0) - COMTibia.get(0), 2) + pow(COM_TibNew.get(1) - COMTibia.get(1), 2));
	dExoShank = sqrt(pow(COM_TibNew.get(0) - COMShank_Exo.get(0), 2) + pow(COM_TibNew.get(1) - COMShank_Exo.get(1), 2));
	dTalus = sqrt(pow(COM_FootNew.get(0) - COMTalus.get(0), 2) + pow(COM_FootNew.get(1) - COMTalus.get(1), 2));
	dExoFoot = sqrt(pow(COM_FootNew.get(0) - COMFoot_Exo.get(0), 2) + pow(COM_FootNew.get(1) - COMTalus.get(1), 2));
	osim_double_adouble IFootTotalz = IFootExo.get(2) + mFoot_Exo * pow(dExoFoot, 2) + ITalus.get(2) + mTalus * pow(dTalus, 2);
	osim_double_adouble IShankTotalz = IShankExo.get(2) + mShank_Exo * pow(dExoShank, 2) + ITibia.get(2) + mTibia * pow(dTibia, 2);

	// resulting inertia
	Inertia ITibiaNew(IShankTotalx, IShankTotaly, IShankTotalz, 0, 0, 0);
	Inertia IFootNew(IFootTotalx, IFootTotaly, IFootTotalz, 0, 0, 0);
	
	// OpenSim model: initialize components
	/// Model
	model = new OpenSim::Model();
	/// Body specifications
	pelvis = new OpenSim::Body("pelvis", 9.127836554216511, Vec3(-0.0655315, 0, 0), Inertia(0.0705367, 0.0705367, 0.0385543, 0, 0, 0));
	femur_l = new OpenSim::Body("femur_l", 7.209107491329666, Vec3(0, - 0.152829, 0), Inertia(0.0838736, 0.0219862, 0.088446, 0, 0, 0));
	femur_r = new OpenSim::Body("femur_r", 7.209107491329666, Vec3(0, -0.152829, 0), Inertia(0.0838736, 0.0219862, 0.088446, 0, 0, 0));
	tibia_l = new OpenSim::Body("tibia_l", mTib_tot, COM_TibNew, ITibiaNew);
	tibia_r = new OpenSim::Body("tibia_r", mTib_tot, COM_TibNew, ITibiaNew);
	talus_l = new OpenSim::Body("talus_l", mFoot_tot, COM_FootNew, IFootNew);
	talus_r = new OpenSim::Body("talus_r", mFoot_tot, COM_FootNew, IFootNew);
	//calcn_l = new OpenSim::Body("calcn_l", 0.9688202167589921, Vec3(0.0913924, 0.0274177, 0), Inertia(0.000906321, 0.00252475, 0.00265422, 0, 0, 0));
	//calcn_r = new OpenSim::Body("calcn_r", 0.9688202167589921, Vec3(0.0913924, 0.0274177, 0), Inertia(0.000906321, 0.00252475, 0.00265422, 0, 0, 0));
	calcn_l = new OpenSim::Body("calcn_l", 0.625241534065349, Vec3(0.066907031412174, 0.033043440839674, 2.064106263057383e-04), Inertia(6.670471633986739e-04, 8.494673219871866e-04, 0.001011813716159, 0, 0, 0));
	calcn_r = new OpenSim::Body("calcn_r", 0.625241534065349, Vec3(0.066907031412174, 0.033043440839674, -2.064106263057383e-04), Inertia(6.670471633986739e-04, 8.494673219871866e-04, 0.001011813716159, 0, 0, 0));
	metatarsi_l = new OpenSim::Body("metatarsi_l", 0.343578682693644, Vec3(0.027675743715505, -0.018380814192906, 0.000077084485412), Inertia(1.834e-04, 6.182e-04, 5.296e-04, 0, 0, 0));
	metatarsi_r = new OpenSim::Body("metatarsi_r", 0.343578682693644, Vec3(0.027675743715505, -0.018380814192906, -0.000077084485412), Inertia(1.834e-04, 6.182e-04, 5.296e-04, 0, 0, 0));
	toes_l = new OpenSim::Body("toes_l", 0.16787716715999804, Vec3(0.0316218, 0.00548355, 0.0159937), Inertia(6.2714132461258e-005, 0.000125428264922516, 6.2714132461258e-005, 0, 0, 0));
	toes_r = new OpenSim::Body("toes_r", 0.16787716715999804, Vec3(0.0316218, 0.00548355, -0.0159937), Inertia(6.2714132461258e-005, 0.000125428264922516, 6.2714132461258e-005, 0, 0, 0));
	torso = new OpenSim::Body("torso", 26.535288186472687, Vec3(-0.0267603, 0.306505, 0), Inertia(1.10388, 0.507805, 1.10388, 0, 0, 0));
	humerus_l = new OpenSim::Body("humerus_l", 1.5753481758205221, Vec3(0, -0.150554, 0), Inertia(0.00775535, 0.00267536, 0.00870516, 0, 0, 0));
	humerus_r = new OpenSim::Body("humerus_r", 1.5753481758205221, Vec3(0, -0.150554, 0), Inertia(0.00775535, 0.00267536, 0.00870516, 0, 0, 0));
	ulna_l = new OpenSim::Body("ulna_l", 0.4708466253448705, Vec3(0, -0.114441, 0), Inertia(0.00206978, 0.000431845, 0.00224518, 0, 0, 0));
	ulna_r = new OpenSim::Body("ulna_r", 0.4708466253448705, Vec3(0, -0.114441, 0), Inertia(0.00206978, 0.000431845, 0.00224518, 0, 0, 0));
	radius_l = new OpenSim::Body("radius_l", 0.4708466253448705, Vec3(0, -0.114441, 0), Inertia(0.00206978, 0.000431845, 0.00224518, 0, 0, 0));
	radius_r = new OpenSim::Body("radius_r", 0.4708466253448705, Vec3(0, -0.114441, 0), Inertia(0.00206978, 0.000431845, 0.00224518, 0, 0, 0));
	hand_l = new OpenSim::Body("hand_l", 0.35458819933379115, Vec3(0, -0.0668239, 0), Inertia(0.000644974415108496, 0.000395516821821017, 0.000968907753638324, 0, 0, 0));
	hand_r = new OpenSim::Body("hand_r", 0.35458819933379115, Vec3(0, -0.0668239, 0), Inertia(0.000644974415108496, 0.000395516821821017, 0.000968907753638324, 0, 0, 0));
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
	/// Hip_l transform
	SpatialTransform st_hip_l;
	st_hip_l[0].setCoordinateNames(OpenSim::Array<std::string>("hip_flexion_l", 1, 1));
	st_hip_l[0].setFunction(new LinearFunction());
	st_hip_l[0].setAxis(Vec3(0, 0, 1));
	st_hip_l[1].setCoordinateNames(OpenSim::Array<std::string>("hip_adduction_l", 1, 1));
	st_hip_l[1].setFunction(new LinearFunction());
	st_hip_l[1].setAxis(Vec3(-1, 0, 0));
	st_hip_l[2].setCoordinateNames(OpenSim::Array<std::string>("hip_rotation_l", 1, 1));
	st_hip_l[2].setFunction(new LinearFunction());
	st_hip_l[2].setAxis(Vec3(0, -1, 0));
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
	/// Knee_l transform
	SpatialTransform st_knee_l;
	st_knee_l[2].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_l", 1, 1));
	st_knee_l[2].setFunction(new LinearFunction());
	st_knee_l[2].setAxis(Vec3(0, 0, 1));
	/// Knee_r transform
	SpatialTransform st_knee_r;
	st_knee_r[2].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_r", 1, 1));
	st_knee_r[2].setFunction(new LinearFunction());
	st_knee_r[2].setAxis(Vec3(0, 0, 1));
	/// Ankle_l transform
	SpatialTransform st_ankle_l;
	st_ankle_l[0].setCoordinateNames(OpenSim::Array<std::string>("ankle_angle_l", 1, 1));
	st_ankle_l[0].setFunction(new LinearFunction());
	st_ankle_l[0].setAxis(Vec3(0.10501355, 0.17402245, 0.97912632));
	/// Ankle_r transform
	SpatialTransform st_ankle_r;
	st_ankle_r[0].setCoordinateNames(OpenSim::Array<std::string>("ankle_angle_r", 1, 1));
	st_ankle_r[0].setFunction(new LinearFunction());
	st_ankle_r[0].setAxis(Vec3(-0.10501355, -0.17402245, 0.97912632));
	/// Subtalar_l transform
	SpatialTransform st_subtalar_l;
	st_subtalar_l[0].setCoordinateNames(OpenSim::Array<std::string>("subtalar_angle_l", 1, 1));
	st_subtalar_l[0].setFunction(new LinearFunction());
	st_subtalar_l[0].setAxis(Vec3(-0.78717961, -0.60474746, -0.12094949));
	/// Subtalar_r transform
	SpatialTransform st_subtalar_r;
	st_subtalar_r[0].setCoordinateNames(OpenSim::Array<std::string>("subtalar_angle_r", 1, 1));
	st_subtalar_r[0].setFunction(new LinearFunction());
	st_subtalar_r[0].setAxis(Vec3(0.78717961, 0.60474746, -0.12094949));
	/// Back transform
	SpatialTransform st_back;
	st_back[0].setCoordinateNames(OpenSim::Array<std::string>("lumbar_extension", 1, 1));
	st_back[0].setFunction(new LinearFunction());
	st_back[0].setAxis(Vec3(0, 0, 1));
	st_back[1].setCoordinateNames(OpenSim::Array<std::string>("lumbar_bending", 1, 1));
	st_back[1].setFunction(new LinearFunction());
	st_back[1].setAxis(Vec3(1, 0, 0));
	st_back[2].setCoordinateNames(OpenSim::Array<std::string>("lumbar_rotation", 1, 1));
	st_back[2].setFunction(new LinearFunction());
	st_back[2].setAxis(Vec3(0, 1, 0));
	/// Shoulder_l transform
	SpatialTransform st_sho_l;
	st_sho_l[0].setCoordinateNames(OpenSim::Array<std::string>("arm_flex_l", 1, 1));
	st_sho_l[0].setFunction(new LinearFunction());
	st_sho_l[0].setAxis(Vec3(0, 0, 1));
	st_sho_l[1].setCoordinateNames(OpenSim::Array<std::string>("arm_add_l", 1, 1));
	st_sho_l[1].setFunction(new LinearFunction());
	st_sho_l[1].setAxis(Vec3(-1, 0, 0));
	st_sho_l[2].setCoordinateNames(OpenSim::Array<std::string>("arm_rot_l", 1, 1));
	st_sho_l[2].setFunction(new LinearFunction());
	st_sho_l[2].setAxis(Vec3(0, -1, 0));
	/// Shoulder_r transform
	SpatialTransform st_sho_r;
	st_sho_r[0].setCoordinateNames(OpenSim::Array<std::string>("arm_flex_r", 1, 1));
	st_sho_r[0].setFunction(new LinearFunction());
	st_sho_r[0].setAxis(Vec3(0, 0, 1));
	st_sho_r[1].setCoordinateNames(OpenSim::Array<std::string>("arm_add_r", 1, 1));
	st_sho_r[1].setFunction(new LinearFunction());
	st_sho_r[1].setAxis(Vec3(1, 0, 0));
	st_sho_r[2].setCoordinateNames(OpenSim::Array<std::string>("arm_rot_r", 1, 1));
	st_sho_r[2].setFunction(new LinearFunction());
	st_sho_r[2].setAxis(Vec3(0, 1, 0));
	/// Elbow_l transform
	SpatialTransform st_elb_l;
	st_elb_l[0].setCoordinateNames(OpenSim::Array<std::string>("elbow_flex_l", 1, 1));
	st_elb_l[0].setFunction(new LinearFunction());
	st_elb_l[0].setAxis(Vec3(-0.22604696, -0.022269, 0.97386183));
	/// Elbow_r transform
	SpatialTransform st_elb_r;
	st_elb_r[0].setCoordinateNames(OpenSim::Array<std::string>("elbow_flex_r", 1, 1));
	st_elb_r[0].setFunction(new LinearFunction());
	st_elb_r[0].setAxis(Vec3(0.22604696, 0.022269, 0.97386183));
	/// Radioulnar_l transform
	SpatialTransform st_radioulnar_l;
	st_radioulnar_l[0].setCoordinateNames(OpenSim::Array<std::string>("pro_sup_l", 1, 1));
	st_radioulnar_l[0].setFunction(new LinearFunction());
	st_radioulnar_l[0].setAxis(Vec3(-0.05639803, -0.99840646, 0.001952));
	/// Radioulnar_r transform
	SpatialTransform st_radioulnar_r;
	st_radioulnar_r[0].setCoordinateNames(OpenSim::Array<std::string>("pro_sup_r", 1, 1));
	st_radioulnar_r[0].setFunction(new LinearFunction());
	st_radioulnar_r[0].setAxis(Vec3(0.05639803, 0.99840646, 0.001952));
	/// Joint specifications
	ground_pelvis = new CustomJoint("ground_pelvis", model->getGround(), Vec3(0), Vec3(0), *pelvis, Vec3(0), Vec3(0), st_ground_pelvis);
	hip_l = new CustomJoint("hip_l", *pelvis, Vec3(-0.065531461457899232, -0.061267748151662485, -0.079019111539209888), Vec3(0), *femur_l, Vec3(0), Vec3(0), st_hip_l);
	hip_r = new CustomJoint("hip_r", *pelvis, Vec3(-0.065531461457899232, -0.061267748151662485, 0.079019111539209888), Vec3(0), *femur_r, Vec3(0), Vec3(0), st_hip_r);
	knee_l = new CustomJoint("knee_l", *femur_l, Vec3(-0.0040453268563842293, -0.35583861269464306, 0), Vec3(0), *tibia_l, Vec3(0), Vec3(0), st_knee_l);
	knee_r = new CustomJoint("knee_r", *femur_r, Vec3(-0.0040453268563842293, -0.35583861269464306, 0), Vec3(0), *tibia_r, Vec3(0), Vec3(0), st_knee_r);
	ankle_l = new CustomJoint("ankle_l", *tibia_l, Vec3(0, -0.41085882914747662, 0), Vec3(0), *talus_l, Vec3(0), Vec3(0), st_ankle_l);
	ankle_r = new CustomJoint("ankle_r", *tibia_r, Vec3(0, -0.41085882914747662, 0), Vec3(0), *talus_r, Vec3(0), Vec3(0), st_ankle_r);
	subtalar_l = new CustomJoint("subtalar_l", *talus_l, Vec3(-0.044572100000000003, -0.038339100000000001, -0.0072382799999999997), Vec3(0), *calcn_l, Vec3(0), Vec3(0), st_subtalar_l);
	subtalar_r = new CustomJoint("subtalar_r", *talus_r, Vec3(-0.044572100000000003, -0.038339100000000001, 0.0072382799999999997), Vec3(0), *calcn_r, Vec3(0), Vec3(0), st_subtalar_r);
	//mtp_l = new PinJoint("mtp_l", *calcn_l, Vec3(0.163409678774199, -0.00182784875586352, -0.000987038328166303), Vec3(0), *toes_l, Vec3(0), Vec3(0));
	//mtp_r = new PinJoint("mtp_r", *calcn_r, Vec3(0.163409678774199, -0.00182784875586352, 0.000987038328166303), Vec3(0), *toes_r, Vec3(0), Vec3(0));
	tmt_l = new PinJoint("tmt_l", *calcn_l, Vec3(0.108274919614407, 0.035560839360244, -0.000452708769279), Vec3(0), *metatarsi_l, Vec3(0), Vec3(0));
	tmt_r = new PinJoint("tmt_r", *calcn_r, Vec3(0.108274919614407, 0.035560839360244, 0.000452708769279), Vec3(0), *metatarsi_r, Vec3(0), Vec3(0));
	mtp_l = new PinJoint("mtp_l", *metatarsi_l, Vec3(0.055134759159792, -0.037388688116107, -0.000534329558887), Vec3(0), *toes_l, Vec3(0), Vec3(0));
	mtp_r = new PinJoint("mtp_r", *metatarsi_r, Vec3(0.055134759159792, -0.037388688116107, 0.000534329558887), Vec3(0), *toes_r, Vec3(0), Vec3(0));
	back = new CustomJoint("back", *pelvis, Vec3(-0.093338312405799553, 0.075541935477359282, 0), Vec3(0), *torso, Vec3(0), Vec3(0), st_back);
	shoulder_l = new CustomJoint("shoulder_l", *torso, Vec3(0.0029380855493219699, 0.37148545291063978, -0.1583124585136152), Vec3(0), *humerus_l, Vec3(0), Vec3(0), st_sho_l);
	shoulder_r = new CustomJoint("shoulder_r", *torso, Vec3(0.0029380855493219699, 0.37148545291063978, 0.1583124585136152), Vec3(0), *humerus_r, Vec3(0), Vec3(0), st_sho_r);
	elbow_l = new CustomJoint("elbow_l", *humerus_l, Vec3(0.012029609243253778, -0.26200145307574269, 0.008781487380376601), Vec3(0), *ulna_l, Vec3(0), Vec3(0), st_elb_l);
	elbow_r = new CustomJoint("elbow_r", *humerus_r, Vec3(0.012029609243253778, -0.26200145307574269, -0.008781487380376601), Vec3(0), *ulna_r, Vec3(0), Vec3(0), st_elb_r);
	radioulnar_l = new CustomJoint("radioulnar_l", *ulna_l, Vec3(-0.006387412252866027, -0.012350385822167704, -0.024766276816626718), Vec3(0), *radius_l, Vec3(0), Vec3(0), st_radioulnar_l);
	radioulnar_r = new CustomJoint("radioulnar_r", *ulna_r, Vec3(-0.006387412252866027, -0.012350385822167704, 0.024766276816626718), Vec3(0), *radius_r, Vec3(0), Vec3(0), st_radioulnar_r);
	radius_hand_l = new WeldJoint("radius_hand_l", *radius_l, Vec3(-0.0083529157504388159, -0.2239357691274558, -0.012922902963153949), Vec3(0), *hand_l, Vec3(0), Vec3(0));
	radius_hand_r = new WeldJoint("radius_hand_r", *radius_r, Vec3(-0.0083529157504388159, -0.2239357691274558, 0.012922902963153949), Vec3(0), *hand_r, Vec3(0), Vec3(0));
	/// Add bodies and joints to model
	model->addBody(pelvis);		    model->addJoint(ground_pelvis);
	model->addBody(femur_l);		model->addJoint(hip_l);
	model->addBody(femur_r);		model->addJoint(hip_r);
	model->addBody(tibia_l);		model->addJoint(knee_l);
	model->addBody(tibia_r);		model->addJoint(knee_r);
	model->addBody(talus_l);		model->addJoint(ankle_l);
	model->addBody(talus_r);		model->addJoint(ankle_r);
	model->addBody(calcn_l);		model->addJoint(subtalar_l);
	model->addBody(calcn_r);		model->addJoint(subtalar_r);
	model->addBody(metatarsi_l);	model->addJoint(tmt_l);
	model->addBody(metatarsi_r);	model->addJoint(tmt_r);
	model->addBody(toes_l);		    model->addJoint(mtp_l);
	model->addBody(toes_r);		    model->addJoint(mtp_r);
	model->addBody(torso);          model->addJoint(back);
	model->addBody(humerus_l);      model->addJoint(shoulder_l);
	model->addBody(humerus_r);      model->addJoint(shoulder_r);
	model->addBody(ulna_l);         model->addJoint(elbow_l);
	model->addBody(ulna_r);         model->addJoint(elbow_r);
	model->addBody(radius_l);       model->addJoint(radioulnar_l);
	model->addBody(radius_r);       model->addJoint(radioulnar_r);
	model->addBody(hand_l);         model->addJoint(radius_hand_l);
	model->addBody(hand_r);         model->addJoint(radius_hand_r);
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
	//Vec3 locSphere_3_r(0.165, -0.01, 0.021183);
	//Vec3 locSphere_4_r(0.165, -0.01, -0.01);
	Vec3 locSphere_3_r(0.0567, -0.0456, 0.020730);
	Vec3 locSphere_4_r(0.0567, -0.0456, -0.01);
	Vec3 locSphere_5_r(0.053154, -0.01, -0.0034173);
	Vec3 locSphere_6_r(1.7381e-06, -0.01, 0.022294);
	Vec3 locSphere_1_l(locSphere_1_r[0], locSphere_1_r[1], -locSphere_1_r[2]);
	Vec3 locSphere_2_l(locSphere_2_r[0], locSphere_2_r[1], -locSphere_2_r[2]);
	Vec3 locSphere_3_l(locSphere_3_r[0], locSphere_3_r[1], -locSphere_3_r[2]);
	Vec3 locSphere_4_l(locSphere_4_r[0], locSphere_4_r[1], -locSphere_4_r[2]);
	Vec3 locSphere_5_l(locSphere_5_r[0], locSphere_5_r[1], -locSphere_5_r[2]);
	Vec3 locSphere_6_l(locSphere_6_r[0], locSphere_6_r[1], -locSphere_6_r[2]);
	/// Left foot contact shere specifications
	HC_1_l = new HuntCrossleyForce_smooth("sphere_1_l", "calcn_l", locSphere_1_l, radiusSphere_s1,
		stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction, transitionVelocity, normal, offset);
	HC_2_l = new HuntCrossleyForce_smooth("sphere_2_l", "calcn_l", locSphere_2_l, radiusSphere_s2,
		stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction, transitionVelocity, normal, offset);
	HC_3_l = new HuntCrossleyForce_smooth("sphere_3_l", "metatarsi_l", locSphere_3_l, radiusSphere_s3,
		stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction, transitionVelocity, normal, offset);
	HC_4_l = new HuntCrossleyForce_smooth("sphere_4_l", "metatarsi_l", locSphere_4_l, radiusSphere_s4,
		stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction, transitionVelocity, normal, offset);
	HC_5_l = new HuntCrossleyForce_smooth("sphere_5_l", "toes_l", locSphere_5_l, radiusSphere_s5,
		stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction, transitionVelocity, normal, offset);
	HC_6_l = new HuntCrossleyForce_smooth("sphere_6_l", "toes_l", locSphere_6_l, radiusSphere_s6,
		stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction, transitionVelocity, normal, offset);
	/// Add left foot contact spheres to model
	model->addComponent(HC_1_l);
	HC_1_l->connectSocket_body_sphere(*calcn_l);
	model->addComponent(HC_2_l);
	HC_2_l->connectSocket_body_sphere(*calcn_l);
	model->addComponent(HC_3_l);
	HC_3_l->connectSocket_body_sphere(*metatarsi_l);
	model->addComponent(HC_4_l);
	HC_4_l->connectSocket_body_sphere(*metatarsi_l);
	model->addComponent(HC_5_l);
	HC_5_l->connectSocket_body_sphere(*toes_l);
	model->addComponent(HC_6_l);
	HC_6_l->connectSocket_body_sphere(*toes_l);
	/// Right foot contact shere specifications
	HC_1_r = new HuntCrossleyForce_smooth("sphere_1_r", "calcn_r", locSphere_1_r, radiusSphere_s1,
		stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction, transitionVelocity, normal, offset);
	HC_2_r = new HuntCrossleyForce_smooth("sphere_2_r", "calcn_r", locSphere_2_r, radiusSphere_s2,
		stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction, transitionVelocity, normal, offset);
	HC_3_r = new HuntCrossleyForce_smooth("sphere_3_r", "metatarsi_r", locSphere_3_r, radiusSphere_s3,
		stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction, transitionVelocity, normal, offset);
	HC_4_r = new HuntCrossleyForce_smooth("sphere_4_r", "metatarsi_r", locSphere_4_r, radiusSphere_s4,
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
	HC_3_r->connectSocket_body_sphere(*metatarsi_r);
	model->addComponent(HC_4_r);
	HC_4_r->connectSocket_body_sphere(*metatarsi_r);
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
	T ua[ndof + 2]; /// joint accelerations (Qdotdots) - controls
	Vector QsUs(NX + 4); /// joint positions (Qs) and velocities (Us) - states

						 // Assign inputs to model variables
						 /// States
	for (int i = 0; i < NX; ++i) QsUs[i] = x[i];
	/// pro_sup dofs are locked so Qs and Qdots are hard coded (0)
	QsUs[NX] = 1.51;
	QsUs[NX + 1] = 0;
	QsUs[NX + 2] = 1.51;
	QsUs[NX + 3] = 0;
	/// Controls
	T ut[ndof + 2];
	for (int i = 0; i < ndof; ++i) ut[i] = u[i];
	/// pro_sup dofs are locked so Qdotdots are hard coded (0)
	/// Need to have a temporary vector to add 0s to the vector before
	/// adjusting for the index difference between OpenSim and Simbody.
	ut[ndof] = 0;
	ut[ndof + 1] = 0;
	/// OpenSim and Simbody have different state orders so we need to adjust
	auto indicesOSInSimbody = getIndicesOSInSimbody(*model);
	for (int i = 0; i < ndofr; ++i) ua[i] = ut[indicesOSInSimbody[i]];

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
	int nmetatarsi_r = model->getBodySet().get("metatarsi_r").getMobilizedBodyIndex();
	int ntoes_r = model->getBodySet().get("toes_r").getMobilizedBodyIndex();
	appliedBodyForces[ncalcn_r] = appliedBodyForces[ncalcn_r] + GRF_1_r + GRF_2_r;
	appliedBodyForces[nmetatarsi_r] = appliedBodyForces[nmetatarsi_r] + GRF_3_r + GRF_4_r;
	appliedBodyForces[ntoes_r] = appliedBodyForces[ntoes_r] + GRF_5_r + GRF_6_r;
	/// Left foot
	Array<osim_double_adouble> Force_values_1_l = HC_1_l->getRecordValues(*state);
	Array<osim_double_adouble> Force_values_2_l = HC_2_l->getRecordValues(*state);
	Array<osim_double_adouble> Force_values_3_l = HC_3_l->getRecordValues(*state);
	Array<osim_double_adouble> Force_values_4_l = HC_4_l->getRecordValues(*state);
	Array<osim_double_adouble> Force_values_5_l = HC_5_l->getRecordValues(*state);
	Array<osim_double_adouble> Force_values_6_l = HC_6_l->getRecordValues(*state);
	SpatialVec GRF_1_l;
	GRF_1_l[0] = Vec3(Force_values_1_l[9], Force_values_1_l[10], Force_values_1_l[11]);
	GRF_1_l[1] = Vec3(Force_values_1_l[6], Force_values_1_l[7], Force_values_1_l[8]);
	SpatialVec GRF_2_l;
	GRF_2_l[0] = Vec3(Force_values_2_l[9], Force_values_2_l[10], Force_values_2_l[11]);
	GRF_2_l[1] = Vec3(Force_values_2_l[6], Force_values_2_l[7], Force_values_2_l[8]);
	SpatialVec GRF_3_l;
	GRF_3_l[0] = Vec3(Force_values_3_l[9], Force_values_3_l[10], Force_values_3_l[11]);
	GRF_3_l[1] = Vec3(Force_values_3_l[6], Force_values_3_l[7], Force_values_3_l[8]);
	SpatialVec GRF_4_l;
	GRF_4_l[0] = Vec3(Force_values_4_l[9], Force_values_4_l[10], Force_values_4_l[11]);
	GRF_4_l[1] = Vec3(Force_values_4_l[6], Force_values_4_l[7], Force_values_4_l[8]);
	SpatialVec GRF_5_l;
	GRF_5_l[0] = Vec3(Force_values_5_l[9], Force_values_5_l[10], Force_values_5_l[11]);
	GRF_5_l[1] = Vec3(Force_values_5_l[6], Force_values_5_l[7], Force_values_5_l[8]);
	SpatialVec GRF_6_l;
	GRF_6_l[0] = Vec3(Force_values_6_l[9], Force_values_6_l[10], Force_values_6_l[11]);
	GRF_6_l[1] = Vec3(Force_values_6_l[6], Force_values_6_l[7], Force_values_6_l[8]);
	int ncalcn_l = model->getBodySet().get("calcn_l").getMobilizedBodyIndex();
	int nmetatarsi_l = model->getBodySet().get("metatarsi_l").getMobilizedBodyIndex();
	int ntoes_l = model->getBodySet().get("toes_l").getMobilizedBodyIndex();
	appliedBodyForces[ncalcn_l] = appliedBodyForces[ncalcn_l] + GRF_1_l + GRF_2_l;
	appliedBodyForces[nmetatarsi_l] = appliedBodyForces[nmetatarsi_l] + GRF_3_l + GRF_4_l;
	appliedBodyForces[ntoes_l] = appliedBodyForces[ntoes_l] + GRF_5_l + GRF_6_l;
	/// knownUdot
	Vector knownUdot(ndofr);
	knownUdot.setToZero();
	for (int i = 0; i < ndofr; ++i) knownUdot[i] = ua[i];

	/// apply the exoskeleton torque
	Vec3 ExoTorque_l = Vec3(0);
	Vec3 ExoTorque_r = Vec3(0);
	ExoTorque_l[2] = u[ndof];
	ExoTorque_r[2] = u[ndof + 1];
	ExoActuation(*model, *tibia_l, *metatarsi_l, *state, ExoTorque_l, appliedBodyForces);
	ExoActuation(*model, *tibia_r, *metatarsi_r, *state, ExoTorque_r, appliedBodyForces);
	
	/// get exoskeleton motor velocity
	osim_double_adouble ExoVel_l = ExoAngularVelocity(*model, *tibia_l, *metatarsi_l, *state, ExoTorque_l);
	osim_double_adouble ExoVel_r = ExoAngularVelocity(*model, *tibia_r, *metatarsi_r, *state, ExoTorque_r);

	/// Calculate residual forces
	Vector residualMobilityForces(ndofr);
	residualMobilityForces.setToZero();
	model->getMatterSubsystem().calcResidualForceIgnoringConstraints(*state,
		appliedMobilityForces, appliedBodyForces, knownUdot,
		residualMobilityForces);

	// Extract several joint origins to set constraints in problem
	Vec3 calcn_or_l = calcn_l->getPositionInGround(*state);
	Vec3 calcn_or_r = calcn_r->getPositionInGround(*state);
	Vec3 femur_or_l = femur_l->getPositionInGround(*state);
	Vec3 femur_or_r = femur_r->getPositionInGround(*state);
	Vec3 hand_or_l = hand_l->getPositionInGround(*state);
	Vec3 hand_or_r = hand_r->getPositionInGround(*state);
	Vec3 tibia_or_l = tibia_l->getPositionInGround(*state);
	Vec3 tibia_or_r = tibia_r->getPositionInGround(*state);
	Vec3 toes_or_l = toes_l->getPositionInGround(*state);
	Vec3 toes_or_r = toes_r->getPositionInGround(*state);

	// Residual forces in OpenSim order
	T res_os[ndofr];
	/// OpenSim and Simbody have different state orders so we need to adjust
	auto indicesSimbodyInOS = getIndicesSimbodyInOS(*model);
	for (int i = 0; i < ndofr; ++i) res_os[i] =
		value<T>(residualMobilityForces[indicesSimbodyInOS[i]]);
	// Extract results
	int nc = 3;
	/// Residual forces
	/// We do want to extract the pro_sup torques (last two -> till NU)
	for (int i = 0; i < ndof; ++i) res[0][i] = res_os[i];
	/// Joint origins
	res[0][ndof] = value<T>(calcn_or_r[0]);   /// calcn_or_r_x
	res[0][ndof + 1] = value<T>(calcn_or_r[2]);   /// calcn_or_r_z
	res[0][ndof + 2] = value<T>(calcn_or_l[0]);   /// calcn_or_l_x
	res[0][ndof + 3] = value<T>(calcn_or_l[2]);   /// calcn_or_l_x
	res[0][ndof + 4] = value<T>(femur_or_r[0]);   /// femur_or_r_x
	res[0][ndof + 5] = value<T>(femur_or_r[2]);   /// femur_or_r_z
	res[0][ndof + 6] = value<T>(femur_or_l[0]);   /// femur_or_l_x
	res[0][ndof + 7] = value<T>(femur_or_l[2]);   /// femur_or_l_z
	res[0][ndof + 8] = value<T>(hand_or_r[0]);    /// hand_or_r_x
	res[0][ndof + 9] = value<T>(hand_or_r[2]);    /// hand_or_r_z
	res[0][ndof + 10] = value<T>(hand_or_l[0]);   /// hand_or_l_x
	res[0][ndof + 11] = value<T>(hand_or_l[2]);   /// hand_or_l_z
	res[0][ndof + 12] = value<T>(tibia_or_r[0]);  /// tibia_or_r_x
	res[0][ndof + 13] = value<T>(tibia_or_r[2]);  /// tibia_or_r_z
	res[0][ndof + 14] = value<T>(tibia_or_l[0]);  /// tibia_or_l_x
	res[0][ndof + 15] = value<T>(tibia_or_l[2]);  /// tibia_or_l_z
	res[0][ndof + 16] = value<T>(toes_or_r[0]);  /// tibia_or_r_x
	res[0][ndof + 17] = value<T>(toes_or_r[2]);  /// tibia_or_r_z
	res[0][ndof + 18] = value<T>(toes_or_l[0]);  /// tibia_or_l_x
	res[0][ndof + 19] = value<T>(toes_or_l[2]);  /// tibia_or_l_z
	/// Exoskeleton motor velocities
	res[0][ndof + 20] = value<T>(ExoVel_r);
	res[0][ndof + 21] = value<T>(ExoVel_l);

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