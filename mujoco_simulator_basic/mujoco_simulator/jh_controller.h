#pragma once
#include <mujoco.h>
#include <Eigen/dense>
#include <iostream>

using namespace Eigen;

extern bool Control_F;
extern bool Control_T;
extern bool Control_R;
extern bool Control_init;
extern double run_command_time;
extern bool showdata;

extern char data_print_data[1000];
extern char data_print_title[1000];








class jh_controller
{
public:
	jh_controller();
	~jh_controller();



	class part {
	public:
		part();
		void name(mjModel *m, mjData *d, char *part_name);
		void refresh();
		void SetContact(Vector3d contact_vector);
		void SetContact(double x, double y, double z);
		void SetPoint(Vector3d contact_vector);
		void SetPoint(double x, double y, double z);

		MatrixXd Rotm;

		MatrixXd Jac_point;

		MatrixXd Jac;

		MatrixXd Jac_COM;
		MatrixXd Jac_COM_p;
		MatrixXd Jac_COM_r;

		MatrixXd Jac_Contact;

		double Mass;

		MatrixXd inertia;

		/*! @Cartesian position of body frame
		*/
		Vector3d xpos;

		/*! @Cartesian position of body COM
		*/
		Vector3d xipos;

		Matrix3d di;
		Matrix3d i_rot;


	private:
		int part_id;
		mjData *d_int;
		mjModel *m_int;

	};








	/*! @initialize jh_controller
	*  initialize variables of controller
	*/
	void init(mjModel *m, mjData *d);
	

	/*! @jh_controller start
	*  current mode is torque controller
	*  Calculated torque to d->ctrl
	*/
	void simulation(mjModel *m, mjData *d);

	







private:

	//Mode bool
	


	static Matrix3d Rotation_base_frame;
	static MatrixXd Rotation2g;

	
	
	part* pR;
	part* pL;
	part* pUR;
	part* pUL;
	part* pBase;

	//Global variables of simulation

	VectorXd torque;
	mjtNum* torque_mj;

	Matrix3d Rot_init;
	Vector3d x_init;
	VectorXd torque_eigen;
	Vector3d desired_x;


	bool initializer;






















	

};


//mujoco simulation functions


/*! @Convert mujoco matrix to eigen matrix.
*  mat_mj = mujoco matrix variable
*  nrow = number of rows of mujoco matrix
*  ncol = number of columns of mujoco matrix
*  ncol = 1 if mat_mj is vector
*/
MatrixXd mj2eigen(mjtNum *mat_mj, int nrow, int ncol);


/*! @Convert mujoco matrix to eigen matrix.
*  mat_mj = mujoco matrix variable
*  start_point = startaddress of mjtnum
*  nrow = number of rows of mujoco matrix
*  ncol = number of columns of mujoco matrix
*  ncol = 1 if mat_mj is vector
*/
MatrixXd mj2eigen(mjtNum *num, int start_point, int nrow, int ncol);

/*! @Allocate mjtNum and insert Eigen Matrix to mjtNum.
*  mat = eigen matrix
*  d = mujoco data d
*/
mjtNum *eigen2mj(MatrixXd mat, mjData *d);

/*! @Return Skew symmatric matrix.
*
*
*/
MatrixXd Skm(Vector3d x);

/*! @Simple Trajectory generation function
*  rT = current time
*  rT_0 = Start time
*  rT_f = Final time
*  rx_0 = Start value
*  rx_dot_0 = Derivative of start value
*  rx_f = Final value
*  rx_dot_f = Derivative of Final value
*  before rT_0, return value is rx_0
*  after rT_f, return value is rx_f
*/
double Cubic(double rT, double rT_0, double rT_f, double rx_0, double rx_dot_0, double rx_f, double rx_dot_f);



Matrix3d Rotation_xyz(double rad_x, double rad_y, double rad_z);

Vector3d rotationMatrixToEulerAngles(Matrix3d &R);













