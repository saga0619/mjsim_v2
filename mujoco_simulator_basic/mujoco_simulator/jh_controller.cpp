#include "jh_controller.h"
using namespace std;


Matrix3d jh_controller::Rotation_base_frame;
MatrixXd jh_controller::Rotation2g;

jh_controller::jh_controller()
{
	std::cout << "Controller Loaded" << std::endl;
}


jh_controller::~jh_controller()
{
}



/** Mujoco Simulation Controller initialize.
* 각 시뮬 변수들 초기화 파트
*
*/
void jh_controller::init(mjModel *m, mjData *d) {
	std::cout << "Controller initializing" << std::endl;
	//Controller value initialize;
	torque.resize(m->nu);
	torque.setZero();
	torque_mj = mj_stackAlloc(d, m->nu);
	mju_zero(torque_mj, m->nu);
	Rotation2g.setIdentity(m->nv, m->nv);


	//Part allocation 
	pR = new part[6];
	pL = new part[6];
	pUR = new part[7];
	pUL = new part[7];
	pBase = new part[3];
	

	//Part name setting 
	char part_buff[20];
	for (int i = 0; i < 6; i++) {
	sprintf_s(part_buff, "R%d", i);
	pR[i].name(m, d, part_buff);	
	sprintf_s(part_buff, "L%d", i);
	pL[i].name(m, d, part_buff);
	}
	for (int i = 0; i < 7; i++) {
		sprintf_s(part_buff, "UR%d", i);
		pUR[i].name(m, d, part_buff);
		sprintf_s(part_buff, "UL%d", i);
		pUL[i].name(m, d, part_buff);
	}
	pBase[0].name(m, d, "Pelvis");
	pBase[1].name(m, d, "Waist");
	pBase[2].name(m, d, "Body");



	//Part Contact Setting
	pR[5].SetContact(0, 0, -0.1);
	pL[5].SetContact(0, 0, -0.1);




	torque_eigen.resize(6);
	torque_eigen.setZero();




	strcpy_s(data_print_title, "JH Controller\nData Print Test");
	std::cout << "Controller Initialize Complete" << std::endl;
}


void jh_controller::simulation(mjModel *m, mjData *d) {
	mjMARKSTACK
		//controller starts		

	/*
	mjtNum* jac_p = mj_stackAlloc(d, 6 * m->nv);
	mjtNum* jac_r = mj_stackAlloc(d, 3 * 6);

	mjtNum* qpos = mj_stackAlloc(d, 6);



	Vector3d pos_com_link[6];
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 3; j++) {
			pos_com_link[i](j) = d->xipos[6 + 3 * i + j];
		}
	}
	mjtNum* j_com_p_link[6];
	mjtNum* j_com_r_link[6];
	for (int i = 0; i < 6; i++) {
		j_com_p_link[i] = mj_stackAlloc(d, 3 * 6);
		j_com_r_link[i] = mj_stackAlloc(d, 3 * 6);
	}

	for (int i = 0; i < 6; i++) {
		mj_jacBodyCom(m, d, j_com_p_link[i], j_com_r_link[i], 2 + i);
	}

	double bodymass[6];
	for (int i = 0; i < 6; i++) {
		bodymass[i] = m->body_mass[i + 2];
	}

	Vector3d gravity;
	gravity << 0, 0, -9.81;
	MatrixXd j_com[6];




	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 6; k++) {
				j_com[i].resize(6, 6);


				j_com[i](j, k) = j_com_p_link[i][6 * j + k];
				j_com[i](3 + j, k) = j_com_r_link[i][6 * j + k];
			}
		}
	}


	VectorXd torque_gravity(6);
	torque_gravity.setZero();

	for (int i = 0; i < 6; i++) {
		torque_gravity = torque_gravity + j_com[i].block<3, 6>(0, 0).transpose()*bodymass[i] * gravity;
	}

	//End effector Position
	Matrix3d Rot_ee;
	Vector3d pos_ee_local;
	Vector3d pos_link5;
	Vector3d pos_ee;

	Rot_ee << d->xmat[63], d->xmat[64], d->xmat[65], d->xmat[66], d->xmat[67], d->xmat[68], d->xmat[69], d->xmat[70], d->xmat[71];
	pos_ee_local << 0, 0.176, 0;
	pos_link5 << d->xpos[21], d->xpos[22], d->xpos[23];
	pos_ee = pos_link5 + Rot_ee*pos_ee_local;

	mjtNum ee_pos_mj[3];
	ee_pos_mj[0] = pos_ee(0);
	ee_pos_mj[1] = pos_ee(1);
	ee_pos_mj[2] = pos_ee(2);

	mj_jac(m, d, jac_p, jac_r, ee_pos_mj, ee_id);


	//Eigen Jacobian
	MatrixXd J_ee(6, 6);

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 6; j++) {
			J_ee(i, j) = jac_p[i * 6 + j];
			J_ee(3 + i, j) = jac_r[i * 6 + j];
		}
	}

	//J_ee divided Positional, Rotational
	MatrixXd J_ee_p(3, 6), J_ee_r(3, 6);


	J_ee_p = J_ee.block<3, 6>(0, 0);
	J_ee_r = J_ee.block<3, 6>(3, 0);


	Matrix3d Rotation_EE_desired;


	//lambda matrix & A matrix
	mjtNum* _qM = mj_stackAlloc(d, 6 * 6);
	mj_fullM(m, _qM, d->qM);
	MatrixXd A_matrix(6, 6);
	A_matrix.resize(6, 6);
	A_matrix.setZero();
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++) {
			A_matrix(i, j) = _qM[i * 6 + j];
		}
	}
	MatrixXd lambda, lambda_inv;
	lambda_inv = J_ee*A_matrix.inverse() *J_ee.transpose();
	lambda = lambda_inv.inverse();

	// end effector velocity, Twist calc
	Vector3d x_dot, xa_dot;
	VectorXd qdot(6);
	for (int i = 0; i < 6; i++) {
		qdot(i) = d->qvel[i];
	}

	x_dot = J_ee_p*qdot;
	xa_dot = J_ee_r*qdot;



	MatrixXd lambda_f, lambda_f_inv;
	MatrixXd J_task(5, 6);
	J_task.block<2, 6>(0, 0) = J_ee.block<2, 6>(0, 0);
	J_task.block<3, 6>(2, 0) = J_ee.block<3, 6>(3, 0);

	lambda_f_inv = J_task*A_matrix.inverse() *J_task.transpose();
	lambda_f = lambda_f_inv.inverse();




	//Get F_star
	VectorXd F_star(6);
	VectorXd torque(6);
	torque.setZero();
	F_star.setZero();



	*/

	for (static bool first = true; first; first = false)
	{
		cout << "Simulation Start" << endl;
		// some code that executes only once
	}



	//1. Base frame Rotation matrix Setting // 
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			Rotation_base_frame(i, j) = d->xmat[mj_name2id(m, mjOBJ_BODY, "Pelvis") * 9 + i*3+ j];
		}
	}	
	Rotation2g.block(3, 3, 3, 3) = Rotation_base_frame.transpose();

	//------------------------------------------//





	//2. Part data refresh //
	for (int i = 0; i < 6; i++) {
		pR[i].refresh();
		pL[i].refresh();
	}
	for (int i = 0; i < 7; i++) {
		pUR[i].refresh();
		pUL[i].refresh();
	}
	for (int i = 0; i < 3; i++) {
		pBase[i].refresh();
	}
	//------------------------------------------//
	


	//3. Mass matrix setting // 
	mjtNum* _qM = mj_stackAlloc(d, m->nv * m->nv);
	mj_fullM(m, _qM, d->qM);
	MatrixXd A_matrix;
	A_matrix.setZero(m->nv, m->nv);
	for (int i = 0; i < m->nv; i++) {
		for (int j = 0; j < m->nv; j++) {
			A_matrix(i, j) = _qM[i * 6 + j];
		}
	}

	A_matrix = Rotation2g.transpose() *A_matrix * Rotation2g;

	MatrixXd A1[6],A2[7],A3[3];
	
	
	



	for (int i = 0; i < 6; i++) {
		A1[i].setZero(m->nv, m->nv);
		A1[i] = pR[i].Mass*(pR[i].Jac_COM.block(0, 0, 3, m->nv)).transpose()*pR[i].Jac_COM.block(0, 0, 3, m->nv)
			+ pR[i].Jac_COM.block(3, 0, 3, m->nv).transpose()*pR[i].inertia*pR[i].Jac_COM.block(3, 0, 3, m->nv)
			+ pL[i].Mass*(pL[i].Jac_COM.block(0, 0, 3, m->nv)).transpose()*pL[i].Jac_COM.block(0, 0, 3, m->nv)
			+ pL[i].Jac_COM.block(3, 0, 3, m->nv).transpose()*pL[i].inertia*pL[i].Jac_COM.block(3, 0, 3, m->nv);

	}
	for (int i = 0; i < 7; i++) {
		A2[i].setZero(m->nv, m->nv);
		A2[i] = pUR[i].Mass*(pUR[i].Jac_COM.block(0, 0, 3, m->nv)).transpose()*pUR[i].Jac_COM.block(0, 0, 3, m->nv)
			+ pUR[i].Jac_COM.block(3, 0, 3, m->nv).transpose()*pUR[i].inertia*pUR[i].Jac_COM.block(3, 0, 3, m->nv)
			+ pUL[i].Mass*(pUL[i].Jac_COM.block(0, 0, 3, m->nv)).transpose()*pUL[i].Jac_COM.block(0, 0, 3, m->nv)
			+ pUL[i].Jac_COM.block(3, 0, 3, m->nv).transpose()*pUL[i].inertia*pUL[i].Jac_COM.block(3, 0, 3, m->nv);
	}	
	for (int i = 0; i < 3; i++) {
		A3[i].setZero(m->nv, m->nv);
		A3[i] = pBase[i].Mass*(pBase[i].Jac_COM.block(0, 0, 3, m->nv)).transpose()*pBase[i].Jac_COM.block(0, 0, 3, m->nv)
			+ pBase[i].Jac_COM.block(3, 0, 3, m->nv).transpose()*pBase[i].inertia*pBase[i].Jac_COM.block(3, 0, 3, m->nv);
	}

	MatrixXd A;
	A.setZero(m->nv, m->nv);
	A = A1[0] + A1[1] + A1[2] + A1[3] + A1[4] + A1[5]
		+ A2[0] + A2[1] + A2[2] + A2[3] + A2[4] + A2[5] + A2[6]
		+ A3[0] + A3[1] + A3[2];
	A = A*0.5;



	A_matrix = A;


	//------------------------------------------//
	





	VectorXd qvel_eigen;
	qvel_eigen.resize(m->nv);
	qvel_eigen = mj2eigen(d->qvel, m->nv, 1);
	VectorXd mystery(6);
	mystery = pR[5].Jac * qvel_eigen;








	//data_print_data -> 시뮬레이션 창에 띄우는 데이터. 
	sprintf_s(data_print_data, "%8.3f%8.3f%8.3f\n%8.3f%8.3f%8.3f\n%8.3f%8.3f%8.3f\n%8.3f%8.3f%8.3f\n%8.3f", mystery(0), mystery(1), mystery(2), mystery(3), mystery(4), mystery(5), d->qvel[0], d->qvel[1], d->qvel[2], d->qvel[3], d->qvel[4], d->qvel[5], A.determinant());
	
	Vector3d Grav_ref;
	Grav_ref << 0, 0, -9.81;
	VectorXd G;
	G.setZero(m->nv);

	for (int i = 0; i < 6; i++) {
		G = G + pR[i].Mass*pR[i].Jac_COM.block(0, 0, 3, m->nv).transpose()*Grav_ref
		 + pL[i].Mass*pL[i].Jac_COM.block(0, 0, 3, m->nv).transpose()*Grav_ref;
	}
	for (int i = 0; i < 7; i++) {
		G = G + pUR[i].Mass*pUR[i].Jac_COM.block(0, 0, 3, m->nv).transpose()*Grav_ref
		+ pUL[i].Mass*pUL[i].Jac_COM.block(0, 0, 3, m->nv).transpose()*Grav_ref;
	}
	for (int i = 0; i < 3; i++) {
		G = G + pBase[i].Mass*pBase[i].Jac_COM.block(0, 0, 3, m->nv).transpose()*Grav_ref;
	}

	MatrixXd J_g;
	J_g.setZero(m->nu, m->nv);
	J_g.block(0, 6, m->nu, m->nu).setIdentity();
	


	MatrixXd J_C, J_C_INV_T;
	J_C.setZero(12, m->nv);
	J_C.block(0, 0, 6, m->nv) = pR[5].Jac_Contact;
	J_C.block(6, 0, 6, m->nv) = pL[5].Jac_Contact;


	MatrixXd Lambda_c;
	Lambda_c=(J_C*A_matrix*(J_C.transpose())).inverse();

	J_C_INV_T = Lambda_c*J_C*A_matrix.inverse();


	

	MatrixXd N_C;
	N_C.setIdentity(m->nv, m->nv);
	N_C =N_C-J_C.transpose()*J_C_INV_T;



	VectorXd torque_grav(m->nu);
	
	

	torque_grav = (J_g*A_matrix.inverse()*N_C*J_g.transpose()).inverse()*J_g*A_matrix.inverse()*N_C*G;

	if(d->time>1.0)
	{
		for (static bool first = true; first; first = false)
		{
			

			// some code that executes only once
		}
		
		// some code that executes only once
	}


	//Control + R, F, T 조합으로 키 명령 
	if (Control_R) {
		if (Control_init) {
			Control_init = !Control_init;
		}
		torque = torque_grav;
		if (d->time > run_command_time + 1.0) {
			Control_R = false;
		}
	}
	else if (Control_F) {
		if (Control_init) {
			Control_init = !Control_init;

			cout << "gravcheck" << endl;
			cout << "grav_matrix :" << endl;
			cout << G << endl;
			cout << " J_g : " << endl;
			cout << "A :" << endl;
			cout << A << endl;
			cout << "A_inv :" << endl;
			cout << A.inverse() << endl;
			cout << J_g << endl;
			cout << " J_C" << endl;
			cout << J_C << endl;
			cout << "J_C_inv_T" << endl;
			cout << J_C_INV_T << endl;
			cout << "lambda_c" << endl;
			cout << Lambda_c << endl;
			cout << "N_C" << endl;
			cout << N_C << endl;
			cout << "torque grav" << endl;
			cout << torque_grav;
			cout << "bonus" << endl;
			cout << (J_g*A_matrix.inverse()*N_C*J_g.transpose()).inverse()*J_g*A_matrix.inverse()*N_C << endl;


		//	cout << "A_matrix" << endl << A_matrix << endl;
			//cout << "A" << endl << A << endl;






			Control_F = !Control_F;





		}
	}
	else if (Control_T) {

		if (Control_init) {
			Control_init = !Control_init;

			std::cout << pR[5].Jac << std::endl << std::endl;
			std::cout << pR[5].xpos << std::endl;
			std::cout << pBase[0].xpos << std::endl << std::endl;
			std::cout << pR[5].xpos - pBase[0].xpos << std::endl;
			cout << "phat :: " << endl;
			std::cout << Skm(pR[5].xpos - pBase[0].xpos) << endl;


			cout << "pur jac" << endl;
			cout << pUR[6].Jac << endl;
			cout << pUR[6].Mass << endl;
			cout << "urhat " << endl;
			cout << Skm(pUR[6].xpos - pBase[0].xpos) << endl;

			std::cout << " target rotMatrix ::" << std::endl;
			std::cout << pR[5].Rotm << std::endl;
			std::cout << "target rotmatx inv:: " << std::endl;
			std::cout << pR[5].Rotm.inverse() << std::endl;
			std::cout << "pelvis rotm::" << std::endl;
			std::cout << pBase[0].Rotm << std::endl;
			std::cout << "Pelvis jac :: " << std::endl;
			std::cout << pBase[0].Jac << std::endl << std::endl;
			std::cout << "pelvis velocity" << std::endl;
			std::cout << mystery << std::endl << endl;
			std::cout << "rotation2g" << endl;
			cout << Rotation2g << endl;



			for (int i = 0; i < 6; i++) {
				A1[i].setZero(m->nv, m->nv);
				A1[i] = pR[i].Mass*(pR[i].Jac_COM.block(0, 0, 3, m->nv)).transpose()*pR[i].Jac_COM.block(0, 0, 3, m->nv)
					+ pR[i].Jac_COM.block(3, 0, 3, m->nv).transpose()*pR[i].inertia*pR[i].Jac_COM.block(3, 0, 3, m->nv)
					+ pL[i].Mass*(pL[i].Jac_COM.block(0, 0, 3, m->nv)).transpose()*pL[i].Jac_COM.block(0, 0, 3, m->nv)
					+ pL[i].Jac_COM.block(3, 0, 3, m->nv).transpose()*pL[i].inertia*pL[i].Jac_COM.block(3, 0, 3, m->nv);



				cout << endl << "A1[" << i << "] : " << endl << A1[i] << endl << endl << endl;
			}
			for (int i = 0; i < 7; i++) {
				A2[i].setZero(m->nv, m->nv);
				A2[i] = pUR[i].Mass*(pUR[i].Jac_COM.block(0, 0, 3, m->nv)).transpose()*pUR[i].Jac_COM.block(0, 0, 3, m->nv)
					+ pUR[i].Jac_COM.block(3, 0, 3, m->nv).transpose()*pUR[i].inertia*pUR[i].Jac_COM.block(3, 0, 3, m->nv)
					+ pUL[i].Mass*(pUL[i].Jac_COM.block(0, 0, 3, m->nv)).transpose()*pUL[i].Jac_COM.block(0, 0, 3, m->nv)
					+ pUL[i].Jac_COM.block(3, 0, 3, m->nv).transpose()*pUL[i].inertia*pUL[i].Jac_COM.block(3, 0, 3, m->nv);
				cout << endl << "A2[" << i << "] : " << endl << A2[i] << endl << endl << endl;
			}
			for (int i = 0; i < 3; i++) {
				A3[i].setZero(m->nv, m->nv);
				A3[i] = pBase[i].Mass*(pBase[i].Jac_COM.block(0, 0, 3, m->nv)).transpose()*pBase[i].Jac_COM.block(0, 0, 3, m->nv)
					+ pBase[i].Jac_COM.block(3, 0, 3, m->nv).transpose()*pBase[i].inertia*pBase[i].Jac_COM.block(3, 0, 3, m->nv);
				cout << endl << "A3[" << i << "] : " << endl << A3[i] << endl << endl << endl;
			}



			
			cout << "//-----------------------------detail test--------------------------------//" << endl;
			int h;
			cin >> h;
			A1[h].setZero(m->nv, m->nv);
			A1[h] = pR[h].Mass*(pR[h].Jac_COM.block(0, 0, 3, m->nv)).transpose()*pR[h].Jac_COM.block(0, 0, 3, m->nv)
				+ pR[h].Jac_COM.block(3, 0, 3, m->nv).transpose()*pR[h].inertia*pR[h].Jac_COM.block(3, 0, 3, m->nv)
				+ pL[h].Mass*(pL[h].Jac_COM.block(0, 0, 3, m->nv)).transpose()*pL[h].Jac_COM.block(0, 0, 3, m->nv)
				+ pL[h].Jac_COM.block(3, 0, 3, m->nv).transpose()*pL[h].inertia*pL[h].Jac_COM.block(3, 0, 3, m->nv);

			
			cout << endl << "R  Mass : " <<pR[h].Mass << endl;
			cout << endl << "all jacobian :" << endl;
			cout << pR[h].Jac_COM << endl;
			cout << endl << "R jac pos :" << endl;
			cout << pR[h].Jac_COM.block(0, 0, 3, m->nv) << endl;
			cout << endl << "R jac rot :" << endl;
			cout << pR[h].Jac_COM.block(3, 0, 3, m->nv) << endl;
			cout << endl << "R inertia :" << endl;
			cout << pR[h].inertia << endl;


			cout << endl << "A1[" << h << "] : " << endl << A1[h] << endl << endl << endl;

			cout << "//-----------------------------END-------------------------------- //" << endl;
		


		}
		Control_T = !Control_T;


	}

	torque = torque_grav;


	//Torque to mujoco
	torque_mj = eigen2mj(torque, d);
	mju_copy(d->ctrl,torque_mj, m->nu);	
	mjFREESTACK
}